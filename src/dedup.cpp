/*
		This program is free software; you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program; if not, write to the Free Software
		Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

		Original implementation by Donald Graft.
		Optimizations and improvements by Klaus Post.
		Variable framerate modifications by Loren Merritt.
*/

#include "internal.h"
#include "info.h"

#define BLKSIZE 32
#define NORM (235*BLKSIZE*BLKSIZE)
#define MAX_COPIES 20
#define DUPVERSION "2.20 beta 1"
#define MYVERSION "0.17"
#define VERSION_PRINTF "DeDup %s by Loren Merritt, based on Dup %s by Donald Graft/Klaus Post, Copyright 2004\n", MYVERSION, DUPVERSION
#define SHOW_NEIGHBORS 10 // show the status of this many frames in each direction

struct FRAMEINFO
{
	unsigned int frame_no;
	PVideoFrame frame;
	unsigned int pitch, pitchY, pitchUV;
	unsigned char *frame_ptr;
	unsigned char *frame_ptrY;
	unsigned char *frame_ptrU;
	unsigned char *frame_ptrV;
	double metric;
	int highest_x;
	int highest_y;
};

struct FRAMEMETRIC
{
	float metric;
	int highest_x;
	int highest_y;
};

enum DROP
{
	DROPPED,
	KEPT_MAXDROPS,
	KEPT_MAXCOPIES,
	KEPT_THRESH
};
static const char* DROP_NAMES = "_dCK";

enum DECWHICH
{
	DECW_FIRST=0,
	DECW_MIDDLE_DOWN=1,
	DECW_MIDDLE_UP=2,
	DECW_LAST=3,
	DECW_BLEND=4,
};

class Dup : public GenericVideoFilter
{
public:
	// In Dedup mode
	Dup(PClip _child, float _threshold, float _threshold2, int _range2,
		float _trigger2, bool _show, bool _dec, int _maxcopies, int _maxdrops,
		int _decwhich, const char* _logfile, const char* _timefile,
		const char* _timeinfile, const char* _ovrfile, const char* _debugfile,
		IScriptEnvironment* _env) :
		GenericVideoFilter(_child), threshold(_threshold), threshold2(_threshold2),
		range2(_range2), trigger2(_trigger2), show(_show), dec(_dec),
		maxcopies(_maxcopies), maxdrops(_maxdrops), decwhich(_decwhich),
                env(_env), pass(2)
	{
		if (!vi.IsYUY2() && !vi.IsYV12())
			env->ThrowError("DeDup: requires YUY2 or YV12 source");
		if (maxcopies > 20)
			env->ThrowError("DeDup: maxcopies must be <= 20");
		if (maxdrops > maxcopies)
			env->ThrowError("DeDup: maxdrops must be <= maxcopies");
		if (threshold < 0.0 || threshold > 100.0 || threshold2 < 0.0 || threshold2 > 100.0)
			env->ThrowError("DeDup: threshold out of range (0.0-100.0)");
		if (trigger2 <= threshold2 || trigger2 <= threshold)
			env->ThrowError("DeDup: trigger2 must be > threshold");
		if (_logfile)
		{
			logfile = fopen(_logfile, "r");
			if (logfile == NULL)
				env->ThrowError("DeDup: failed to open logfile");
		}
		else
			logfile = NULL;
		if (_timefile)
		{
			timefile = fopen(_timefile, "w");
			if (timefile == NULL)
				env->ThrowError("DeDup: failed to open timefile");
		}
		else
			timefile = NULL;
		if (_timeinfile)
		{
			timeinfile = fopen(_timeinfile, "r");
			if (timeinfile == NULL)
				env->ThrowError("DeDup: failed to open timeinfile");
		}
		else
			timeinfile = NULL;
		if (_ovrfile)
		{
			ovrfile = fopen(_ovrfile, "r");
			if (ovrfile == NULL)
				env->ThrowError("DeDup: failed to open ovrfile");
		}
		else
			ovrfile = NULL;
		if (_debugfile)
		{
			debugfile = fopen(_debugfile, "w");
			if (debugfile == NULL)
				env->ThrowError("DeDup: failed to open debugfile");
		}
		else
			debugfile = NULL;

		if (env->GetCPUFlags() & CPUF_INTEGER_SSE) have_isse = true;
		else have_isse = false;
		copyframe = env->NewVideoFrame(vi);
		cache_out.frame_no = -1;

		xblocks = (vi.width+BLKSIZE-1) / BLKSIZE;
		yblocks = (vi.height+BLKSIZE-1) / BLKSIZE;

		num_iframes = vi.num_frames;
		sum = NULL;

		/* For safety in case someone came in without doing it. */
		__asm emms;

		LoadFirstPass();
	}
	// In metric collection mode
	Dup(PClip _child, bool _chroma, const char* _logfile, int _search, IScriptEnvironment* _env) :
	    GenericVideoFilter(_child), chroma(_chroma), range2(_search), env(_env), pass(1)
	{
		if (!vi.IsYUY2() && !vi.IsYV12())
			env->ThrowError("DeDup: requires YUY2 or YV12 source");
		if (_logfile)
		{
			logfile = fopen(_logfile, "w");
			if (logfile == NULL)
				env->ThrowError("DeDup: failed to open logfile");
			fprintf(logfile, VERSION_PRINTF);
		}
		else
			env->ThrowError("DeDup: option 'log' required");
		if (_search < 1)
			env->ThrowError("DeDup: 'search' must be >= 1");
		//// Will enable when I get around to caching
		//if (_search > MAX_COPIES)
		//	env->ThrowError("DeDup: 'search' must be <= %d", MAX_COPIES);

		if (env->GetCPUFlags() & CPUF_INTEGER_SSE) have_isse = true;
		else have_isse = false;

		xblocks = (vi.width+BLKSIZE-1) / BLKSIZE;
		yblocks = (vi.height+BLKSIZE-1) / BLKSIZE;

		num_iframes = vi.num_frames;
		sum = new unsigned int [xblocks * yblocks];
		metrics_done = new bool [num_iframes];
		if (sum == NULL || metrics_done == NULL)
			env->ThrowError("DeDup: cannot allocate needed memory");
		memset(metrics_done, 0, num_iframes * sizeof(*metrics_done));
		metrics = NULL;
		keep = NULL;
		mapend = mapstart = mapinv = NULL;

		/* For safety in case someone came in without doing it. */
		__asm emms;
	}
	~Dup()
	{
		if (sum)          delete [] sum;
		if (metrics)      delete [] metrics;
		if (metrics_done) delete [] metrics_done;
		if (keep)         delete [] keep;
		if (mapend)       delete [] mapend;
		if (mapstart)     delete [] mapstart;
		if (mapinv)       delete [] mapinv;
	}
	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* _env);
	void isse_scenechange(const BYTE* c_plane, const BYTE* tplane, int height, int width, int pitch,  int t_pitch, int* blk_values);
	void isse_scenechange_16(const BYTE* c_plane, const BYTE* tplane, int height, int width, int pitch,  int t_pitch, int* blk_values);
	void mmx_average_planes(BYTE* dst_plane, const BYTE** src_planes, int width_mod8, int planes, int div);

private:
	//// Options
	float threshold, threshold2;
	int range2;
	float trigger2;
	bool show, dec, chroma;
	int decwhich;      // which of a run of dups to keep
	int maxcopies;     // max consecutive frames to merge (blend / copy)
	int maxdrops;      // max consecutive frames to drop
	int pass;          // 1 => collect metrics, 2 => decimate
	FILE* logfile;     // load/save stats for 2 passes
	FILE* timeinfile;  // load matroska timecodes (vfr input)
	FILE* timefile;    // save matroska timecodes
	FILE* ovrfile;     // override calculated dups (or just about any other option)
	FILE* debugfile;   // report decisions here

	//// State
	bool have_isse;
	struct FRAMEINFO cache[MAX_COPIES+1];
	struct FRAMEINFO cache_out;
	int cache_count;
	PVideoFrame copyframe;
	IScriptEnvironment* env;

	int xblocks, yblocks;
	unsigned int *sum, highest_sum;

	FRAMEMETRIC* metrics; // metrics[i] = diff(i, i+1)
	bool* metrics_done;   // which metrics have we calculated yet?
	char* keep;           // keep[i] = (0 => drop, 1 => output due to minframerate, 2 => keep)
	unsigned int* mapend; // mapend[out_frame] = in_frame (end of run)
	unsigned int* mapstart; // start of run
	unsigned int* mapinv; // mapinv[in_frame] = out_frame
	unsigned int num_iframes, num_oframes; // total frames before and after dropping dups

	int row_size, row_sizeY, row_sizeUV;
	int height, heightY, heightUV;

	//// Methods
	void LoadFirstPass();
	void CalcMetric(int n);
	PVideoFrame ConstructFrame(int n);
	void CompareFrames(FRAMEINFO& frame1, FRAMEINFO& frame0);
	FRAMEINFO LoadFrame(int n);
	void CopyFrame(PVideoFrame& to, PVideoFrame& from);
	void DrawBox(PVideoFrame& frame, int box_x, int box_y, bool crossp);

	inline void DrawString(PVideoFrame& dst, int x, int y, const char *s);
	void ThrowError(const char* msg);
	template<class T> void ThrowError(const char* msg, T arg);
};

void Dup::DrawString(PVideoFrame& dst, int x, int y, const char *s)
{
	if (vi.IsYUY2()) ::DrawStringYUY2(dst, x, y, s);
	else ::DrawString(dst, x, y, s);
}

void Dup::ThrowError(const char* msg)
{
	env->ThrowError(msg);
}

template<class T>
void Dup::ThrowError(const char* msg, T arg)
{
	// unsafe: I would use snprintf, but MSVC doesn't have it.
	// no, I don't ever free buf
	char* buf = new char[strlen(msg) + 20];
	sprintf(buf, msg, arg);
	env->ThrowError(buf);
}

PVideoFrame __stdcall Dup::GetFrame(int n, IScriptEnvironment* _env)
{
	env = _env;
	if(pass == 1)
	{
		CalcMetric(n);
		return child->GetFrame(n, env);
	}
	else
	{
		return ConstructFrame(n);
	}
}

void Dup::LoadFirstPass()
{
	float metric;
	unsigned int i, j;
	int highest_x, highest_y;
	unsigned int frame_no, frame_next;
	char buf[201];

	float* thresholds = new float [num_iframes];
	float* threshold2s = new float [num_iframes];
	float* trigger2s = new float [num_iframes];
	int* range2s = new int [num_iframes];
	int* maxdropss = new int [num_iframes];
	int* maxcopiess = new int [num_iframes];

	metrics = new FRAMEMETRIC [num_iframes];
	metrics_done = new bool [num_iframes];
	keep = new char [num_iframes];
	mapend = new unsigned int [num_iframes];
	mapstart = new unsigned int [num_iframes];
	mapinv = new unsigned int [num_iframes];

	if(metrics == NULL || keep == NULL || mapend == NULL || mapstart == NULL
			|| mapinv == NULL || metrics_done == NULL || thresholds == NULL
			|| threshold2s == NULL || trigger2s == NULL || range2s == NULL
			|| maxdropss == NULL || maxcopiess == NULL)
		env->ThrowError("DeDup: cannot allocate needed memory");

	memset(metrics_done, 0, num_iframes * sizeof(bool));
	memset(keep, 0, num_iframes * sizeof(char));
	for(i=0; i<num_iframes; i++)
	{
		thresholds[i] = threshold;
		threshold2s[i] = threshold2;
		trigger2s[i] = trigger2;
		range2s[i] = range2;
		maxdropss[i] = maxdrops;
		maxcopiess[i] = maxcopies;
	}

	// read stats from first pass
	if(logfile != NULL)
	{
		while(fgets(buf, 200, logfile))
		{
			if(!sscanf(buf, "frm %u: diff from frm %u = %f%% at (%d,%d)", &frame_no, &frame_next, &metric, &highest_x, &highest_y))
				continue;
			if(frame_no >= num_iframes)
				continue;
			if(frame_next != frame_no+1)
				// the data might be useful, but I don't know what to do with it
				continue;
			metrics[frame_no].metric = metric;
			metrics[frame_no].highest_x = highest_x;
			metrics[frame_no].highest_y = highest_y;
			metrics_done[frame_no] = true;
		}
		fclose(logfile);
	}

	// read override file
	if(ovrfile != NULL)
	{
		char* pos;
		int linenum = 0;
		while(fgets(buf, 200, ovrfile))
		{
			linenum++;
			pos = buf + strspn(buf, " \t");
			if(pos[0] == '#' || pos[0] == '\0')
				continue;

			unsigned int frame0, frame1;
			if(sscanf(pos, " %u-%u", &frame0, &frame1) == 2)
			{
				if(frame0 > frame1)
					ThrowError("override section of negative length at line %d", linenum);
			}
			else if(sscanf(buf, " %u", &frame0))
				frame1 = frame0;
			else
				ThrowError("can't parse override at line %d", linenum);

			pos += strcspn(pos, " \t");
			pos += strspn(pos, " \t");
			if(pos[0] == '\0')
				ThrowError("can't parse override at line %d", linenum);

			while(pos[0] != '\0')
			{
				if(pos[0] == '#')
					break;
				if(sscanf(pos, "threshold=%f", &thresholds[frame0])
				   || sscanf(pos, "thresh=%f", &thresholds[frame0]))
					for(i=frame0+1; i<=frame1 && i < num_iframes; i++)
					{
						thresholds[i] = thresholds[frame0];
						if(thresholds[i] > threshold2s[i])
							threshold2s[i] = thresholds[i];
					}
				else if(sscanf(pos, "threshold2=%f", &threshold2s[frame0])
				        || sscanf(pos, "thresh2=%f", &threshold2s[frame0]))
					for(i=frame0+1; i<=frame1 && i < num_iframes; i++)
					{
						threshold2s[i] = threshold2s[frame0];
						if(range2s[i] <= 0)
							range2s[i] = range2;
						if(thresholds[i] > threshold2s[i])
							thresholds[i] = threshold2s[i];
					}
				else if(sscanf(pos, "trigger2=%f", &trigger2s[frame0]))
					for(i=frame0+1; i<=frame1 && i < num_iframes; i++)
					{
						trigger2s[i] = trigger2s[frame0];
						if(range2s[i] <= 0)
							range2s[i] = range2;
					}
				else if(sscanf(pos, "range2=%d", &range2s[frame0]))
				{
					if(range2s[frame0] < 0)
						ThrowError("range2 < 0 at line %d", linenum);
					for(i=frame0+1; i<=frame1 && i < num_iframes; i++)
						range2s[i] = range2s[frame0];
				}
				else if(sscanf(pos, "maxdrops=%d", &maxdropss[frame0]))
				{
					if(maxdropss[frame0] < 0)
						ThrowError("maxdrops < 0 at line %d", linenum);
					for(i=frame0+1; i<=frame1 && i < num_iframes; i++)
						maxdropss[i] = maxdropss[frame0];
				}
				else if(sscanf(pos, "maxcopies=%d", &maxcopiess[frame0]))
				{
					if(maxcopiess[frame0] < 0)
						ThrowError("maxcopies < 0 at line %d", linenum);
					for(i=frame0+1; i<=frame1 && i < num_iframes; i++)
						maxcopiess[i] = maxcopiess[frame0];
				}
				// word of the form /[kd]+/ forces keep or drop
				// FIXME: doesn't check for non-kd letters
				else if(pos[0] == 'k' || pos[0] == 'd')
					for(j=0, i=frame0; i<=frame1 && i < num_iframes; i++, j++)
					{
						if(pos[j] != 'k' && pos[j] != 'd')
							j = 0;
						if(pos[j] == 'k')
						{
							thresholds[i] = threshold2s[i] = -1.0f;
							range2s[i] = 0;
						}
						else
						{
							thresholds[i] = threshold2s[i] = 101.0f;
							range2s[i] = 0;
						}
						metrics_done[i] = true;
					}
				else
					ThrowError("unrecognized override at line %d", linenum);

				pos += strcspn(pos, " \t");
				pos += strspn(pos, " \t");
			} // End while(pos[0] != '\0')
		} // End while(fgets(buf, 200, ovrfile))
		fclose(ovrfile);
	} // End if(ovrfile)

	// warn of any holes in the first pass not filled by overrides
	for(i=0; i<num_iframes; i++)
		if(!metrics_done[i])
		{
			if(logfile == NULL && ovrfile == NULL)
				env->ThrowError("DeDup: no first pass logfile, nor ovrfile");
			else if(!logfile)
				env->ThrowError("DeDup: no first pass logfile, and not all frames overridden");
			else
				ThrowError("DeDup: first pass missed frame %d", i);
		}
	delete [] metrics_done;
	metrics_done = NULL;
	logfile = NULL;
	ovrfile = NULL;

	// decide which frames to drop
	num_oframes = 0;
	int numdropped = 0;
	int start = 0;
	for(i=0; i<num_iframes; i++)
	{
		numdropped++;

		// dual threshold, for high motion / low fps
		float& thresh = thresholds[i];
		float thresh2 = threshold2s ? threshold2s[i] : threshold2;
		float trig2 = trigger2s ? trigger2s[i] : trigger2;
		int r2 = range2s ? range2s[i] : range2;
		int maxcopy = maxcopiess ? maxcopiess[i] : maxcopies;
		int maxdrop = maxdropss ? maxdropss[i] : maxdrops;

		if(r2 > 0)
		{
			float neighbor0 = 0., neighbor1 = 0.;
			for(j = (int)i<r2 ? 0 : i-r2; j<i; j++)
				if(metrics[j].metric > neighbor0)
					neighbor0 = metrics[j].metric;
			for(j = i+r2 > num_iframes-1 ? num_iframes-1 : i+r2; j>i; j--)
				if(metrics[j].metric > neighbor1)
					neighbor1 = metrics[j].metric;
			if(neighbor1 < neighbor0)
				neighbor0 = neighbor1;
			if(neighbor0 > thresh2)
			{
				float ratio = (neighbor0 - threshold2s[i]) / (trig2 - thresh2);
				ratio = ratio < 0 ? 0 : ratio > 1 ? 1 : ratio;
				thresh = thresh2 * ratio + thresh * (1-ratio);
			}
		}

		int status = DROPPED;
		if(metrics[i].metric >= thresh)
			status = KEPT_THRESH;
		else if(numdropped >= maxcopy)
			status = KEPT_MAXCOPIES;

		if(status != DROPPED)
		{
			int j;
			int numsplits = 1 + (numdropped-1) / maxdrops;
			for(j=0; j<numsplits; j++)
			{
				mapend[num_oframes] = i;
				mapstart[num_oframes] = start;
				num_oframes++;
			}
			for(j=1; j<numsplits; j++)
			{
				int k = i - numdropped + (j * numdropped) / numsplits;
				keep[k] = KEPT_MAXDROPS;
			}
			keep[i] = status;

			numdropped = 0;
			start = i+1;
		}
	} // End for(i=0; i<num_iframes; i++)
	if(dec)
		vi.num_frames = num_oframes;

	// calc i<->o mappings ignoring maxdrops
	for(i=j=0; i<num_iframes; i++)
	{
		mapinv[i] = j;
		if(keep[i])
			j++;
	}

	// calc timestamps
	double* timestamps = new double [num_iframes+1]; // beginning of the frame, in ms
	if(timeinfile)
	{
		unsigned int frame0, frame1;
		double globalfps = -1, rangefps;
		unsigned int i = 0, linenum = 1;
		double time = 0.;
		fgets(buf, 200, timeinfile);
		if(strncmp(buf, "# timecode format v1", 20) == 0)
		{
			while(fgets(buf, 200, timeinfile))
			{
				linenum++;
				if(buf[0] == '#' || buf[strspn(buf, " \t\r\n")] == '\0')
					continue;
				if(sscanf(buf, "assume %lf", &globalfps)
				|| sscanf(buf, "Assume %lf", &globalfps))
					continue;
				if(globalfps <= 0)
					env->ThrowError("DeDup: timeinfile: no 'assume fps' line");

				if(sscanf(buf, "%u,%u,%lf", &frame0, &frame1, &rangefps) < 3)
					ThrowError("DeDup: timeinfile: can't parse line %u", linenum);
				if(frame1 < frame0 || rangefps <= 0)
					ThrowError("DeDup: timeinfile: inconsistent data at line %u", linenum);
				for(; i < frame0; i++)
				{
					timestamps[i] = time;
					time += 1000. / globalfps;
				}
				for(; i <= frame1; i++)
				{
					timestamps[i] = time;
					time += 1000. / rangefps;
				}
			}
			for(; i < num_iframes; i++)
			{
				timestamps[i] = time;
				time += 1000. / globalfps;
			}
		}
		else if(strncmp(buf, "# timecode format v2", 20) == 0)
		{
			for(i = 0; i < num_iframes; i++)
			{
				if(!fgets(buf, 200, timeinfile))
					env->ThrowError("DeDup: timeinfile doesn't contain enough timestamps");
				if(buf[0] == '#' || buf[strspn(buf, " \t\r\n")] == '\0')
					continue;
				if(sscanf(buf, "%lf", &timestamps[i]) < 1)
					env->ThrowError("DeDup: timeinfile not a valid Matroska timecode v2 file");
			}
		}
	} // End if(timeinfile)
	else
	{
		double mspf = 1000. * double(vi.fps_denominator) / vi.fps_numerator;
		for(i=0; i<num_iframes; i++)
			timestamps[i] = i * mspf;
	}

	if(timefile != NULL)
	{
		// print timestamps
		fprintf(timefile, "# timecode format v2\n");
		fprintf(timefile, "%.6lf\n", 0.);
		for(i=0; i<num_iframes-1; i++)
			if(keep[i])
				fprintf(timefile, "%.6lf\n", timestamps[i+1]);
		fclose(timefile);
		timefile = NULL;
	}

	if(debugfile)
	{
		fprintf(debugfile, VERSION_PRINTF);
		fprintf(debugfile, "\ninframe (outframe) @ time: kept?   metric <> threshold\n");
		for(i=0; i<num_iframes; i++)
			fprintf(debugfile, "%d (%d) @ %.2fs: %d   %7.4f%% %c %.4f%%\n",
					i, mapinv[i], timestamps[i] / 1000., keep[i], metrics[i].metric,
					keep[i] == KEPT_THRESH ? '>' : '<', thresholds[i]);

		fclose(debugfile);
		debugfile = NULL;
	}

	delete [] thresholds;
	delete [] threshold2s;
	delete [] trigger2s;
	delete [] range2s;
	delete [] maxdropss;
	delete [] maxcopiess;
} // End Dup::LoadFirstPass()

void Dup::CalcMetric(int n)
{
	if(metrics_done[n])
		return;

	// FIXME: save frames (don't rely on Avisynth's cache for large lookahead)
	cache[0] = LoadFrame(n);
	for(int skip=1; skip<=range2; skip++)
	{
		if(n+skip < (int)num_iframes)
		{
			// compare frame n to n+skip
			cache[0].metric = 0.0;
			cache[1] = LoadFrame(n+skip);
			cache[1].metric = 0.0;
			CompareFrames(cache[0], cache[1]);
		}
		else
		{
			// always keep the last frame
			cache[0].metric = 100.0;
			cache[0].highest_x = 0;
			cache[0].highest_y = 0;
		}

		fprintf(logfile, "frm %d: diff from frm %d = %2.4f%% at (%d,%d)\n",
			n, n+skip, cache[0].metric, cache[0].highest_x, cache[0].highest_y);
	}
	fflush(logfile);
	metrics_done[n] = 1;
}

PVideoFrame Dup::ConstructFrame(int n)
{
	char buf[80];
	int offset_remainX = (vi.width&(~(BLKSIZE-1)));  // Offset into the frame (pixels)
	int offset_remainY = (vi.height/BLKSIZE)*BLKSIZE;  // yposition the remaining pixels start (lines)
	int remainX = vi.width&(BLKSIZE-1) + offset_remainX;      // Where do the remaining pixels end? (pixels)
	int remainY = (vi.height&(BLKSIZE-1)) + offset_remainY;   // Where do the remaining pixels end? (lines)
	int i;
	int n0; // start of the run of duplicates
	int n1; // end of the run of duplicates
	int ni; // dec ? n1 : input frame number
	int nk; // frame kept

	if (dec)
	{
		n0 = mapstart[n];
		ni = n1 = mapend[n];
		/* Modify copyframe to be a blend of all the frames in the
		   string of duplicates. Skip if there's only one frame */
		if (decwhich==DECW_BLEND && n0 < n1)
		{
			int i;
			int dpitchY, dpitchUV, row_sizeY2, heightY2;
			BYTE* dst_planeY;
			BYTE* dst_planeU;
			BYTE* dst_planeV;
			const BYTE** src_planesY=0;
			const BYTE** src_planesU=0;
			const BYTE** src_planesV=0;
			int* src_pitchY=0;
			int* src_pitchUV=0;
			int planesY=0;
			int planesUV=0;
			cache_count = n1-n0+1;

			if (vi.IsYUY2())
			{
				dst_planeY = copyframe->GetWritePtr();
				src_planesY = new const BYTE*[cache_count];
				src_pitchY = new int[cache_count];
					  dpitchY = copyframe->GetPitch();
				row_sizeY2=row_size;
				heightY2=height;
			}
			else
			{
				dpitchY = copyframe->GetPitch(PLANAR_Y);
				dpitchUV = copyframe->GetPitch(PLANAR_U);
				dst_planeY = copyframe->GetWritePtr(PLANAR_Y);
				dst_planeU = copyframe->GetWritePtr(PLANAR_U);
				dst_planeV = copyframe->GetWritePtr(PLANAR_V);
				src_planesY = new const BYTE*[cache_count];
				src_planesU = new const BYTE*[cache_count];
				src_planesV = new const BYTE*[cache_count];
				src_pitchY = new int[cache_count];
				src_pitchUV = new int[cache_count];
				row_sizeY2=row_sizeY;
				heightY2=heightY;
			}
			for (i = 0; i < cache_count; i++)
			{
				cache[i] = LoadFrame(n0+i);
				if (vi.IsYUY2())
				{
					 src_planesY[planesY] = cache[i].frame_ptr;
					 src_pitchY[planesY] = cache[i].pitch;
					 planesY++;
				}
				else  //YV12
				{
					src_planesY[planesY] = cache[i].frame_ptrY;
					src_pitchY[planesY] = cache[i].pitchY;
				    planesY++;

					src_planesU[planesUV] = cache[i].frame_ptrU;
					src_planesV[planesUV] = cache[i].frame_ptrV;
					src_pitchUV[planesUV] = cache[i].pitchUV;
					planesUV++;
				}
			}  // End for i
			// Blend Y
			if (planesY)
			{
				int c_div=32768/planesY;
				for (int j=0;j<heightY2;j++)
				{
					mmx_average_planes(dst_planeY,src_planesY,row_sizeY2, planesY-1, c_div);
					dst_planeY+=dpitchY;
					for (int i=0;i<planesY;i++)
						src_planesY[i]+=src_pitchY[i];
				}
			} // End if planesY
			if (planesUV)
			{
				int c_div=32768/planesUV;
				for (int j=0;j<heightUV;j++)
				{
					mmx_average_planes(dst_planeU,src_planesU,row_sizeUV, planesUV-1, c_div);
					mmx_average_planes(dst_planeV,src_planesV,row_sizeUV, planesUV-1, c_div);
					dst_planeU+=dpitchUV;
					dst_planeV+=dpitchUV;
					for (int i=0;i<planesUV;i++)
					{
						src_planesU[i]+=src_pitchUV[i];
						src_planesV[i]+=src_pitchUV[i];
					}
				}
			} // End if planesUV
			delete[] src_planesY;
			delete[] src_pitchY;
			if (src_planesU) delete[] src_planesU;
			if (src_planesV) delete[] src_planesV;
			if (src_pitchUV) delete[] src_pitchUV;
			nk = n1;
		} // End if blend
		else
		{
			nk = (decwhich == DECW_FIRST) ? n0
			   : (decwhich == DECW_LAST)  ? n1
			   : (decwhich == DECW_MIDDLE_DOWN) ? n0 + (n1-n0)/2
			   : /*(decwhich == DECW_MIDDLE_UP) ?*/ n0 + (n1-n0+1)/2;
			if(cache_out.frame_no != nk)
				cache_out = LoadFrame(nk);
			copyframe = cache_out.frame;

		}
	} // End if dec
	else
	{
		ni = nk = n;
		n0 = mapstart[mapinv[ni]];
		n1 = mapend[mapinv[ni]];
		if(cache_out.frame_no != nk)
			cache_out = LoadFrame(nk);
		copyframe = cache_out.frame;
	}

	if (show)
	{
		/* Generate text overlay. */
		env->MakeWritable(&copyframe);
		sprintf(buf, "DeDup %s", MYVERSION);
		DrawString(copyframe, 0, 0, buf);
		sprintf(buf, "Copyright 2004 Loren Merritt/Donald Graft/Klaus Post");
		DrawString(copyframe, 0, 1, buf);
		sprintf(buf, "frm %d: diff from frm %d = %2.2f%%", ni, ni+1, metrics[ni].metric);
		DrawString(copyframe, 0, 2, buf);
		if (decwhich==DECW_BLEND && n1 > n0)
			sprintf(buf, "ofrm %d: blending %d through %d", n, n0, n1);
		else
			sprintf(buf, "ofrm %d: using ifrm %d", n, nk);
		DrawString(copyframe, 0, 3, buf);

		for(i = -SHOW_NEIGHBORS; i <= SHOW_NEIGHBORS; i++)
			buf[i+SHOW_NEIGHBORS] = (n+i >= 0 && n+i < (int)num_iframes) ? DROP_NAMES[keep[n+i]] : ' ';
		buf[2*SHOW_NEIGHBORS+1] = '\0';
		DrawString(copyframe, 0, 4, buf);
		DrawString(copyframe, SHOW_NEIGHBORS, 5, "^");
	}

	return copyframe;
}

void Dup::CompareFrames(FRAMEINFO& frame1, FRAMEINFO& frame0)
{
	const unsigned char *srcp, *src0p;
	int highest_x, highest_y;
	int i, j, x, y;
	int offset_remainX = (vi.width&(~(BLKSIZE-1)));  // Offset into the frame (pixels)
	int offset_remainY = (vi.height/BLKSIZE)*BLKSIZE;  // yposition the remaining pixels start (lines)
	int remainX = vi.width&(BLKSIZE-1) + offset_remainX;      // Where do the remaining pixels end? (pixels)
	int remainY = (vi.height&(BLKSIZE-1)) + offset_remainY;   // Where do the remaining pixels end? (lines)

	/* Clear the block sums. */
	for (i = 0; i < yblocks; i++)
		for (j = 0; j < xblocks; j++)
			sum[i*xblocks+j] = 0;
	/* Do the comparison. */
	if (vi.IsYUY2())
	{
		src0p = frame0.frame_ptr;
		srcp = frame1.frame_ptr;

		for (y = 0; y < height; y++)
		{
			for (x = 0; x < row_size;)
			{
				sum[(y/BLKSIZE)*xblocks + x/(2*BLKSIZE)] += abs((int)srcp[x] - (int)src0p[x]);
				chroma ? x++ : x+=2;
			}
			srcp += frame1.pitch;
			src0p += frame0.pitch;
		}
	}
	else
	{
		src0p = frame0.frame_ptrY;
		srcp = frame1.frame_ptrY;
		if (have_isse)
		{
			isse_scenechange(srcp, src0p, heightY, row_sizeY, frame1.pitchY,  frame0.pitchY,(int*)sum);
			// Right remaining
			for (y = 0; y < remainY; y++)
			{
				for (x = offset_remainX; x < remainX; x++)
					sum[(y/BLKSIZE)*xblocks + x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
				srcp += frame1.pitchY;
				src0p += frame0.pitchY;
			}
			// Bottom remaining
			src0p = frame0.frame_ptrY+ (frame0.pitchY*offset_remainY);
			srcp = frame1.frame_ptrY + (frame1.pitchY*offset_remainY);
			for (y = offset_remainY; y < heightY; y++)
			{
				for (x = 0; x < row_sizeY; x++)
					sum[(y/BLKSIZE)*xblocks + x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
				srcp += frame1.pitchY;
				src0p += frame0.pitchY;
			}
		}
		else
		{
			for (y = 0; y < heightY; y++)
			{
				for (x = 0; x < row_sizeY; x++)
					sum[(y/BLKSIZE)*xblocks + x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
				srcp += frame1.pitchY;
				src0p += frame0.pitchY;
			}
		}
		if (chroma)
		{
			src0p = frame0.frame_ptrU;
			srcp = frame1.frame_ptrU;
			if (have_isse)
			{
				isse_scenechange_16(srcp, src0p, heightUV, row_sizeUV, frame1.pitchUV,  frame0.pitchUV,(int*)sum);
				// Right remaining
				for (y = 0; y < (remainY>>1); y++)
				{
					for (x = (offset_remainX>>1); x < row_sizeUV; x++)
						sum[2*(y/BLKSIZE)*xblocks + 2*x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
					srcp += frame1.pitchUV;
					src0p += frame0.pitchUV;
				}
				// Bottom remaining
				src0p = frame0.frame_ptrU+ ((frame0.pitchUV*offset_remainY)>>1);
				srcp = frame1.frame_ptrU + ((frame1.pitchUV*offset_remainY)>>1);
				for (y = (offset_remainY>>1); y < heightUV; y++)
				{
					for (x = 0; x < row_sizeUV; x++)
						sum[2*(y/BLKSIZE)*xblocks + 2*x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
					srcp += frame1.pitchUV;
					src0p += frame0.pitchUV;
				}
			}
			else
			{
				for (y = 0; y < heightUV; y++)
				{
					for (x = 0; x < row_sizeUV; x++)
						sum[(2*y/BLKSIZE)*xblocks + 2*x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
					srcp += frame1.pitchUV;
					src0p += frame0.pitchUV;
				}
			}

			src0p = frame0.frame_ptrV;
			srcp = frame1.frame_ptrV;
			if (have_isse)
			{
				isse_scenechange_16(srcp, src0p, heightUV, row_sizeUV, frame1.pitchUV,  frame0.pitchUV,(int*)sum);
				// Right remaining
				for (y = 0; y < (remainY>>1); y++)
				{
					for (x = (offset_remainX>>1); x < (remainX>>1); x++)
						sum[2*(y/BLKSIZE)*xblocks + 2*x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
					srcp += frame1.pitchUV;
					src0p += frame0.pitchUV;
				}
				// Bottom remaining
				src0p = frame0.frame_ptrV+ ((frame0.pitchUV*offset_remainY)>>1);
				srcp = frame1.frame_ptrV + ((frame1.pitchUV*offset_remainY)>>1);
				for (y = (offset_remainY>>1); y < (heightUV); y++)
				{
					for (x = 0; x < row_sizeUV; x++)
						sum[2*(y/BLKSIZE)*xblocks + 2*x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
					srcp += frame1.pitchUV;
					src0p += frame0.pitchUV;
				}
			}
			else
			{
				for (y = 0; y < heightUV; y++)
				{
					for (x = 0; x < row_sizeUV; x++)
						sum[(2*y/BLKSIZE)*xblocks + 2*x/BLKSIZE] += abs((int)srcp[x] - (int)src0p[x]);
					srcp += frame1.pitchUV;
					src0p += frame0.pitchUV;
				}
			}
		}
	}

	/* Now find the 32x32 block that has the greatest difference. */
	highest_sum = 0;
	highest_x = highest_y = 0;
	for (i = 0; i < yblocks; i++)
	{
		for (j = 0; j < xblocks; j++)
		{
			if (sum[i * xblocks + j] > highest_sum)
			{
				highest_sum = sum[i * xblocks + j];
				if (vi.IsYUY2()) highest_x = j * BLKSIZE * 2;
				else highest_x = j * BLKSIZE;
				highest_y = i * BLKSIZE;
			}
		}
	}
	/* Calculate the percentage difference for the block and store the results. */
	frame1.metric = (highest_sum * 100.0) / NORM;
	frame1.highest_x = highest_x;
	frame1.highest_y = highest_y;
}

FRAMEINFO Dup::LoadFrame(int n)
{
	FRAMEINFO out;
	out.frame_no = n;
	out.frame = child->GetFrame(n, env);
	if (vi.IsYUY2())
	{
		out.frame_ptr = (unsigned char *) out.frame->GetReadPtr();
		out.pitch = out.frame->GetPitch();
		row_size = out.frame->GetRowSize();
		height = out.frame->GetHeight();
	}
	else
	{
		out.frame_ptrY = (unsigned char *) out.frame->GetReadPtr(PLANAR_Y);
		out.frame_ptrU = (unsigned char *) out.frame->GetReadPtr(PLANAR_U);
		out.frame_ptrV = (unsigned char *) out.frame->GetReadPtr(PLANAR_V);
		out.pitchY = out.frame->GetPitch(PLANAR_Y);
		row_sizeY = out.frame->GetRowSize(PLANAR_Y);
		heightY = out.frame->GetHeight(PLANAR_Y);
		out.pitchUV = out.frame->GetPitch(PLANAR_U);
		row_sizeUV = out.frame->GetRowSize(PLANAR_U);
		heightUV = out.frame->GetHeight(PLANAR_U);
	}
	return out;
}

void Dup::CopyFrame(PVideoFrame& to, PVideoFrame& from)
{
	if (vi.IsYUY2())
	{
		env->BitBlt(to->GetWritePtr(), to->GetPitch(),
					from->GetReadPtr(), from->GetPitch(),
					from->GetRowSize(), from->GetHeight());
	}
	else
	{
		env->BitBlt(to->GetWritePtr(PLANAR_Y),
					to->GetPitch(PLANAR_Y),
					from->GetReadPtr(PLANAR_Y),
					from->GetPitch(PLANAR_Y),
					from->GetRowSize(PLANAR_Y),
					from->GetHeight(PLANAR_Y));
		env->BitBlt(to->GetWritePtr(PLANAR_U),
					to->GetPitch(PLANAR_U),
					from->GetReadPtr(PLANAR_U),
					from->GetPitch(PLANAR_U),
					from->GetRowSize(PLANAR_U),
					from->GetHeight(PLANAR_U));
		env->BitBlt(to->GetWritePtr(PLANAR_V),
					to->GetPitch(PLANAR_V),
					from->GetReadPtr(PLANAR_V),
					from->GetPitch(PLANAR_V),
					from->GetRowSize(PLANAR_V),
					from->GetHeight(PLANAR_V));
	}
}

void Dup::DrawBox(PVideoFrame& frame, int box_x, int box_y, bool crossp)
{
	int x, y;
	int xlim, ylim, xtmp, ytmp;
	unsigned char *dstp;
	if (vi.IsYUY2())
	{
		int pitch = frame->GetPitch();
		dstp = frame->GetWritePtr();
		xlim = box_x + 2*BLKSIZE;
		if (xlim > row_size) xlim = row_size;
		ylim = box_y + BLKSIZE;
		if (ylim > height) ylim = height;
		for (y = box_y; y < ylim; y++)
		{
			(dstp + y * pitch)[box_x] = 235;
			xtmp = box_x+2*(BLKSIZE - 1);
			if (xtmp < row_size)
				(dstp + y * pitch)[xtmp] = 235;
		}
		for (x = box_x; x < xlim; x+=4)
		{
			(dstp + (box_y) * pitch)[x] = 235;
			(dstp + (box_y) * pitch)[x+2] = 235;
			(dstp + (box_y) * pitch)[x+1] = 128;
			(dstp + (box_y) * pitch)[x+3] = 128;
			ytmp = box_y + BLKSIZE - 1;
			if (ytmp < height)
			{
				(dstp + ytmp * pitch)[x] = 235;
				(dstp + ytmp * pitch)[x+2] = 235;
				(dstp + ytmp * pitch)[x+1] = 128;
				(dstp + ytmp * pitch)[x+3] = 128;
			}
		}
		if (crossp)
		{
			for (y = box_y, x = 0; y < ylim; y++, x++)
			{
				xtmp = box_x+2*x;
				if (xtmp < row_size)
					(dstp + y * pitch)[box_x+2*x] = 235;
				xtmp = box_x+2*(BLKSIZE - 1 - x);
				if (xtmp < row_size)
					(dstp + y * pitch)[box_x+2*(BLKSIZE - 1 - x)] = 235;
			}
		}
	}
	else
	{
		int pitchY = frame->GetPitch(PLANAR_Y);
		dstp = frame->GetWritePtr(PLANAR_Y);
		xlim = box_x + BLKSIZE;
		if (xlim > row_sizeY) xlim = row_sizeY;
		ylim = box_y + BLKSIZE;
		if (ylim > heightY) ylim = heightY;
		for (y = box_y; y < ylim; y++)
		{
			(dstp + y * pitchY)[box_x] = 235;
			xtmp = box_x + (BLKSIZE - 1);
			if (xtmp < row_sizeY)
				(dstp + y * pitchY)[xtmp] = 235;
		}
		for (x = box_x; x < xlim; x++)
		{
			(dstp + (box_y) * pitchY)[x] = 235;
			ytmp = box_y + BLKSIZE - 1;
			if (ytmp < heightY)
				(dstp + ytmp * pitchY)[x] = 235;
		}
		if (crossp)
		{
			for (y = box_y, x = 0; y < ylim; y++, x++)
			{
				xtmp = box_x+x;
				if (xtmp < row_sizeY)
					(dstp + y * pitchY)[box_x+x] = 235;
				xtmp = box_x+(BLKSIZE - 1 - x);
				if (xtmp < row_sizeY)
					(dstp + y * pitchY)[box_x+(BLKSIZE - 1 - x)] = 235;
			}
		}
	}
}


AVSValue __cdecl Create_DeDup(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	float threshold = 0.3f;
	float threshold2 = 0.5f;
	int range2 = 0;
	float trigger2 = 5.0f;
	bool show = false;
	bool dec = true;
	int  maxcopies = 12;
	int  maxdrops = 1;
	int  decwhich = DECW_MIDDLE_UP;
	char* logfile = NULL;
	char* timefile = NULL;
	char* timeinfile = NULL;
	char* ovrfile = NULL;
	char* debugfile = NULL;

	return new Dup(args[0].AsClip(),
	(float)args[1].AsFloat(threshold),  // threshold for duplicate declaration
	(float)args[2].AsFloat(threshold2), // threshold when a pattern is detected
		args[3].AsInt(range2),          // range to search for decimate pattern
	(float)args[4].AsFloat(trigger2),   // what counts as a pattern
		args[5].AsBool(show),           // show biggest difference area
		args[6].AsBool(dec),            // decimate
		args[7].AsInt(maxcopies),       // max successive copies to emit
		args[8].AsInt(maxdrops),        // max successive frames to drop
		args[9].AsInt(decwhich),        // which of the duplicates do we keep?
		args[10].AsString(logfile),     // save metrics
		args[11].AsString(timefile),    // save timecodes
		args[12].AsString(timeinfile),  // load timecodes
		args[13].AsString(ovrfile),     // override per-frame settings
		args[14].AsString(debugfile),   // output decisions
		env);
}

AVSValue __cdecl Create_DupMetric(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	bool chroma = true;
	char* logfile = NULL;
	int  search = 1;

	return new Dup(args[0].AsClip(),
		args[1].AsBool(chroma),         // use chroma in differencing
		args[2].AsString(logfile),      // save metrics
		args[3].AsInt(search),          // compare frames up to this far apart
		env);
}

 /***
  * Accumulated differences of two planes.
  *
  * This routine is for testing luma planes.
  * The accumulated differences for each 32x32 box is written directly to blk_values.
  * Boxes not fitting within mod32 width sizes are filled with '0'.
  * (c) 2002, Donald Graft (algorithm)
  * (c) 2003, Klaus Post (ISSE code)
  ***/


void Dup::isse_scenechange(const BYTE* c_plane, const BYTE* tplane, int height, int width, int pitch, int t_pitch, int* blk_values) {
  __declspec(align(8)) static __int64 full = 0xffffffffffffffffi64;
  int wp=(width/BLKSIZE)*BLKSIZE;
  int hp=(height/BLKSIZE)*BLKSIZE;
  int pad_blk=(wp-width!=0);

  int y=0;
 __asm {
    mov esi, c_plane
    mov edi, tplane
    mov ebx,0
    jmp yloopover
    align 16
yloop:
    mov eax,[pad_blk]
    cmp eax,0
    je no_pad
    mov eax,blk_values
    mov [eax],0
    add eax,4
    mov blk_values,eax
no_pad:

    mov ebx, [y]
    mov edx, pitch    //copy pitch
    mov ecx, t_pitch    //copy pitch
    add ebx, 32
    shl edx,5
    shl ecx,5
    add edi,ecx     // add pitch to both planes
    add esi,edx
    mov y, ebx
yloopover:
    cmp ebx,[hp]
    jge endframe
    xor ebx, ebx  // X pos.
    align 16
xloop:
    cmp ebx,[wp]
    jge yloop
    mov eax,ebx       // Width (esi)
    mov ecx,ebx       // Width (edi)

    pxor mm6,mm6   // We maintain two sums, for better pairablility
    pxor mm7,mm7
    mov edx, 32
y_loop_inner:
    movq mm0,[esi+eax]
     movq mm2,[esi+eax+8]
    movq mm1,[edi+ecx]
     movq mm3,[edi+ecx+8]
    psadbw mm0,mm1    // Sum of absolute difference
     psadbw mm2,mm3
    paddd mm6,mm0     // Add...
     paddd mm7,mm2
    movq mm0,[esi+eax+16]
     movq mm2,[esi+eax+24]
    movq mm1,[edi+ecx+16]
     movq mm3,[edi+ecx+24]
    psadbw mm0,mm1
     psadbw mm2,mm3
    paddd mm6,mm0
     paddd mm7,mm2

    add eax,pitch
    add ecx,t_pitch

    dec edx
    jnz y_loop_inner

    mov eax,blk_values
    paddd mm6,mm7
    movd [eax],mm6
    add eax,4
    mov blk_values,eax

    add ebx,32
    jmp xloop

endframe:
    emms
  }
}

 /***
  * Accumulated differences of two planes.
  *
  * This routine is for testing chroma planes.
  * The accumulated differences for each 16x16 box is ADDED to the current values.
  * (c) 2002, Donald Graft (algorithm)
  * (c) 2003, Klaus Post (ISSE code)
  ***/

void Dup::isse_scenechange_16(const BYTE* c_plane, const BYTE* tplane, int height, int width, int pitch, int t_pitch, int* blk_values) {
  __declspec(align(8)) static __int64 full = 0xffffffffffffffffi64;
  int wp=(width/16)*16;
  int hp=(height/16)*16;
  int y=0;
  int pad_blk=(wp-width!=0);
 __asm {
    mov esi, c_plane
    mov edi, tplane
    mov ebx,0
    jmp yloopover
    align 16
yloop:
    mov eax,[pad_blk]
    cmp eax,0
    je no_pad
    mov eax,blk_values
    mov [eax],0
    add eax,4
    mov blk_values,eax
no_pad:
    mov ebx, [y]
    mov edx, pitch    //copy pitch
    mov ecx, t_pitch    //copy pitch
    add ebx, 16;
    shl edx,4
    shl ecx,4
    add edi,ecx     // add pitch to both planes
    add esi,edx
    mov y, ebx
yloopover:
    cmp ebx,[hp]
    jge endframe
    xor ebx, ebx  // X pos.
    align 16
xloop:
    cmp ebx,[wp]
    jge yloop
    mov eax,ebx       // Width (esi)
    mov ecx,ebx       // Width (edi)
    pxor mm6,mm6   // We maintain two sums, for better pairablility
    pxor mm7,mm7
    mov edx, 16
y_loop_inner:
    movq mm0,[esi+eax]
     movq mm2,[esi+eax+8]
    movq mm1,[edi+ecx]
     movq mm3,[edi+ecx+8]
    psadbw mm0,mm1    // Sum of absolute difference
     psadbw mm2,mm3
    paddd mm6,mm0     // Add...
     paddd mm7,mm2

    add eax,pitch
    add ecx,t_pitch

    dec edx
    jnz y_loop_inner

    mov eax,blk_values
    movd mm5,[eax]
    paddd mm6,mm7
    paddd mm6,mm5
    movd [eax],mm6
    add eax,4
    mov blk_values,eax

    add ebx,16
    jmp xloop


    jmp xloop
endframe:
    emms
  }
}

 /**
  * Blends one line of several frames equally weighed.
  *
  * An array of pointers is delivered as source frames.
  *
  * "planes" is the number of planes that should be blended.
  * "div" (divisor) should be multiplied by 32768.
  * (c) 2003, Klaus Post
  */

void Dup::mmx_average_planes(BYTE* dst_plane, const BYTE** src_planes, int width_mod8, int planes, int div) {
  __declspec(align(8)) static __int64 low_ffff = 0x000000000000ffffi64;

  __int64 div64 = (__int64)(div) | ((__int64)(div)<<16) | ((__int64)(div)<<32) | ((__int64)(div)<<48);
  div>>=1;
  __int64 add64 = (__int64)(div) | ((__int64)(div)<<32);

  if (planes<=0) return;
  __asm {
    mov esi,dst_plane;
    xor eax,eax          // EAX will be plane offset (all planes).
    align 16
testplane:
    cmp eax, [width_mod8]
    jge outloop

    movq mm0,[esi+eax]  // Load current frame pixels
     pxor mm2,mm2        // Clear mm2
    movq mm6,mm0
     movq mm7,mm0
    mov edi,[src_planes];  // Adress of planeP array is now in edi
    mov ebx,[planes]   // How many planes (this will be our counter)
    punpcklbw mm6,mm2    // mm0 = lower 4 pixels
     punpckhbw mm7,mm2     // mm1 = upper 4 pixels
    lea edi,[edi+ebx*4]

    align 16
kernel_loop:
    mov edx,[edi]
    movq mm4,[edx+eax]      // Load 8 pixels from test plane
     pxor mm1,mm1
    movq mm5,mm4
    punpcklbw mm4,mm1         // mm4 = lower pixels
     punpckhbw mm5,mm1        // mm5 = upper pixels
    paddusw mm6,mm4
     paddusw mm7,mm5

    sub edi,4
    dec ebx
    jnz kernel_loop
     // Multiply (or in reality divides) added values
    movq mm4,[add64]
    pxor mm5,mm5
     movq mm0,mm6
    movq mm1,mm6
     punpcklwd mm0,mm5         // low,low
    movq mm6,[div64]
    punpckhwd mm1,mm5         // low,high
     movq mm2,mm7
    pmaddwd mm0,mm6
     punpcklwd mm2,mm5         // high,low
     movq mm3,mm7
     paddd mm0,mm4
    pmaddwd mm1,mm6
     punpckhwd mm3,mm5         // high,high
     psrld mm0,15
     paddd mm1,mm4
    pmaddwd mm2,mm6
     packssdw mm0, mm0
     psrld mm1,15
     paddd mm2,mm4
    pmaddwd mm3,mm6
     packssdw mm1, mm1
     psrld mm2,15
     paddd mm3,mm4
    psrld mm3,15
     packssdw mm2, mm2
    packssdw mm3, mm3
     packuswb mm0,mm5
    packuswb mm1,mm5
     packuswb mm2,mm5
    packuswb mm3,mm5
     movq mm4, [low_ffff]
    pand mm0, mm4;
     pand mm1, mm4;
    pand mm2, mm4;
     pand mm3, mm4;
    psllq mm1, 16
    psllq mm2, 32
     por mm0,mm1
    psllq mm3, 48
    por mm2,mm3
    por mm0,mm2
    movq [esi+eax],mm0

    add eax,8   // Next 8 pixels
    jmp testplane
outloop:
    emms
  }

}


extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit2(IScriptEnvironment* env)
{
    env->AddFunction("DeDup", "c[threshold]f[threshold2]f[range2]i[trigger2]f[show]b[dec]b[maxcopies]i[maxdrops]i[decwhich]i[log]s[times]s[timesin]s[ovr]s[debug]s", Create_DeDup, 0);
    env->AddFunction("DupMC", "c[chroma]b[log]s[search]i", Create_DupMetric, 0);
    return 0;
}

