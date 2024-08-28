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

		The author can be contacted at:
		Donald Graft
		neuron2@attbi.com.
*/

#include "internal.h"
#include "info.h"

#define BLKSIZE 32
#define NORM (235*BLKSIZE*BLKSIZE)
#define MAX_COPIES 20
#define DUPVERSION "2.20 beta 1"
#define MYVERSION "0.1"
#define VERSION_PRINTF "DeDup %s by Loren Merritt, based on Dup %s by Donald Graft/Klaus Post, Copyright 2004\n", MYVERSION, DUPVERSION

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

class Dup : public GenericVideoFilter
{
public:
	// In Dedup mode
    Dup(PClip _child, float _threshold, bool _show, bool _dec, int _maxcopies,
		int _maxdrops, bool _blend, const char* _logfile, const char* _timefile,
		const char* _ovrfile, const char* _debugfile, IScriptEnvironment* _env) :
	    GenericVideoFilter(_child), threshold(_threshold), show(_show), dec(_dec),
		maxcopies(_maxcopies), maxdrops(_maxdrops), blend(_blend), env(_env), pass(2)
	{
		if (!vi.IsYUY2() && !vi.IsYV12())
			env->ThrowError("DeDup: requires YUY2 or YV12 source");
		if (maxcopies > 20)
			env->ThrowError("DeDup: maxcopies must be <= 20");
		if (threshold < 0.0 || threshold > 100.0)
			env->ThrowError("DeDup: threshold out of range (0.0-100.0)");
		//if (blend && !dec)
		//	env->ThrowError("DeDup: blend=true requires dec=true");
		if (_logfile)
		{
			logfile = fopen(_logfile, "r");
			if (logfile == NULL)
				env->ThrowError("DeDup: failed to open logfile");
		}
		else
			env->ThrowError("DeDup: option 'log' required");
		if (_timefile)
		{
			timefile = fopen(_timefile, "w");
			if (timefile == NULL)
				env->ThrowError("DeDup: failed to open timefile");
		}
		else
			env->ThrowError("DeDup: option 'times' required");
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

		xblocks = (vi.width+BLKSIZE-1) / BLKSIZE;
		yblocks = (vi.height+BLKSIZE-1) / BLKSIZE;

		num_iframes = vi.num_frames;
		sum = NULL;

		/* For safety in case someone came in without doing it. */
		__asm emms;

		LoadFirstPass();
	}
	// In metric collection mode
    Dup(PClip _child, bool _chroma, const char* _logfile, IScriptEnvironment* _env) :
	    GenericVideoFilter(_child), chroma(_chroma), env(_env), pass(1)
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
		thresholds = NULL;

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
		if (thresholds)   delete [] thresholds;
	}
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* _env);
    void isse_scenechange(const BYTE* c_plane, const BYTE* tplane, int height, int width, int pitch,  int t_pitch, int* blk_values);
    void isse_scenechange_16(const BYTE* c_plane, const BYTE* tplane, int height, int width, int pitch,  int t_pitch, int* blk_values);
    void mmx_average_planes(BYTE* dst_plane, const BYTE** src_planes, int width_mod8, int planes, int div);

private:
	//// Options
	float threshold;
	bool show, dec, chroma, blend;
	int maxcopies;     // max consecutive frames to merge (blend / copy)
	int maxdrops;      // max consecutive frames to drop
	bool have_isse;
	int pass;          // 1 => collect metrics, 2 => decimate
	FILE* logfile;     // load/save stats for 2 passes
	FILE* timefile;    // save matroska timecodes
	FILE* ovrfile;     // override calculated dups (or just about any other option)
	FILE* debugfile;   // report decisions here

	//// Per-frame options
	float* thresholds;

	//// State
	struct FRAMEINFO cache[MAX_COPIES+1];
	int cache_count;
	PVideoFrame copyframe;
	IScriptEnvironment* env;

	int xblocks, yblocks;
	unsigned int *sum, highest_sum;

	FRAMEMETRIC* metrics; // metrics[i] = diff(i, i+1)
	bool* metrics_done;   // which metrics have we calculated yet?
	char* keep;           // keep[i] = (0 => drop, 1 => output due to minframerate, 2 => keep)
	unsigned int* mapend; // mapend[out_frame] = in_frame (end of blend)
	unsigned int* mapstart; // start of blend
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

	inline void DrawString(PVideoFrame& dst, int x, int y, const char *s)
	{
		if (vi.IsYUY2()) ::DrawStringYUY2(dst, x, y, s);
		else ::DrawString(dst, x, y, s);
	}
};

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

	metrics = new FRAMEMETRIC [num_iframes];
	metrics_done = new bool [num_iframes];
	keep = new char [num_iframes];
	mapend = new unsigned int [num_iframes];
	mapstart = new unsigned int [num_iframes];
	mapinv = new unsigned int [num_iframes];
	thresholds = new float [num_iframes];

	if(metrics == NULL || keep == NULL || mapend == NULL || mapstart == NULL
			|| mapinv == NULL || metrics_done == NULL || thresholds == NULL)
		env->ThrowError("Dup: cannot allocate needed memory");

	memset(metrics_done, 0, num_iframes * sizeof(*metrics_done));
	memset(keep, 0, num_iframes * sizeof(*keep));
	for(i=0; i<num_iframes; i++)
		thresholds[i] = threshold;

	// read stats from first pass
	char buf[201];
	while(fgets(buf, 200, logfile))
	{
		if(!sscanf(buf, "frm %u: diff from frm %u = %f%% at (%d,%d)", &frame_no, &frame_next, &metric, &highest_x, &highest_y))
			continue;
		if(frame_next != frame_no+1)
			continue;
		metrics[frame_no].metric = metric;
		metrics[frame_no].highest_x = highest_x;
		metrics[frame_no].highest_y = highest_y;
		metrics_done[frame_no] = true;
	}
	fclose(logfile);
	logfile = NULL;

	for(i=0; i<num_iframes; i++)
		if(!metrics_done[i])
			// should print frame number?
			env->ThrowError("DeDup: incomplete first pass");
	delete [] metrics_done;
	metrics_done = NULL;

	// read override file
	if(ovrfile != NULL)
	{
		char* pos;
		int len;

		while(fgets(buf, 200, ovrfile))
		{
			if(buf[0] == '#')
				continue;

			unsigned int frame0, frame1;
			if(sscanf(buf, " %u-%u", &frame0, &frame1) == 2)
				;// do nothing
			else if(sscanf(buf, " %u", &frame0))
				frame1 = frame0;
			else
				//invalid line
				continue;

			pos = buf + strspn(buf, " \t");
			pos += strcspn(pos, " \t");
			pos += strspn(pos, " \t");
			while(*pos != '\0')
			{
				if(sscanf(pos, "threshold=%f", &thresholds[frame0]))
					for(i=frame0+1; i<=frame1; i++)
						thresholds[i] = thresholds[frame0];
				// word of the form /[kd]+/ forces keep or drop
				/*
				else if((len = strspn(pos, "kd")) > 0 && len == strcspn(pos, " \t"))
					for(i=frame0; i<=frame1; i++)
						thresholds[i] = (pos[(i-frame0)%len] == 'k') ? -1.0f : 101.0f;
				*/
				// FIXME: doesn't check for non-kd words
				else if(pos[0] == 'k' || pos[0] == 'd')
					for(j=0, i=frame0; i<=frame1; i++, j++)
					{
						if(pos[j] != 'k' && pos[j] != 'd')
							j = 0;
						if(pos[j] == 'k')
							thresholds[i] = -1.0f;
						else
							thresholds[i] = 101.0f;
					}


				pos += strcspn(pos, " \t");
				pos += strspn(pos, " \t");
			}
		}
		fclose(ovrfile);
		ovrfile = NULL;
	}

	// decide which frames to drop
	// TODO: replace the original dup algorithm with my dual threshold / neighborhood
	num_oframes = 0;
	int numdropped = 0;
	for(i=0; i<num_iframes; i++)
	{
		mapend[num_oframes] = i;
		numdropped++;
		if(metrics[i].metric > thresholds[i])
		{
			keep[i] = KEPT_THRESH;
			num_oframes++;
			numdropped = 0;
		} else if(numdropped > maxcopies) {
			keep[i] = KEPT_MAXCOPIES;
			num_oframes++;
			numdropped = 0;
		} else if(numdropped % maxdrops == 0) {
			keep[i] = KEPT_MAXDROPS;
			num_oframes++;
		}
	}
	if(dec)
		vi.num_frames = num_oframes;

	// correct the sources of frames kept due to maxdrops
	int blendstart = 0;
	for(i=0; i<num_oframes; i++)
	{
		mapstart[i] = blendstart;
		if(keep[mapend[i]] && keep[mapend[i]] != KEPT_MAXDROPS)
			blendstart = mapend[i]+1;
	}
	for(i=num_oframes-2; i>0; i--)
		if(keep[mapend[i]] == KEPT_MAXDROPS)
			mapend[i] = mapend[i+1];
	for(i=j=0; i<num_oframes; i++)
		for(; j<=mapend[i]; j++)
			mapinv[j] = i;

	// calc timestamps
	double mspf = 1000. * double(vi.fps_denominator) / vi.fps_numerator;
	fprintf(timefile, "# timecode format v2\n");
	fprintf(timefile, "%.6f\n", 0.0f);
	for(i=1; i<num_iframes; i++)
		if(keep[i-1])
			fprintf(timefile, "%.6lf\n", i * mspf);

	fclose(timefile);
	timefile = NULL;

	if(debugfile)
	{
		fprintf(debugfile, VERSION_PRINTF);
		fprintf(debugfile, "\ninframe (outframe): kept?   metric <> threshold\n");
		for(i=0; i<num_iframes; i++)
			fprintf(debugfile, "%d (%d): %d   %.4f%% %c %.4f%%\n", i, mapinv[i], keep[i], metrics[i].metric, keep[i] == KEPT_THRESH ? '>' : '<', thresholds[i]);

		fclose(debugfile);
		debugfile = NULL;
	}
}

void Dup::CalcMetric(int n)
{
	if(metrics_done[n])
		return;

	if(n < num_iframes-1)
	{
		// compare frame n to n+1
		cache[0] = LoadFrame(n);
		cache[0].metric = 0.0;
		cache[1] = LoadFrame(n+1);
		cache[1].metric = 0.0;
		CompareFrames(cache[0], cache[1]);
	}
	else
	{
		// always keep the last frame
		cache[0] = LoadFrame(n);
		cache[0].metric = 100.0;
		cache[0].highest_x = 0;
		cache[0].highest_y = 0;
	}

	fprintf(logfile, "frm %d: diff from frm %d = %2.4f%% at (%d,%d)\n",
		n, n+1, cache[0].metric, cache[0].highest_x, cache[0].highest_y);
}

PVideoFrame Dup::ConstructFrame(int n)
{
	char buf[80];
	int offset_remainX = (vi.width&(~(BLKSIZE-1)));  // Offset into the frame (pixels)
	int offset_remainY = (vi.height/BLKSIZE)*BLKSIZE;  // yposition the remaining pixels start (lines)
	int remainX = vi.width&(BLKSIZE-1) + offset_remainX;      // Where do the remaining pixels end? (pixels)
	int remainY = (vi.height&(BLKSIZE-1)) + offset_remainY;   // Where do the remaining pixels end? (lines)

	int n0, n1, ni;

	PVideoFrame showframe;

	if (dec)
	{
		n0 = mapstart[n]; // beginning of the run of duplicates
		ni = n1 = mapend[n];     // end of the run of duplicates
		/* If blend=true, modify copyframe to be a blend of all the frames in the
		   string of duplicates. Skip if there's only one frame */
		if (blend && n0 < n1)
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

		} // End if blend
		else
		{
			copyframe = child->GetFrame(n1, env);
		}
	} // End if dec
	else
	{
		ni = n;
		n0 = mapstart[mapinv[ni]];
		n1 = mapend[mapinv[ni]];
		copyframe = child->GetFrame(n, env);
	}

	if (show)
	{
		/* Generate show data overlay. */
		sprintf(buf, "DeDup %s", MYVERSION);
		DrawString(copyframe, 0, 0, buf);
		sprintf(buf, "Copyright 2004 Loren Merritt/Donald Graft/Klaus Post");
		DrawString(copyframe, 0, 1, buf);
		sprintf(buf, "frm %d: diff from frm %d = %2.2f%%", ni, ni+1, metrics[ni].metric);
		DrawString(copyframe, 0, 3, buf);
		if (!dec)
			DrawBox(copyframe, metrics[ni].highest_x, metrics[ni].highest_y, (metrics[ni].metric < threshold));
		if (blend && n1 > n0)
			sprintf(buf, "Blending %d through %d", n0, n1);
		else
			sprintf(buf, "Using frm %d", n1);
		DrawString(copyframe, 0, 4, buf);
	}

	/* Return the appropriate frame. */
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
	float threshold = 1;
	bool show = false;
	bool dec = true;
	int  maxcopies = 12;
	int  maxdrops = 3;
	bool blend = false;
	char* logfile = NULL;
	char* timefile = NULL;
	char* ovrfile = NULL;
	char* debugfile = NULL;

    return new Dup(args[0].AsClip(),
		args[1].AsFloat(threshold),		// threshold for duplicate declaration
		args[2].AsBool(show),			// show biggest difference area
		args[3].AsBool(dec),			// decimate
		args[4].AsInt(maxcopies),		// max successive copies to emit
		args[5].AsInt(maxdrops),		// max successive frames to drop
		args[6].AsBool(blend),			// blend the duplicates
		args[7].AsString(logfile),		// save metrics
		args[8].AsString(timefile),		// save timecodes
		args[9].AsString(ovrfile),		// save timecodes
		args[10].AsString(debugfile),	// save timecodes
		env);
}

AVSValue __cdecl Create_DupMetric(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	bool chroma = true;
	char* logfile = NULL;

    return new Dup(args[0].AsClip(),
		args[1].AsBool(chroma),			// use chroma in differencing
		args[2].AsString(logfile),		// save metrics
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
    env->AddFunction("DeDup", "c[threshold]f[show]b[dec]b[maxcopies]i[maxdrops]i[blend]b[log]s[times]s[ovr]s[debug]s", Create_DeDup, 0);
	env->AddFunction("DupMC", "c[chroma]b[log]s", Create_DupMetric, 0);
    return 0;
}
