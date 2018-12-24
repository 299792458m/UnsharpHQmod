/*
	UnsharpHQ ver 0.5
	Copyright (c) 2011-2014, Lucas De Lio (lucasdelio@gmail.com).

	This file is part of UnsharpHQ.

    UnsharpHQ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    UnsharpHQ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with UnsharpHQ.  If not, see <http://www.gnu.org/licenses/>.
*/
    
#include "windows.h"
#include "avisynth.h"
#include "math.h"

#define VERDE1 100
#define VERDE2 100


#define uc unsigned char

#define threshold 20
#define sharp_str 4
#define smoot_str .5

#include "utiles.h"

class UnsharpHQ : public GenericVideoFilter {   
	float MIN,MAX;
	int thres;
	bool show;
	int mode;

public:
	UnsharpHQ(PClip _child,int , float ,float ,bool ,int, IScriptEnvironment* env);
	~UnsharpHQ();
	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
};



	UnsharpHQ::UnsharpHQ(PClip _child,int _thres,float _MAX,float _MIN,bool _show,int _mode,IScriptEnvironment* env) :
	GenericVideoFilter(_child) {
		thres = _thres;
		MIN = _MIN;
		MAX = _MAX;
		show = _show;
		mode = _mode;
		
		if( !vi.IsYV12() )
			env->ThrowError("UnSharpHQ supports YV12 only");

		if( thres<0 || thres>99 || MIN<0 || MIN>4 || MAX<0 || MAX>99)
			env->ThrowError("UnsharpHQ(THRESHOLD=20, SHARPSTR=4.0 , SMOOTH=0.5 , SHOW=FALSE)  \n "
			"        USAGE:\n"
			"      THRES [ 0.0 , 99.0]\n"
			"   SHARPSTR [ 0.0 , 99.0]\n"
			" DENOISESTR [ 0.0 ,  4.0]\n"
		);
		if( mode >=2 )
			env->ThrowError("MODE is 0 or 1");

	}

	UnsharpHQ::~UnsharpHQ() {}


PVideoFrame __stdcall UnsharpHQ::GetFrame(int n, IScriptEnvironment* env) {

	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst = env->NewVideoFrame(vi);

	if (vi.IsYV12()) {
		const unsigned char* srcp = src->GetReadPtr(PLANAR_Y);

		unsigned char * dstp = dst->GetWritePtr(PLANAR_Y);
		int src_pitch = src->GetPitch(PLANAR_Y);
		int dst_pitch = dst->GetPitch(PLANAR_Y);
		int rowsize = dst->GetRowSize(PLANAR_Y);  
		int height = dst->GetHeight(PLANAR_Y);

		if (mode==0) FieldUnSharpY(dstp, srcp, rowsize,height,src_pitch,dst_pitch,(float)thres,-(MIN+1.0f),MAX+1.0f,show,mode);
		else FieldUnSharpY(dstp, srcp, rowsize,height,src_pitch,dst_pitch,(float)thres,-MIN,MAX,show,mode);
	
		//- - - - - - - - U - - - - - - - - - -  
		
		srcp = src->GetReadPtr(PLANAR_U);
		dstp = dst->GetWritePtr(PLANAR_U);
		src_pitch = src->GetPitch(PLANAR_U);
		dst_pitch = dst->GetPitch(PLANAR_U);
		rowsize = dst->GetRowSize(PLANAR_U);  
		height = dst->GetHeight(PLANAR_U);
		
		if (!show)
			Fieldcopy(dstp, srcp, rowsize,height, dst_pitch, src_pitch);
		else
			FieldNull(dstp,rowsize,height, dst_pitch,VERDE1);
			
		//- - - - - - - - V - - - - - - - - - -  
		srcp = src->GetReadPtr(PLANAR_V);
		dstp = dst->GetWritePtr(PLANAR_V);

		if (!show)
			Fieldcopy(dstp, srcp, rowsize,height, dst_pitch, src_pitch);
		else
			FieldNull(dstp,rowsize,height, dst_pitch,VERDE2);

	}	
	return dst;
}



AVSValue __cdecl Create_UnsharpHQ(AVSValue args, void* user_data, IScriptEnvironment* env) {
	return new UnsharpHQ(args[0].AsClip()
		,args[1].AsInt(threshold)
		,(float)args[2].AsFloat(sharp_str)
		,(float)args[3].AsFloat(smoot_str)
		,args[4].AsBool(false)
		,args[5].AsInt(0)
	, env);
}

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit2(IScriptEnvironment* env) {
    env->AddFunction("UnsharpHQ", "c[THRESHOLD]i[SHARPSTR]f[SMOOTH]f[SHOW]b[MODE]i", Create_UnsharpHQ, 0);

    // c - Video Clip
    // i - Integer number
    // f - Float number
    // s - String
    // b - boolean

    return "`UnsharpHQ' UnsharpHQ plugin";
 }
