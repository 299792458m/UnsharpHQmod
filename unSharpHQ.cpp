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

template<int bits_per_pixel>
void FieldUnSharpYHBD(unsigned char* dstp0, const unsigned char* srcp, int rowsize, int height, int src_pitch, int dst_pitch,
	float thres, float MIN, float MAX, bool show)
{
	constexpr int max_pixel_size = (1 << bits_per_pixel);
	constexpr int max_pixel_value = max_pixel_size - 1;
	constexpr int cnvsft8toHBD = bits_per_pixel - 8;
	//__declspec(align(16)) int32_t funcion[max_pixel_size];
	int32_t* funcion= (int32_t*)_aligned_malloc(max_pixel_size * sizeof(int32_t), 16);	//�T�C�Y���傫��(16bit�ł�256kB�ɂȂ�)���߃q�[�v��

	//if (show) {
	//	MIN = -4;
	//	MAX = 60;
	//}

	int ocho = 8;	//1<<PRECISION�����A�Ӗ��s�� �s�v�ł́H

	//MIN = -(SMOOTH+1.0)�ŁAMAX = SHAPSTR+1.0�ɂȂ��Ă�

	double giro, atanMIN, atanMAX, A, S;
	float fmin, fmax;	//MIN/MAX��1���g�����ꂽ�l���n����遨�Ăяo������ύX���Ċg�����Ȃ��悤�ɂ��� �K�v�ɉ����ā��Œ���


	//i==giro��1�ɂȂ�悤��mode1���I�t�Z�b�g���Ă�����ۂ����Ă邪�A�`��i=giro��atan(0)�̌`�ɂȂ��Ă�(smooth�������܂�l������ĂȂ�)
	fmin = MIN * ocho / 4 - 0.99999f;	//-1.0�ɂ���Ɖ��Z���x��-1.0001�Ƃ��ɂȂ�̂�
	fmax = MAX * ocho / 4;

	giro = thres;
	atanMIN = atan((0 - giro) / curva);			//0�Ȃ̂Ŋ��ʂł������Ă�����������
	atanMAX = atan((255 - giro) / curva);		//funcion��(i-giro)/curv�ł��邱�Ƃ��l����ƁA������(i-giro)��������
	A = (fmax - fmin) / (atanMAX - atanMIN);		//(atanMIN,atanMAX)��(MIN,MAX)�̌X��
	S = fmin - A * atanMIN - 1.0;					//(atanMIN,atanMAX)��(MIN,MAX)�̃I�t�Z�b�g

	for (int ii = 0; ii < bits_per_pixel; ii++) {				//i=giro��funcion��1�ɂȂ�悤�ɌJ��Ԃ����Z ���x�͂悭�m��񂪂����񎟎�������Ȃ炱��Ȃ���
		giro = thres - tan(-S / A) * curva;
		atanMIN = atan((0 - giro) / curva);
		atanMAX = atan((255 - giro) / curva);
		A = (fmax - fmin) / (atanMAX - atanMIN);
		S = fmin - A * atanMIN - 1.0;
	}
	for (int i = 0; i < max_pixel_size; i++) {
		funcion[i] = (__int32)(A * atan(((double)i/(1<< cnvsft8toHBD) - giro) / curva) + S + 1.0);	//������HBD���Z����
		//test[i]  = ( A*atan( (i -giro) / curva  ) + S +1.0);
	}


	const uint16_t* row0 = (uint16_t*)srcp;		//���ꂪ���S���ǂ����͒m��Ȃ��E�E�E
	const uint16_t* row1 = (uint16_t*)(srcp + src_pitch);
	const uint16_t* row2 = (uint16_t*)(srcp + src_pitch * 2);
	const uint16_t* row3 = (uint16_t*)(srcp + src_pitch * 3);
	const uint16_t* row4 = (uint16_t*)(srcp + src_pitch * 4);

	uint16_t* dstp = (uint16_t*)dstp0;	//16bit�փL���X�g

	memcpy((BYTE*)dstp, (BYTE*)row0, rowsize);		//1,2�s�ڂ̓\�[�X���R�s�[
	dstp = dstp + dst_pitch / sizeof(uint16_t);		//pitch��byte�P�ʂȂ̂�
	memcpy((BYTE*)dstp, (BYTE*)row1, rowsize);
	dstp = dstp + dst_pitch / sizeof(uint16_t);

	int width = rowsize / sizeof(uint16_t);
	int w;
	for (int h = 2; h < height - 2; h++) {
		{
			dstp[0] = row2[0];			//1,2��ڂ̓\�[�X���R�s�[
			dstp[1] = row2[1];
			for (w = 2; w < width - 2; w++) {
				int avrg, diff, A, B;

				avrg = row0[w]
					+ row1[w - 1] + row1[w + 1]
					+ row2[w - 2] + row2[w + 2]
					+ row3[w - 1] + row3[w + 1]
					+ row4[w];
				avrg = avrg / 8;

				diff = abs(avrg - row0[w])
					+ abs(avrg - row1[w - 1]) + abs(avrg - row1[w + 1])
					+ abs(avrg - row2[w - 2]) + abs(avrg - row2[w + 2])
					+ abs(avrg - row3[w - 1]) + abs(avrg - row3[w + 1])
					+ abs(avrg - row4[w]);			 //asm�ł�abs��2�̑傫����-���������ł���Ă�
				diff = min(max(diff, 0), max_pixel_value);

				B = funcion[diff];
				A = B + ocho;		// (1+B/ocho)*ocho

				dstp[w] = min(max((row2[w] * A - avrg * B) / ocho, 0), max_pixel_value);
			}
			dstp[width - 2] = row2[width - 2];			//last-1,last��̓\�[�X���R�s�[
			dstp[width - 1] = row2[width - 1];
		}

		row0 = row1;
		row1 = row2;
		row2 = row3;
		row3 = row4;
		row4 += src_pitch/ sizeof(uint16_t);

		dstp += dst_pitch/ sizeof(uint16_t);
	}
	memcpy((BYTE*)dstp, (BYTE*)row2, rowsize);			//last-1,last�s�ڂ̓\�[�X���R�s�[
	dstp = dstp + dst_pitch / sizeof(uint16_t);
	memcpy((BYTE*)dstp, (BYTE*)row3, rowsize);
	dstp = dstp + dst_pitch / sizeof(uint16_t);

	_aligned_free(funcion);
}

class UnsharpHQ : public GenericVideoFilter {   
	float MIN,MAX;
	int thres;
	bool show;
	int mode;

public:
	UnsharpHQ(PClip _child,int , float ,float ,bool ,int,IScriptEnvironment* env);
	~UnsharpHQ();
	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
};



UnsharpHQ::UnsharpHQ(PClip _child,int _thres,float _MAX,float _MIN,bool _show,int _mode, IScriptEnvironment* env) :
GenericVideoFilter(_child) {
	thres = _thres;
	MIN = _MIN;	//SMOOTH
	MAX = _MAX;	//SHARPSTR
	show = _show;
	mode = _mode;
		
	if( !vi.IsYUV() || !vi.IsPlanar() )
		env->ThrowError("UnSharpHQ: supports only planar YUV(YV12,YUV420..)");	//���`�����l�����R�s�[���ĂȂ��̂�YUVA���_��

	if (vi.BitsPerComponent() > 16)
		env->ThrowError("UnsharpHQ: supports only 8-16bit integer input. not float!");

	if( thres<0 || thres>99 || MIN<0 || MIN>4 || MAX<0 || MAX>99)
		env->ThrowError("UnsharpHQ(THRESHOLD=20, SHARPSTR=4.0 , SMOOTH=0.5 , SHOW=FALSE)  \n "
		"        USAGE:\n"
		"  THRESHOLD [ 0.0 , 99.0]\n"
		"   SHARPSTR [ 0.0 , 99.0]\n"
		"     SMOOTH [ 0.0 ,  4.0]\n"
	);
	if( mode >=2 ) env->ThrowError("MODE is 0 to 1");
}

UnsharpHQ::~UnsharpHQ() {}


PVideoFrame __stdcall UnsharpHQ::GetFrame(int n, IScriptEnvironment* env) {

	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst = env->NewVideoFrame(vi);
	int bits_per_pixel = vi.BitsPerComponent();

	if (vi.IsYUV()) {
		if (bits_per_pixel == 8) {
			const uint8_t* srcp = src->GetReadPtr(PLANAR_Y);

			uint8_t* dstp = dst->GetWritePtr(PLANAR_Y);
			int src_pitch = src->GetPitch(PLANAR_Y);
			int dst_pitch = dst->GetPitch(PLANAR_Y);
			int rowsize = dst->GetRowSize(PLANAR_Y);
			int height = dst->GetHeight(PLANAR_Y);

			FieldUnSharpY(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, show, mode);

			//- - - - - - - - U - - - - - - - - - -  

			srcp = src->GetReadPtr(PLANAR_U);
			dstp = dst->GetWritePtr(PLANAR_U);
			src_pitch = src->GetPitch(PLANAR_U);
			dst_pitch = dst->GetPitch(PLANAR_U);
			rowsize = dst->GetRowSize(PLANAR_U);
			height = dst->GetHeight(PLANAR_U);

			if (!show)
				Fieldcopy(dstp, srcp, rowsize, height, dst_pitch, src_pitch);
			else
				FieldNull(dstp, rowsize, height, dst_pitch, VERDE1);

			//- - - - - - - - V - - - - - - - - - -  
			srcp = src->GetReadPtr(PLANAR_V);
			dstp = dst->GetWritePtr(PLANAR_V);
			src_pitch = src->GetPitch(PLANAR_V);
			dst_pitch = dst->GetPitch(PLANAR_V);
			rowsize = dst->GetRowSize(PLANAR_V);
			height = dst->GetHeight(PLANAR_V);

			if (!show)
				Fieldcopy(dstp, srcp, rowsize, height, dst_pitch, src_pitch);
			else
				FieldNull(dstp, rowsize, height, dst_pitch, VERDE2);

		}
		else {	//HBD
			const uint8_t* srcp = src->GetReadPtr(PLANAR_Y);
			uint8_t* dstp = dst->GetWritePtr(PLANAR_Y);
			int src_pitch = src->GetPitch(PLANAR_Y);
			int dst_pitch = dst->GetPitch(PLANAR_Y);	//���̍s�܂ł�byte��(�p�f�B���O�܂�)
			int rowsize = dst->GetRowSize(PLANAR_Y);	//1�s��byte��
			int height = dst->GetHeight(PLANAR_Y);		//���񂠂邩(pixel��)

			switch (bits_per_pixel) {
			case 10:	FieldUnSharpYHBD<10>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, show); break;	//mode�͍폜
			case 12:	FieldUnSharpYHBD<12>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, show); break;
			case 14:	FieldUnSharpYHBD<14>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, show); break;
			case 16:
			default:	FieldUnSharpYHBD<16>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, show); break;
			}
			
			//- - - - - - - - U - - - - - - - - - -  

			srcp = src->GetReadPtr(PLANAR_U);
			dstp = dst->GetWritePtr(PLANAR_U);
			src_pitch = src->GetPitch(PLANAR_U);
			dst_pitch = dst->GetPitch(PLANAR_U);
			rowsize = dst->GetRowSize(PLANAR_U);
			height = dst->GetHeight(PLANAR_U);

			if (!show)
				Fieldcopy(dstp, srcp, rowsize, height, dst_pitch, src_pitch);
			else
				FieldNull(dstp, rowsize, height, dst_pitch, VERDE1);

			//- - - - - - - - V - - - - - - - - - -  
			srcp = src->GetReadPtr(PLANAR_V);
			dstp = dst->GetWritePtr(PLANAR_V);
			src_pitch = src->GetPitch(PLANAR_V);
			dst_pitch = dst->GetPitch(PLANAR_V);
			rowsize = dst->GetRowSize(PLANAR_V);
			height = dst->GetHeight(PLANAR_V);

			if (!show)
				Fieldcopy(dstp, srcp, rowsize, height, dst_pitch, src_pitch);
			else
				FieldNull(dstp, rowsize, height, dst_pitch, VERDE2);

		}
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

const AVS_Linkage* AVS_linkage = 0;
extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors) {
	AVS_linkage = vectors;
    env->AddFunction("UnsharpHQ", "c[THRESHOLD]i[SHARPSTR]f[SMOOTH]f[SHOW]b[MODE]i", Create_UnsharpHQ, 0);

    // c - Video Clip
    // i - Integer number
    // f - Float number
    // s - String
    // b - boolean

    return "`UnsharpHQ' UnsharpHQ plugin";
 }
