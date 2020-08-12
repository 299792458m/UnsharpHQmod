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
#include <intrin.h>

#define VERDE1 100
#define VERDE2 100


#define uc unsigned char

#define threshold 20
#define sharp_str 4
#define smoot_str .5

#include "utiles.h"

template<int bits_per_pixel>
void FieldUnSharpY_HBD(unsigned char* dstp0, const unsigned char* srcp, int rowsize, int height, int src_pitch, int dst_pitch,
	float thres, float MIN, float MAX, int mode, int opt, IScriptEnvironment* env)
{
	constexpr int max_pixel_size = (1 << bits_per_pixel);
	constexpr int max_pixel_value = max_pixel_size - 1;
	constexpr int cnvsft8toHBD = bits_per_pixel - 8;
	//__declspec(align(16)) int32_t funcion[max_pixel_size];
	int32_t* funcion= (int32_t*)_aligned_malloc(max_pixel_size * sizeof(int32_t), 32);	//サイズが大きい(16bitでは256kBになる)ためヒープへ

	int ocho = 1 << PRECISION;	//意味不明 不要にできるのでは？

	//MIN = -(SMOOTH+1.0)で、MAX = SHAPSTR+1.0になってる

	double giro, atanMIN, atanMAX, A, S;
	float fmin, fmax;	//MIN/MAXは1ずつ拡張された値が渡される→呼び出し元を変更して拡張しないようにした 必要に応じて↓で調整


	//i==giroで1になるようにmode1をオフセットしてそれっぽくしてるが、形はi=giroでatan(0)の形になってる(smooth側があまり考慮されてない)
	fmin = MIN * ocho / 4 - 0.99999f;	//-1.0にすると演算精度上-1.0001とかになるので
	fmax = MAX * ocho / 4;

	giro = thres;
	atanMIN = atan((0 - giro) / curva);			//0なので括弧でくくっても同じだけど
	atanMAX = atan((255 - giro) / curva);		//funcionが(i-giro)/curvであることを考えると、ここも(i-giro)が正しい
	A = (fmax - fmin) / (atanMAX - atanMIN);		//(atanMIN,atanMAX)→(MIN,MAX)の傾き
	S = fmin - A * atanMIN - 1.0;					//(atanMIN,atanMAX)→(MIN,MAX)のオフセット

	for (int ii = 0; ii < bits_per_pixel; ii++) {				//i=giroでfuncionが1になるように繰り返し演算 精度はよく知らんがもし二次収束するならこんなもん
		giro = thres - tan(-S / A) * curva;
		atanMIN = atan((0 - giro) / curva);
		atanMAX = atan((255 - giro) / curva);
		A = (fmax - fmin) / (atanMAX - atanMIN);
		S = fmin - A * atanMIN - 1.0;
	}
	for (int i = 0; i < max_pixel_size; i++) {
		funcion[i] = (__int32)(A * atan(((double)i/(1<< cnvsft8toHBD) - giro) / curva) + S + 1.0);	//ここでHBD換算する
		//test[i]  = ( A*atan( (i -giro) / curva  ) + S +1.0);
	}


	const uint16_t* row0 = (uint16_t*)srcp;		//これが安全かどうかは知らない・・・
	const uint16_t* row1 = (uint16_t*)(srcp + src_pitch);
	const uint16_t* row2 = (uint16_t*)(srcp + src_pitch * 2);
	const uint16_t* row3 = (uint16_t*)(srcp + src_pitch * 3);
	const uint16_t* row4 = (uint16_t*)(srcp + src_pitch * 4);

	uint16_t* dstp = (uint16_t*)dstp0;	//16bitへキャスト

	memcpy((BYTE*)dstp, (BYTE*)row0, rowsize);		//1,2行目はソースをコピー
	dstp = dstp + dst_pitch / sizeof(uint16_t);		//pitchはbyte単位なので
	memcpy((BYTE*)dstp, (BYTE*)row1, rowsize);
	dstp = dstp + dst_pitch / sizeof(uint16_t);

	int width = rowsize / sizeof(uint16_t);
	for (int h = 2; h < height - 2; h++) {
		{
			dstp[0] = row2[0];			//1,2列目はソースをコピー
			dstp[1] = row2[1];
			int w = 2;
			if (opt>=3) {	//AVX2
				for (; w < width - 9; w+=8) {

					//		- - - - -		- - i - -	0
					//		- a b c -		- j - k -	1
					//		- d X e -	=)	l - X - m	2
					//		- f g h -		- n - o -	3
					//		- - - - -		- - p - -	4

					__declspec(align(32)) int32_t DIFs[8];	//max_pixel_valueまで
					__declspec(align(32)) int32_t Bs[8];

					__m256i avrg, diff;
					__m256i j32, k32, n32, o32;
					__m256i i32, p32, l32, m32;

					// AVG -- 33
					auto xmm0 = _mm_loadu_si128((const __m128i*)(row0 + w));		//i
					auto xmm1 = _mm_loadu_si128((const __m128i*)(row4 + w));		//p
					auto xmm2 = _mm_loadu_si128((const __m128i*)(row2 + w - 2));	//l
					auto xmm3 = _mm_loadu_si128((const __m128i*)(row2 + w + 2));	//m

					i32 = _mm256_cvtepu16_epi32(xmm0);			//int16->int32
					p32 = _mm256_cvtepu16_epi32(xmm1);
					l32 = _mm256_cvtepu16_epi32(xmm2);
					m32 = _mm256_cvtepu16_epi32(xmm3);

					xmm0 = _mm_loadu_si128((const __m128i*)(row1 + w - 1));	//j
					xmm1 = _mm_loadu_si128((const __m128i*)(row1 + w + 1));	//k
					xmm2 = _mm_loadu_si128((const __m128i*)(row3 + w - 1));	//n
					xmm3 = _mm_loadu_si128((const __m128i*)(row3 + w + 1));	//o

					j32 = _mm256_cvtepu16_epi32(xmm0);	//int16->int32	_mm256_unpacklo_epi16(xmm0,zero)を2度やるのとどっちが速いか・・・
					k32 = _mm256_cvtepu16_epi32(xmm1);
					n32 = _mm256_cvtepu16_epi32(xmm2);
					o32 = _mm256_cvtepu16_epi32(xmm3);

					avrg = _mm256_add_epi32(i32, p32);
					avrg = _mm256_add_epi32(avrg, l32);
					avrg = _mm256_add_epi32(avrg, m32);
					avrg = _mm256_add_epi32(avrg, j32);
					avrg = _mm256_add_epi32(avrg, k32);
					avrg = _mm256_add_epi32(avrg, n32);
					avrg = _mm256_add_epi32(avrg, o32);

					avrg = _mm256_srli_epi32(avrg, 3);		//>>3(/8個平均)


					// DIF -- 61
					auto abs = _mm256_sub_epi32(i32, avrg);	//i
					diff = _mm256_abs_epi32(abs);			//minmaxでやるよりabsのほうが速い？

					abs = _mm256_sub_epi32(p32, avrg);		//p
					abs = _mm256_abs_epi32(abs);
					diff = _mm256_add_epi32(diff, abs);

					abs = _mm256_sub_epi32(l32, avrg);		//l
					abs = _mm256_abs_epi32(abs);
					diff = _mm256_add_epi32(diff, abs);

					abs = _mm256_sub_epi32(m32, avrg);		//m
					abs = _mm256_abs_epi32(abs);
					diff = _mm256_add_epi32(diff, abs);

					abs = _mm256_sub_epi32(j32, avrg);		//j
					abs = _mm256_abs_epi32(abs);
					diff = _mm256_add_epi32(diff, abs);

					abs = _mm256_sub_epi32(k32, avrg);		//k
					abs = _mm256_abs_epi32(abs);
					diff = _mm256_add_epi32(diff, abs);

					abs = _mm256_sub_epi32(n32, avrg);		//n
					abs = _mm256_abs_epi32(abs);
					diff = _mm256_add_epi32(diff, abs);

					abs = _mm256_sub_epi32(o32, avrg);		//o
					abs = _mm256_abs_epi32(abs);
					diff = _mm256_add_epi32(diff, abs);	//|i-avrg|+|p-avrg|+|l-avrg|+|m-avrg|+|j-avrg|+|k-avrg|+|n-avrg|+|o-avrg|

					diff = _mm256_min_epu32(diff, _mm256_set1_epi32(max_pixel_value));	//max_pixel_valueに飽和 0方向はpackusで実施
					//diff = _mm256_max_epu32(diff, _mm256_setzero_si256()); いらなくね？

					//xmm0 = _mm_packus_epi32(_mm256_extracti128_si256(diff, 0), _mm256_extracti128_si256(diff, 1));	//diff uint16 飽和付
					//_mm_storeu_si128((__m128i*)DIFs, xmm0);
					_mm256_storeu_si256((__m256i*)DIFs, diff);


					// FUNCION 44
					for (int i = 0; i < 8; i++) {
						Bs[i] = funcion[DIFs[i]];	//Bは負になりうるので注意
					}

					__m256i ymm0,ymm1;
					ymm0 = _mm256_set1_epi32(ocho);				//ocho
					auto B = _mm256_loadu_si256((__m256i*)Bs);	//B
					auto A = _mm256_add_epi32(B, ymm0);			//A=B+ocho
					xmm0 = _mm_loadu_si128((const __m128i*)(row2 + w));	//X(中央値)
					auto X = _mm256_cvtepu16_epi32(xmm0);
					ymm0 = _mm256_mullo_epi32(A, X);			//X * A			org*(ocho + funcion[diff]) 結果の下位32bitで十分？
					ymm1 = _mm256_mullo_epi32(avrg, B);			//AVG * B		avr*funcion[diff]
					ymm0 = _mm256_sub_epi32(ymm0, ymm1);		//X*A - AVG*B
					ymm0 = _mm256_max_epi32(ymm0, _mm256_setzero_si256()); //マイナスになる時はここでリミットしないとおかしくなる
					ymm0 = _mm256_srli_epi32(ymm0, PRECISION);	//dst=org + funcion[diff] *(org-avr)/ocho = ( org*(ocho + funcion[diff]) - avr*funcion[diff] )	/ocho
					ymm0 = _mm256_min_epu32(ymm0, _mm256_set1_epi32(max_pixel_value));	//max_pixel_valueに飽和 0方向はpackusで実施
					xmm0 = _mm_packus_epi32(_mm256_extracti128_si256(ymm0, 0), _mm256_extracti128_si256(ymm0, 1));	//int32 -> int16 with unsigned satulation
					_mm_storeu_epi16((__m128i*)(dstp + w), xmm0);
				}
			}
			
			for (; w < width - 2; w++) {
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
					+ abs(avrg - row4[w]);			 //
				diff = min(max(diff, 0), max_pixel_value);

				B = funcion[diff];
				A = B + ocho;		// (1+B/ocho)*ocho
					
				dstp[w] = min(max((row2[w] * A - avrg * B) / ocho, 0), max_pixel_value);
			}

			dstp[width - 2] = row2[width - 2];			//last-1,last列はソースをコピー
			dstp[width - 1] = row2[width - 1];
		}

		row0 = row1;
		row1 = row2;
		row2 = row3;
		row3 = row4;
		row4 += src_pitch/ sizeof(uint16_t);

		dstp += dst_pitch/ sizeof(uint16_t);
	}
	memcpy((BYTE*)dstp, (BYTE*)row2, rowsize);			//last-1,last行目はソースをコピー
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
	int opt;
	int cpuf;

public:
	UnsharpHQ(PClip _child,int , float ,float ,bool ,int, int ,IScriptEnvironment* env);
	~UnsharpHQ();
	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);

};



UnsharpHQ::UnsharpHQ(PClip _child,int _thres,float _MAX,float _MIN,bool _show,int _mode, int _opt, IScriptEnvironment* env) :
GenericVideoFilter(_child) {
	thres = _thres;
	MIN = _MIN;	//SMOOTH
	MAX = _MAX;	//SHARPSTR
	show = _show;
	mode = _mode;
	cpuf = env->GetCPUFlags();
	
	if ((_opt == 0) || (_opt == 3))	//0:auto 1:c 2:sse2 3:AVX2
		opt = (cpuf & CPUF_AVX2) ? 3 : (cpuf & CPUF_SSE2) ? 2 : 1;	//AVX2 or SSE2
	else if (_opt == 2)
		opt = (cpuf & CPUF_SSE2) ? 2 : 1;							//SSE2
	else
		opt = 1;													//C


	if( !vi.IsYUV() || !vi.IsPlanar() )
		env->ThrowError("UnSharpHQ: supports only planar YUV(YV12,YUV420..)");	//αチャンネルもコピーしてないのでYUVAもダメ

	if (vi.BitsPerComponent() > 16)
		env->ThrowError("UnsharpHQ: supports only 8-16bit integer input. not float!");

	if( thres<0 || thres>99 || MIN<0 || MIN>4 || MAX<0 || MAX>99)
		env->ThrowError("UnsharpHQ(THRESHOLD=20, SHARPSTR=4.0 , SMOOTH=0.5 , SHOW=FALSE)  \n "
		"        USAGE:\n"
		"  THRESHOLD [ 0.0 , 99.0]\n"
		"   SHARPSTR [ 0.0 , 99.0]\n"
		"     SMOOTH [ 0.0 ,  4.0]\n"
	);
	//if( mode >=2 ) env->ThrowError("MODE is 0 to 1");
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

			FieldUnSharpY(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, show, mode, opt);

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
			int dst_pitch = dst->GetPitch(PLANAR_Y);	//次の行までのbyte数(パディング含む)
			int rowsize = dst->GetRowSize(PLANAR_Y);	//1行のbyte数
			int height = dst->GetHeight(PLANAR_Y);		//何列あるか(pixel数)

			switch (bits_per_pixel) {
			case 10:	FieldUnSharpY_HBD<10>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, mode, opt, env); break;	//mode,showは削除
			case 12:	FieldUnSharpY_HBD<12>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, mode, opt, env); break;
			case 14:	FieldUnSharpY_HBD<14>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, mode, opt, env); break;
			case 16:
			default:	FieldUnSharpY_HBD<16>(dstp, srcp, rowsize, height, src_pitch, dst_pitch, (float)thres, -MIN, MAX, mode, opt, env); break;
			}
			
			//- - - - - - - - U - - - - - - - - - -  

			srcp = src->GetReadPtr(PLANAR_U);
			dstp = dst->GetWritePtr(PLANAR_U);
			src_pitch = src->GetPitch(PLANAR_U);
			dst_pitch = dst->GetPitch(PLANAR_U);
			rowsize = dst->GetRowSize(PLANAR_U);
			height = dst->GetHeight(PLANAR_U);

			Fieldcopy(dstp, srcp, rowsize, height, dst_pitch, src_pitch);

			//- - - - - - - - V - - - - - - - - - -  
			srcp = src->GetReadPtr(PLANAR_V);
			dstp = dst->GetWritePtr(PLANAR_V);
			src_pitch = src->GetPitch(PLANAR_V);
			dst_pitch = dst->GetPitch(PLANAR_V);
			rowsize = dst->GetRowSize(PLANAR_V);
			height = dst->GetHeight(PLANAR_V);

			Fieldcopy(dstp, srcp, rowsize, height, dst_pitch, src_pitch);

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
		, args[6].AsInt(0)
		, env);
}

const AVS_Linkage* AVS_linkage = 0;
extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors) {
	AVS_linkage = vectors;
    env->AddFunction("UnsharpHQ", "c[THRESHOLD]i[SHARPSTR]f[SMOOTH]f[SHOW]b[MODE]i[OPT]i", Create_UnsharpHQ, 0);

    // c - Video Clip
    // i - Integer number
    // f - Float number
    // s - String
    // b - boolean

    return "`UnsharpHQ' UnsharpHQ plugin";
 }
