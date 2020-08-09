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
    
#include "math.h"

#define PRECISION 3 //osea 4
#define curva 30


void FieldUnSharpY(unsigned char* dstp,const unsigned char* srcp
				  ,int rowsize,int height,int src_pitch,int dst_pitch,
				  float thres,float MIN,float MAX,bool show,int mode)
{

	//__int16 ocho = (__int16) pow( (float)2, (float)PRECISION);
	__int16 ocho = (1<<PRECISION);
	__int16 OCHO[4]; for (int i=0;i<4;i++) OCHO[i]= ocho;
	
	//MIN = -(SMOOTH+1.0)�ŁAMAX = SHAPSTR+1.0�ɂȂ��Ă�

	double giro,atanMIN,atanMAX,A,S;
	__declspec(align(16)) __int16 funcion[256];
	//volatile double test[256];
	float fmin,fmax;	//MIN/MAX��1���g�����ꂽ�l���n����遨�Ăяo������ύX���Ċg�����Ȃ��悤�ɂ��� �K�v�ɉ����ā��Œ���


	if (mode==0){	//org�݊�
		fmin=MIN-1.0f;
		fmax=MAX+1.0f;

		giro = thres;
		atanMIN = atan(0-giro/curva);
		atanMAX = atan(256-giro/curva);
		A= (fmax-fmin) / (atanMAX-atanMIN) ;
		S= fmin - A*atanMIN ;
	
		giro= giro - tan(-S/A)*curva;
		atanMIN = atan(0-giro/curva);
		atanMAX = atan(256-giro/curva);
		A= (fmax-fmin) / (atanMAX-atanMIN) ;
		S= fmin - A*atanMIN ;

		for (int i=0;i<256;i++){
			funcion[i]= (__int16)(((__int16) A*atan( ((float)i -giro) / curva  ) + S))*ocho/4;	//short��S�ȊO�Ɋ|���� �o�O�H �d�݂�ύX�������߈ȑO�ƌ݊���ۂ��߂Ɋ��int�Ő؂�̂ĂĂ�
			//test[i] = (A*atan( ((float)i -giro) / curva  ) + S)*ocho/4;
		}
	}
	//else if (mode==1){	//i==giro��0�ɂȂ�悤�� �o�͂�round���Ă�̂ł�����Ɣ����E�E�E
	//	fmin=MIN*ocho/4-0.49999f;
	//	fmax=MAX*ocho/4;

	//	giro = thres;
	//	atanMIN = atan((0-giro)/curva);			//0�Ȃ̂Ŋ��ʂł������Ă�����������
	//	atanMAX = atan((255-giro)/curva);		//funcion��(i-giro)/curv�ł��邱�Ƃ��l����ƁA������(i-giro)�������� fncion��giro/curv�ɂ���Ȃ琳�������E�E�E ���ƁA255�܂ł���ˁH
	//	A= (fmax-fmin) / (atanMAX-atanMIN) ;		//(atanMIN,atanMAX)��(MIN,MAX)�̌X��
	//	S= fmin - A*atanMIN ;					//(atanMIN,atanMAX)��(MIN,MAX)�̃I�t�Z�b�g

	//	for (int ii=0;ii<8;ii++){				//i=giro��funcion��0�ɂȂ�悤�ɌJ��Ԃ����Z
	//		giro= thres - tan(-S/A)*curva;	
	//		atanMIN = atan((0-giro)/curva);
	//		atanMAX = atan((255-giro)/curva);
	//		A= (fmax-fmin) / (atanMAX-atanMIN) ;
	//		S= fmin - A*atanMIN ;
	//	}
	//	for (int i=0;i<256;i++){
	//		funcion[i]  = (__int16)round( A*atan( (i -giro) / curva  ) + S );	//���������������Ǝv��
	//		//test[i]  = ( A*atan( (i -giro) / curva  ) + S );
	//	}
	//}
	else{	//i==giro��1�ɂȂ�悤��mode1���I�t�Z�b�g���Ă�����ۂ����Ă邪�A�`��i=giro��atan(0)�̌`�ɂȂ��Ă�(smooth�������܂�l������ĂȂ�)
		fmin=MIN*ocho/4-0.99999f;	//-1.0�ɂ���Ɖ��Z���x��-1.0001�Ƃ��ɂȂ�̂�
		fmax=MAX*ocho/4;

		giro = thres;
		atanMIN = atan((0-giro)/curva);			//0�Ȃ̂Ŋ��ʂł������Ă�����������
		atanMAX = atan((255-giro)/curva);		//funcion��(i-giro)/curv�ł��邱�Ƃ��l����ƁA������(i-giro)�������� fncion��giro/curv�ɂ���Ȃ琳�������E�E�E ���ƁA255�܂ł���ˁH
		A= (fmax-fmin) / (atanMAX-atanMIN) ;		//(atanMIN,atanMAX)��(MIN,MAX)�̌X��
		S= fmin - A*atanMIN -1.0;					//(atanMIN,atanMAX)��(MIN,MAX)�̃I�t�Z�b�g

		for (int ii=0;ii<8;ii++){				//i=giro��funcion��1�ɂȂ�悤�ɌJ��Ԃ����Z
			giro= thres - tan(-S/A)*curva;	
			atanMIN = atan((0-giro)/curva);
			atanMAX = atan((255-giro)/curva);
			A= (fmax-fmin) / (atanMAX-atanMIN) ;
			S= fmin - A*atanMIN -1.0;
		}
		for (int i=0;i<256;i++){
			funcion[i]  = (__int16)(A*atan( (i -giro) / curva  ) + S +1.0);
			//test[i]  = ( A*atan( (i -giro) / curva  ) + S +1.0);
		}
	}



	const unsigned char* row0=srcp;
	const unsigned char* row1=srcp+src_pitch;
	const unsigned char* row2=srcp+2*src_pitch;
	const unsigned char* row3=srcp+3*src_pitch;
	const unsigned char* row4=srcp+4*src_pitch;
	
	memcpy( (BYTE*) dstp, (BYTE*) row0, rowsize);			//1,2�s�ڂ̓\�[�X���R�s�[
	dstp = dstp + dst_pitch;
	memcpy( (BYTE*) dstp, (BYTE*) row1, rowsize);
	dstp = dstp + dst_pitch;

	int w;
	uc AVGs[8];
	uc DIFs[8];

	for (int h=2; h < height-2; h++){
#ifdef _WIN64	//�œK���͂��ĂȂ�
		dstp[0] = row2[0];			//1,2��ڂ̓\�[�X���R�s�[
		dstp[1] = row2[1];
		for (w = 2; w < rowsize - 2; w++) {
			int avrg, diff, A, B;

			avrg = row0[w] + row4[w];
			avrg += row1[w + 1] + row2[w + 2] + row3[w + 1];
			avrg += row1[w - 1] + row2[w - 2] + row3[w - 1];
			avrg = avrg / 8;

			diff = abs(avrg - row0[w]) + abs(avrg - row4[w]); //asm�ł�abs��2�̑傫����-���������ł���Ă�
			diff += abs(avrg - row1[w + 1]) + abs(avrg - row2[w + 2]) + abs(avrg - row3[w + 1]);
			diff += abs(avrg - row1[w - 1]) + abs(avrg - row2[w - 2]) + abs(avrg - row3[w - 1]);
			diff = min(max(diff, 0), 255);

			B = funcion[diff];
			A = B + ocho;		// (1+B/ocho)*ocho

			dstp[w] = min(max((row2[w] * A - avrg * B) / ocho, 0), 255);
		}
		dstp[rowsize - 2] = row2[rowsize - 2];			//last-1,last��̓\�[�X���R�s�[
		dstp[rowsize - 1] = row2[rowsize - 1];
#else
		//0, 1 , height-1 , height, los copio
		dstp[0] = row2[0];			//1,2��ڂ̓\�[�X���R�s�[
		dstp[1] = row2[1];
		for (w = 2; w < rowsize - 2; w += 8) {

			//		- - - - -		- - o - -
			//		- o o o -		- o - o -
			//		- o X o -	=)	o - X - o
			//		- o o o -		- o - o -
			//		- - - - -		- - o - -

			//		- - - - -		- - i - -	0
			//		- a b c -		- j - k -	1
			//		- d X e -	=)	l - X - m	2
			//		- f g h -		- n - o -	3
			//		- - - - -		- - p - -	4

			__asm {
				// AVG -- 33
				pxor xmm0, xmm0;
				mov eax, row1;
				mov ebx, w;
				mov ecx, row3;

				movq xmm1, qword ptr[eax + ebx - 1];	//j[0-7]	�l��
				movq xmm2, qword ptr[eax + ebx + 1];	//k
				movq xmm3, qword ptr[ecx + ebx - 1];	//n
				movq xmm4, qword ptr[ecx + ebx + 1];	//o
				punpcklbw xmm1, xmm0;		//byte->word
				punpcklbw xmm2, xmm0;
				punpcklbw xmm3, xmm0;
				punpcklbw xmm4, xmm0;

				paddw xmm1, xmm2;					//j+k+n+o
				paddw xmm1, xmm3;
				paddw xmm1, xmm4;

				mov eax, row0;
				mov ecx, row4;
				mov edx, row2;

				movq xmm2, qword ptr[eax + ebx];		//i
				movq xmm3, qword ptr[ecx + ebx];		//p
				movq xmm4, qword ptr[edx + ebx + 2];	//l
				movq xmm5, qword ptr[edx + ebx - 2];	//m
				punpcklbw xmm2, xmm0;		//byte->word
				punpcklbw xmm3, xmm0;
				punpcklbw xmm4, xmm0;
				punpcklbw xmm5, xmm0;

				paddw xmm1, xmm2;					//j+k+n+o +i+p+l+m
				paddw xmm1, xmm3;
				paddw xmm1, xmm4;
				paddw xmm1, xmm5;

				psraw xmm1, 3;						// >>3(/8����)
				packuswb xmm1, xmm1;		//word->byte
				movq qword ptr[AVGs], xmm1;		//avrg

			// - - - - - - - - - - - - - - - - - - - - - - - - - - -
			// DIF -- 61

				movq mm1, qword ptr[AVGs];
				pxor mm0, mm0;
				//MM2 = PIXEL DE AL LADO ���̃s�N�Z��

				movq mm2, qword ptr[eax + ebx];		//i
				movq mm3, mm1;						//avrg
				movq mm4, mm2;
				pmaxub mm2, mm1;					//mm2=max(i,avrg)
				pminub mm3, mm4;					//mm3=min(i,avrg)
				psubusb mm2, mm3;					// |i-avrg|
				paddusb   mm0, mm2;				//mmm0��zero�Ȃ̂ŁAmov�Ɠ���
				movq mm2, qword ptr[ecx + ebx];		//p
				movq mm3, mm1;						//avrg
				movq mm4, mm2;
				pmaxub mm2, mm1;
				pminub mm3, mm4;
				psubusb mm2, mm3;
				paddusb   mm0, mm2;
				movq mm2, qword ptr[edx + ebx + 2];	//l
				movq mm3, mm1;
				movq mm4, mm2;
				pmaxub mm2, mm1;
				pminub mm3, mm4;
				psubusb mm2, mm3;
				paddusb   mm0, mm2;
				movq mm2, qword ptr[edx + ebx - 2];	//m
				movq mm3, mm1;
				movq mm4, mm2;
				pmaxub mm2, mm1;
				pminub mm3, mm4;
				psubusb mm2, mm3;
				paddusb   mm0, mm2;				//|i-avrg|+|p-avrg|+|l-avrg|+|m-avrg|

				mov eax, row1;
				mov ecx, row3;

				movq mm2, qword ptr[eax + ebx - 1];	//j
				movq mm3, mm1;
				movq mm4, mm2;
				pmaxub mm2, mm1;
				pminub mm3, mm4;
				psubusb mm2, mm3;
				paddusb   mm0, mm2;
				movq mm2, qword ptr[eax + ebx + 1];	//k
				movq mm3, mm1;
				movq mm4, mm2;
				pmaxub mm2, mm1;
				pminub mm3, mm4;
				psubusb mm2, mm3;
				paddusb   mm0, mm2;
				movq mm2, qword ptr[ecx + ebx - 1];	//n
				movq mm3, mm1;
				movq mm4, mm2;
				pmaxub mm2, mm1;
				pminub mm3, mm4;
				psubusb mm2, mm3;
				paddusb   mm0, mm2;
				movq mm2, qword ptr[ecx + ebx + 1];	//o
				movq mm3, mm1;
				movq mm4, mm2;
				pmaxub mm2, mm1;
				pminub mm3, mm4;
				psubusb mm2, mm3;
				paddusb   mm0, mm2;				//|i-avrg|+|p-avrg|+|l-avrg|+|m-avrg|+|j-avrg|+|k-avrg|+|n-avrg|+|o-avrg|	0/255�Ń��~�b�g

				movq qword ptr[DIFs], mm0;
				emms;
			}
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			// FUNCION 44

			__declspec(align(16)) __int16 B[8];
			for (int i = 0; i < 8; i++) {
				B[i] = funcion[DIFs[i]];
			}

			if (show)
				__asm {
				mov edi, dstp;
				//movq xmm1 , qword ptr[B];
				//movhpd xmm1 , qword ptr[B+8];
				movdqu xmm1, B;

				packuswb xmm1, xmm1;
				movq qword ptr[edi + ebx], xmm1;
				emms;
			}
			else
				__asm {
				//ahora xmm1 = B	��
				//      xmm4 = A = B+4			//�����ǂ����̃o�[�W�������炩16�ł͂Ȃ�4�ɂȂ��Ă�(���̃R�����g����������Ή��Z���̃I�[�o�t���[���)
				movddup xmm6, qword ptr[OCHO];	//[0 4]*4*2 double����Ȃ����������ǁA�o���邩��悢
				//pxor xmm0 , xmm0;

				//movq xmm1 , qword ptr[B];
				//movhpd xmm1 , qword ptr[B+8];
				movdqu xmm1, B;					//sse2�g����1�R�}���h�ɂȂ邯�Ǒ��x�͂��܂�ς��Ȃ��E�E�E ���Ȃ݂ɍŋ߂�CPU��align����Ă��a�ł�u�ł��R�X�g�͓����炵���H

				movaps xmm5, xmm1;
				movaps xmm4, xmm1;
				paddw  xmm4, xmm6;				//A=B+ocho


				//PARTE FINAL
				//xmm1 = (medio)
				//xmm2 = (AVGs)
				//xmm4 = A
				//xmm5 = B
				mov eax, row2;
				mov ebx, w;
				mov edi, dstp;

				movq xmm1, qword ptr[eax + ebx];		//X(�����l)
				movq xmm2, qword ptr[AVGs];

				punpcklbw xmm1, xmm0;
				punpcklbw xmm2, xmm0;

				/*
									// - - - - - -		A=B+4�ɂ���΃I�[�o�[�t���[���Ȃ����炱���v��Ȃ��A�Ə����Ă��� ��������float�ɂ��邩
								// ESTA PARTE ES EL REMPLAZO , PERO MAS LARGO
								// PARA EVITAR EL OVERFLOW CON B GRANDE, PERO
								// SI USO A=B+4 EN LUGAR DE 16, SE EVITA
								// SINO TENGO QUE USAR ESTOS EN FLOAT
									movaps xmm7 , xmm1;
									movaps xmm6 , xmm2;
									//xmm1 = Medio * A
									pmullw xmm1, xmm4;
									pmulhw xmm7, xmm4;
									//xmm2 = AVG * B
									pmullw xmm2, xmm5;

									pmulhw xmm6, xmm5;

									movaps xmm0, xmm1;
									punpcklwd xmm1, xmm7;
									punpckhwd xmm0, xmm7;

									movaps xmm4, xmm2;
									punpcklwd xmm2, xmm6;
									punpckhwd xmm4, xmm6;
									//de 32 tengo xmm1,xmm0 y xmm2,xmm4
									//resto
									psubd xmm1,xmm2;
									psubd xmm0,xmm4;
									//divido
									psrad xmm1 , 4;
									psrad xmm0 , 4;

									//de int32 a int16
									packssdw xmm1 , xmm1;
									packssdw xmm0 , xmm0;
									//los uno
									movsd xmm0, xmm1;
									movaps xmm1, xmm0;
				*/

				//xmm1 = Medio * A
				pmullw xmm1, xmm4;
				//xmm2 = AVG * B
				pmullw xmm2, xmm5;
				//xmm1 = Z0 *A - AVG *B (con saturacion)	dst=org + funcion[diff] *(org-avr)/ocho = ( org*(ocho + funcion[diff]) - avr*funcion[diff] )	/ocho
				psubw  xmm1, xmm2;
				//divido por 4 (shift 2 en ints 16) y desempaco a 8
				psraw xmm1, PRECISION;
				// - - - - - - - - - - 

				packuswb xmm1, xmm1;
				movq qword ptr[edi + ebx], xmm1;
				emms;
			}

		}
		//copio los dos ultimos pixels
		dstp[rowsize - 2] = row2[rowsize - 2];			//last-1,last��̓\�[�X���R�s�[
		dstp[rowsize - 1] = row2[rowsize - 1];

#endif

		row0=row1;
		row1=row2;
		row2=row3;
		row3=row4;
		row4+=src_pitch;
		
		dstp+=dst_pitch;
	}
	memcpy((BYTE*) dstp, (BYTE*) row2, rowsize);			//last-1,last�s�ڂ̓\�[�X���R�s�[
	dstp = dstp + dst_pitch;
	memcpy((BYTE*) dstp, (BYTE*) row3, rowsize);
}


int abso(int a){
	if (a<0) return -a;
	return a;
}

void Fieldcopy(void *dest, const void *src, size_t count, 
				int rows, int dst_pitch, int src_pitch)
{
BYTE* pDest = (BYTE*) dest;
BYTE* pSrc = (BYTE*) src;
	
	for (int i=0; i < rows; i++){
		memcpy(pDest, pSrc, count);
		pSrc += src_pitch;
		pDest += dst_pitch;
	}
}
void FieldNull(unsigned char* dstp,int rowsize,int rows,int dst_pitch,uc color){
	for (int h=0; h < rows; h++){
		for (int w=0; w < rowsize; w++){
			dstp[w]=color;
		}
	dstp+=dst_pitch;
	}
}