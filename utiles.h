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
#define VERDE1 100
#define VERDE2 100


#define threshold 20
#define sharp_str 4
#define smoot_str .5


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
void FieldNull(unsigned char* dstp,int rowsize,int rows,int dst_pitch, unsigned char color){
	for (int h=0; h < rows; h++){
		for (int w=0; w < rowsize; w++){
			dstp[w]=color;
		}
	dstp+=dst_pitch;
	}
}