/*============================================================================

This C source file is part of the SoftPosit Posit Arithmetic Package
by S. H. Leong (Cerlane).

Copyright 2017, 2018 A*STAR.  All rights reserved.

This C source file was based on SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3d, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016, 2017 The Regents of the
University of California.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions, and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions, and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the University nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

#include <stdlib.h>

#include "platform.h"
#include "internals.h"

typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

typedef struct {
    int128_t quot;
    int128_t rem;
} int128div_t;

int128div_t int128div(int128_t numer, int128_t denom) {
    int128div_t result;
    result.quot = numer / denom;
    result.rem  = numer % denom;
    return result;
}

posit64_t p64_div( posit64_t pA, posit64_t pB )
{
    union ui64_p64 uA, uB, uZ;
    uint_fast64_t uiA, uiB, fracA, fracB, regA, regime, regB, tmp;
    bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast8_t kA=0;
	int_fast64_t expA;
	uint128_t frac128A, frac128Z, rem;
	int128div_t divresult;

	uA.p = pA;
	uiA = uA.ui;
	uB.p = pB;
	uiB = uB.ui;

	//Zero or infinity
	if ( uiA==0x8000000000000000ULL || uiB==0x8000000000000000ULL || uiB==0){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x8000000000000000ULL;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x8000000000000000ULL;
#endif
		return uZ.p;
	}
	else if (uiA==0){
#ifdef SOFTPOSIT_EXACT

		uZ.ui.v = 0;
		if ( (uiA==0 && uiA.ui.exact) || (uiB==0 && uiB.ui.exact) )
			uZ.ui.exact = 1;
		else
			uZ.ui.exact = 0;
#else
		uZ.ui = 0;
#endif
		return uZ.p;
	}

	signA = signP64UI( uiA );
	signB = signP64UI( uiB );
	signZ = signA ^ signB;
	if(signA) uiA = (-uiA & 0xFFFFFFFFFFFFFFFFULL);
	if(signB) uiB = (-uiB & 0xFFFFFFFFFFFFFFFFULL);
	regSA = signregP64UI(uiA);
	regSB = signregP64UI(uiB);

	tmp = (uiA<<2)&0xFFFFFFFFFFFFFFFFULL;
	if (regSA){
		while (tmp>>63){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFFULL;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>63)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFFULL;
		}
		tmp&=0x7FFFFFFFFFFFFFFFULL;
	}
	expA = tmp>>61; //to get 2 bits
	fracA = ((tmp<<1) | 0x4000000000000000ULL) & 0x7FFFFFFFFFFFFFFFULL;
	frac128A = (uint128_t) fracA << 62;

	tmp = (uiB<<2)&0xFFFFFFFFFFFFFFFFULL;
	if (regSB){
		while (tmp>>63){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFFULL;
		}
	}
	else{
		kA++;
		while (!(tmp>>63)){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFFULL;
		}
		tmp&=0x7FFFFFFFFFFFFFFFULL;
	}
	expA -= tmp>>61;
	fracB = ((tmp<<1) | 0x4000000000000000ULL) & 0x7FFFFFFFFFFFFFFFULL;

	divresult = int128div(frac128A,(uint128_t)fracB);
	frac128Z = divresult.quot;
	rem = divresult.rem;

	if (expA<0){
		expA+=4;
		kA--;
	}
	if (frac128Z!=0){
		rcarry = frac128Z >> 62; // this is the hidden bit (14th bit) , extreme right bit is bit 0
		if (!rcarry){
			if (expA==0){
				kA--;
				expA=3;
			}
			else
				expA--;
			frac128Z<<=1;
		}
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x4000000000000000ULL>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFFFFFFFFFFULL - (0x7FFFFFFFFFFFFFFFULL>>regA);
	}
	if(regA>62){
		//max or min pos. exp and frac does not matter.
		(regSA) ? (uZ.ui= 0x7FFFFFFFFFFFFFFFULL): (uZ.ui=0x1);
	}
	else{
		//remove carry and rcarry bits and shift to correct position
		frac128Z &= 0x3FFFFFFFFFFFFFFFULL;

		fracA = (uint_fast64_t)frac128Z >> (regA+2);

		if (regA<=60){
			bitNPlusOne = (frac128Z >> (regA +1)) & 0x1;
			expA<<= (60-regA);
			if (bitNPlusOne) ( ((1<<(regA+1))-1) & frac128Z ) ? (bitsMore=1) : (bitsMore=0);
		}
		else {
			if (regA==62){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==61){
				bitNPlusOne = expA&0x1;
				expA>>=1; //taken care of by the pack algo
			}
			if (frac128Z>0){
				fracA=0;
				bitsMore =1;
			}

		}
		if (rem) bitsMore =1;

		uZ.ui = packToP64UI(regime, expA, fracA);
		if (bitNPlusOne) uZ.ui += (uZ.ui&1) | bitsMore;
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFFULL;
	return uZ.p;
}

