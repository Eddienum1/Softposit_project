
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

#include "platform.h"
#include "internals.h"

typedef unsigned __int128 uint128_t;

posit64_t p64_mul( posit64_t pA, posit64_t pB ){


	union ui64_p64 uA, uB, uZ;
	uint_fast64_t uiA, uiB;
	uint_fast64_t regA, fracA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast64_t expA;
	int_fast8_t kA=0;
	uint128_t frac128Z;

	uA.p = pA;
	uiA = uA.ui;
	uB.p = pB;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif
	//NaR or Zero
	if ( uiA==0x8000000000000000 || uiB==0x8000000000000000 ){

#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x8000000000000000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x8000000000000000;
#endif
		return uZ.p;
	}
	else if (uiA==0 || uiB==0){
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

	if(signA) uiA = (-uiA & 0xFFFFFFFFFFFFFFFF);
	if(signB) uiB = (-uiB & 0xFFFFFFFFFFFFFFFF);

	regSA = signregP64UI(uiA);
	regSB = signregP64UI(uiB);

	tmp = (uiA<<2)&0xFFFFFFFFFFFFFFFF;
	if (regSA){
		while (tmp>>63){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>63)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
		}
		tmp&=0x7FFFFFFFFFFFFFFF;
	}
	expA = tmp>>61; //to get 2 bits
	fracA = ((tmp<<1) | 0x4000000000000000) & 0x7FFFFFFFFFFFFFFF;

	tmp = (uiB<<2)&0xFFFFFFFFFFFFFFFF;
	if (regSB){
		while (tmp>>63){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
		}
	}
	else{
		kA--;
		while (!(tmp>>63)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
		}
		tmp&=0x7FFFFFFFFFFFFFFF;
	}
	expA += tmp>>61;
	frac128Z = (uint128_t) fracA * (((tmp<<1) | 0x4000000000000000) & 0x7FFFFFFFFFFFFFFF);

	if (expA>3){
		kA++;
		expA&=0x3; // -=4
	}

	rcarry = frac128Z>>125;//3rd bit of frac128Z
	if (rcarry){
		expA++;
		if (expA>3){
			kA ++;
			expA&=0x3;
		}
		frac128Z>>=1;
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x4000000000000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFFFFFFFFFF - (0x7FFFFFFFFFFFFFFF>>regA);
	}


	if(regA>62){
		//max or min pos. exp and frac does not matter.
		(regSA) ? (uZ.ui= 0x7FFFFFFFFFFFFFFF): (uZ.ui=0x1);
	}
	else{
		//remove carry and rcarry bits and shift to correct position (2 bits exp, so + 1 than 16 bits)
		frac128Z = (frac128Z & (((uint128_t)1 << 124) - 1)) >> regA;
		fracA = (uint_fast64_t) (frac128Z>>64);
		if (regA<=60){
			bitNPlusOne |= (0x8000000000000000 & frac128Z);
			expA<<= (60-regA);
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
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}

		}
		//sign is always zero
		uZ.ui = packToP64UI(regime, expA, fracA);
		//n+1 frac bit is 1. Need to check if another bit is 1 too if not round to even
		if (bitNPlusOne){
			if (0x7FFFFFFFFFFFFFFF & frac128Z)  bitsMore=1;
			uZ.ui += (uZ.ui&1) | bitsMore;
		}
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFF;
	return uZ.p;

}

