/*============================================================================

This C source file is part of the SoftPosit Posit Arithmetic Package
by S. H. Leong (Cerlane).

Copyright 2017, 2018 A*STAR.  All rights reserved.

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

#include <math.h>

#include <quadmath.h>

#include "platform.h"
#include "internals.h"

#ifdef SOFTPOSIT_QUAD

void checkQuadExtraP64TwoBits(__float128 f64, __float128 temp, bool *bitsNPlusOne, bool *bitsMore) {
    temp /= 2;
    if (temp <= f64) {
        *bitsNPlusOne = 1;
        f64 -= temp;
    }
    if (f64 > 0)
        *bitsMore = 1;
}

uint_fast64_t convertQuadFractionP64(__float128 f64, uint_fast16_t fracLength, bool *bitNPlusOne, bool *bitsMore) {

    uint_fast64_t frac = 0;

    if (f64 == 0)
        return 0;
    else if (f64 == INFINITY)
        return 0x8000000000000000;  // Most significant bit set

    f64 -= 1; // remove hidden bit

    if (fracLength == 0) {
        checkQuadExtraP64TwoBits(f64, 1.0, bitNPlusOne, bitsMore);
    } else {
        __float128 temp = 1;
        while (true) {
            temp /= 2;
            if (temp <= f64) {
                f64 -= temp;

                fracLength--;
                frac = (frac << 1) + 1; // shift in one

                if (f64 == 0) {
                    frac <<= fracLength; // fill remaining with 0s
                    break;
                }

                if (fracLength == 0) {
                    checkQuadExtraP64TwoBits(f64, temp, bitNPlusOne, bitsMore);
                    break;
                }
            } else {
                frac <<= 1; // shift in zero
                fracLength--;

                if (fracLength == 0) {
                    checkQuadExtraP64TwoBits(f64, temp, bitNPlusOne, bitsMore);
                    break;
                }
            }
        }
    }

    return frac;
}

posit32_t convertQuadToP32(__float128 f32){
	union ui32_p32 uZ;
	bool sign, regS;
	uint_fast32_t reg, frac=0;
	int_fast32_t exp=0;
	bool bitNPlusOne=0, bitsMore=0;

	(f32>=0) ? (sign=0) : (sign=1);

	if (f32 == 0 ){
		uZ.ui = 0;
		return uZ.p;
	}
	else if(f32 == INFINITY || f32 == -INFINITY || f32 == NAN){
		uZ.ui = 0x80000000;
		return uZ.p;
	}
	else if (f32 == 1) {
		uZ.ui = 0x40000000;
		return uZ.p;
	}
	else if (f32 == -1){
		uZ.ui = 0xC0000000;
		return uZ.p;
	}
	else if (f32 >= 1.329227995784916e+36){
		//maxpos
		uZ.ui = 0x7FFFFFFF;
		return uZ.p;
	}
	else if (f32 <= -1.329227995784916e+36){
		// -maxpos
		uZ.ui = 0x80000001;
		return uZ.p;
	}
	else if(f32 <= 7.52316384526264e-37 && !sign){
		//minpos
		uZ.ui = 0x1;
		return uZ.p;
	}
	else if(f32 >= -7.52316384526264e-37 && sign){
		//-minpos
		uZ.ui = 0xFFFFFFFF;
		return uZ.p;
	}
	else if (f32>1 || f32<-1){
		if (sign){
			//Make negative numbers positive for easier computation
			f32 = -f32;
		}

		regS = 1;
		reg = 1; //because k = m-1; so need to add back 1
		// minpos
		if (f32 <= 7.52316384526264e-37){
			uZ.ui = 1;
		}
		else{
			//regime
			while (f32>=16){
				f32 *=0.0625;  // f32/=16;
				reg++;
			}

			while (f32>=2){
				f32*=0.5;
				exp++;
			}

			int fracLength = 28-reg;
			if (fracLength<0){
				//remove hidden bit

				//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
				if(reg==29){
					bitNPlusOne = exp&0x1;
					exp>>=1; //taken care of by the pack algo
				}
				else{//reg=30
					bitNPlusOne=exp>>1;
					bitsMore=exp&0x1;
					exp=0;
				}
				if (f32!=1){
					bitsMore =1;
					frac=0;
				}
			}
			else
				frac = convertQuadFractionP32 (f32, fracLength, &bitNPlusOne, &bitsMore);

			if (reg>30 ){
				(regS) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
			}
			//rounding off fraction bits
			else{
				uint_fast32_t regime = 1;
				if (regS) regime = ( (1<<reg)-1 ) <<1;
				if (reg<=28)  exp<<= (28-reg);
				uZ.ui = ((uint32_t) (regime) << (30-reg)) + ((uint32_t) exp ) + ((uint32_t)(frac));
				uZ.ui += (bitNPlusOne & (uZ.ui&1)) | ( bitNPlusOne & bitsMore);
			}
			if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
			/*if(sign)
					uZ.ui = -uZ.ui;*/
		}
	}
	else if (f32 < 1 || f32 > -1 ){


		if (sign){
			//Make negative numbers positive for easier computation
			f32 = -f32;
		}
		regS = 0;
		reg = 0;

		while (f32<1){
			f32 *= 16;
			reg++;
		}

		while (f32>=2){
			f32*=0.5;
			exp++;
		}

		//only possible combination for reg=15 to reach here is 7FFF (maxpos) and FFFF (-minpos)
		//but since it should be caught on top, so no need to handle
		int fracLength = 28-reg;
		if (fracLength<0){
			//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
			if(reg==29){
				bitNPlusOne = exp&0x1;
				exp>>=1; //taken care of by the pack algo
			}
			else{//reg=30
				bitNPlusOne=exp>>1;
				bitsMore=exp&0x1;
				exp=0;
			}
			if (f32!=1){//because of hidden bit
				bitsMore =1;
				frac=0;
			}
		}
		else
			frac = convertQuadFractionP32 (f32, fracLength, &bitNPlusOne, &bitsMore);

		if (reg>30 ){
			(regS) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
		}
		//rounding off fraction bits
		else{
			uint_fast32_t regime = 1;
			if (regS) regime = ( (1<<reg)-1 ) <<1;
			if (reg<=28)  exp<<= (28-reg);
			uZ.ui = ((uint32_t) (regime) << (30-reg)) + ((uint32_t) exp ) + ((uint32_t)(frac));
			uZ.ui += (bitNPlusOne & (uZ.ui&1)) | ( bitNPlusOne & bitsMore);
		}
		if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	}
	else {
		//NaR - for NaN, INF and all other combinations
		uZ.ui = 0x80000000;
	}
	return uZ.p;

}

posit_2_t convertQuadToPX2(__float128 f32, int x){

	union ui32_pX2 uZ;
	bool sign, regS;
	uint_fast32_t reg, frac=0;
	int_fast32_t exp=0;
	bool bitNPlusOne=0, bitsMore=0;

	(f32>=0) ? (sign=0) : (sign=1);

	if (f32 == 0 ){
		uZ.ui = 0;
		return uZ.p;
	}
	else if(f32 == INFINITY || f32 == -INFINITY || f32 == NAN){
		uZ.ui = 0x80000000;
		return uZ.p;
	}
	else if (f32 == 1) {
		uZ.ui = 0x40000000;
		return uZ.p;
	}
	else if (f32 == -1){
		uZ.ui = 0xC0000000;
		return uZ.p;
	}
	/*else if (f32 >= 1.329227995784916e+36){
		//maxpos
		uZ.ui = 0x7FFFFFFF;
		return uZ.p;
	}
	else if (f32 <= -1.329227995784916e+36){
		// -maxpos
		uZ.ui = 0x80000001;
		return uZ.p;
	}
	else if(f32 <= 7.52316384526264e-37 && !sign){
		//minpos
		uZ.ui = 0x1;
		return uZ.p;
	}
	else if(f32 >= -7.52316384526264e-37 && sign){
		//-minpos
		uZ.ui = 0xFFFFFFFF;
		return uZ.p;
	}*/
	else if (f32>1 || f32<-1){
		if (sign){
			//Make negative numbers positive for easier computation
			f32 = -f32;
		}

		regS = 1;
		reg = 1; //because k = m-1; so need to add back 1
		// minpos
		if (x==32 && f32 <= 7.52316384526264e-37){
			uZ.ui = 1;
		}
		else{
			//regime
			while (f32>=16){
				f32 *=0.0625;  // f32/=16;
				reg++;
			}
			while (f32>=2){
				f32*=0.5;
				exp++;
			}

			int fracLength = x-4-reg;
			if (fracLength<0){
				//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
				if(reg==x-3){
					bitNPlusOne = exp&0x1;
					//exp>>=1; //taken care of by the pack algo
					exp&=0x2;
				}
				else{//reg=30
					bitNPlusOne=exp>>1;
					bitsMore=exp&0x1;
					exp=0;
				}
				if (f32!=1){//because of hidden bit
					bitsMore =1;
					frac=0;
				}
			}
			else
				frac = convertQuadFractionP32 (f32, fracLength, &bitNPlusOne, &bitsMore);


			if (reg>(x-2) ){
				uZ.ui=(regS) ? (0x7FFFFFFF & ((int32_t)0x80000000>>(x-1)) ): (0x1 << (32-x));
			}
			//rounding off fraction bits
			else{

				uint_fast32_t regime = 1;
				if (regS) regime = ( (1<<reg)-1 ) <<1;
				if (reg<=28)  exp<<= (28-reg);
				uZ.ui = ((uint32_t) (regime) << (30-reg)) + ((uint32_t) exp ) + ((uint32_t)(frac<<(32-x)));
				//minpos
				if (uZ.ui==0 && frac>0){
					uZ.ui = 0x1 << (32-x);
				}
				if (bitNPlusOne)
					uZ.ui +=  ( ((uZ.ui>>(32-x)) & 0x1) | bitsMore ) << (32-x);
			}
			if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;

		}
	}
	else if (f32 < 1 || f32 > -1 ){
		if (sign){
			//Make negative numbers positive for easier computation
			f32 = -f32;
		}
		regS = 0;
		reg = 0;

		//regime
		while (f32<1){
			f32 *= 16;
			reg++;
		}

		while (f32>=2){
			f32*=0.5;
			exp++;
		}


		int fracLength = x-4-reg;
		if (fracLength<0){
			//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
			if(reg==x-3){
				bitNPlusOne = exp&0x1;
				//exp>>=1; //taken care of by the pack algo
				exp&=0x2;
			}
			else{//reg=30
				bitNPlusOne=exp>>1;
				bitsMore=exp&0x1;
				exp=0;
			}
			if (f32!=1){//because of hidden bit
				bitsMore =1;
				frac=0;
			}
		}
		else
			frac = convertQuadFractionP32 (f32, fracLength, &bitNPlusOne, &bitsMore);

		if (reg>(x-2) ){
			uZ.ui=(regS) ? (0x7FFFFFFF & ((int32_t)0x80000000>>(x-1)) ): (0x1 << (32-x));
		}
		//rounding off fraction bits
		else{

			uint_fast32_t regime = 1;
			if (regS) regime = ( (1<<reg)-1 ) <<1;

			if (reg<=28)  exp<<= (28-reg);

			uZ.ui = ((uint32_t) (regime) << (30-reg)) + ((uint32_t) exp ) + ((uint32_t)(frac<<(32-x)));
			//minpos
			if (uZ.ui==0 && frac>0){
				uZ.ui = 0x1 << (32-x);
			}

			if (bitNPlusOne){
				uZ.ui += ( ((uZ.ui>>(32-x)) & 0x1) | bitsMore ) << (32-x);
			}


		}
		if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;

	}
	else {
		//NaR - for NaN, INF and all other combinations
		uZ.ui = 0x80000000;
	}
	return uZ.p;
}

#endif

void checkExtraP64TwoBits(__float128 f64, __float128 temp, bool * bitsNPlusOne, bool * bitsMore ){
	temp /= 2;
	if (temp<=f64){
		*bitsNPlusOne = 1;
		f64-=temp;
	}
	if (f64>0)
		*bitsMore = 1;
}
uint_fast64_t convertFractionP64(__float128 f64, uint_fast32_t fracLength, bool * bitsNPlusOne, bool * bitsMore ){

	uint_fast64_t frac=0;

	if(f64==0) return 0;
	else if(f64==INFINITY) return 0x8000000000000000;

	f64 -= 1; //remove hidden bit
	if (fracLength==0)
		checkExtraP64TwoBits(f64, 1.0, bitsNPlusOne, bitsMore);
	else{
		__float128 temp = 1;
		while (true){
			temp /= 2;
			if (temp<=f64){
				f64-=temp;
				fracLength--;
				frac = (frac<<1) + 1; //shift in one

				if (f64==0){
					frac <<= (uint_fast32_t)fracLength;
					break;
				}

				if (fracLength == 0){
					checkExtraP64TwoBits(f64, temp, bitsNPlusOne, bitsMore);
					break;
				}
			}
			else{

				frac <<= 1; //shift in a zero
				fracLength--;
				if (fracLength == 0){
					checkExtraP64TwoBits(f64, temp, bitsNPlusOne, bitsMore);
					break;
				}
			}
		}
	}

	return frac;
}

posit64_t convertDoubleToP64(double f64){

	union ui64_p64 uZ;
	bool sign, regS;
	uint_fast64_t reg, frac=0;
	int_fast64_t exp=0;
	bool bitNPlusOne=0, bitsMore=0;

	(f64>=0) ? (sign=0) : (sign=1);

	if (f64 == 0 ){
		uZ.ui = 0;
		return uZ.p;
	}
	else if(f64 == INFINITY || f64 == -INFINITY || f64 == NAN){
		uZ.ui = 0x8000000000000000;
		return uZ.p;
	}
	else if (f64 == 1) {
		uZ.ui = 0x4000000000000000;
		return uZ.p;
	}
	else if (f64 == -1){
		uZ.ui = 0xC000000000000000;
		return uZ.p;
	}
	else if (f64 >= 4.523128485832664e+74){ //2^248(63 bits are all 1, k = 62, (2^4)^62 = 2^248)
		//maxpos
		uZ.ui = 0x7FFFFFFFFFFFFFFF;
		return uZ.p;
	}
	else if (f64 <= -4.523128485832664e+74){ //-2^248
		// -maxpos
		uZ.ui = 0x8000000000000001;
		return uZ.p;
	}
	else if(f64 <= 2.210859150104178e-75 && !sign){ //2^-248(62 bits are all 0 and 1 bit 1, k = -62, (2^4)^-62 = 2^-248)
		//minpos
		uZ.ui = 0x1;
		return uZ.p;
	}
	else if(f64 >= -2.210859150104178e-75 && sign){ //-2^-248
		//-minpos
		uZ.ui = 0xFFFFFFFFFFFFFFFF;
		return uZ.p;
	}
	else if (f64>1 || f64<-1){
		if (sign){
			//Make negative numbers positive for easier computation
			f64 = -f64;
		}

		regS = 1;
		reg = 1; //because k = m-1; so need to add back 1
		// minpos
		if (f64 <= 2.210859150104178e-75){
			uZ.ui = 1;
		}
		else{
			//regime
			while (f64>=16){
				f64 *=0.0625;  // f64/=16, useed = 2^4;
				reg++;
			}
			while (f64>=2){
				f64 *= 0.5;
				exp++;
			}

			int16_t fracLength = 60-reg;

			if (fracLength<0){
				//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
				if(reg==61){
					bitNPlusOne = exp&0x1;
					exp>>=1; //taken care of by the pack algo
				}
				else{//reg=62
					bitNPlusOne=exp>>1;
					bitsMore=exp&0x1;
					exp=0;
				}
				if (f64!=1){//because of hidden bit
					bitsMore =1;
					frac=0;
				}
			}
			else
				frac = convertFractionP64 (f64, fracLength, &bitNPlusOne, &bitsMore);


			if (reg>62 ){
				(regS) ? (uZ.ui= 0x7FFFFFFFFFFFFFFF): (uZ.ui=0x1);
			}
			//rounding off fraction bits
			else{

				uint_fast64_t regime = 1;
				if (regS) regime = ( (1<<reg)-1 ) <<1;
				if (reg<=60)  exp<<= (60-reg);
				uZ.ui = ((uint64_t) (regime) << (62-reg)) + ((uint64_t) exp ) + ((uint64_t)(frac));
				uZ.ui += (bitNPlusOne & (uZ.ui&1)) | ( bitNPlusOne & bitsMore);
			}
			if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFF;

		}
	}
	else if (f64 < 1 || f64 > -1 ){
		if (sign){
			//Make negative numbers positive for easier computation
			f64 = -f64;
		}
		regS = 0;
		reg = 0;

		//regime
		while (f64<1){
			f64 *= 16;
			reg++;
		}

		while (f64>=2){
			f64*=0.5;
			exp++;
		}


		//only possible combination for reg=15 to reach here is 7FFF (maxpos) and FFFF (-minpos)
		//but since it should be caught on top, so no need to handle
		int_fast16_t fracLength =60-reg;
		if (fracLength<0){
			//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
			if(reg==61){
				bitNPlusOne = exp&0x1;
				exp>>=1; //taken care of by the pack algo
			}
			else{//reg=62
				bitNPlusOne=exp>>1;
				bitsMore=exp&0x1;
				exp=0;
			}
			if (f64!=1){//because of hidden bit
				bitsMore =1;
				frac=0;
			}
		}
		else
			frac = convertFractionP64 (f64, fracLength, &bitNPlusOne, &bitsMore);


		if (reg>62 ){
			(regS) ? (uZ.ui= 0x7FFFFFFFFFFFFFFF): (uZ.ui=0x1);
		}
		//rounding off fraction bits
		else{

			uint_fast64_t regime = 1;
			if (regS) regime = ( (1<<reg)-1 ) <<1;
			if (reg<=60)  exp<<= (60-reg);
			uZ.ui = ((uint64_t) (regime) << (62-reg)) + ((uint64_t) exp ) + ((uint64_t)(frac));
			uZ.ui += (bitNPlusOne & (uZ.ui&1)) | ( bitNPlusOne & bitsMore);
		}
		if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFF;

	}
	else {
		//NaR - for NaN, INF and all other combinations
		uZ.ui = 0x8000000000000000;
	}
	return uZ.p;
}


posit64_t convertFloatToP64(float a){
	return convertDoubleToP64((double) a );
}

/*
posit_2_t convertDoubleToPX2(double f32, int x){

	union ui32_pX2 uZ;
	bool sign, regS;
	uint_fast32_t reg, frac=0;
	int_fast32_t exp=0;
	bool bitNPlusOne=0, bitsMore=0;

	(f32>=0) ? (sign=0) : (sign=1);

	if (f32 == 0 ){
		uZ.ui = 0;
		return uZ.p;
	}
	else if(f32 == INFINITY || f32 == -INFINITY || f32 == NAN){
		uZ.ui = 0x80000000;
		return uZ.p;
	}
	else if (f32 == 1) {
		uZ.ui = 0x40000000;
		return uZ.p;
	}
	else if (f32 == -1){
		uZ.ui = 0xC0000000;
		return uZ.p;
	}
	/*else if (f32 >= 1.329227995784916e+36){
		//maxpos
		uZ.ui = 0x7FFFFFFF;
		return uZ.p;
	}
	else if (f32 <= -1.329227995784916e+36){
		// -maxpos
		uZ.ui = 0x80000001;
		return uZ.p;
	}
	else if(f32 <= 7.52316384526264e-37 && !sign){
		//minpos
		uZ.ui = 0x1;
		return uZ.p;
	}
	else if(f32 >= -7.52316384526264e-37 && sign){
		//-minpos
		uZ.ui = 0xFFFFFFFF;
		return uZ.p;
	}(there is a * / here)
	else if (f32>1 || f32<-1){
		if (sign){
			//Make negative numbers positive for easier computation
			f32 = -f32;
		}

		regS = 1;
		reg = 1; //because k = m-1; so need to add back 1
		// minpos
		if (x==32 && f32 <= 7.52316384526264e-37){
			uZ.ui = 1;
		}
		else{
			//regime
			while (f32>=16){
				f32 *=0.0625;  // f32/=16;
				reg++;
			}
			while (f32>=2){
				f32*=0.5;
				exp++;
			}

			int fracLength = x-4-reg;
			if (fracLength<0){
				//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
				if(reg==x-3){
					bitNPlusOne = exp&0x1;
					//exp>>=1; //taken care of by the pack algo
					exp&=0x2;
				}
				else{//reg=30
					bitNPlusOne=exp>>1;
					bitsMore=exp&0x1;
					exp=0;
				}
				if (f32!=1){//because of hidden bit
					bitsMore =1;
					frac=0;
				}
			}
			else
				frac = convertFractionP32 (f32, fracLength, &bitNPlusOne, &bitsMore);


			if (reg>(x-2) ){
				uZ.ui=(regS) ? (0x7FFFFFFF & ((int32_t)0x80000000>>(x-1)) ): (0x1 << (32-x));
			}
			//rounding off fraction bits
			else{
				uint_fast32_t regime = 1;
				if (regS) regime = ( (1<<reg)-1 ) <<1;

				if (x==32 && reg==29) exp>>=1;
				else if (reg<=28)  exp<<= (28-reg);

				uZ.ui = ((uint32_t) (regime) << (30-reg)) + ((uint32_t) exp ) + ((uint32_t)(frac<<(32-x)));
				//minpos
				if (uZ.ui==0 && frac>0){
					uZ.ui = 0x1 << (32-x);
				}
				if (bitNPlusOne)
					uZ.ui +=  ( ((uZ.ui>>(32-x)) & 0x1) | bitsMore ) << (32-x);
			}
			if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;

		}
	}
	else if (f32 < 1 || f32 > -1 ){
		if (sign){
			//Make negative numbers positive for easier computation
			f32 = -f32;
		}
		regS = 0;
		reg = 0;

		//regime
		while (f32<1){
			f32 *= 16;
			reg++;
		}

		while (f32>=2){
			f32*=0.5;
			exp++;
		}

		int fracLength = x-4-reg;
		if (fracLength<0){
			//in both cases, reg=29 and 30, e is n+1 bit and frac are sticky bits
			if(reg==x-3){
				bitNPlusOne = exp&0x1;
				//exp>>=1; //taken care of by the pack algo
				exp&=0x2;
			}
			else{//reg=30
				bitNPlusOne=exp>>1;
				bitsMore=exp&0x1;
				exp=0;
			}
			if (f32!=1){//because of hidden bit
				bitsMore =1;
				frac=0;
			}
		}
		else
			frac = convertFractionP32 (f32, fracLength, &bitNPlusOne, &bitsMore);

		if (reg>(x-2) ){
			uZ.ui=(regS) ? (0x7FFFFFFF & ((int32_t)0x80000000>>(x-1)) ): (0x1 << (32-x));
		}
		//rounding off fraction bits
		else{

			uint_fast32_t regime = 1;
			if (regS) regime = ( (1<<reg)-1 ) <<1;

			if (x==32 && reg==29) exp>>=1;
			else if (reg<=28)  exp<<= (28-reg);

			uZ.ui = ((uint32_t) (regime) << (30-reg)) + ((uint32_t) exp ) + ((uint32_t)(frac<<(32-x)));
			//minpos
			if (uZ.ui==0 && frac>0){
				uZ.ui = 0x1 << (32-x);
			}

			if (bitNPlusOne){
				uZ.ui += ( ((uZ.ui>>(32-x)) & 0x1) | bitsMore ) << (32-x);
			}
		}
		if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;

	}
	else {
		//NaR - for NaN, INF and all other combinations
		uZ.ui = 0x80000000;
	}
	return uZ.p;
}
*/




