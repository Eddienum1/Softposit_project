#include "platform.h"
#include "internals.h"
/*
posit64_t softposit_mulAddP64(
    uint_fast64_t uiA, uint_fast64_t uiB, uint_fast64_t uiC, uint_fast64_t op) {

    union ui64_p64 uZ;
    uint_fast64_t regA, regZ, fracA, fracZ, regime, tmp;
    bool signA, signB, signC, signZ, regSA, regSB, regSC, regSZ, bitNPlusOne = 0, bitsMore = 0, rcarry;
    int_fast32_t expA, expC, expZ;
    int_fast16_t kA = 0, kC = 0, kZ = 0, shiftRight;
    __uint128_t frac128C = 0, frac128Z = 0;

    if (uiA == 0x8000000000000000 || uiB == 0x8000000000000000 || uiC == 0x8000000000000000) {
        uZ.ui = 0x8000000000000000;
        return uZ.p;
    }
    else if (uiA == 0 || uiB == 0) {
        if (op == softposit_mulAdd_subC)
			uZ.ui = -uiC;
		else
			uZ.ui = uiC;
		return uZ.p;
    }

    signA = signP64UI(uiA);
    signB = signP64UI(uiB);
    signC = signP64UI(uiC);
    signZ = signA ^ signB;

    if (signA) uiA = (-uiA & 0xFFFFFFFFFFFFFFFF);
    if (signB) uiB = (-uiB & 0xFFFFFFFFFFFFFFFF);
    if (signC) uiC = (-uiC & 0xFFFFFFFFFFFFFFFF);

    regSA = signregP64UI(uiA);
    regSB = signregP64UI(uiB);
    regSC = signregP64UI(uiC);

    tmp = (uiA << 2) & 0xFFFFFFFFFFFFFFFF;
    if (regSA) {
        while (tmp >> 63) {
            kA++;
            tmp = (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
        }
    } else {
        kA = -1;
        while (!(tmp >> 63)) {
            kA--;
            tmp = (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
        }
        tmp &= 0x7FFFFFFFFFFFFFFF;
    }
    expA = tmp >> 61;
    fracA = ((tmp << 3) | 0x8000000000000000) & 0xFFFFFFFFFFFFFFFF;

    tmp = (uiB << 2) & 0xFFFFFFFFFFFFFFFF;
    if (regSB) {
        while (tmp >> 63) {
            kA++;
            tmp = (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
        }
    } else {
        kA--;
        while (!(tmp >> 63)) {
            kA--;
            tmp = (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
        }
        tmp &= 0x7FFFFFFFFFFFFFFF;
    }
    expA += tmp >> 61;

    frac128Z = (__uint128_t)fracA * (((tmp << 2) | 0x8000000000000000) & 0xFFFFFFFFFFFFFFFF);

    if (expA > 3) {
        kA++;
        expA &= 0x3;
    }

    rcarry = frac128Z >> 127;
    if (rcarry) {
        expA++;
        if (expA > 3) {
            kA++;
            expA &= 0x3;
        }
        frac128Z >>= 1;
    }

    if (uiC != 0) {
        tmp = (uiC<<2) & 0xFFFFFFFFFFFFFFFF;
        if (regSC) {
            while (tmp >> 63) {
                kC++;
                tmp = (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
            }
        } else {
            kC = -1;
            while (!(tmp >> 63)) {
                kC--;
                tmp = (tmp<<1) & 0xFFFFFFFFFFFFFFFF;
            }
            tmp &= 0x7FFFFFFFFFFFFFFF;
        }
        expC = tmp >> 61;
        frac128C = ((__uint128_t)((tmp << 1) | 0x4000000000000000)) << 64;

        shiftRight = ((kA - kC) << 2) + (expA - expC);
        if (shiftRight < 0) {
            if (shiftRight <= -127) {
                bitsMore = 1;
                frac128Z = 0;
            } else if ((frac128Z << (128 + shiftRight)) != 0) bitsMore = 1;

            if (signZ == signC)
                frac128Z = frac128C + (frac128Z >> -shiftRight);
            else {
                frac128Z = frac128C - (frac128Z >> -shiftRight);
                signZ = signC;
                if (bitsMore) frac128Z -= 1;
            }
            kZ = kC;
            expZ = expC;
        } else if (shiftRight > 0) {
            if (shiftRight >= 127) {
                bitsMore = 1;
                frac128C = 0;
            } else if ((frac128C << (128 - shiftRight)) != 0) bitsMore = 1;

            if (signZ == signC)
                frac128Z += (frac128C >> shiftRight);
            else {
                frac128Z -= (frac128C >> shiftRight);
                if (bitsMore) frac128Z -= 1;
            }
            kZ = kA;
            expZ = expA;
        } else {
            if (frac128C == frac128Z && signZ != signC) {
                uZ.ui = 0;
                return uZ.p;
            } else {
                if (signZ == signC)
                    frac128Z += frac128C;
                else {
                    if (frac128Z < frac128C) {
                        frac128Z = frac128C - frac128Z;
                        signZ = signC;
                    } else {
                        frac128Z -= frac128C;
                    }
                }
            }
            kZ = kA;
            expZ = expA;
        }

        rcarry = (frac128Z >> 127);
        if (rcarry) {
            expZ++;
            if (expZ > 3) {
                kZ++;
                expZ &= 0x3;
            }
            frac128Z = (frac128Z >> 1) & 0x7FFFFFFFFFFFFFFF;
        } else {
            if (frac128Z != 0) {
                while ((frac128Z >> 123) == 0) {
                    kZ--;
                    frac128Z <<= 4;
                }
                while ((frac128Z >> 126) == 0) {
                    expZ--;
                    frac128Z <<= 1;
                    if (expZ < 0) {
                        kZ--;
                        expZ = 3;
                    }
                }
            }
        }
    } else {
        kZ = kA;
        expZ = expA;
    }

    if (kZ < 0) {
        regZ = -kZ;
        regSZ = 0;
        regime = 0x4000000000000000 >> regZ;
    } else {
        regZ = kZ + 1;
        regSZ = 1;
        regime = 0x7FFFFFFFFFFFFFFF - (0x7FFFFFFFFFFFFFFF >> regZ);
    }

    if (regZ > 62) {
        (regSZ) ? (uZ.ui= 0x7FFFFFFFFFFFFFFF): (uZ.ui=0x1);
    } 
    else {

        if (regZ <= 60){
            frac128Z &= 0x3FFFFFFFFFFFFFFF;
            fracZ = (uint_fast64_t)(frac128Z >> (regZ + 66));
            bitNPlusOne |= (0x2000000000000000 & (frac128Z >>regZ ) ) ;
            expZ <<= (60 - regZ);
        }
        else {
            if (regZ == 62) {
                bitNPlusOne = expZ & 0x2;
                bitsMore = expZ & 0x1;
                expZ = 0;
            } else if (regZ == 61) {
                bitNPlusOne = expZ & 0x1;
                expZ >>= 1;
            }
            if (fracZ > 0) {
                fracZ = 0;
                bitsMore = 1;
            }
        }

        uZ.ui = packToP64UI(regime, expZ, fracZ);
        if (bitNPlusOne) {
            if ((frac128Z << (64 - regZ)) != 0) bitsMore = 1;
            uZ.ui += (uZ.ui & 1) | bitsMore;
        }
    }

    if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFF;
    return uZ.p;
}
*/
posit64_t softposit_mulAddP64(
    uint_fast64_t uiA, uint_fast64_t uiB, uint_fast64_t uiC, uint_fast64_t op ){

    union ui64_p64 uZ;
    uint_fast64_t regA, regZ, fracA, fracZ, regime, tmp;
    bool signA, signB, signC, signZ, regSA, regSB, regSC, regSZ;
    bool bitNPlusOne=0, bitsMore=0, rcarry;
    int_fast64_t expA, expC, expZ;
    int_fast32_t kA=0, kC=0, kZ=0, shiftRight;
    __uint128_t frac128C, frac128Z;

    //NaR
    if ( uiA==0x8000000000000000ULL || uiB==0x8000000000000000ULL || uiC==0x8000000000000000ULL ){
        uZ.ui = 0x8000000000000000ULL;
        return uZ.p;
    }
    else if (uiA==0 || uiB==0){
        if (op == softposit_mulAdd_subC)
            uZ.ui = -uiC;
        else
            uZ.ui = uiC;
        return uZ.p;
    }

    signA = signP64UI( uiA );
    signB = signP64UI( uiB );
    signC = signP64UI( uiC );
    signZ = signA ^ signB;

    if(signA) uiA = (-uiA & 0xFFFFFFFFFFFFFFFFULL);
    if(signB) uiB = (-uiB & 0xFFFFFFFFFFFFFFFFULL);
    if(signC) uiC = (-uiC & 0xFFFFFFFFFFFFFFFFULL);

    regSA = signregP64UI(uiA);
    regSB = signregP64UI(uiB);
    regSC = signregP64UI(uiC);

    // Decode uiA
    tmp = (uiA<<2);
    if (regSA){
        while (tmp>>63){
            kA++;
            tmp <<= 1;
        }
    } else {
        kA=-1;
        while (!(tmp>>63)){
            kA--;
            tmp <<= 1;
        }
        tmp &= 0x7FFFFFFFFFFFFFFFULL;
    }
    expA = tmp>>61;  // es=2
    fracA = ((tmp<<2) | 0x8000000000000000ULL);

    // Decode uiB, accumulate exp+k
    tmp = (uiB<<2);
    if (regSB){
        while (tmp>>63){
            kA++;
            tmp <<= 1;
        }
    } else {
        kA--;
        while (!(tmp>>63)){
            kA--;
            tmp <<= 1;
        }
        tmp &= 0x7FFFFFFFFFFFFFFFULL;
    }
    expA += tmp>>61;
    frac128Z = (__uint128_t) fracA * (((tmp<<2) | 0x8000000000000000ULL));

    if (expA>3){
        kA++;
        expA&=0x3;
    }

    rcarry = (frac128Z>>127); // first bit
    if (rcarry){
        expA++;
        if (expA>3){
            kA++;
            expA&=0x3;
        }
        frac128Z >>= 1;
    }

    if (uiC!=0){
        tmp = (uiC<<2);
        if (regSC){
            while (tmp>>63){
                kC++;
                tmp <<= 1;
            }
        } else {
            kC=-1;
            while (!(tmp>>63)){
                kC--;
                tmp <<= 1;
            }
            tmp &= 0x7FFFFFFFFFFFFFFFULL;
        }
        expC = tmp>>61;
        frac128C = (((__uint128_t)(tmp<<1) | 0x4000000000000000ULL) & 0x7FFFFFFFFFFFFFFFULL)<<64;
        shiftRight = ((kA-kC)<<2) + (expA-expC);

        if (shiftRight<0){ // |uiC| > |Prod|
            if (shiftRight<=-127){
                bitsMore = 1;
                frac128Z = 0;
            }
            else if ((frac128Z<<(128+shiftRight))!=0) bitsMore = 1;
            if (signZ==signC)
                frac128Z = frac128C + (frac128Z>>-shiftRight);
            else {
                frac128Z = frac128C - (frac128Z>>-shiftRight);
                signZ=signC;
                if (bitsMore) frac128Z-=1;
            }
            kZ = kC;
            expZ = expC;

        } else if (shiftRight>0){
            if(shiftRight>=127) {
                bitsMore = 1;
                frac128C = 0;
            }
            else if ((frac128C<<(128-shiftRight))!=0) bitsMore = 1;
            if (signZ==signC)
                frac128Z = frac128Z + (frac128C>>shiftRight);
            else{
                frac128Z = frac128Z - (frac128C>>shiftRight);
                if (bitsMore) frac128Z-=1;
            }
            kZ = kA;
            expZ = expA;

        } else {
            if(frac128C==frac128Z && signZ!=signC){
                uZ.ui = 0;
                return uZ.p;
            } else {
                if (signZ==signC)
                    frac128Z += frac128C;
                else {
                    if (frac128Z<frac128C){
                        frac128Z = frac128C - frac128Z;
                        signZ = signC;
                    }
                    else{
                        frac128Z -= frac128C;
                    }
                }
            }
            kZ = kA;
            expZ = expA;
        }

        rcarry = (frac128Z>>127);
        if(rcarry){
            expZ++;
            if (expZ>3){
                kZ++;
                expZ&=0x3;
            }
            frac128Z = (frac128Z>>1) & (((__uint128_t)1 << 127) - 1);
        }
        else {
            if (frac128Z!=0){
                while((frac128Z>>123)==0){
                    kZ--;
                    frac128Z<<=4;
                }
                while((frac128Z>>126)==0){
                    expZ--;
                    frac128Z<<=1;
                    if (expZ<0){
                        kZ--;
                        expZ=3;
                    }
                }
            }
        }
    }
    else{
        kZ = kA;
        expZ=expA;
    }

    if(kZ<0){
        regZ = -kZ;
        regSZ = 0;
        regime = 0x4000000000000000ULL>>regZ;
    }
    else{
        regZ = kZ+1;
        regSZ=1;
        regime = 0x7FFFFFFFFFFFFFFFULL - (0x7FFFFFFFFFFFFFFFULL>>regZ);
    }

    if(regZ>62){
        (regSZ) ? (uZ.ui= 0x7FFFFFFFFFFFFFFFULL): (uZ.ui=0x1);
    }
    else{
        if (regZ<=60){
            frac128Z &= (((__uint128_t)1 << 126) - 1);
            fracZ = (uint_fast64_t)(frac128Z >> (regZ + 66));
            bitNPlusOne |= (0x2000000000000000ULL & (frac128Z >> regZ));
            expZ <<= (60-regZ);
        }
        else {
            if (regZ==62){
                bitNPlusOne = expZ&0x2;
                bitsMore = (expZ&0x1);
                expZ = 0;
            }
            else if (regZ==61){
                bitNPlusOne = expZ&0x1;
                expZ>>=1;
            }
            if (fracZ>0){
                fracZ=0;
                bitsMore =1;
            }
        }
        uZ.ui = packToP64UI(regime, expZ, fracZ);

        if (bitNPlusOne){
            if ( (frac128Z<<(64-regZ)) ) bitsMore =1;
            uZ.ui += (uZ.ui&1) | bitsMore;
        }
    }
    if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFFULL;
    return uZ.p;
}
