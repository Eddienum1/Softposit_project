#include "platform.h"
#include "internals.h"

posit64_t softposit_subMagsP64(uint_fast64_t uiA, uint_fast64_t uiB) {

    uint_fast16_t regA, regB;
    uint_fast64_t fracA=0, regime, tmp;
    __uint128_t frac128A=0, frac128B=0;
    bool sign, regSA, regSB, ecarry=0, bitNPlusOne=0, bitsMore=0;
    int_fast8_t kA=0;
    int_fast32_t expA=0;
    int_fast16_t shiftRight;
    union ui64_p64 uZ;

    sign = signP64UI(uiA);
    if (sign)
        uiA = -uiA & 0xFFFFFFFFFFFFFFFF;
    else
        uiB = -uiB & 0xFFFFFFFFFFFFFFFF;

    if (uiA == uiB) {
        uZ.ui = 0;
        return uZ.p;
    }

    if ((int_fast64_t)uiA < (int_fast64_t)uiB) {
        uiA ^= uiB;
        uiB ^= uiA;
        uiA ^= uiB;
        sign = !sign;
    }

    regSA = signregP64UI(uiA);
    regSB = signregP64UI(uiB);

    tmp = (uiA << 2) & 0xFFFFFFFFFFFFFFFF;
    if (regSA) {
        while (tmp >> 63) {
            kA++;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFF;
        }
    } else {
        kA = -1;
        while (!(tmp >> 63)) {
            kA--;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFF;
        }
        tmp &= 0x7FFFFFFFFFFFFFFF;
    }

    expA = tmp >> 61;
    frac128A = ((__uint128_t)(0x8000000000000000ULL | (tmp << 3)) << 64);
    shiftRight = kA;

    tmp = (uiB << 2) & 0xFFFFFFFFFFFFFFFF;
    if (regSB) {
        while (tmp >> 63) {
            shiftRight--;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFF;
        }
    } else {
        shiftRight++;
        while (!(tmp >> 63)) {
            shiftRight++;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFF;
        }
        tmp &= 0x7FFFFFFFFFFFFFFF;
    }

    frac128B = ((__uint128_t)(0x8000000000000000ULL | (tmp << 3)) << 64);

    shiftRight = (shiftRight << 2) + expA - (tmp >> 61);

    if (shiftRight > 127) {
        uZ.ui = uiA;
        if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFF;
        return uZ.p;
    } else {
        frac128B >>= shiftRight;
    }

    frac128A -= frac128B;

    while ((frac128A >> 115) == 0) {
        kA--;
        frac128A <<= 4;
    }

    ecarry = (frac128A & ((__uint128_t)1 << 127));
    while (!ecarry) {
        if (expA == 0) {
            kA--;
            expA = 3;
        } else {
            expA--;
        }
        frac128A <<= 1;
        ecarry = (frac128A & ((__uint128_t)1 << 127));
    }

    if (kA < 0) {
        regA = -kA;
        regSA = 0;
        regime = 0x4000000000000000ULL >> regA;
    } else {
        regA = kA + 1;
        regSA = 1;
        regime = 0x7FFFFFFFFFFFFFFFULL - (0x7FFFFFFFFFFFFFFFULL >> regA);
    }

    if (regA > 62) {
        uZ.ui = regSA ? 0x7FFFFFFFFFFFFFFFULL : 0x1;
    } else {
        frac128A = (frac128A & (((__uint128_t)1 << 127) - 1)) >> (regA + 2);
        fracA = (uint_fast64_t)(frac128A >> 64);

        if (regA <= 60) {
            bitNPlusOne = (frac128A >> 63) & 1;
            expA <<= (60 - regA);
        } else {
            if (regA == 62) {
                bitNPlusOne = expA & 0x2;
                bitsMore = expA & 0x1;
                expA = 0;
            } else if (regA == 61) {
                bitNPlusOne = expA & 0x1;
                expA >>= 1;
            }
            if (fracA > 0) {
                fracA = 0;
                bitsMore = 1;
            }
        }

        uZ.ui = packToP64UI(regime, expA, fracA);

        if (bitNPlusOne) {
            if ((frac128A & 0x7FFFFFFFFFFFFFFFULL) != 0) bitsMore = 1;
            uZ.ui += (uZ.ui & 1) | bitsMore;
        }
    }

    if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFF;
    return uZ.p;
}
