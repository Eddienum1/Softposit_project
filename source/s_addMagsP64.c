#include "platform.h"
#include "internals.h"

posit64_t softposit_addMagsP64(uint_fast64_t uiA, uint_fast64_t uiB) {

    uint_fast16_t regA;
    uint_fast64_t frac64A = 0, frac64B = 0, fracA = 0, regime, tmp;
    bool sign, regSA, regSB, rcarry = 0, bitNPlusOne = 0, bitsMore = 0;
    int_fast8_t kA = 0;
    int_fast32_t expA;
    int_fast16_t shiftRight;
    union ui64_p64 uZ;

    sign = signP64UI(uiA);
    if (sign) {
        uiA = -uiA & 0xFFFFFFFFFFFFFFFFULL;
        uiB = -uiB & 0xFFFFFFFFFFFFFFFFULL;
    }

    if ((int_fast64_t)uiA < (int_fast64_t)uiB) {
        uint_fast64_t temp = uiA;
        uiA = uiB;
        uiB = temp;
    }

    regSA = signregP64UI(uiA);
    regSB = signregP64UI(uiB);

    tmp = (uiA << 2);
    if (regSA) {
        while (tmp >> 63) {
            kA++;
            tmp = (tmp << 1);
        }
    } else {
        kA = -1;
        while (!(tmp >> 63)) {
            kA--;
            tmp = (tmp << 1);
        }
        tmp &= 0x7FFFFFFFFFFFFFFFULL;
    }

    expA = tmp >> 61;  // 3 bits exponent
    frac64A = ((0x8000000000000000ULL | (tmp << 2)) & 0xFFFFFFFFFFFFFFFFULL) << 32;
    shiftRight = kA;

    tmp = (uiB << 2);
    if (regSB) {
        while (tmp >> 63) {
            shiftRight--;
            tmp = (tmp << 1);
        }
    } else {
        shiftRight++;
        while (!(tmp >> 63)) {
            shiftRight++;
            tmp = (tmp << 1);
        }
        tmp &= 0x7FFFFFFFFFFFFFFFULL;
    }

    frac64B = ((0x8000000000000000ULL | (tmp << 2)) & 0xFFFFFFFFFFFFFFFFULL) << 32;
    shiftRight = (shiftRight << 3) + expA - (tmp >> 61);  // 3-bit exp for P64

    if (shiftRight > 63) {
        frac64B = 0;
    } else {
        frac64B >>= shiftRight;
    }

    frac64A += frac64B;
    rcarry = frac64A & 0x80000000000000000000000000000000ULL;

    if (rcarry) {
        expA++;
        if (expA > 7) {
            kA++;
            expA &= 0x7;
        }
        frac64A >>= 1;
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

    if (regA > 60) {
        uZ.ui = regSA ? 0x7FFFFFFFFFFFFFFFULL : 0x1ULL;
    } else {
        frac64A = (frac64A & 0x7FFFFFFFFFFFFFFFULL) >> (regA + 3);  // remove hidden and exp bits
        fracA = frac64A >> 32;

        if (regA <= 58) {
            bitNPlusOne = (frac64A & 0x80000000ULL);
            expA <<= (58 - regA);
        } else {
            if (regA == 60) {
                bitNPlusOne = expA & 0x4;
                bitsMore = (expA & 0x3);
                expA = 0;
            } else if (regA == 59) {
                bitNPlusOne = expA & 0x2;
                bitsMore = (expA & 0x1);
                expA >>= 1;
            } else if (regA == 58) {
                bitNPlusOne = expA & 0x1;
                expA >>= 2;
            }
            if (fracA > 0) {
                fracA = 0;
                bitsMore = 1;
            }
        }

        uZ.ui = packToP64UI(regime, expA, fracA);
        if (bitNPlusOne) {
            if (frac64A & 0x7FFFFFFFULL) bitsMore = 1;
            uZ.ui += (uZ.ui & 1) | bitsMore;
        }
    }

    if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFFULL;
    return uZ.p;
}
