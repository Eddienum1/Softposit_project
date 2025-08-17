#include "platform.h"
#include "internals.h"

posit32_t p64_to_p32(posit64_t pA) {

    union ui64_p64 uA;
    union ui32_p32 uZ;
    uint_fast64_t uiA, tmp = 0, exp_frac64A;
    uint_fast32_t regime, exp_frac = 0;
    bool sign, regSA, bitsMore = 0, bitNPlusOne = 0;
    int_fast16_t kA = 0, regA;

    uA.p = pA;
    uiA = uA.ui;

    if (uiA == 0x8000000000000000ULL || uiA == 0) {
        uZ.ui = (uint32_t)(uiA >> 32);
        return uZ.p;
    }

    sign = signP64UI(uiA);
    if (sign) uiA = -uiA & 0xFFFFFFFFFFFFFFFFULL;

    if (uiA > 0x7FFE000000000000ULL) {
        uZ.ui = 0x7FFFFFFF;
    } else if (uiA < 0x0001000000000000ULL) {
        uZ.ui = 0x1; 
    } else {
        regSA = signregP64UI(uiA);

        tmp = (uiA << 2) & 0xFFFFFFFFFFFFFFFFULL;
        if (regSA) {
            while (tmp >> 63) {
                kA++;
                tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFFULL;
            }
        } else {
            kA = -1;
            while (!(tmp >> 63)) {
                kA--;
                tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFFULL;
            }
            tmp &= 0x7FFFFFFFFFFFFFFFULL;
        }

        exp_frac64A = tmp << 1;

        if (kA < 0) {
            regA = (-kA) << 1;
            if (exp_frac64A & 0x8000000000000000ULL) regA--;
            exp_frac64A = (exp_frac64A << 1) & 0xFFFFFFFFFFFFFFFFULL;
            regSA = 0;
            regime = 0x40000000 >> regA;
        } else {
            regA = (kA << 1) + 1;
            if (exp_frac64A & 0x8000000000000000ULL) regA++;
            exp_frac64A = (exp_frac64A << 1) & 0xFFFFFFFFFFFFFFFFULL;
            regSA = 1;
            regime = 0x7FFFFFFF - (0x7FFFFFFF >> regA);
        }

        if ((exp_frac64A >> (33 + regA)) & 0x1) bitNPlusOne = 1;
        if (regA < 30) exp_frac = (uint32_t)(exp_frac64A >> (34 + regA));

        uZ.ui = regime + exp_frac;

        if (bitNPlusOne) {
            if ((exp_frac64A << (31 - regA)) & 0xFFFFFFFFFFFFFFFFULL) bitsMore = 1;
            uZ.ui += (bitNPlusOne & (uZ.ui & 1)) | (bitNPlusOne & bitsMore);
        }
    }

    if (sign) uZ.ui = (-uZ.ui) & 0xFFFFFFFF;
    return uZ.p;
}
