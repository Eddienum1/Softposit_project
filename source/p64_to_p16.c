#include "platform.h"
#include "internals.h"

posit16_t p64_to_p16(posit64_t pA) {

    union ui64_p64 uA;
    union ui16_p16 uZ;
    uint_fast64_t uiA, tmp = 0, exp_frac64A = 0;
    uint_fast16_t regime = 0;
    bool sign, regSA, bitsMore = 0, bitNPlusOne = 0;
    int_fast16_t kA = 0, regA;

    uA.p = pA;
    uiA = uA.ui;

    if (uiA == 0x8000000000000000ULL || uiA == 0) {
        uZ.ui = (uint16_t)((uiA >> 48) & 0xFFFF);
        return uZ.p;
    }

    sign = signP64UI(uiA);
    if (sign) uiA = -uiA & 0xFFFFFFFFFFFFFFFFULL;

    if (uiA > 0x7FC0000000000000ULL) {
        uZ.ui = 0x7FFF;
    } else if (uiA < 0x0040000000000000ULL) {
        uZ.ui = 0x0001;
    } else {
        regSA = signregP64UI(uiA);

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

        exp_frac64A = tmp;

        if (kA < 0) {
            regA = ((-kA) << 1) - (exp_frac64A >> 62);
            if (regA == 0) regA = 1;
            regSA = 0;
            regime = (regA > 14) ? 0x0001 : (0x4000 >> regA);
        } else {
            regA = ((kA << 1) + (exp_frac64A >> 62) + 1);
            regSA = 1;
            regime = 0x7FFF - (0x7FFF >> regA);
        }

        exp_frac64A = (exp_frac64A << 2); 

        if (regA > 13) {
            uZ.ui = regime;
        } else {
            uZ.ui = regime | (exp_frac64A >> (regA + 49));
        }

        if ((exp_frac64A >> (regA + 48)) & 0x1) {
            bitsMore = exp_frac64A & ((1ULL << (regA + 48)) - 1);
            uZ.ui += ((uZ.ui & 1) | bitsMore);
        }
    }

    if (sign) uZ.ui = -uZ.ui & 0xFFFF;
    return uZ.p;
}
