#include "platform.h"
#include "internals.h"

uint_fast64_t pX2_to_ui64(posit_2_t pA) {
    posit64_t p64 = {.v = pA.v};
    return p64_to_ui64(p64);
}

uint_fast64_t p64_to_ui64(posit64_t pA) {

    union ui64_p64 uA;
    uint_fast64_t uiA, mask, iZ, tmp;
    uint_fast32_t scale = 0;
    bool bitLast, bitNPlusOne, bitsMore;

    uA.p = pA;
    uiA = uA.ui;


    if (uiA >= 0x8000000000000000ULL)
        return 0;


    if (uiA <= 0x3800000000000000ULL) {
        return 0;
    }
    else if (uiA < 0x4400000000000000ULL) {
        return 1;
    }
    else if (uiA <= 0x4A00000000000000ULL) {
        return 2;
    }
    else if (uiA > 0x7FFFFFFFFFFFFFFFULL) {
        return 0xFFFFFFFFFFFFFFFFULL;
    }

    uiA -= 0x4000000000000000ULL;


    while (0x2000000000000000ULL & uiA) {
        scale += 4;
        uiA = (uiA - 0x2000000000000000ULL) << 1;
    }


    uiA <<= 1;


    if (0x2000000000000000ULL & uiA) scale += 2;
    if (0x1000000000000000ULL & uiA) scale += 1;


    iZ = ((uiA | 0x1000000000000000ULL) & 0x1FFFFFFFFFFFFFFFULL) << 34;

    if (scale < 62) {
        mask = 0x4000000000000000ULL >> scale;

        bitLast = (iZ & mask);
        mask >>= 1;
        tmp = (iZ & mask);
        bitNPlusOne = tmp;
        iZ ^= tmp;
        tmp = iZ & (mask - 1);
        iZ ^= tmp;

        if (bitNPlusOne) {
            if (bitLast | tmp) iZ += (mask << 1);
        }

        iZ = iZ >> (62 - scale);
    }
    else if (scale > 62) {
        iZ = iZ << (scale - 62);
    }

    return iZ;
}
