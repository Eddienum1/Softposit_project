#include "platform.h"
#include "internals.h"

uint_fast32_t p64_to_ui32(posit64_t pA) {

    union ui64_p64 uA;
    uint_fast64_t uiA, iZ64, tmp, mask;
    uint_fast32_t iZ, scale = 0;
    bool bitLast, bitNPlusOne;

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
    else if (uiA > 0x7FBFFFFFFFFFFFFFULL) {
        return 0xFFFFFFFF;
    }
    else {
        uiA -= 0x4000000000000000ULL;
        while (0x2000000000000000ULL & uiA) {
            scale += 4;
            uiA = (uiA - 0x2000000000000000ULL) << 1;
        }
        uiA <<= 1;

        if (0x2000000000000000ULL & uiA) scale += 2;
        if (0x1000000000000000ULL & uiA) scale += 1;

        iZ64 = ((uiA | 0x1000000000000000ULL) & 0x1FFFFFFFFFFFFFFFULL) << 34;

        mask = 0x4000000000000000ULL >> scale;
        bitLast = (iZ64 & mask);
        mask >>= 1;
        tmp = (iZ64 & mask);
        bitNPlusOne = tmp;
        iZ64 ^= tmp;
        tmp = iZ64 & (mask - 1);
        iZ64 ^= tmp;

        if (bitNPlusOne) {
            if (bitLast | tmp) iZ64 += (mask << 1);
        }

        iZ = (uint64_t)iZ64 >> (62 - scale);
    }

    return iZ;
}
