#include "platform.h"
#include "internals.h"

uint_fast32_t p64_to_ui32(posit64_t pA) {

    union ui64_p64 uA;
    uint_fast64_t uiA;
    uint_fast32_t iZ, scale = 0;
    bool bitLast, bitNPlusOne;
    __uint128_t iZ128, mask, tmp;

    uA.p = pA;
    uiA = uA.ui;

    if (uiA >= 0x8000000000000000)
        return 0;

    if (uiA <= 0x3800000000000000) {
        return 0;
    }
    else if (uiA < 0x4400000000000000) {
        return 1;
    }
    else if (uiA <= 0x4A00000000000000) {
        return 2;
    }
    else if (uiA > 0x7FBFFFFFFFFFFFFF) {
        return 0xFFFFFFFF;
    }
    else {
        uiA -= 0x4000000000000000;
        while (0x2000000000000000 & uiA) {
            scale += 4;
            uiA = (uiA - 0x2000000000000000) << 1;
        }
        uiA <<= 1;

        if (0x2000000000000000 & uiA) scale += 2;
        if (0x1000000000000000 & uiA) scale += 1;

        iZ128 = ((__uint128_t)(uiA | 0x1000000000000000) & 0x1FFFFFFFFFFFFFFF) << 66;

    mask = ((__uint128_t)0x4000000000000000 << 64) >> scale; 

    bitLast = (iZ128 & mask);
    mask >>= 1;
    tmp = (iZ128 & mask);
    bitNPlusOne = tmp;
    iZ128 ^= tmp;
    tmp = iZ128 & (mask - 1);
    iZ128 ^= tmp;

    if (bitNPlusOne) {
        if (bitLast | tmp) iZ128 += (mask << 1);
    }

    iZ = (uint64_t)(iZ128 >> (126 - scale));
    }

    return iZ;
}
