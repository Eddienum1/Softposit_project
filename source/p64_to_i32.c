#include "platform.h"
#include "internals.h"

int_fast32_t p64_to_i32(posit64_t pA) {
    union ui64_p64 uA;
    __uint128_t iZ128, mask, tmp;
    int_fast32_t iZ;
    uint_fast64_t uiA;
    uint_fast32_t scale = 0;
    bool bitLast, bitNPlusOne, bitsMore, sign;

    uA.p = pA;
    uiA = uA.ui;

    if (uiA == 0x8000000000000000) return 0;  // NaR

    sign = uiA >> 63;
    if (sign) uiA = -uiA & 0xFFFFFFFFFFFFFFFF;

    // Rounding cases
    if (uiA <= 0x3800000000000000) return 0;            
    else if (uiA < 0x4400000000000000) return 1;        
    else if (uiA <= 0x4A00000000000000) return 2;       
    else if (uiA > 0x7FAFFFFFFFFFFFFF)              
        return (sign) ? (-2147483648) : (2147483647);

    // General case
    uiA -= 0x4000000000000000;
    while (0x2000000000000000 & uiA) {
        scale += 4;
        uiA = (uiA - 0x2000000000000000) << 1;
    }
    uiA <<= 1;  // Skip regime stop bit

    if (uiA & 0x2000000000000000) scale += 2;
    if (uiA & 0x1000000000000000) scale += 1;

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

    return sign ? -iZ : iZ;
}
