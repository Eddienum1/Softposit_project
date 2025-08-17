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

    if (uiA == 0x8000000000000000ULL) return 0;  // NaR

    sign = uiA >> 63;
    if (sign) uiA = -uiA & 0xFFFFFFFFFFFFFFFFULL;

    // Rounding cases
    if (uiA <= 0x7000000000000000ULL) return 0;            
    else if (uiA < 0x8800000000000000ULL) return 1;        
    else if (uiA <= 0x9200000000000000ULL) return 2;       
    else if (uiA > 0xFFBFFFFFFFFFFFFFULL)              
        return sign ? INT32_MIN : INT32_MAX;

    // General case
    uiA -= 0x8000000000000000ULL;
    while (0x4000000000000000ULL & uiA) {
        scale += 4;
        uiA = (uiA - 0x4000000000000000ULL) << 1;
    }
    uiA <<= 1;  // Skip regime stop bit

    if (uiA & 0x4000000000000000ULL) scale += 2;
    if (uiA & 0x2000000000000000ULL) scale += 1;

    iZ128 = ((__uint128_t)(uiA | 0x2000000000000000ULL) & 0x3FFFFFFFFFFFFFFFULL) << 66;

    mask = (__uint128_t)1 << (127 - scale);

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

    iZ = (uint64_t)(iZ128 >> (125 - scale));

    return sign ? -iZ : iZ;
}
