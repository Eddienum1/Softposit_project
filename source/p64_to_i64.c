#include "platform.h"
#include "internals.h"

int_fast64_t p64_to_i64(posit64_t pA) {
    union ui64_p64 uA;
    __uint128_t iZ128, mask, tmp;
    int_fast64_t iZ;
    uint_fast64_t uiA;
    uint_fast32_t scale = 0;
    bool bitLast, bitNPlusOne, bitsMore, sign;

    uA.p = pA;
    uiA = uA.ui;

    if (uiA == 0x8000000000000000ULL) return 0;  // NaR

    sign = uiA >> 63;
    if (sign) uiA = -uiA & 0xFFFFFFFFFFFFFFFFULL;

    if (uiA <= 0x7000000000000000ULL) return 0;              
    else if (uiA < 0x8800000000000000ULL) return 1;         
    else if (uiA <= 0x9200000000000000ULL) return 2;        
    else if (uiA > 0xFFFDFFFFFFFFFFFFULL)                  
        return sign ? INT64_MIN : INT64_MAX;

    // Decode regime
    uiA -= 0x8000000000000000ULL;
    while (uiA & 0x4000000000000000ULL) {
        scale += 4;
        uiA = (uiA - 0x4000000000000000ULL) << 1;
    }

    // Skip regime stop bit
    uiA <<= 1;

    // Decode exponent (2 bits for es=2)
    if (uiA & 0x4000000000000000ULL) scale += 2;
    if (uiA & 0x2000000000000000ULL) scale += 1;

    // Left-justify fraction into iZ128
    iZ128 = (__uint128_t)((uiA | 0x2000000000000000ULL) & 0x3FFFFFFFFFFFFFFFULL) << 66;

    if (scale < 127) {
        mask = (__uint128_t)1 << (127 - scale);

        bitLast = (iZ128 & mask);
        mask >>= 1;
        tmp = (iZ128 & mask);
        bitNPlusOne = tmp;
        iZ128 ^= tmp;

        tmp = iZ128 & (mask - 1);
        bitsMore = tmp;
        iZ128 ^= tmp;

        if (bitNPlusOne && (bitLast | bitsMore)) {
            iZ128 += (mask << 1);
        }

        iZ = (int_fast64_t)(iZ128 >> (126 - scale)); // right-justify to int64
    } else if (scale <= 126) {
        iZ = (int_fast64_t)(iZ128 >> (126 - scale));
    } else {
        iZ = (int_fast64_t)(iZ128 << (scale - 126));
    }

    return sign ? -iZ : iZ;
}
