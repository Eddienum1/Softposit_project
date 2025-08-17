#include <stdint.h>

#include "platform.h"
#include "internals.h"

posit64_t ui32_to_p64(uint32_t a) {
    int_fast8_t k, log2 = 31; 
    union ui64_p64 uZ;
    uint_fast64_t uiA;
    uint_fast64_t expA, mask = 0x80000000ULL, fracA;

    if (a > 4294966271U) {
        uiA = 0x7FC0000000000000ULL;
    }
    else if (a < 0x2) {
        uiA = ((uint_fast64_t)a) << 62;
    }
    else {
        fracA = a;
        while (!(fracA & mask)) {
            log2--;
            fracA <<= 1;
        }

        k = log2 >> 2;

        expA = ((uint_fast64_t)(log2 & 0x3)) << (59 - k);

        fracA ^= mask;

        uiA = (0x7FFFFFFFFFFFFFFFULL ^ (0x3FFFFFFFFFFFFFFFULL >> k)) | expA | (fracA >> (k + 5));

        mask = 0x10ULL << k;

        if (mask & fracA) {
            if (((mask - 1) & fracA) | ((mask << 1) & fracA)) {
                uiA++;
            }
        }
    }

    uZ.ui = uiA;
    return uZ.p;
}
