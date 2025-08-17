#include <stdint.h>
#include "platform.h"
#include "internals.h"

posit64_t ui64_to_p64(uint64_t a) {
    int_fast8_t k, log2 = 63;  
    union ui64_p64 uZ;
    uint_fast64_t uiA;
    uint_fast64_t fracA, mask = 0x8000000000000000ULL;
    uint_fast64_t expA;

   
    if (a > 18445618173802708991ULL) { 
        uiA = 0x7FFFFFFFFFFFFFFFULL;    
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
