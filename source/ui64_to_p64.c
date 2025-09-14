#include <stdint.h>
#include "platform.h"
#include "internals.h"

posit64_t ui64_to_p64(uint64_t a) {
    int_fast8_t k, log2 = 63;  
    union ui64_p64 uZ;
    uint_fast64_t uiA;
    uint_fast64_t fracA, mask = 0x8000000000000000;
    uint_fast64_t expA;

   
    if (a > 18445618173802708991ULL) { 
        uiA = 0x7FFFC00000000000;
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

        k = (log2 >> 2);

		expA = (uint_fast64_t)(log2 & 0x3) << (59 - k);

		fracA = (fracA ^ mask);

		uiA = ((uint_fast64_t)0x7FFFFFFFFFFFFFFF ^ ((uint_fast64_t)0x3FFFFFFFFFFFFFFF >> k ))
          | (uint_fast64_t)expA
          | (uint_fast64_t)(fracA >> (k + 4));
	
		mask = 0x4000000000000000 << k;  //bitNPlusOne

		if (mask & fracA) {
			if (((mask - 1) & fracA) | ((mask << 1) & fracA)) uiA++;
		}
    }

    uZ.ui = uiA;
    return uZ.p;
}
