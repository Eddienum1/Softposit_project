#include "platform.h"
#include "internals.h"

posit64_t p64_sub(posit64_t a, posit64_t b) {

    union ui64_p64 uA, uB, uZ;
    uint_fast64_t uiA, uiB;

    uA.p = a;
    uiA = uA.ui;
    uB.p = b;
    uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
    uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

    // NaR (Not a Real)
    if (uiA == 0x8000000000000000ULL || uiB == 0x8000000000000000ULL) {
#ifdef SOFTPOSIT_EXACT
        uZ.ui.v = 0x8000000000000000ULL;
        uZ.ui.exact = 0;
#else
        uZ.ui = 0x8000000000000000ULL;
#endif
        return uZ.p;
    }

    // Zero
    else if (uiA == 0 || uiB == 0) {
#ifdef SOFTPOSIT_EXACT
        uZ.ui.v = (uiA | -uiB);
        uZ.ui.exact = 0;
#else
        uZ.ui = (uiA | -uiB);
#endif
        return uZ.p;
    }

    // Different signs: add magnitudes
    if ((uiA ^ uiB) >> 63)
        return softposit_addMagsP64(uiA, (-uiB & 0xFFFFFFFFFFFFFFFFULL));
    else
        return softposit_subMagsP64(uiA, (-uiB & 0xFFFFFFFFFFFFFFFFULL));
}
