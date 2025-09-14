#include "platform.h"
#include "internals.h"

posit64_t softposit_addMagsP64(uint_fast64_t uiA, uint_fast64_t uiB) {

    uint_fast16_t regA, regB;
    unsigned __int128 frac128A = 0, frac128B = 0; // 128-bit intermediate
    uint_fast64_t fracA = 0, regime, tmp;
    bool sign, regSA, regSB, rcarry = 0, bitNPlusOne = 0, bitsMore = 0;
    int_fast32_t kA = 0; // widen kA
    int_fast32_t expA, expB;
    int_fast32_t shiftRight;
    union ui64_p64 uZ;

    sign = signP64UI(uiA);
    if (sign) {
        uiA = -uiA & 0xFFFFFFFFFFFFFFFFULL;
        uiB = -uiB & 0xFFFFFFFFFFFFFFFFULL;
    }

    if ((int_fast64_t)uiA < (int_fast64_t)uiB) {
        uiA ^= uiB;
        uiB ^= uiA;
        uiA ^= uiB;
    }

    regSA = signregP64UI(uiA);
    regSB = signregP64UI(uiB);

    // decode A
    tmp = (uiA << 2) & 0xFFFFFFFFFFFFFFFFULL;
    if (regSA) {
        while (tmp >> 63) {
            kA++;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFFULL;
        }
    } else {
        kA = -1;
        while (!(tmp >> 63)) {
            kA--;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFFULL;
        }
        tmp &= 0x7FFFFFFFFFFFFFFFULL;
    }

    expA = (int_fast32_t)((tmp >> 62) & 0x3);  // es = 2
    // build 128-bit significand for A: hidden bit (1<<62) | (tmp<<1 & mask62), then <<64
    frac128A = ( (unsigned __int128)0x4000000000000000ULL | ((unsigned __int128)tmp << 1 & (unsigned __int128)0x7FFFFFFFFFFFFFFFULL) );
    frac128A <<= 64;

    shiftRight = kA;

    // decode B
    tmp = (uiB << 2) & 0xFFFFFFFFFFFFFFFFULL;
    if (regSB) {
        while (tmp >> 63) {
            shiftRight--;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFFULL;
        }
    } else {
        shiftRight++;
        while (!(tmp >> 63)) {
            shiftRight++;
            tmp = (tmp << 1) & 0xFFFFFFFFFFFFFFFFULL;
        }
        tmp &= 0x7FFFFFFFFFFFFFFFULL;
    }

    expB = (int_fast32_t)((tmp >> 62) & 0x3);
    frac128B = ( (unsigned __int128)0x4000000000000000ULL | ((unsigned __int128)tmp << 1 & (unsigned __int128)0x7FFFFFFFFFFFFFFFULL) );
    frac128B <<= 64;

    // compute shift in bits: 4*k + expA - expB
    shiftRight = (shiftRight << 2) + expA - expB;

    // shift frac128B right with sticky
    if (shiftRight > 127) {
        if (frac128B != 0) bitsMore = 1;
        frac128B = 0;
    } else if (shiftRight > 0) {
        unsigned __int128 lowmask = ( ((unsigned __int128)1 << shiftRight) - 1 );
        unsigned __int128 sticky = frac128B & lowmask;
        frac128B >>= shiftRight;
        if (sticky) bitsMore = 1;
    }
    // else shiftRight <= 0 : no right shift (rare given uiA>=uiB)

    // add aligned significands
    frac128A += frac128B;

    // carry out normalization (top bit at position 127)
    rcarry = ( (frac128A >> 127) & 1 );
    if (rcarry) {
        expA++;
        if (expA > 3) {
            kA++;
            expA &= 0x3;
        }
        // any low bit lost by >>1 becomes sticky
        if ((frac128A & 1) != 0) bitsMore = 1;
        frac128A >>= 1;
    }

    if (kA < 0) {
        regA = (uint_fast16_t)(-kA);
        regSA = 0;
        regime = 0x4000000000000000ULL >> regA;
    } else {
        regA = (uint_fast16_t)(kA + 1);
        regSA = 1;
        regime = 0x7FFFFFFFFFFFFFFFULL - (0x7FFFFFFFFFFFFFFFULL >> regA);
    }

    if (regA > 62) {
        uZ.ui = regSA ? 0x7FFFFFFFFFFFFFFFULL : 0x1ULL;
    } else {
        // remove hidden + exp bits: drop (regA + es) bits from top of frac128A
        unsigned __int128 mask127 = ( ((unsigned __int128)1 << 127) - 1 );
        unsigned __int128 fracField = (frac128A & mask127) >> (regA + 2); // drop regime + es

        // high 64 bits become fracA (like p32 used >>32)
        fracA = (uint_fast64_t)(fracField >> 64);

        if (regA <= 60) {
            // bitNPlusOne is the next bit below the top 64 bits
            bitNPlusOne = ( (fracField >> 63) & 0x1 );
            expA <<= (60 - regA);
        } else {
            if (regA == 62) {
                bitNPlusOne = (expA & 0x2) != 0;
                bitsMore |= (expA & 0x1) != 0;
                expA = 0;
            } else if (regA == 61) {
                bitNPlusOne = (expA & 0x1) != 0;
                expA >>= 1;
            }
            if (fracA > 0) { fracA = 0; bitsMore = 1; }
        }

        uZ.ui = packToP64UI(regime, expA, fracA);

        if (bitNPlusOne) {
            if ((uint64_t)(fracField & 0xFFFFFFFFFFFFFFFFULL) != 0) bitsMore = 1;
            uZ.ui += (uZ.ui & 1) | (bitsMore ? 1ULL : 0ULL);
        }
    }

    if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFFFFFFFFFFULL;
    return uZ.p;
}

