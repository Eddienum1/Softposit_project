#include "platform.h"
#include "internals.h"

extern const uint_fast64_t softposit_approxRecipSqrt0[];
extern const uint_fast64_t softposit_approxRecipSqrt1[];

posit64_t p64_sqrt(posit64_t pA) {
    union ui64_p64 uA;
    uint_fast64_t uiA, uiZ, fracA, expA, index, r0, shift, expZ, mask;
    __uint128_t eSqrR0, fracZ, recipSqrt, sigma0, sqrSigma0, negRem, shiftedFracZ;
    int_fast32_t eps, shiftZ;

    uA.p = pA;
    uiA = uA.ui;

    // Handle NaR and negative values
    if (uiA & 0x8000000000000000ULL) {
        uA.ui = 0x8000000000000000ULL;
        return uA.p;
    } else if (!uiA) {
        return uA.p;
    }

    // Decode regime and exponent
    if (uiA & 0x4000000000000000ULL) {
        shiftZ = -2;
        while (uiA & 0x4000000000000000ULL) {
            shiftZ += 2;
            uiA = (uiA << 1) & 0xFFFFFFFFFFFFFFFFULL;
        }
    } else {
        shiftZ = 0;
        while (!(uiA & 0x4000000000000000ULL)) {
            shiftZ -= 2;
            uiA = (uiA << 1) & 0xFFFFFFFFFFFFFFFFULL;
        }
    }

    uiA &= 0x3FFFFFFFFFFFFFFFULL;
    expA = (uiA >> 60);
    shiftZ += (expA >> 1);
    expA = (0x1 ^ (expA & 0x1));
    uiA &= 0x0FFFFFFFFFFFFFFFULL;
    fracA = (uiA | 0x1000000000000000ULL);  // add hidden bit (adjust based on ES)

    // Lookup 1/sqrt approximation
    index = ((fracA >> 56) & 0xE) + expA;
    eps = (int_fast32_t)((fracA >> 40) & 0xFFFF);
    r0 = softposit_approxRecipSqrt0[index]
         - ((softposit_approxRecipSqrt1[index] * eps) >> 20);

    // Newton-Raphson refinement
    eSqrR0 = (__uint128_t)r0 * r0;
    if (!expA) eSqrR0 <<= 1;
    sigma0 = ((((__uint128_t)1 << 64) - ((eSqrR0 * fracA) >> 20)) & 0xFFFFFFFFFFFFFFFFULL);
    recipSqrt = ((__uint128_t)r0 << 20) + (((__uint128_t)r0 * sigma0) >> 21);

    sqrSigma0 = (sigma0 * sigma0) >> 35;
    recipSqrt += (((recipSqrt + (recipSqrt >> 2) - ((uint_fast64_t)r0 << 19)) * sqrSigma0) >> 46);

    fracZ = ((__uint128_t)fracA * recipSqrt) >> 63;
    if (expA) fracZ = (fracZ >> 1);

    // Assemble regime
    expZ = shiftZ & 0x3;
    if (shiftZ < 0) {
        shift = (-1 - shiftZ) >> 2;
        uiZ = 0x2000000000000000ULL >> shift;
    } else {
        shift = shiftZ >> 2;
        uiZ = 0x7FFFFFFFFFFFFFFFULL - (0x3FFFFFFFFFFFFFFFULL >> shift);
    }

    // Final refinement
    fracZ++;
    if (!(fracZ & 0xF)) {
        shiftedFracZ = fracZ >> 1;
        negRem = (shiftedFracZ * shiftedFracZ) & 0x3FFFFFFFFFFFFFFFFULL;
        if (negRem & 0x20000000000000000ULL) {
            fracZ |= 1;
        } else {
            if (negRem) fracZ--;
        }
    }

    // Round to nearest
    fracZ &= 0xFFFFFFFFFFFFFFFFULL;
    mask = (1ULL << (4 + shift));
    if (mask & fracZ) {
        if (((mask - 1) & fracZ) | ((mask << 1) & fracZ)) fracZ += (mask << 1);
    }

    // Assemble final result
    uA.ui = uiZ | (expZ << (59 - shift)) | (uint_fast64_t)(fracZ >> (5 + shift));
    return uA.p;
}
