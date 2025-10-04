/*
#include "platform.h"
#include "internals.h"

extern const uint_fast64_t softposit_approxRecipSqrt0[];
extern const uint_fast64_t softposit_approxRecipSqrt1[];

posit64_t p64_sqrt(posit64_t pA) {
    union ui64_p64 uA;
    uint_fast64_t index, r0, shift, fracA, expZ, expA;
    uint_fast64_t mask, uiA, uiZ;
    __uint128_t eSqrR0, fracZ, recipSqrt, sigma0, sqrSigma0, negRem, shiftedFracZ;
    int_fast32_t eps, shiftZ;

    uA.p = pA;
    uiA = uA.ui;

    // Handle NaR and negative values
    if (uiA & 0x8000000000000000) {
        uA.ui = 0x8000000000000000;
        return uA.p;
    } 
    
    else if (!uiA) {
        return uA.p;
    }

    // Decode regime and exponent
    if (uiA & 0x4000000000000000) {
        shiftZ = -2;
        while (uiA & 0x4000000000000000) {
            shiftZ += 2;
            uiA = (uiA << 1) & 0xFFFFFFFFFFFFFFFF;
        }
    } else {
        shiftZ = 0;
        while (!(uiA & 0x4000000000000000)) {
            shiftZ -= 2;
            uiA = (uiA << 1) & 0xFFFFFFFFFFFFFFFF;
        }
    }

    uiA &= 0x3FFFFFFFFFFFFFFF;
    expA = (uiA >> 60);
    shiftZ += (expA >> 1);
    expA = (0x1 ^ (expA & 0x1));
    uiA &= 0x0FFFFFFFFFFFFFFF;
    fracA = (uiA | 0x1000000000000000);  // add hidden bit (adjust based on ES)

    // Lookup 1/sqrt approximation
    index = ((fracA >> 56) & 0xE) + expA;
    eps = (int_fast32_t)((fracA >> 41) & 0xFFFFFFFF);
    r0 = softposit_approxRecipSqrt0[index]
         - ((softposit_approxRecipSqrt1[index] * eps) >> 20);

    // Newton-Raphson refinement
    eSqrR0 = (__uint128_t)r0 * r0;
    if (!expA) eSqrR0 <<= 1;
    sigma0 = ((((__uint128_t)1 << 64) - ((eSqrR0 * fracA) >> 52)) & 0xFFFFFFFFFFFFFFFF);
    recipSqrt = ((__uint128_t)r0 << 20) + (((__uint128_t)r0 * sigma0) >> 21);

    sqrSigma0 = (sigma0 * sigma0) >> 67;
    recipSqrt += (((recipSqrt + (recipSqrt >> 2) - ((uint_fast64_t)r0 << 19)) * sqrSigma0) >> 46);

    fracZ = ((__uint128_t)fracA * recipSqrt) >> 63;
    if (expA) fracZ = (fracZ >> 1);

    // Assemble regime
    expZ = shiftZ & 0x3;
    if (shiftZ < 0) {
        shift = (-1 - shiftZ) >> 2;
        uiZ = 0x2000000000000000 >> shift;
    } else {
        shift = shiftZ >> 2;
        uiZ = 0x7FFFFFFFFFFFFFFF - (0x3FFFFFFFFFFFFFFF >> shift);
    }

    // Final refinement
    fracZ++;
    if (!(fracZ & 0xF)) {
        shiftedFracZ = fracZ >> 1;
        negRem = (shiftedFracZ * shiftedFracZ) & 0x1FFFFFFFFFFFFFFF;
        if (negRem & 0x1000000000000000) {
            fracZ |= 1;
        } else {
            if (negRem) fracZ--;
        }
    }

    // Round to nearest
    fracZ &= 0xFFFFFFFFFFFFFFFF;
    mask = (1ULL << (4 + shift));
    if (mask & fracZ) {
        if (((mask - 1) & fracZ) | ((mask << 1) & fracZ)) fracZ += (mask << 1);
    }

    // Assemble final result
    uA.ui = uiZ | (expZ << (59 - shift)) | (uint_fast64_t)(fracZ >> (5 + shift));
    return uA.p;
}

*/

#include <math.h>
#include "platform.h"
#include "internals.h"


posit64_t p64_sqrt( posit64_t pA ) {
    union ui64_p64 uA;
    uA.p = pA;
    uint_fast64_t uiA = uA.ui;

    if (uiA & 0x8000000000000000ULL) {
        uA.ui = 0x8000000000000000ULL;
        return uA.p;
    }

    if (uiA == 0ULL) {
        return uA.p;
    }

    double dv = convertP64ToDouble(pA);
    if (dv < 0.0) {
        uA.ui = 0x8000000000000000ULL;
        return uA.p;
    }
    double res = sqrt(dv);

    return convertDoubleToP64(res);
}