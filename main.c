#include "softposit.h"

int main() {

    posit64_t pA, pB, pZ;
    double dA, dB, dZ;

    pA = castP64(0x50000664656ADC55);
    pB = castP64(0x5681656461331664);

    //Testing posit64 convert to double
    
    dA = convertP64ToDouble(pA);
    dB = convertP64ToDouble(pB);
    
    printf("pA convert to double is %.15f\n", dA);
    printf("pB convert to double is %.15f\n", dB);

    printf("\n");

    //Testing double convert to posit64

    pA = convertDoubleToP64(dA);
    pB = convertDoubleToP64(dB);

    printf("dA convert to posit64 is ");
    printHex(pA.v);

    printf("dB convert to posit64 is ");
    printHex(pB.v);

    printf("\n");

    //Testing division

    pZ = p64_div(pA, pB);
    dZ = convertP64ToDouble(pZ);
    printf("pA / pB = %.15f\n", dZ);

    printf("\n");

    //Testing mul

    pZ = p64_mul(pA, pB);
    dZ = convertP64ToDouble(pZ);
    printf("pA * pB = %.15f\n", dZ);

    printf("\n");

    //Testing comparison

    if(p64_le(pA, pB)){
        if(p64_eq(pA, pB)){
            printf("pA is equal to pB\n");
        }
        else printf("pA is smaller than pB\n");
    }
    else printf("pA is greater than pB\n");

    printf("\n");

    //Testing round to integer

    pA = p64_roundToInt(pA);
    pB = p64_roundToInt(pB);

    dA = convertP64ToDouble(pA);
    dB = convertP64ToDouble(pB);

    printf("pA round to int is %g\n", dA);
    printf("pA round to int is %g\n", dB);

    printf("\n");
    return 0;
}
