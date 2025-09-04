#include "softposit.h"

int main() {

    posit64_t pA, pB, pC, pZ;
    double dA, dB, dC, dZ;

    pA = castP64(0x50616151ABCD5955);
    pB = castP64(0x5681656461331664);
    pC = castP64(0x1646AB3315DE255F);

    //Testing posit64 convert to double
    
    dA = convertP64ToDouble(pA);
    dB = convertP64ToDouble(pB);
    dC = convertP64ToDouble(pC);
    
    printf("pA convert to double is %.15f\n", dA);
    printf("pB convert to double is %.15f\n", dB);
    printf("pC convert to double is %.15f\n", dC);

    printf("\n");

    //Testing double convert to posit64

    pA = convertDoubleToP64(dA);
    pB = convertDoubleToP64(dB);

    printf("dA convert to posit64 is ");
    printHex(pA.v);

    printf("dB convert to posit64 is ");
    printHex(pB.v);

    printf("\n");

    //Reassign to avoid errors
    pA = castP64(0x50616151ABCD5955);
    pB = castP64(0x5681656461331664);
    pC = castP64(0x1646AB3315DE255F);

    //Testing division

    pZ = p64_div(pA, pB);
    dZ = convertP64ToDouble(pZ);
    printf("pA / pB in Hex: ");
    printHex(pZ.v);
    printf("pA / pB in double: %.15f\n", dZ);

    printf("\n");

    //Testing mul

    pZ = p64_mul(pA, pB);
    dZ = convertP64ToDouble(pZ);
    printf("pA * pB in Hex: ");
    printHex(pZ.v);
    printf("pA * pB in double: %.15f\n", dZ);

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

    /*
    //Testing sqrt

    pZ = p64_sqrt(pA);
    dZ = convertP64ToDouble(pZ);
    printf("sqrt(pA) in Hex: ");
    printHex(pZ.v);
    printf("sqrt(pA) in double: %.15f\n", dZ);

    printf("\n");

    pZ = p64_sqrt(pB);
    dZ = convertP64ToDouble(pZ);
    printf("sqrt(pB) in Hex: ");
    printHex(pZ.v);
    printf("sqrt(pB) in double: %.15f\n", dZ);

    printf("\n");

    //Testing add

    pZ = p64_add(pA, pB);
    dZ = convertP64ToDouble(pZ);
    printf("pA + pB in Hex: ");
    printHex(pZ.v);
    printf("pA + pB in double: %.15f\n", dZ);

    printf("\n");

    //Testing sub

    pZ = p64_sub(pA, pB);
    dZ = convertP64ToDouble(pZ);
    printf("pA - pB in Hex: ");
    printHex(pZ.v);
    printf("pA - pB in double: %.15f\n", dZ);

    printf("\n");

    //Testing muladd

    pZ = p64_mulAdd(pA, pB, pC);
    dZ = convertP64ToDouble(pZ);
    printf("pA * pB + pC in Hex: ");
    printHex(pZ.v);
    printf("pA * pB + pC in double: %.15f\n", dZ);

    printf("\n");

    */

    //Testing int32 to posit64

    int32_t iA = 13545696;
    pZ = i32_to_p64(iA);
    dZ = convertP64ToDouble(pZ);
    printf("iA convert to posit64 in Hex: ");
    printHex(pZ.v);
    printf("iA convert to posit64 in double: %.15f\n", dZ);

    printf("\n");

    //Testing int64 to posit64

    int64_t iB = INT64_MIN;
    pZ = i64_to_p64(iB);
    dZ = convertP64ToDouble(pZ);
    printf("iB convert to posit64 in Hex: ");
    printHex(pZ.v);
    printf("iB convert to posit64 in double: %.15f\n", dZ);

    printf("\n");

    return 0;
}
