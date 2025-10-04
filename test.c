#include "softposit.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *fp;
    posit64_t pA, pB, pC, pZ;
    double dZ;
    unsigned long long hex;
    unsigned long long buffer[200];
    int count = 0, i, numTests, line = 0;

    setvbuf(stdout, NULL, _IONBF, 0);

    pC = castP64(0x1646AB3315DE255F);

    fp = fopen("posit64_output.txt", "r");
    if (fp == NULL) {
        perror("Failed to open file");
        return 1;
    }

    while (count < 200) {
        int ret = fscanf(fp, "%llx", &hex);
        line++;
        if (ret == 1) {
            buffer[count++] = hex;
        } else if (ret == EOF) {
            break;
        } else {
            printf("Line %d: failed to read hex, skipping\n", line);
            char linebuf[128];
            fgets(linebuf, sizeof(linebuf), fp);
        }
    }
    fclose(fp);

    if (count < 2) {
        printf("Not enough data in the file, at least two hex values needed\n");
        return 1;
    }

    printf("Read %d hex values from the file\n", count);

    numTests = count / 2;
    if (numTests > 50) numTests = 50;

    for (i = 0; i < numTests; i++) {
        pA.v = buffer[2*i];
        pB.v = buffer[2*i+1];

        printf("===== Testing #%d =====\n", i+1);
        printf("pA in Hex: 0x%016llx\n", (unsigned long long)pA.v);
        printf("pB in Hex: 0x%016llx\n", (unsigned long long)pB.v);
        printf("pC in Hex: 0x%016llx\n", (unsigned long long)pC.v);

        printf("\n");
        // p64_add
        pZ = p64_add(pA, pB);
        dZ = convertP64ToDouble(pZ);
        printf("pA + pB: 0x%016llx  %.15f\n", (unsigned long long)pZ.v, dZ);

        // p64_sub
        pZ = p64_sub(pA, pC);
        dZ = convertP64ToDouble(pZ);
        printf("pA - pC: 0x%016llx  %.15f\n", (unsigned long long)pZ.v, dZ);

        // p64_mul
        pZ = p64_mul(pB, pC);
        dZ = convertP64ToDouble(pZ);
        printf("pB * pC: 0x%016llx  %.15f\n", (unsigned long long)pZ.v, dZ);

        // p64_div
        pZ = p64_div(pA, pC);
        dZ = convertP64ToDouble(pZ);
        printf("pA / pC: 0x%016llx  %.15f\n", (unsigned long long)pZ.v, dZ);

        // p64_sqrt
        pZ = p64_sqrt(pA);
        dZ = convertP64ToDouble(pZ);
        printf("sqrt(pA): 0x%016llx  %.15f\n", (unsigned long long)pZ.v, dZ);

        // p64_mulAdd (pA * pB + pC)
        pZ = p64_mulAdd(pA, pB, pC);
        dZ = convertP64ToDouble(pZ);
        printf("pA * pB + pC: 0x%016llx  %.15f\n", (unsigned long long)pZ.v, dZ);

        printf("\n");

        // p64_to_i32
        int32_t i32_val = p64_to_i32(pA);
        printf("pA to int32: %d\n", i32_val);

        // p64_to_i64
        int64_t i64_val = p64_to_i64(pA);
        printf("pA to int64: %ld\n", i64_val);

        // p64_to_ui32
        uint32_t ui32_val = p64_to_ui32(pA);
        printf("pA to uint32: %u\n", ui32_val);

        // p64_to_ui64
        uint64_t ui64_val = p64_to_ui64(pA);
        printf("pA to uint64: %lu\n", ui64_val);

        printf("\n");

        //Testing comparison
        printf("pA = %.15f pB = %.15f\n", convertP64ToDouble(pA), convertP64ToDouble(pB));
        if(p64_le(pA, pB)){
            if(p64_eq(pA, pB)){
                printf("pA is equal to pB\n");
            }
            else printf("pA is smaller than pB\n");
        }
        else printf("pA is greater than pB\n");

        printf("\n");

        //round to int
        pZ = p64_roundToInt(pA);
        dZ = convertP64ToDouble(pZ);
        printf("pA round to int is %g\n", dZ);

        printf("\n");

    }
    printf("===== End =====\n");
    printf("\n");
    int32_t test_i32 = -12345;
    int64_t test_i64 = -9876543210LL;
    uint32_t test_ui32 = 12345U;
    uint64_t test_ui64 = 9876543210ULL;

    posit64_t from_i32 = i32_to_p64(test_i32);
    posit64_t from_i64 = i64_to_p64(test_i64);
    posit64_t from_ui32 = ui32_to_p64(test_ui32);
    posit64_t from_ui64 = ui64_to_p64(test_ui64);

    printf("==== Testing data types to Posit64 ====\n");
    printf("i32_to_p64(-12345): 0x%016llx  %.15f\n", (unsigned long long)from_i32.v, convertP64ToDouble(from_i32));
    printf("i64_to_p64(-9876543210): 0x%016llx  %.15f\n", (unsigned long long)from_i64.v, convertP64ToDouble(from_i64));
    printf("ui32_to_p64(12345): 0x%016llx  %.15f\n", (unsigned long long)from_ui32.v, convertP64ToDouble(from_ui32));
    printf("ui64_to_p64(9876543210): 0x%016llx  %.15f\n", (unsigned long long)from_ui64.v, convertP64ToDouble(from_ui64));

    return 0;
}
