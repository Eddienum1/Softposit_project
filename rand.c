#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "softposit.h"


double randomDouble(double min, double max) {
    unsigned long long r = ((unsigned long long)rand() << 31) ^ rand();
    double frac = (double)r / (double)ULLONG_MAX;
    return min + frac * (max - min);
}

int main(void) {
    FILE *fpDouble, *fpPosit;
    int i;

    fpDouble = fopen("double_output.txt", "w");
    fpPosit  = fopen("posit64_output.txt", "w");

    if (fpDouble == NULL || fpPosit == NULL) {
        perror("Failed to open file");
        return 1;
    }

    srand((unsigned)time(NULL));

    // --- Output 100 random double ---
    for (i = 0; i < 100; i++) {
        double d = randomDouble(0, (double)INT_MAX);
        fprintf(fpDouble, "%.15f\n", d);
    }

    // --- Output 100 random posit64 ---
    for (i = 0; i < 100; i++) {
        double dRand = randomDouble(0, (double)INT_MAX);
        posit64_t p = convertDoubleToP64(dRand);
        fprintf(fpPosit, "%016llx\n", (unsigned long long)p.v);
    }

    fclose(fpDouble);
    fclose(fpPosit);

    printf("double_output.txt and posit64_output.txt have been generated\n");
    return 0;
}
