#include "softposit.h"
#include <math.h>
#include <inttypes.h>
#include <string.h>

uint64_t double_to_u64(double d){
    uint64_t u;
    memcpy(&u, &d, sizeof u);
    return u;
}

int main(){
    //Compare the ULP(Unit in the Last Place)
    //It shows the advantage of posit number at the precision of the small number
    
    //2^(-28) <= d1 < 2^28 is better than double
    //d1 >= 2^32 or d1 < 2^(-32) is worse than double
    double value = 1.0;
    for(int i = 1; i <= 28; i++) value *= 2; 

    //Double ULP
    double d1 = value;
    double d2 = nextafter(d1, d1 + 1.0);//The next nubmer that can display from d1 to d1 + 1
    double diff = d2 - d1;
    uint64_t bits;

    bits = double_to_u64(d1);
    printf("d1 in Hex is 0x%016" PRIx64 "\n", bits);
    bits = double_to_u64(d2);
    printf("d2 in Hex is 0x%016" PRIx64 "\n", bits);
    printf("\n");

    //Posit64 ULP
    posit64_t p1, p2;
    union ui64_p64 u;
    p1 = convertDoubleToP64(value);

    u.p = p1;
    u.ui += 1;
    p2 = u.p;

    printf("p1 in Hex is ");
    printHex(p1.v);

    printf("p2 in Hex is ");
    printHex(p2.v);
 
    double pd1 = convertP64ToDouble(p1);
    double pd2 = convertP64ToDouble(p2);
    
    printf("\n");

    //The diff is smaller, the higher precision around that value
    //The diff is larger, the lower precision around that value
    printf("Double at value: d1 = %.20f\n", d1);
    printf("Double next number: d2 = %.20f\n", d2);
    printf("Double diff = %.20e\n", diff);//should be a very small number

    printf("\n");

    printf("Posit64 at value: p1 = %.20f\n", pd1);
    printf("Posit64 next number: p2 = %.20f\n", pd2);
    printf("Posit64 diff = %.20e\n", pd2 - pd1);//should be 0

    printf("\n");

    //It shows that we do not need three different bit pattern to represent infinity, -infinity, and Nan in posit system
    posit64_t pA;
    double dA;
    dA = INFINITY;
    bits = double_to_u64(dA);
    pA = convertDoubleToP64(dA);
    printf("double infinity is 0x%016" PRIx64 "\n", bits);
    printf("Posit64 infinity is ");
    printHex(pA.v);
    printf("\n");

    dA = -INFINITY;
    bits = double_to_u64(dA);
    pA = convertDoubleToP64(dA);
    printf("double negative infinity is 0x%016" PRIx64 "\n", bits);
    printf("Posit64 negative infinity is ");
    printHex(pA.v);
    printf("\n");

    dA = NAN;
    bits = double_to_u64(dA);
    pA = convertDoubleToP64(dA);
    printf("double NAN is 0x%016" PRIx64 "\n", bits);
    printf("Posit64 NAN is ");
    printHex(pA.v);
    printf("\n");

    return 0;
}
