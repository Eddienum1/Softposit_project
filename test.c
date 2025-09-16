#include "softposit.h"
#include <math.h>

int main() {
/*
    int N = 20;
    double approx_log = (1.0 - pow(10.0, -N)) / 9.0; // ~1/9
    double exact = exp(approx_log);

    float prod_f = 1.0f;
    for (int n = 1; n <= N; n++) {
        prod_f *= (1.0f + powf(10.0f, -n));
    }

    double prod_d = 1.0;
    for (int n = 1; n <= N; n++) {
        prod_d *= (1.0 + pow(10.0, -n));
    }

    posit64_t prod_p = convertDoubleToP64(1.0);
    for (int n = N; n >= 1; n--) {
        posit64_t term = convertDoubleToP64(1.0 + pow(10.0, -n));
        prod_p = p64_mul(prod_p, term);
    }

    double result = convertP64ToDouble(prod_p);

    printf("The exact value is %.15f\n", exact);
    printf("The result of float is %.15f\n", prod_f);
    printf("The result of double is %.15f\n", prod_d);
    printf("The result of posit64 is %.15f\n", result);
*/

    //Compare the ULP(Unit in the Last Place)
    //It shows the advantage of posit number at the precision of the extremely small number
    
    //Double ULP
    double d1 = 1.0;
    double d2 = nextafter(d1, 2.0);//The next nubmer that can display from 1.0 to 2.0
    double double_diff = d2 - d1;

    //Posit64 ULP
    posit64_t p1, p2;
    union ui64_p64 u;
    p1 = convertDoubleToP64(1.0);

    u.p = p1;
    u.ui += 1;
    p2 = u.p;

    printf("pd1 in Hex is ");
    printHex(p1.v);

    printf("pd2 in Hex is ");
    printHex(p2.v);
 
    double pd1 = convertP64ToDouble(p1);
    double pd2 = convertP64ToDouble(p2);

    printf("\n");

    //The output shows that double can represent the small difference (ULP) near 1.0
    //However, the difference between two adjacent posit64 numbers near 1.0 is smaller than double's precision
    //So when converting to double, this tiny difference cannot be displayed
    printf("Double @1.0: d1 = %.17f\n", d1);
    printf("Double candidate: d2 = %.17f\n", d2);
    printf("Double diff = %.20e\n", double_diff);//should be a very small number

    printf("\n");

    printf("Posit64 @1.0: p1 = %.17f\n", pd1);
    printf("Posit64 candidate: p2 = %.17f\n", pd2);
    printf("Posit64 diff = %.20e\n", pd2 - pd1);//should be 0

    return 0;
}
