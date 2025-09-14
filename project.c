#include "softposit.h"
#include <math.h>

int main() {

    int option = 0, temp = 0; //function choice and temp for flag
    bool flag = false; //do it again or not

    //at most 3 kinds of parameter and an answer
    posit32_t pA1, pB1, pC1, pZ1;
    posit64_t pA2, pB2, pC2, pZ2;
    double dA, dB, dC, dZ;
    unsigned long long bits; //set bits to posit64
    int32_t iA; //for convert i32
    int64_t iB; //for convert i64
    uint32_t uiA;
    uint64_t uiB;

    do{
        printf("This is a simple calculator\n");
        printf("These are the function you can do\n");
        printf("1.add\n2.sub\n3.mul\n4.div\n5.muladd\n6.sqrt\n7.comparison\n8.round to int\n9.convert data type\n");
        printf("Please choose your action: ");

        if(scanf("%d", &option) != 1){
            printf("Invalid input!\n");
            return 1;
        }

        switch(option){

            //addition
            case 1:
                printf("Enter two numbers: ");
                if(scanf("%lf %lf", &dA, &dB) != 2){
                    printf("Invalid input");
                    return 1;
                }
                
                dZ = dA + dB;
                printf("%f + %f = %f in double data type\n", dA, dB, dZ);

                pA1 = convertDoubleToP32(dA);
                pB1 = convertDoubleToP32(dB);
                pZ1 = p32_add(pA1, pB1);
                dZ = convertP32ToDouble(pZ1);
                printf("%f + %f = %f in posit32 data type\n", dA, dB, dZ);

                pA2 = convertDoubleToP64(dA);
                pB2 = convertDoubleToP64(dB);
                pZ2 = p64_add(pA2, pB2);
                dZ = convertP64ToDouble(pZ2);
                printf("%f + %f = %f in posit64 data type\n", dA, dB, dZ);

                break;
            
            //subtraction
            case 2:
                printf("Enter two numbers: ");
                if(scanf("%lf %lf", &dA, &dB) != 2){
                    printf("Invalid input");
                    return 1;
                }
    
                dZ = dA - dB;
                printf("%f - %f = %f in double data type\n", dA, dB, dZ);

                pA1 = convertDoubleToP32(dA);
                pB1 = convertDoubleToP32(dB);
                pZ1 = p32_sub(pA1, pB1);
                dZ = convertP32ToDouble(pZ1);
                printf("%f - %f = %f in posit32 data type\n", dA, dB, dZ);

                pA2 = convertDoubleToP64(dA);
                pB2 = convertDoubleToP64(dB);
                pZ2 = p64_sub(pA2, pB2);
                dZ = convertP64ToDouble(pZ2);
                printf("%f - %f = %f in posit64 data type\n", dA, dB, dZ);

                break;
            
            //multiplication
            case 3:
                printf("Enter two numbers: ");
                if(scanf("%lf %lf", &dA, &dB) != 2){
                    printf("Invalid input!\n");
                    return 1;
                }

                dZ = dA * dB;
                printf("%f * %f = %f in double data type\n", dA, dB, dZ);

                pA1 = convertDoubleToP32(dA);
                pB1 = convertDoubleToP32(dB);
                pZ1 = p32_mul(pA1, pB1);
                dZ = convertP32ToDouble(pZ1);
                printf("%f * %f = %f in posit32 data type\n", dA, dB, dZ);

                pA2 = convertDoubleToP64(dA);
                pB2 = convertDoubleToP64(dB);
                pZ2 = p64_mul(pA2, pB2);
                dZ = convertP64ToDouble(pZ2);
                printf("%f * %f = %f in posit64 data type\n", dA, dB, dZ);

                break;

            //division
            case 4:
                printf("Enter two numbers: ");
                if(scanf("%lf %lf", &dA, &dB) != 2){
                    printf("Invalid input!\n");
                    return 1;
                }
                
                dZ = dA / dB;
                printf("%f / %f = %f in double data type\n", dA, dB, dZ);

                pA1 = convertDoubleToP32(dA);
                pB1 = convertDoubleToP32(dB);
                pZ1 = p32_div(pA1, pB1);
                dZ = convertP32ToDouble(pZ1);
                printf("%f / %f = %f in posit32 data type\n", dA, dB, dZ);

                pA2 = convertDoubleToP64(dA);
                pB2 = convertDoubleToP64(dB);
                pZ2 = p64_div(pA2, pB2);
                dZ = convertP64ToDouble(pZ2);
                printf("%f / %f = %f in posit64 data type\n", dA, dB, dZ);

                break;
            
            case 5:
                printf("Enter three numbers: ");
                if(scanf("%lf %lf %lf", &dA, &dB, &dC) != 3){
                    printf("Invalid input!\n");
                    return 1;
                }

                dZ = dA * dB;
                dZ = dZ + dC;
                printf("%f * %f + %f = %f in double data type\n", dA, dB, dC, dZ);

                pA1 = convertDoubleToP32(dA);
                pB1 = convertDoubleToP32(dB);
                pC1 = convertDoubleToP32(dC);
                pZ1 = p32_mulAdd(pA1, pB1, pC1);
                dZ = convertP32ToDouble(pZ1);
                printf("%f * %f + %f = %f in posit32 data type\n", dA, dB, dC, dZ);

                pA2 = convertDoubleToP64(dA);
                pB2 = convertDoubleToP64(dB);
                pC2 = convertDoubleToP64(dC);
                pZ2 = p64_mulAdd(pA2, pB2, pC2);
                dZ = convertP64ToDouble(pZ2);
                printf("%f * %f + %f = %f in posit64 data type\n", dA, dB, dC, dZ);

                break;
            
            case 6:
            /*
                printf("Enter one number: ");
                if(scanf("%lf", &dA) != 1){
                    printf("Invalid input!\n");
                }

                dZ = sqrt(dA);
                printf("The square root of %f is %f in double data type\n", dA, dZ);

                pA1 = convertDoubleToP32(dA);
                pZ1 = p32_sqrt(pA1);
                dZ = convertP32ToDouble(pZ1);
                printf("The square root of %f is %f in posit32 data type\n", dA, dZ);

                pA2 = convertDoubleToP64(dA);
                pZ2 = p64_sqrt(pA2);
                dZ = convertP64ToDouble(pZ2);
                printf("The square root of %f is %f in posit32 data type\n", dA, dZ);

                break;
            */

            //comparison
            case 7:
                printf("Enter two numbers: ");
                if(scanf("%lf %lf", &dA, &dB) != 2){
                    printf("Invalid input!\n");
                    return 1;
                }
                
                if(dA <= dB){
                    if(dA == dB) printf("%f is equal to %f in double data type\n", dA, dB);
                    else printf("%f is less than %f in double data type\n", dA, dB);
                }
                else printf("%f is greater than %f in double data type\n", dA, dB);

                pA1 = convertDoubleToP32(dA);
                pB1 = convertDoubleToP32(dB);
                if(p32_le(pA1, pB1)){
                    if(p32_eq(pA1, pB1)){
                        printf("%f is equal to %f in posit32 data type\n", dA, dB);
                    }
                    else printf("%f is less than %f in posit32 data type\n", dA, dB);
                }
                else printf("%f is greater than %f in posit32 data type\n", dA, dB);

                pA2 = convertDoubleToP64(dA);
                pB2 = convertDoubleToP64(dB);
                if(p64_le(pA2, pB2)){
                    if(p64_eq(pA2, pB2)){
                        printf("%f is equal to %f in posit64 data type\n", dA, dB);
                    }
                    else printf("%f is less than %f in posit64 data type\n", dA, dB);
                }
                else printf("%f is greater than %f in posit64 data type\n", dA, dB);

                break;

            //round to int
            case 8:
                printf("Enter one number: ");
                if(scanf("%lf", &dA) != 1){
                    printf("Invalid input!\n");
                }

                dZ = rint(dA);
                printf("%f round to int is %f in double data type\n", dA, dZ);

                pA1 = convertDoubleToP32(dA);
                pZ1 = p32_roundToInt(pA1);
                dZ = convertP32ToDouble(pZ1);
                printf("%f round to int is %f in posit32 data type\n", dA, dZ);

                pA2 = convertDoubleToP64(dA);
                pZ2 = p64_roundToInt(pA2);
                dZ = convertP64ToDouble(pZ2);
                printf("%f round to int is %f in posit64 data type\n", dA, dZ);

                break;

            //change data type
            case 9:
                printf("1.double convert to posit64\n"
                    "2.posit64 convert to double\n"
                    "3.int32 convert to posit64\n"
                    "4.int64 convert to posit64\n"
                    "5.posit64 convert to i32\n"
                    "6.posit64 convert to i64\n"
                    "7.posit64 convert to ui32\n"
                    "8.posit64 convert to ui64\n"
                    "9.ui32 convert to posit64\n"
                    "10.ui64 convert to posit64\n");

                printf("Select the function to convert: ");

                if(scanf("%d", &option) != 1){
                    printf("Invalid input!\n");
                    return 1;
                }

                switch(option){
                    case 1:
                        printf("Enter one number(in double): ");
                        if(scanf("%lf", &dA) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }
                        
                        pA1 = convertDoubleToP32(dA);
                        pA2 = convertDoubleToP64(dA);

                        printf("%f convert to posit32 is ", dA);
                        printHex(pA1.v);
                        printf("%f convert to posit64 is ", dA);
                        printHex(pA2.v);

                        break;

                    case 2:
                        printf("Enter the bit pattern(in Hex): ");
                        if(scanf("%llx", &bits) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }

                        pA2.v = bits;
                        dA = convertP64ToDouble(pA2);
                        printHex(pA2.v);
                        printf("convert to double is %f\n", dA);
                        break;

                    case 3:
                        printf("Enter one number(in int32): ");
                        if(scanf("%d", &iA) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }

                        pA1 = i32_to_p32(iA);
                        pA2 = i32_to_p64(iA);

                        printf("%d convert to posit32 is ", iA);
                        printHex(pA1.v);
                        printf("%d convert to posit64 is ", iA);
                        printHex(pA2.v);

                        break;

                    case 4:
                        printf("Enter one number(in int64): ");
                        if(scanf("%ld", &iB) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }
                    
                        pA1 = i64_to_p32(iB);
                        pA2 = i64_to_p64(iB);

                        printf("%ld convert to posit32 is ", iB);
                        printHex(pA1.v);
                        printf("%ld convert to posit64 is ", iB);
                        printHex(pA2.v);

                        break;

                    case 5:
                        printf("Enter the bit pattern(in Hex): ");
                        if(scanf("%llx", &bits) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }
                        pA2.v = bits;
                        iA = p64_to_i32(pA2);
                        printHex(pA2.v);
                        printf("convert to i32 is %d\n", iA);

                        break;
                    
                    case 6:
                        printf("Enter the bit pattern(in Hex): ");
                        if(scanf("%llx", &bits) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }
                        pA2.v = bits;
                        iA = p64_to_i64(pA2);
                        printHex(pA2.v);
                        printf("convert to i64 is %d\n", iA);

                        break;
                    
                    case 7:
                        printf("Enter the bit pattern(in Hex): ");
                        if(scanf("%llx", &bits) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }
                        pA2.v = bits;
                        uiA = p64_to_ui32(pA2);
                        printHex(pA2.v);
                        printf("convert to ui32 is %d\n", uiA);

                        break;

                    case 8:
                        printf("Enter the bit pattern(in Hex): ");
                        if(scanf("%llx", &bits) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }
                        pA2.v = bits;
                        uiB = p64_to_ui64(pA2);
                        printHex(pA2.v);
                        printf("convert to ui64 is %ld\n", uiB);

                    case 9:
                        printf("Enter one number(in uint32): ");
                        if(scanf("%d", &uiA) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }

                        pA1 = ui32_to_p32(uiA);
                        pA2 = ui32_to_p64(uiA);

                        printf("%d convert to posit32 is ", uiA);
                        printHex(pA1.v);
                        printf("%d convert to posit64 is ", uiA);
                        printHex(pA2.v);

                        break;

                    case 10:
                         printf("Enter one number(in uint64): ");
                        if(scanf("%ld", &uiB) != 1){
                            printf("Invalid input!\n");
                            return 1;
                        }

                        pA1 = ui64_to_p32(uiB);
                        pA2 = ui64_to_p64(uiB);

                        printf("%ld convert to posit32 is ", uiB);
                        printHex(pA1.v);
                        printf("%ld convert to posit64 is ", uiB);
                        printHex(pA2.v);

                        break;

                    default:

                        break;

                }

            default:

                break;
        }

        printf("Do you want to do another round? (0 = No, 1 = Yes): ");
        if (scanf("%d", &temp) != 1){
            printf("Invalid input!\n");
            return 1;
        }
        flag = (temp != 0);

    }while(flag == true);

    return 0;
}
