## Softposit_project

### GroupMembers: 莊順閎 王政崴

## 1. Introduction

Our goal is to build the 64-bit softposit system with 2 exponent bits.

Let the user can use some function(ex. add, sub...)to calculate the number or change the data type and make sure the answer is correct.

## 2. Build

Type the instruction in the terminal

```
cd SoftPosit/build/Linux-x86_64-GCC
make -j6 all
```

## 3. Link

If your source code is for example "main.c" and you want to create an executable "main".
The "main.c" should place at the same directory.

NOTICE: You have to go back to Softposit/ 

```
gcc -O2 -Isource/include -o main main.c build/Linux-x86_64-GCC/softposit.a -lm
./main
```

## 4. Function

advantage.c is for showing the range that the precision of 64 bit posit number is better than IEEE double

main.c is for testing the specific number and the result

project.c is a small calculator for comparing the difference between double, 32 bit posit number, and 64 bit posit number

rand.c is for generating random double number and 64 bit posit number

test.c is for testing the random number and the result
