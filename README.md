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

```
gcc -O2 -Isource/include -o main main.c build/Linux-x86_64-GCC/softposit.a -lm
./main
```


