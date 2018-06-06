#include <stdio.h>
#include <math.h>

using namespace std;

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Class definitions
#include "Vector.h"
#include "Tensor.h"

// Class functions
#include "Vector.c"
#include "Tensor.c"

// Prototypes
Tensor Dyadic_Product(const Vector V1,const Vector V2);

////////////////////////////////////////////////////////////////////////////////
// Function definitions

Tensor Dyadic_Product(const Vector V1,const Vector V2) {
  Tensor T;

  /* Assign the elements of our dyadic product using nested for loop. Note that
     We only cycle through the columns. Let S denote the dyadic product of V1
     and v2. The jth column of S is equal to V1*V2[j] (scale the vector V1 by
     the jth component of V2)
  */
  for(int i = 0; i < 3; i++) {
    T(i,0) = V1[i]*V2[0];
    T(i,1) = V1[i]*V2[1];
    T(i,2) = V2[i]*V2[2];
  } //   for(int i = 0; i < 3; i++) {

  return T;
} // Tensor Dyatic_Product(const Vector V2,const Vector V2) {

int main(int argc, char *argv[]) {

  Vector V1 = {1,2,3};
  V1(2) = 5;
  V1[1] = 4;
  V1.Print();
  printf("\n");

  Vector V2 = {1,2,3};
  V2 += V2;
  V2 += {1,2,3};
  V2.Print();
  printf("\n");

  Tensor T0 = Dyadic_Product(V1,V2);
  T0.Print();
  printf("\n");

  Tensor T1(1,2,3,
            4,5,6,
            7,8,9);
  Tensor T2 = {1,4,7,
            2,5,8,
            3,6,9};

  Tensor T3 = T1*T2;
  T3.Print();
  printf("\n");

  Tensor T4 = T1;
  T4 += T1;
  T4 += {1,2,3,4,5,6,7,8,9};
  T4.Print();
  printf("\n");

  Tensor T5 = {1,5,8,
               12,6,3,
               6,9,11};
  Tensor T6 = T5.Inverse();
  T6.Print();

  return 0;
} // int main() {
