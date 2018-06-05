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
Tensor Dyadic_Product(Vector V1, Vector V2);

////////////////////////////////////////////////////////////////////////////////
// Function definitions

Tensor Dyatic_Product(Vector V1, Vector V2) {
  Tensor S;

  /* Assign the elements of our dyadic product using nested for loop. Note that
     We only cycle through the columns. Let S denote the dyadic product of V1
     and v2. The jth column of S is equal to V1*V2[j] (scale the vector V1 by
     the jth component of V2)
  */
  for(int j = 0; j < 3; j++) {
    S.T[0*3 + j] = V1.V[0]*V2.V[j];
    S.T[1*3 + j] = V1.V[1]*V2.V[j];
    S.T[2*3 + j] = V1.V[2]*V2.V[j];
  } //   for(int j = 0; j < 3; j++) {

  return S;
} // Tensor Dyatic_Product(Vector V2, Vector V2) {

int main(int argc, char *argv[]) {

  Vector V1 = {1,2,3};
  V1(2) = 5;
  V1[1] = 4;
  V1.Print();

  Tensor T1(1,2,3,
            4,5,6,
            7,8,9);
  Tensor T2 = {1,4,7,
            2,5,8,
            3,6,9};

  Tensor T3 = T1*T2;

  T3.Print();

  return 0;
} // int main() {
