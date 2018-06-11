#include <stdio.h>
#include <math.h>

using namespace std;

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Class definitions
#include "Vector.h"
#include "Tensor.h"
#include "Tests.h"

// Prototypes
Tensor Dyadic_Product(const Vector & V1,const Vector & V2);

// Class functions
#include "Vector.c"
#include "Tensor.c"
#include "Tests.c"

////////////////////////////////////////////////////////////////////////////////
// Function definitions

Tensor Dyadic_Product(const Vector & V1,const Vector & V2) {
  Tensor T;

  /* Assign the elements of our dyadic product using nested for loop. Note that
     We only cycle through the columns. Let S denote the dyadic product of V1
     and v2. The jth column of S is equal to V1*V2[j] (scale the vector V1 by
     the jth component of V2)
  */
  for(int i = 0; i < 3; i++) {
    T(i,0) = V1[i]*V2[0];
    T(i,1) = V1[i]*V2[1];
    T(i,2) = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++) {

  return T;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {

int main() {

  Vector_Tests();
  Tensor_Tests();

  return 0;
} // int main() {
