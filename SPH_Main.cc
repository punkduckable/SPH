#include <stdio.h>
#include <math.h>
#include <time.h>

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Definitions
#define PI 3.1415926535897932384626

// Header files
#include "Classes.h"
#include "Vector.h"
#include "Tensor.h"
#include "Particle.h"
#include "List.h"
#include "Tests.h"

// Prototypes
Tensor Dyadic_Product(const Vector & V1,const Vector & V2);

// Source files
#include "Vector.c"
#include "Tensor.c"
#include "Particle.c"
#include "List.c"
#include "Tests.c"

////////////////////////////////////////////////////////////////////////////////
// Initialize static members of particle class

double Particle::density = 1;
double Particle::E = 1;
double Particle::alpha = 50;                    // What was used in the paper
Vector M_Set;
Vector Particle::M = M_Set;
double Particle::mu = 1;
double Particle::mu0 = 1;
double Particle::k1 = 1;
double Particle::k2 = 1;
double Particle::h = 4;
double Particle::Shape_Function_Amp = 15./(PI*h*h*h*h*h*h);

////////////////////////////////////////////////////////////////////////////////
// Function definitions

Tensor Dyadic_Product(const Vector & V1,const Vector & V2) {
  Tensor T;

  /* Assign the elements of our dyadic product using nested for loop. Note that
     We only cycle through the columns. Let S denote the dyadic product of V1
     and v2. The jth column of S is equal to V1*V2[j] (scale the vector V1 by
     the jth component of V2)
  */

  /* Unrolled loop (runs slower than 1 rolled loop for some reason)
  T[3*0 + 0] = V1[0]*V2[0];            // i = 0, j = 0
  T[3*0 + 1] = V1[0]*V2[1];            // i = 0, j = 1
  T[3*0 + 2] = V1[0]*V2[2];            // i = 0, j = 2

  T[3*1 + 0] = V1[1]*V2[0];            // i = 1, j = 0
  T[3*1 + 1] = V1[1]*V2[1];            // i = 1, j = 1
  T[3*1 + 2] = V1[1]*V2[2];            // i = 1, j = 2

  T[3*2 + 0] = V1[2]*V2[0];            // i = 2, j = 0
  T[3*2 + 1] = V1[2]*V2[1];            // i = 2, j = 1
  T[3*2 + 2] = V1[2]*V2[2];            // i = 2, j = 2
  */

  // Old loop. (works better than 9 statements with O2 optimization)
  for(int i = 0; i < 3; i++) {
    T[3*i + 0] = V1[i]*V2[0];
    T[3*i + 1] = V1[i]*V2[1];
    T[3*i + 2] = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++)

  return T;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {

int main() {
  // Run Vector, Tensor, List tests.
  //Vector_Tests();
  //Tensor_Tests();
  //List_Tests();
  Particle_Tests();
  //Timing_Tests();

  Vector V1{1,2,3}, V2{92.392,-203.29, 5.2039};
  Tensor T1, T2;
  Tensor A = {1, -20, 39,
              6 ,2.293, -32.3020,
              .20392, .592, -.0001993};
  T1 = Dyadic_Product(V1, (A.Inverse()*V2));
  T2 = Dyadic_Product(V1, V2)*Transpose(Inverse(A));

  T1.Print();
  T2.Print();

  return 0;
} // int main() {
