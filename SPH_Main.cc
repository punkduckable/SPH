#include <stdio.h>
#include <math.h>
#include <time.h>

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Definitions
#define PI 3.1415926535897932384626

// Header files
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
double Particle::h = 2;
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
  for(int i = 0; i < 3; i++) {
    T(i,0) = V1[i]*V2[0];
    T(i,1) = V1[i]*V2[1];
    T(i,2) = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++) {

  return T;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {

int main() {
  // Run Vector, Tensor, List tests.
  //Vector_Tests();
  Tensor_Tests();
  //List_Tests();
  //Particle_Tests();
  Timing_Tests();

  return 0;
} // int main() {
