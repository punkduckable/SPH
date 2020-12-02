#include "Vector.h"
#include "Tensor/Tensor.h"
#include "Errors.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* File description:

In this file, I define methods for Vector objects. Theese methods are
designed to make vector objects work just like mathematical vectors.  */

////////////////////////////////////////////////////////////////////////////////
// Constructors, Destructor

Vector::Vector(void) {
  #ifdef OPERATION_COUNT
    OP_Count::V_Default_Constructor++;             // Increment operator count
  #endif
} // Vector::Vector(void) {



Vector::Vector(const double v0, const double v1, const double v2) {
  // Initialize components of vector using supplied components
  (*this).Ar[0] = v0;
  (*this).Ar[1] = v1;
  (*this).Ar[2] = v2;

  #ifdef OPERATION_COUNT
    OP_Count::V_Component_Constructor++;           // Increment operator count
  #endif
} // Vector::Vector(const double v0, const double v1, const double v2) {



Vector::Vector(const Vector & V_In) {
  // Initialize components of vector using supplied components
  (*this).Ar[0] = V_In.Ar[0];
  (*this).Ar[1] = V_In.Ar[1];
  (*this).Ar[2] = V_In.Ar[2];

  #ifdef OPERATION_COUNT
    OP_Count::V_Copy_Constructor++;                // Increment operator count
  #endif
} // Vector::Vector(const Vector & V_In) {



Vector::~Vector(void) { }





////////////////////////////////////////////////////////////////////////////////
// Vector equality

Vector & Vector::operator=(const double V_In[3]) {
  // Assign components of vector to V_In array
  (*this).Ar[0] = V_In[0];
  (*this).Ar[1] = V_In[1];
  (*this).Ar[2] = V_In[2];

  #ifdef OPERATION_COUNT
    OP_Count::V_Equality++;                        // Increment operator count
  #endif

  // Return this Vector
  return *this;
} // Vector & Vector::operator=(const double V_In[3]) {



Vector & Vector::operator=(const Vector & V_In) {
  // Assign components of V using V_In.
  (*this).Ar[0] = V_In.Ar[0];
  (*this).Ar[1] = V_In.Ar[1];
  (*this).Ar[2] = V_In.Ar[2];

  #ifdef OPERATION_COUNT
    OP_Count::V_Equality++;                        // Increment operator count
  #endif

  // Return this vector
  return *this;
} // Vector & Vector::operator=(const Vector V_In) {





////////////////////////////////////////////////////////////////////////////////
// Simple arithmetic operators

Vector Vector::operator+(const Vector & V_In) const {
  // Declare a sum Vector. This will be used to store the sum
  Vector Sum;

  Sum.Ar[0] = (*this).Ar[0] + V_In.Ar[0];
  Sum.Ar[1] = (*this).Ar[1] + V_In.Ar[1];
  Sum.Ar[2] = (*this).Ar[2] + V_In.Ar[2];

  #ifdef OPERATION_COUNT
    OP_Count::V_V_Addition++;                      // Increment operator count
  #endif

  return Sum;
} // Vector Vector::operator+(const Vector V) const {



Vector Vector::operator-(const Vector & V_In) const{
  // Declare a Diff vector. This will be used to store the difference.
  Vector Diff;

  Diff.Ar[0] = (*this).Ar[0] - V_In.Ar[0];
  Diff.Ar[1] = (*this).Ar[1] - V_In.Ar[1];
  Diff.Ar[2] = (*this).Ar[2] - V_In.Ar[2];

  #ifdef OPERATION_COUNT
    OP_Count::V_V_Subtraction++;                   // Increment operator count
  #endif

  return Diff;
} // Vector Vector::operator-(const Vector V) const {



Vector Vector::operator*(const double c) const {
  // Declare product vector
  Vector Prod;

  // Scale components of Prod by c.
  Prod.Ar[0] = (*this).Ar[0]*c;
  Prod.Ar[1] = (*this).Ar[1]*c;
  Prod.Ar[2] = (*this).Ar[2]*c;

  #ifdef OPERATION_COUNT
    OP_Count::V_S_Multiplication++;                // Increment operator count
  #endif

  return Prod;
} // Vector Vector::operator*(const double c) const {



Vector Vector::operator/(const double c) const {
  /* Check for divide by zero
  Note: we keep this as an exception rather than an assertion because dividing
  by zero can happen in particular simulations (rather than just in buggy code). */
  if(c == 0) {
    throw Divide_By_Zero("Divide by Zero exception: Thrown by Vector::operator/\n"
                         "You tried to divide a vector by zero. Bad!\n");
  } // if(c == 0) {

  #ifdef OPERATION_COUNT
    OP_Count::V_S_Multiplication++;                // This prevents double counting. /s uses *s
    OP_Count::V_S_Division++;                      // Increment operator count
  #endif

  return (*this)*(1./c);
} // Vector Vector::operator/(const double c) const {





////////////////////////////////////////////////////////////////////////////////
// Compound arithmetic operators

Vector & Vector::operator+=(const Vector & V_In) {
  (*this).Ar[0] += V_In.Ar[0];
  (*this).Ar[1] += V_In.Ar[1];
  (*this).Ar[2] += V_In.Ar[2];

  #ifdef OPERATION_COUNT
    OP_Count::Compound_V_V_Addition++;             // Increment operator count
  #endif

  // Return this vector
  return (*this);
} // Vector & Vector:operator+=(const Vector & V_In) {



Vector & Vector::operator-=(const Vector & V_In) {
  (*this).Ar[0] -= V_In.Ar[0];
  (*this).Ar[1] -= V_In.Ar[1];
  (*this).Ar[2] -= V_In.Ar[2];

  #ifdef OPERATION_COUNT
    OP_Count::Compound_V_V_Subtraction++;        // Increment operator count
  #endif

  // Return this vector
  return (*this);
} // Vector & Vector::operator-=(const Vector & V_In) {



Vector & Vector::operator*=(const double c) {
  (*this).Ar[0] *= c;
  (*this).Ar[1] *= c;
  (*this).Ar[2] *= c;

  #ifdef OPERATION_COUNT
    OP_Count::Compound_V_S_Multiplication++;       // Increment operator count
  #endif

  // Return this vector (now scalled by c)
  return (*this);
} // Vector & Vector::operator*=(const double c) {





////////////////////////////////////////////////////////////////////////////////
// Component access: (), []

double & Vector::operator[](const unsigned index) {
  /* Check if index is > 3.
  Vectors only have 3 components (with indicies 0, 1, 2) */
  assert(index < 3);

  return (*this).Ar[index];
} // double & Vector::operator[](const unsigned index) {

double & Vector::operator()(const unsigned index) { return (*this)[index]; }



double Vector::operator[](const unsigned index) const {
  /* Check if index is > 3.
  Vectors only have 3 components (with indicies 0, 1, 2) */
  assert(index < 3);

  return (*this).Ar[index];
} // double Vector::operator()(const unsigned index) const {

double Vector::operator()(const unsigned index) const { return (*this)[index]; }





////////////////////////////////////////////////////////////////////////////////
// Boolean operators

bool Vector::operator==(const Vector & V_In) const {
  /* Check that the distance between components is <= Vector Epsilon.
  If it is then return false. */
  for(unsigned i = 0; i < 3; i++) {
    double d = (*this).Ar[i] - V_In.Ar[i];
    if( d < -Vector_Epsilon || d > Vector_Epsilon) { return false; }
  } // for(unsigned i = 0; i < 3; i++) {

  // If all components are sufficiently close, return true.
  return true;
} // bool Vector::operator==(const Vector & V_In) const {



bool Vector::operator!=(const Vector & V_In) const {
  // Return negation of (*this) == V_In.
  return !((*this) == V_In);
} // bool Vector::operator!=(const Vector & V_In) const {





////////////////////////////////////////////////////////////////////////////////
// Other methods

void Vector::Print(void) const {
  printf("< %9.2e, %9.2e, %9.2e>\n", (*this).Ar[0], (*this).Ar[1], (*this).Ar[2]);
} // void Print(void) const {



double Vector::Magnitude(void) const {
  #ifdef OPERATION_COUNT
    OP_Count::V_Magnitude++;                       // Increment operator count
  #endif

  return sqrt((*this).Ar[0]*(*this).Ar[0] +
              (*this).Ar[1]*(*this).Ar[1] +
              (*this).Ar[2]*(*this).Ar[2]);
} // double Vector::Magnitude(void) const {



double Vector::Max_Component(void) const {
  double Max = (*this).Ar[0];
  if((*this).Ar[1] > Max) { Max = (*this).Ar[1]; }
  if((*this).Ar[2] > Max) { Max = (*this).Ar[2]; }

  return Max;
} // double Vetor::Max_Component(void) const {



const double* Vector::Get_Ar(void) const {
  /* This function returns the address of the vector's internal array.

  This can be used to bypass the operator access methods and, thereby, improve
  runtime. However, it is extremely risky (no checks at all). Only use this if
  you know what you're doing. */
  return (*this).Ar;
} // const double* Vector::Get_Ar(void) const {





////////////////////////////////////////////////////////////////////////////////
// Functions of a single vector

Vector operator*(double c, const Vector & V_In) { return V_In*c; }

double Max_Component(const Vector & V_In) { return V_In.Max_Component(); }

double Magnitude(const Vector & V_In)  { return V_In.Magnitude(); }

void Print(const Vector & V_In) { V_In.Print(); }






////////////////////////////////////////////////////////////////////////////////
// Functions of multiple vectors

double Dot_Product(const Vector & V1, const Vector & V2) {
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
} // double Dot_Product(const Vector & V1, const Vector & V2) {



Vector Cross_Product(const Vector & V1, const Vector & V2) {
  Vector V1_x_V2;

  /* We assign the elements of the cross product V1 x V2 to the corresponding
  elements of V1_x_V2. */
  V1_x_V2[0] = V1[1]*V2[2] - V1[2]*V2[1];
  V1_x_V2[1] = V1[2]*V2[0] - V1[0]*V2[2];
  V1_x_V2[2] = V1[0]*V2[1] - V1[1]*V2[0];

  return V1_x_V2;
} // Vector Cross_Product(const Vector & V1, const Vector & V2) {



Tensor Dyadic_Product(const Vector & V1,const Vector & V2) {
  Tensor S;

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

  // Old loop. (works better than fully unrolled when using O2 optimization)
  for(int i = 0; i < 3; i++) {
    S[3*i + 0] = V1[i]*V2[0];
    S[3*i + 1] = V1[i]*V2[1];
    S[3*i + 2] = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++)

  #ifdef OPERATION_COUNT
    OP_Count::Dyadic_Product++;                    // Increment operator count
  #endif

  return S;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {g
