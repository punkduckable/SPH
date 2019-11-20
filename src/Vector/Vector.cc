#if !defined(VECTOR_SOURCE)
#define VECTOR_SOURCE

#include "Vector.h"

/* In this file, I define methods for Vector objects. Theese methods are
designed to make vector objects work just like mathematical vectors. In general,
I use [] to access components of the vectors. However, () and [] are defined in
the same way/give the same result (V[1]= V(1)).
*/

////////////////////////////////////////////////////////////////////////////////
// Constructors, Destructor

Vector::Vector(void) {
  //OP_Count::V_Default_Constructor++;             // Increment operator count (See SPH Diagnostics)
  //printf("Vector default constructor\n");
} // Vector::Vector(void) {



Vector::Vector(const double v0, const double v1, const double v2) {
  // Initialize components of vector using supplied components
  V[0] = v0;
  V[1] = v1;
  V[2] = v2;

  //OP_Count::V_Component_Constructor++;           // Increment operator count (See SPH Diagnostics)
  //printf("Vector component constructor\n");
} // Vector::Vector(const double v0, const double v1, const double v2) {



Vector::Vector(const Vector & V_In) {
  // Initialize components of vector using supplied components
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];

  //OP_Count::V_Copy_Constructor++;                // Increment operator count (See SPH Diagnostics)
  //printf("Vector copy constructor\n");
} // Vector::Vector(const Vector & V_In) {



Vector::~Vector(void) {
  //printf("Vector destroyed\n");
}



////////////////////////////////////////////////////////////////////////////////
// Vector equality

Vector & Vector::operator=(const double V_In[3]) {
  // Assign components of vector to V_In array
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];

  //OP_Count::V_Equality++;                        // Increment operator count (See SPH Diagnostics)

  // Return this Vector
  return *this;
} // Vector & Vector::operator=(const double V_In[3]) {



Vector & Vector::operator=(const Vector & V_In) {
  // Assign components of V using V_In.
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];

  //OP_Count::V_Equality++;                        // Increment operator count (See SPH Diagnostics)

  // Return this vector
  return *this;
} // Vector & Vector::operator=(const Vector V_In) {



////////////////////////////////////////////////////////////////////////////////
// Simple arithmetic operators

Vector Vector::operator+(const Vector & V_In) const {
  // Declare a sum Vector. This will be used to store the sum
  Vector Sum;

  /*  Add the components of the two vectors, store the results in the Sum Vector.

      Notice that I choose to write out the three components here rather than
      using a for loop. I did this because a for loop incurs overhead.
      This overhead comes from declaring the increment variable, running the
      test condition on each iteration and then performing the update to the
      iteration variable. Since there are only three components, it makes more
      sense to write out the component updates explicitly rather than using a
      for loop since this method will have less overhead/will run faster.

      I choose to use a similar optimization on other vector functions
  */
  Sum[0] = V[0] + V_In[0];
  Sum[1] = V[1] + V_In[1];
  Sum[2] = V[2] + V_In[2];

  //OP_Count::V_V_Addition++;                      // Increment operator count (See SPH Diagnostics)

  return Sum;
} // Vector Vector::operator+(const Vector V) const {



Vector Vector::operator-(const Vector & V_In) const{
  // Declare a Diff vector. This will be used to store the difference.
  Vector Diff;

  Diff[0] = V[0] - V_In[0];
  Diff[1] = V[1] - V_In[1];
  Diff[2] = V[2] - V_In[2];

  //OP_Count::V_V_Subtraction++;                   // Increment operator count (See SPH Diagnostics)

  return Diff;
} // Vector Vector::operator-(const Vector V) const {



Vector Vector::operator*(const double c) const {
  // Declare product vector
  Vector Prod;

  // Scale components of Prod by c.
  Prod[0] = V[0]*c;
  Prod[1] = V[1]*c;
  Prod[2] = V[2]*c;

  //OP_Count::V_S_Multiplication++;                // Increment operator count (See SPH Diagnostics)

  return Prod;
} // Vector Vector::operator*(const double c) const {



Vector Vector::operator/(const double c) const {
  // Declare product vector
  Vector Quotient;

  // Check for divide by zero
  if(c == 0) {
    printf("Vector /: Divide by zero error!\n");
    return Quotient;
  }

  // Divide components of V by c
  Quotient[0] = V[0]/c;
  Quotient[1] = V[1]/c;
  Quotient[2] = V[2]/c;

  //OP_Count::V_S_Division++;                      // Increment operator count (See SPH Diagnostics)

  return Quotient;
} // Vector Vector::operator/(const double c) const {



////////////////////////////////////////////////////////////////////////////////
// Compound arithmetic operators

Vector & Vector::operator+=(const Vector & V_In) {
  V[0] = V[0] + V_In[0];
  V[1] = V[1] + V_In[1];
  V[2] = V[2] + V_In[2];

  //OP_Count::Compound_V_V_Addition++;             // Increment operator count (See SPH Diagnostics)

  // Return this vector
  return *this;
} // Vector & Vector:operator+=(const Vector & V_In) {



Vector & Vector::operator+=(const double V_In[3]) {
  V[0] = V[0] + V_In[0];
  V[1] = V[1] + V_In[1];
  V[2] = V[2] + V_In[2];

  //OP_Count::Compound_V_V_Addition++;             // Increment operator count (See SPH Diagnostics)

  // Return this vector
  return *this;
} // Vector & Vector::operator+=(const double V_In[3]) {



Vector & Vector::operator-=(const Vector & V_In) {
  V[0] = V[0] - V_In[0];
  V[1] = V[1] - V_In[1];
  V[2] = V[2] - V_In[2];

  //OP_Count::Compound_V_V_Subtraction++;        // Increment operator count (See SPH Diagnostics)

  // Return this vector
  return *this;
} // Vector & Vector::operator-=(const double V_In[3]) {



Vector & Vector::operator*=(const double c) {
  // Scale the components of V by c.
  V[0] = V[0]*c;
  V[1] = V[1]*c;
  V[2] = V[2]*c;

  //OP_Count::Compound_V_S_Multiplication++;       // Increment operator count (See SPH Diagnostics)

  // Return this vector (now scalled by c)
  return *this;
} // Vector & Vector::operator*=(const double c) {



////////////////////////////////////////////////////////////////////////////////
// Component access: (), []

double & Vector::operator()(const uByte index) {
  /* Check if index is > 3. Note that this function only accepts an unsigned
  integer input. Thus, there is no possibility of negative numbers. Therefore,
  we only need to check that the index is < 3. */
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double & Vector::operator()(const uByte index) {



double Vector::operator()(const uByte index) const {
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double Vector::operator()(const uByte index) const {



////////////////////////////////////////////////////////////////////////////////
// Other methods

void Vector::Print(void) const {
  printf("< %9.2e, %9.2e, %9.2e>\n",V[0], V[1], V[2]);
} // void Print(void) const {



double Vector::Magnitude(void) const {
  return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
} // double Magnitude(void) const {



double Max_Component(const Vector & V_In) {
  if(V_In[0] > V_In[1] && V_In[0] > V_In[2])
    return V_In[0];
  else if(V_In[1] > V_In[2])
    return V_In[1];
  else
    return V_In[2];
} // double Max_Component(const Vector V_In) {



////////////////////////////////////////////////////////////////////////////////
// Friend methods

Vector operator*(double c, const Vector & V_In) {
  return V_In*c;
} //Vector operator*(double c, const Vector & V_In) {



double Magnitude(const Vector & V_In)  {
  //OP_Count::V_Magnitude++;                       // Increment operator count (See SPH Diagnostics)
  return V_In.Magnitude();
} // double Magnitude(const Vector V_In)  {



double Vector_Dot_Product(const Vector & V1, const Vector & V2) {
  //OP_Count::V_Dot_Product++;                     // Increment operator count (See SPH Diagnostics)
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}



void Print(const Vector & V_In) {
  V_In.Print();
} // void Print(const Vector & V_In) {



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

  // Old loop. (works better than 9 statements with O2 optimization)
  for(int i = 0; i < 3; i++) {
    S[3*i + 0] = V1[i]*V2[0];
    S[3*i + 1] = V1[i]*V2[1];
    S[3*i + 2] = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++)

  //OP_Count::Dyadic_Product++;                    // Increment operator count (See SPH Diagnostics)

  return S;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {g

#endif
