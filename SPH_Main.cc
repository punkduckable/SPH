#include <stdio.h>
#include <math.h>

using namespace std;

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Classes
class Vector {
  public:
    double V[3];                    // Stores the three components of the Vector

    Vector(void);                                  // Default Constructor
    Vector(const double V_in[3]);                  // Vector based constructor
    Vector(const double v0,
           const double v1,
           const double v2);                       // Component based constructor

    Vector operator+(const Vector V_In) const;     // Addition overload (so we can add vectors)
    Vector operator-(const Vector V_In) const;     // Subtraction overload (so we can subtract vectors)
    Vector operator=(const double V_In[3]);        // Initialize a vector to an array
    Vector operator=(const Vector V_In);           // Initialize a vector to another vector!
    double& operator()(const uByte index);
    double operator()(const uByte index) const;
    double& operator[](const uByte index);
    double operator[](const uByte index) const;

    void Print(void) const;               // Print vector components
}; // class Vector {

class Tensor {
  public:
    double T[9];

    Tensor(void);                                 // Default constructor
    Tensor(double t11, double t12, double t13,
           double t21, double t22, double t23,
           double t31, double t32, double t33);   // Component constructor

    Tensor operator+(const Tensor S_In) const;
    Tensor operator*(const Tensor S_In) const;
    Vector operator*(const Vector V_In) const;
    Tensor operator=(const double S_In[9]);
    Tensor operator=(const Tensor S_In);
    double& operator()(const uByte row, const uByte col);
    double operator()(const uByte row, const uByte col) const;
    void Print(void);
}; // class Tensor {

// Prototypes
Tensor Dyadic_Product(Vector V1, Vector V2);

////////////////////////////////////////////////////////////////////////////////
// Vector method definitions
Vector::Vector(void) {
  // Initialize components of vector to zero (no input supplied, assume zero)
  V[0] = 0;
  V[1] = 0;
  V[2] = 0;
} // Vector::Vector(void) {

Vector::Vector(const double V_in[3]) {
  // Initialize components of vector using supplied array
  V[0] = V_in[0];
  V[1] = V_in[1];
  V[2] = V_in[2];
} // Vector::Vector(const double v_in[3]) {

Vector::Vector(const double v0, const double v1, const double v2) {
  // Initialize components of vector using supplied components
  V[0] = v0;
  V[1] = v1;
  V[2] = v2;
} // Vector::Vector(const double v0, const double v1, const double v2) {

Vector Vector::operator+(const Vector V_In) const {
  // Declare a sum Vector. This will be used to store the sum
  Vector Sum;

  /*  Add the components of the two vectors, store the results in the Sum Vector.

      Notice that I choose to write out the three components here rather than
      using a for loop. I did this because a for loop incures overhead.
      This overhead comes from declaring the incremenet variable, running the
      test condition on each iteration and then performing the update to the
      iteration variable. Since there are only three components, it makes more
      sense to write out the component updates explicitly rather than using a
      for loop since this method will have less overhead/will run faster.
  */
  Sum.V[0] = V[0] + V_In.V[0];
  Sum.V[1] = V[1] + V_In.V[1];
  Sum.V[2] = V[2] + V_In.V[2];

  return Sum;
} // Vector Vector::operator+(const Vector V) const {

Vector Vector::operator-(const Vector V_In) const{
  // Declare a Diff vector. This will be used to store the difference.
  Vector Diff;

  /* Subtract componnets of the two vectors, store in the Diff Vector.

     Notice that I choose to write out the three components here rather than
     using a for loop. I did this because a for loop incures overhead.
     This overhead comes from declaring the incremenet variable, running the
     test condition on each iteration and then performing the update to the
     iteration variable. Since there are only three components, it makes more
     sense to write out the component updates explicitly rather than using a
     for loop since this method will have less overhead/will run faster.
  */
  Diff.V[0] = V[0] - V_In.V[0];
  Diff.V[1] = V[1] - V_In.V[1];
  Diff.V[2] = V[2] - V_In.V[2];

  return Diff;
} // Vector Vector::operator-(const Vector V) const {

Vector Vector::operator=(const double V_In[3]) {
  // Assign components of vector to V_In array
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];

  // Return this Vector
  return *this;
}

Vector Vector::operator=(const Vector V_In) {
  // Assign components of V using V_In.
  V[0] = V_In.V[0];
  V[0] = V_In.V[0];
  V[0] = V_In.V[0];

  // Return this vector
  return *this;
} // Vector Vector::operator=(const Vector V_In) {

double& Vector::operator()(const uByte index) {
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double& Vector::operator()(const uByte index) {

  double Vector::operator()(const uByte index) const {
    if(index >= 3)
      printf("Index out of bounds");

    return V[index];
  } // double Vector::operator()(const uByte index) const {

double& Vector::operator[](const uByte index) {
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double& Vector::operaotr[](const uByte index) {

double Vector::operator[](const uByte index) const {
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double Vector::operaotr[](const uByte index) const {

void Vector::Print(void) const {
  printf("< %4.2f, %4.2f, %4.2f>\n",V[0], V[1], V[2]);
} // void Print(void) const {

////////////////////////////////////////////////////////////////////////////////
// Tensor method definitions
Tensor::Tensor(void) {
  // Defualt constructor: Set 9 components to zero.
  T[0] = T[1] = T[2] =
  T[3] = T[4] = T[5] =
  T[6] = T[7] = T[8] = 0;
} // Tensor::Tensor(void) {

Tensor::Tensor(double t11, double t12, double t13,
               double t21, double t22, double t23,
               double t31, double t32, double t33) {

    // Set the 9 individual components of T using inputs.
    T[0] = t11;
    T[1] = t12;
    T[2] = t13;
    T[3] = t21;
    T[4] = t22;
    T[5] = t23;
    T[6] = t31;
    T[7] = t32;
    T[8] = t33;
} // Tensor:Tensor(double t11,.... double t33) {

Tensor Tensor::operator+(const Tensor S_In) const {
  // Declare some vector to store the sum
  Tensor Sum;

  /* Compute the 9 elements of the sum.
     I choose to use a for loop here, despite the extra overhead, because
     writing out 9 individual lines was too hard to maintain (lots of room for
     typos, hard to change how the function operates). I am hopeful that the
     compuler willl unroll the loop for me.
  */
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++) {
      Sum(i,j) = T[3*i + j] + S_In(i,j);
    }
  } // for(int i = 0; i < 9; i++){

  return Sum;
} // Tensor Tensor::operator+(const Tensor S_In) const {

Tensor Tensor::operator*(const Tensor S_In) const{
  // Declare product tensor
  Tensor Prod;

  // Use nested for loops to calculate the 9 elements of the sum. This uses
  // standard matrix matrix multiplication.
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        Prod(i,j) += T[3*i + k]*S_In(k,j);
      } // for(int k = 0; k < 3; k++) {
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  return Prod;
} // Tensor Tensor::operator*(const Tensor S_In) const{

Vector Tensor::operator*(const Vector V_In) const {
  // Declare product vector (matrix vector product is a vector)
  Vector Prod;

  // Calculate componnets of Prod using a for loop. This uses standard matrix
  // vector multiplication
  for(int i = 0; i < 3; i++) {
    Prod.V[i] = T[3*i]*V_In.V[0] +
                T[3*i+1]*V_In.V[1] +
                T[3*i+2]*V_In.V[2];
  } //   for(int i = 0; i < 3; i++) {

  return Prod;
} // Vector Tensor::operator*(const Vector V_In) const {

Tensor Tensor::operator=(const double V_In[9]) {
  for(int i = 0; i < 9; i++) {
    for(int j = 0; j < 9; j++) {
      T[i] = V_In[i];
    }
  } //   for(int i = 0; i < 9; i++) {

  return *this;
} // Tensor Tensor::operator=(const double V_In[9]) {

Tensor Tensor::operator=(const Tensor S_In) {
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i + j] = S_In(i,j);
    }
  }

  return *this;
} // Tensor Tensor::operator=(const Tensor S_In) {

double& Tensor::operator()(const uByte row, const uByte col) {
  if(row >= 3 || col >=3)
      printf("Index out of bounds\n");

  return T[3*row + col];
} // double& Tensor::operator()(const uByte row, const uByte col) {

double Tensor::operator()( const uByte row, const uByte col) const {
  if(row >= 3 || col >=3)
      printf("Index out of bounds\n");

  return T[3*row + col];
} // double Tensor::operator()(const uByte row, const ubyte col) const {

void Tensor::Print(void) {
  for(int i = 0; i < 3; i++) {
    printf("| %6.2f %6.2f %6.2f |\n",T[i*3], T[i*3+1], T[i*3+2]);
  }
} // void Tensor::Print(void) {

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
