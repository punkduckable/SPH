#include <stdio.h>
#include <math.h>

using namespace std;

// Type definitions
typedef signed char Byte;

// Classes
class Vector {
  public:
    double V[3];                    // The three components of the Vector

    Vector(void);                   // Default Constructor
    Vector(double V_in[3]);         // Vector based constructor
    Vector(double v0,
           double v1,
           double v2);              // Component based constructor

    Vector operator+(Vector V_In);     // Addition overload (so we can add vectors)
    Vector operator-(Vector V_In);     // Subtraction overload (so we can subtract vectors)
    Vector operator=(double V_In[3]);
    double& operator()(Byte index);
    double& operator[](Byte index);
    void Print(void);               // Print vector components
}; // class Vector {

class Tensor {
  public:
    double T[9];

    Tensor(void);                                 // Default constructor
    Tensor(double t11, double t12, double t13,
           double t21, double t22, double t23,
           double t31, double t32, double t33);   // Component constructor

    Tensor operator+(Tensor S);
    Tensor operator*(Tensor S);
    Vector operator*(Vector V);
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

Vector::Vector(double V_in[3]) {
  // Initialize components of vector using supplied array
  V[0] = V_in[0];
  V[1] = V_in[1];
  V[2] = V_in[2];
} // Vector::Vector(double v_in[3]) {

Vector::Vector(double v0, double v1, double v2) {
  // Initialize components of vector using supplied components
  V[0] = v0;
  V[1] = v1;
  V[2] = v2;
} // Vector::Vector(double v0, double v1, double v2) {

Vector Vector::operator+(Vector V_In) {
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
} // Vector Vector::operator+(Vector V) {

Vector Vector::operator-(Vector V_In) {
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
} // Vector Vector::operator-(Vector V) {

Vector Vector::operator=(double V_In[3]) {
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];

  return *this;
}

double& Vector::operator()(Byte index) {
  return V[index];
} // double& Vector::operator()(Byte index) {

double& Vector::operator[](Byte index) {
  return V[index];
} // double& Vector::operaotr[](Byte index) {

void Vector::Print(void) {
  printf("< %4.2f, %4.2f, %4.2f>\n",V[0], V[1], V[2]);
} // void Print(void) {

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

Tensor Tensor::operator+(Tensor S) {
  // Declare some vector to store the sum
  Tensor Sum;

  /* Compute the 9 elements of the sum.
     I choose to use a for loop here, despite the extra overhead, because
     writing out 9 individual lines was too hard to maintain (lots of room for
     typos, hard to change how the function operates). I am hopeful that the
     compuler willl unroll the loop for me.
  */
  for(int i = 0; i < 9; i++){
    Sum.T[i] = T[i] + S.T[i];
  } // for(int i = 0; i < 9; i++){

  return Sum;
} // Tensor Tensor::operator+(Tensor S) {

Tensor Tensor::operator*(Tensor S) {
  // Declare product tensor
  Tensor Prod;

  // Use nested for loops to calculate the 9 elements of the sum. This uses
  // standard matrix matrix multiplication.
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        Prod.T[3*i+j] += T[3*i+k]*S.T[3*k+j];
      } // for(int k = 0; k < 3; k++) {
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  return Prod;
} // Tensor Tensor::operator*(Tensor S) {

Vector Tensor::operator*(Vector V_In) {
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
} // Vector Tensor::operator*(Vector V) {

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
  Tensor T2(1,4,7,
            2,5,8,
            3,6,9);

  Tensor T3 = T1*T2;

  T3.Print();


  return 0;
} // int main() {
