#if !defined(TENSOR_HEADER)
#define TENSOR_HEADER

/*  File Description:

This is the header file for my Tensor class. Most of the methods for this
class define operations with Tensors, Vectors, and Scalars. The intent of this
is to make Tensor objects act just like second order tensors. As such, things
like vector Tensor Tensor products, Tensor Vector products, Scalar
multiplication, Tensor equality, Tensor Addition, Tensor Inverse, etc... have
been defined. Each method has a comment explaining what that method does. T, T1,
and T2 denote 2nd order Tensor objects, S[9] denotes a tensor stored as a 9 element
array, V denotes a Vector object, c denotes a scalar. */

#include "Vector/Vector.h"
#include "Classes.h"
#include "Errors.h"

// Used for ^ method (to define ^T and stuff).
const unsigned T = 2;

class Tensor {
  private:
    double Ar[9];                                  // Holds the 9 components of the Tensor

  public:
    // Constructors, destructor
    Tensor(void);                                 // Default constructor
    Tensor(double t11, double t12, double t13,
           double t21, double t22, double t23,
           double t31, double t32, double t33);   // Component constructor
    Tensor(const Tensor & T_In);                  // Copy Constructor

    ~Tensor(void);                                // Destructor



    // Tensor equality
    Tensor & operator=(const Tensor & T_In);      // Tensor equaltiy (defines T1 = T2)



    // Simple arithmetic operators
    Tensor operator+(const Tensor & T_In) const;  // Tensor addition (defines T1 + T2)
    Tensor operator-(const Tensor & T_In) const;  // Tensor subtraction (defines T1 - T2);
    Tensor operator*(const Tensor & T_In) const;  // Tensor-Tensor multiplication (Defines T1*T2)
    Vector operator*(const Vector & V_In) const;  // Tensor-Vector multiplication (defines T*V)
    Tensor operator*(const double c) const;       // Scalar multiplication (defines T*c)
    Tensor operator/(const double c) const;       // Scalar division (defines T/c)



    // Compound arithmetic operators
    Tensor & operator+=(const Tensor & T_In);     // Compound Tensor addition (defines T1 += T2)
    Tensor & operator-=(const Tensor & T_In);     // Compound tensor subtraction (defines T1 -= T2)
    Tensor & operator*=(const Tensor & T_In);     // Compound Tensor-Tensor multiplication (defines T1 *= T2)



    // Component access (public)
    double & operator()(const unsigned row,
                       const unsigned col);       // Write to a component of the tensor (defines T(i,j) = ....)
    double & operator[](const unsigned index);    // Write to a component of the tensor
    double operator()(const unsigned row,
                      const unsigned col) const;  // Read a component of the tensor (defines ... = T(i,j)) (runs checks, safer)
    double operator[](const unsigned index) const;// Read a component of the tensor (no checks, faster)



    // Inverse methods
    Tensor operator^(const int exp);
    Tensor Inverse(void) const;                   // Tensor Inverse. Returns T^(-1)



    // Other methods
    double Determinant(void) const;               // Tensor Determinant. Returns Det(T)
    Tensor Transpose(void) const;                 // Tensor Transpose. Returns T^T  (Tranpose of T)
    const Vector Eigenvalues(const char Mode = 'A') const; // Returns eigenvalues of the Tensor
    void Print(void) const;                       // Print tensor components
}; // class Tensor {

// Tensor methods that don't belong in the Tensor class.
Tensor operator*(double c,
                 const Tensor & T_In); // Scalar multiplication (defines c*T)
Tensor Inverse(const Tensor & T_In);
double Determinant(const Tensor & T_In);
Tensor Transpose(const Tensor & T_In);
const Vector Eigenvalues(const Tensor & S_In,
                         const char Mode = 'A');
void Print(const Tensor & T_In);



double Tensor_Dot_Product(const Tensor & T1,
                                 const Tensor & T2);
#endif
