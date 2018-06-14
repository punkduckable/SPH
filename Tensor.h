#if !defined(_TENSOR_HEADER)
#define _TENSOR_HEADER

/*  This is the header file for my Tensor class. Most of the methods for this
class define operations with Tensors, Vectors, and Scalars. The intent of this
is to make Tensor objects act just like second order tensors. As such, things
like vector Tensor Tensor products, Tensor Vector products, Scalar
multiplication, Tensor equality, Tensor Addition, Tensor Inverse, etc... have
been defined. Each method has a comment explaining what that method does. T, T1,
and T2 denote 2nd order Tensor objects, T[9] denotes a tensor stored as a 9 element
array, V denotes a Vector object, c denotes a scalar. */

class Tensor {
  private:
    double T[9];                                  // Holds the 9 components of the Tensor

  public:
    Tensor(void);                                 // Default constructor
    Tensor(double t11, double t12, double t13,
           double t21, double t22, double t23,
           double t31, double t32, double t33);   // Component constructor
    Tensor(const Tensor & T_In);                  // Copy Constructor

    ~Tensor(void);                                // Destructor

    Tensor operator+(const Tensor & T_In) const;  // Tensor addition (defines T1 + T2)
    Tensor operator-(const Tensor & T_In) const;  // Tensor subtraction (defines T1 - T2);
    Tensor operator*(const Tensor & T_In) const;  // Tensor-Tensor multiplication (Defines T1*T2)
    Vector operator*(const Vector & V_In) const;  // Tensor-Vector multiplication (defines T*V)
    Tensor operator*(const double c) const;       // Scalar multiplication (defines T*c)
    Tensor operator/(const double c) const;       // Scalar division (defines T/c)

    Tensor & operator+=(const Tensor & T_In);     // Compound Tensor addition (defines T1 += T2)
    Tensor & operator+=(const double T_In[9]);    // Compound tensor addition (defines T1 += T2[9])
    Tensor & operator*=(const double c);          // Compound scalar multiplication (defines T *= c)
    Tensor & operator*=(const Tensor & T_In);     // Compound Tensor-Tensor multiplication (defines T1 *= T2)

    Tensor & operator=(const double T_In[9]);     // Tensor equality (Defines T1 = T2[9])
    Tensor & operator=(const Tensor & T_In);      // Tensor equaltiy (defines T1 = T2)

    double & operator()(const uByte row,
                      const uByte col);           // Component access (defines T(i,j))
    double operator()(const uByte row,
                      const uByte col) const;

    Tensor Inverse(void) const;                   // Tensor Inverse. Returns T^(-1)
    double Determinant(void) const;               // Tensor Determinant. Returns Det(T)
    Tensor Transpose(void) const;                 // Tensor Transpose. Returns T^T (Tranpose of T)
    void Print(void) const;                       // Print tensor components

    friend Tensor operator*(double c,
                            const Tensor & T_In); // Scalar multiplication (defines c*T)
    friend Tensor Inverse(const Tensor & T_In);
    friend double Determinant(const Tensor & T_In);
    friend Tensor Transpose(const Tensor & T_In);
    friend double Tensor_Dot_Product(const Tensor & T1, const Tensor & T2);
    friend void Print(const Tensor & T_In);
}; // class Tensor {

// Tensor functions that don't belong in the tensor class
Tensor Inverse(const Tensor & T_In);
double Determinant(const Tensor & T_In);
Tensor Transpose(const Tensor & T_In);
double Tensor_Dot_Product(const Tensor & T1, const Tensor & T2);
void Print(const Tensor & T_In);

#endif
