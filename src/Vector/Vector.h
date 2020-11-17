#if !defined(VECTOR_HEADER)
#define VECTOR_HEADER
//#define VECTOR_MONITOR

/*  File description:

This is the header file for my vector class. Most of the methods for this
class define operations with vectors and scalars. The intent of this is to make
vector objects act just like mathematical vectors. As such, things like vector
addition, scalar multiplication, vector equality, etc.. have been defined.
Each method has a comment explaining what that method does. V, V1, and V2 denote
Vector objects, V[3] denotes a vector stored as an array, and c denotes a scalar
constant. */

#include "Classes.h"

/* Smallest allowed difference for equality
Getting exact equality (even when you should) is difficult with floating point
operations (because of roundoffs). Thus, we say that two floating point numbers
are "equal" if they are sufficiently close to one another. Vector_Epsilon is
the largest allowed distance between two "equal" floats.

This is used by the == nd != operators */
const double Vector_Epsilon = 1e-16;

class Vector {
  private:
    double V[3];                                           // Holds the three components of the Vector

  public:
    // Constructors, destructor
    Vector(void);                                          // Default Constructor
    Vector(const double v0,
           const double v1,
           const double v2);                               // Component based constructor
    Vector(const Vector & V_In);                           // Copy constructor

    ~Vector(void);                                         // Destructor


    // Vector equality
    Vector & operator=(const double V_In[3]);              // Vector equality (defines V1 = V2[3])
    Vector & operator=(const Vector & V_In);               // Vector equalitty (defines V1 = V2)


    // Operator overloading
    Vector operator+(const Vector & V_In) const;           // Vector addition (defines V1 + V2)
    Vector operator-(const Vector & V_In) const;           // Vector Subtraction (defines V1 - V2)
    Vector operator*(const double c) const;                // Scalar multiplication (defines V*c)
    Vector operator/(const double c) const;                // Scalar Divide (defines V/c).

    Vector & operator+=(const Vector & V_In);              // Compound Vector addition (Defines V1 += V2)
    Vector & operator-=(const Vector & V_In);              // Compounted Vector subtraction (Defines V1 -= V2)
    Vector & operator*=(const double c);                   // Compound Scalar multiplication (Defines V *= c)

    double & operator()(const unsigned index);                // Write to a component of the vector
    double & operator[](const unsigned index);                // Write to a component of the vector
    double operator()(const unsigned index) const;            // Read a component of the vector
    double operator[](const unsigned index) const;            // Read a component of the vector

    bool operator==(const Vector & V_In) const;
    bool operator!=(const Vector & V_In) const;


    // Other methods
    double Magnitude(void) const;                          // Returns magnitude of vector
    double Max_Component(void) const;                      // Returns maximum component of vector.
    void Print(void) const;                                // Print vector components
}; // class Vector {


// Functions of a single vector
Vector operator*(double c,                      // Scalar multiplication (defines c*V)
                 const Vector & V_In);
double Max_Component(const Vector & V_In);      // Returns maximum component of vector.
double Magnitude(const Vector & V_In);
void Print(const Vector & V_In);


// Functions of multiple vectors
double Dot_Product(const Vector & V1,
                   const Vector & V2);
Vector Cross_Product(const Vector & V1,
                     const Vector & V2);
Tensor Dyadic_Product(const Vector & V1,
                      const Vector & V2);


#endif
