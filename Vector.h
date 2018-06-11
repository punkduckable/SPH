#if !defined(_VECTOR_HEADER)
#define _VECTOR_HEADER

/*  This is the header file for my vector class. Most of the methods for this
class define operations with vectors and scalars. The intent of this is to make
vector objects act just like mathematical vectors. As such, things like vector
addition, scalar multiplication, vector equality, etc.. have been defined.
Each method has a comment explaining what that method does. V, V1, and V2 denote
Vector objects, V[3] denotes a vector stored as an array, and c denotes a scalar
constant. */

class Vector {
  private:
    double V[3];                                     // Holds the three components of the Vector

  public:
    Vector(void);                                    // Default Constructor
    Vector(const double v0,
           const double v1,
           const double v2);                         // Component based constructor
    Vector(const Vector & V_In);                     // Copy constructor

    ~Vector(void);                                   // Destructor

    Vector operator+(const Vector & V_In) const;     // Vector addition (defines V1 + V2)
    Vector operator-(const Vector & V_In) const;     // Vector Subtraction (defines V1 - V2)
    Vector operator*(const double c) const;          // Scalar multiplication (defines V*c)
    Vector operator/(const double c) const;          // Scalar Divide (defines V/c).

    Vector & operator+=(const Vector & V_In);        // Compound Vector addition (Defines V1 += V2)
    Vector & operator+=(const double V_In[3]);       // Compound Vector addition with a 3 element array (Defines V1 += V2[3])
    Vector & operator*=(const double c);             // Compound Scalar multiplication (defines V *= c)

    Vector & operator=(const double V_In[3]);        // Vector equality (defines V1 = V2[3])
    Vector & operator=(const Vector & V_In);         // Vector equalitty (defines V1 = V2)

    double & operator()(const uByte index);          // () component access (defines V(n))
    double operator()(const uByte index) const;
    double & operator[](const uByte index);          // [] component access (defines V[n])
    double operator[](const uByte index) const;

    void Print(void) const;                          // Print vector components
    double Magnitude(void) const;

    friend Vector operator*(double c, const Vector & V_In);  // Scalar multiplication (defines c*V)
    friend double Magnitude(const Vector V_In);
    friend double Vector_Dot_Product(const Vector V1, const Vector V2);
}; // class Vector {

#endif
