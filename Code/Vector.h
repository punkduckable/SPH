#if !defined(VECTOR_HEADER)
#define VECTOR_HEADER

/*  This is the header file for my vector class. Most of the methods for this
class define operations with vectors and scalars. The intent of this is to make
vector objects act just like mathematical vectors. As such, things like vector
addition, scalar multiplication, vector equality, etc.. have been defined.
Each method has a comment explaining what that method does. V, V1, and V2 denote
Vector objects, V[3] denotes a vector stored as an array, and c denotes a scalar
constant. */

class Vector {
  friend class Tensor;

  private:
    double V[3];                                           // Holds the three components of the Vector

    double & operator[](const uByte index) {               // Write to a component of vector (no checks, faster)
      return V[index];
    };                                                     // Read a componnet of the vector (no checks, faster)
    double operator[](const uByte index) const {
      return V[index];
    };

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
    Vector & operator+=(const double V_In[3]);             // Compound Vector addition with a 3 element array (Defines V1 += V2[3])
    Vector & operator-=(const Vector & V_In);
    Vector & operator*=(const double c);                   // Compound Scalar multiplication (defines V *= c)

    double & operator()(const uByte index);                // Write to a component of the vector (runs checks, safer)
    double operator()(const uByte index) const;            // Read a component of the vector (runs checks, safer)

    // Magnitude method
    double Magnitude(void) const;                          // Returns magnitude of vector

    // Other methods

    // Friends
    friend const Vector Eigenvalues(const Tensor & S_In,
                                    const char Mode);
    friend double Max_Component(const Vector & V_In);      // Returns maximum component of vector.
    friend Vector operator*(double c,                      // Scalar multiplication (defines c*V)
                            const Vector & V_In);
    friend double Magnitude(const Vector & V_In);
    friend double Vector_Dot_Product(const Vector & V1,
                                     const Vector & V2);
    friend void Print(const Vector & V_In);
    friend Tensor Dyadic_Product(const Vector & V1,
                                 const Vector & V2);

    // Temporary friends (should remove)
    friend void Particle_Tests(void);
    friend void Particle_Debugger::Export_Pariticle_Forces(const unsigned int Num_Particles,
                                                           const Particle * Particles);
    friend void VTK_File::Export_Pariticle_Positions(const unsigned int Num_Particles,
                                                     const Particle * Particles);

    // Printing functions
    void Print(void) const;                                // Print vector components
}; // class Vector {

// Vector functions that don't belong in the vector class
double Magnitue(const Vector & V_In);
double Vector_Dot_Product(const Vector & V1,
                          const Vector & V2);
void Print(const Vector & V_In);

#endif
