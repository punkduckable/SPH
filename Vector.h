#if !defined(_VECTOR_HEADER)
#define _VECTOR_HEADER

class Vector {
  public:
    double V[3];                    // Stores the three components of the Vector

    Vector(void);                                  // Default Constructor
    Vector(const double v0,
           const double v1,
           const double v2);                       // Component based constructor

    Vector operator+(const Vector V_In) const;     // Addition overload (so we can add vectors)
    Vector operator-(const Vector V_In) const;     // Subtraction overload (so we can subtract vectors)

    Vector operator+=(const Vector V_In);
    Vector operator+=(const double V_In[3]);

    Vector operator=(const double V_In[3]);        // Initialize a vector to an array
    Vector operator=(const Vector V_In);           // Initialize a vector to another vector!

    double& operator()(const uByte index);
    double operator()(const uByte index) const;
    double& operator[](const uByte index);
    double operator[](const uByte index) const;

    void Print(void) const;               // Print vector components
}; // class Vector {

#endif
