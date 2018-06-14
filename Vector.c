#if !defined(_VECTOR_SOURCE)
#define _VECTOR_SOURCE

/* In this file, I define methods for Vector objects. Theese methods are
designed to make vector objects work just like mathematical vectors. In general,
I use [] to access components of the vectors. However, () and [] are defined in
the same way/give the same result (V[1]= V(1)).
*/

////////////////////////////////////////////////////////////////////////////////
// Constructors, Destructor

Vector::Vector(void) {
  //printf("Vector default constructor\n");

  // Initialize components of vector to zero (no input supplied, assume zero)
  V[0] = 0;
  V[1] = 0;
  V[2] = 0;
} // Vector::Vector(void) {

Vector::Vector(const double v0, const double v1, const double v2) {
  //printf("Vector component constructor\n");

  // Initialize components of vector using supplied components
  V[0] = v0;
  V[1] = v1;
  V[2] = v2;
} // Vector::Vector(const double v0, const double v1, const double v2) {

Vector::Vector(const Vector & V_In) {
  //printf("Vector copy constructor\n");

  // Initialize components of vector using supplied components
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];
} // Vector::Vector(const Vector & V_In) {

Vector::~Vector(void) {
  //printf("Vector destroyed\n");
}



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

  return Sum;
} // Vector Vector::operator+(const Vector V) const {

Vector Vector::operator-(const Vector & V_In) const{
  // Declare a Diff vector. This will be used to store the difference.
  Vector Diff;

  Diff[0] = V[0] - V_In[0];
  Diff[1] = V[1] - V_In[1];
  Diff[2] = V[2] - V_In[2];

  return Diff;
} // Vector Vector::operator-(const Vector V) const {

Vector Vector::operator*(const double c) const {
  // Declare product vector
  Vector Prod;

  // Scale components of Prod by c.
  Prod[0] = V[0]*c;
  Prod[1] = V[1]*c;
  Prod[2] = V[2]*c;

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

  return Quotient;
} // Vector Vector::operator/(const double c) const {



////////////////////////////////////////////////////////////////////////////////
// Compound arithmetic operators

Vector & Vector::operator+=(const Vector & V_In) {
  V[0] = V[0] + V_In[0];
  V[1] = V[1] + V_In[1];
  V[2] = V[2] + V_In[2];

  // Return this vector
  return *this;
} // Vector & Vector:operator+=(const Vector & V_In) {

Vector & Vector::operator+=(const double V_In[3]) {
  V[0] = V[0] + V_In[0];
  V[1] = V[1] + V_In[1];
  V[2] = V[2] + V_In[2];

  // Return this vector
  return *this;
} // Vector & Vector::operator+=(const double V_In[3]) {

Vector & Vector::operator*=(const double c) {
  // Scale the components of V by c.
  V[0] = V[0]*c;
  V[1] = V[1]*c;
  V[2] = V[2]*c;

  // Return this vector (now scalled by c)
  return *this;
} // Vector & Vector::operator*=(const double c) {



////////////////////////////////////////////////////////////////////////////////
// Vector equality

Vector & Vector::operator=(const double V_In[3]) {
  // Assign components of vector to V_In array
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];

  // Return this Vector
  return *this;
} // Vector & Vector::operator=(const double V_In[3]) {

Vector & Vector::operator=(const Vector & V_In) {
  // Assign components of V using V_In.
  V[0] = V_In[0];
  V[1] = V_In[1];
  V[2] = V_In[2];

  // Return this vector
  return *this;
} // Vector & Vector::operator=(const Vector V_In) {



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

double & Vector::operator[](const uByte index) {
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double & Vector::operator[](const uByte index) {

double Vector::operator[](const uByte index) const {
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double Vector::operator[](const uByte index) const {



////////////////////////////////////////////////////////////////////////////////
// Other methods

void Vector::Print(void) const {
  printf("< %4.2f, %4.2f, %4.2f>\n",V[0], V[1], V[2]);
} // void Print(void) const {

double Vector::Magnitude(void) const {
  return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
} // double Magnitude(void) const {



////////////////////////////////////////////////////////////////////////////////
// Friend methods

Vector operator*(double c, const Vector & V_In) {
  return V_In*c;
} //Vector operator*(double c, const Vector & V_In) {

double Magnitude(const Vector & V_In)  {
  return V_In.Magnitude();
} // double Magnitude(const Vector V_In)  {

double Vector_Dot_Product(const Vector & V1, const Vector & V2) {
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

void Print(const Vector & V_In) {
  V_In.Print();
} // void Print(const Vector & V_In) {

#endif
