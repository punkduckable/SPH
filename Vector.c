#if !defined(_VECTOR_SOURCE)
#define _VECTOR_SOURCE

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
      using a for loop. I did this because a for loop incurs overhead.
      This overhead comes from declaring the increment variable, running the
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

  /* Subtract components of the two vectors, store in the Diff Vector.

     Notice that I choose to write out the three components here rather than
     using a for loop. I did this because a for loop incurs overhead.
     This overhead comes from declaring the increment variable, running the
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
} // double& Vector::operator[](const uByte index) {

double Vector::operator[](const uByte index) const {
  if(index >= 3)
    printf("Index out of bounds");

  return V[index];
} // double Vector::operator[](const uByte index) const {

void Vector::Print(void) const {
  printf("< %4.2f, %4.2f, %4.2f>\n",V[0], V[1], V[2]);
} // void Print(void) const {

#endif
