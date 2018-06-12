#if !defined(_TENSOR_SOURCE)
#define _TENSOR_SOURCE

/* In this file, I define methods for tensor objects. It should be noted that
tensor data is stored in a 9 element array rather than a 3x3 matrix. This is
done because single dimension arrays tend to perform better and are, in my
personal opinion, easier to work with. As a result, successive elements in a
column are 3 elements away from each other. Successive elements in a row are one
element away from each other. Therefore, T[3*i+j] is the same as T(i,j). */

////////////////////////////////////////////////////////////////////////////////
// Constuctors and destructor

Tensor::Tensor(void) {
  printf("Tensor default constructor\n");

  // Default constructor: Set 9 components to zero.
  T[0] = T[1] = T[2] =
  T[3] = T[4] = T[5] =
  T[6] = T[7] = T[8] = 0;
} // Tensor::Tensor(void) {

Tensor::Tensor(double t11, double t12, double t13,
               double t21, double t22, double t23,
               double t31, double t32, double t33) {
    printf("Tensor component constrctor\n");

    // Set the 9 individual components of T using inputs.
    T[0] = t11; T[1] = t12; T[2] = t13; // Row 1
    T[3] = t21; T[4] = t22; T[5] = t23; // Row 2
    T[6] = t31; T[7] = t32; T[8] = t33; // Row 3
} // Tensor:Tensor(double t11,.... double t33) {

Tensor::Tensor(const Tensor & T_In) {
  printf("Tensor Copy contructor\n");

  for(int i = 1; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i + j] = T_In(i,j);
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 1; i < 3; i++) {
} // Tensor::Tensor(const Tensor & T_In) {

Tensor::~Tensor(void) {
  printf("Tensor destroyed\n");
}


////////////////////////////////////////////////////////////////////////////////
// Simple arithmetic operators

Tensor Tensor::operator+(const Tensor & T_In) const {
  // Declare some Tensor to store the sum
  Tensor Sum;

  /* Compute the 9 elements of the sum.
     I choose to use a for loop here, despite the extra overhead, because
     writing out 9 individual lines was too hard to maintain (lots of room for
     typos, hard to change how the function operates). I am hopeful that the
     compiler will unroll the loop for me.
  */
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++) {
      Sum(i,j) = T[3*i + j] + T_In(i,j);
    }
  } // for(int i = 0; i < 9; i++){

  return Sum;
} // Tensor Tensor::operator+(const Tensor & T_In) const {

Tensor Tensor::operator-(const Tensor & T_In) const {
  // Declare some Tensor to store the difference
  Tensor Diff;

  /* Compute the 9 elements of the sum.
     I choose to use a for loop here, despite the extra overhead, because
     writing out 9 individual lines was too hard to maintain (lots of room for
     typos, hard to change how the function operates). I am hopeful that the
     compiler will unroll the loop for me.
  */
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++) {
      Diff(i,j) = T[3*i + j] - T_In(i,j);
    }
  } // for(int i = 0; i < 9; i++){

  return Diff;
} // Tensor Tensor::operator-(const Tensor & T_In) const {

Tensor Tensor::operator*(const Tensor & T_In) const{
  // Declare product tensor
  Tensor Prod;

  /* Calcualate Tensor Tensor product using nested for loops. Store each element
   Of the product in the Tensor 'Prod'.

   Let A and B be tensors. The (i,j) element of AB is the dot product of the
   ith row of A and the jth column of B. This nromally means that we'd need
   three nested loops to calculate a Tensor Tensor product; one for the rows,
   one for the columns, and one for the the dot product for each element.
   However, since the Tensors we are working with will always be of dimension
   3x3, the dot product has a predictable form. Therefore, rather than include
   a third nested loop (which would add more overhead), we write out the 3 terms
   of the dot product explicitly.

   It should be noted that this small change will significantly reduce overhead.
   The reason is that the 3rd nested loop would run 9 times (once for each
   element) this means that the 3 iterations of the 3rd nested loop would run
   9 times each, leading to a total of 27 equivalent iterations. By contrast,
   the j loop only runs 3 times leading to 9 equivalent iterations and the i
   loop only runs once for 3 equivelent iterations. Therefore, the innermost
   loop is by far the most expensive. Removing it should have an appreciable
   impact on performance. */
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
        Prod(i,j) = T[3*i + 0]*T_In(0,j) +
                    T[3*i + 1]*T_In(1,j) +
                    T[3*i + 2]*T_In(2,j);
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  return Prod;
} // Tensor Tensor::operator*(const Tensor & T_In) const{

Vector Tensor::operator*(const Vector & V_In) const {
  // Declare product vector (matrix vector product is a vector)
  Vector Prod;

  // Calculate components of vector tensor product. Store results Prod.
  for(int i = 0; i < 3; i++) {
    Prod[i] =   T[3*i + 0]*V_In[0] +
                T[3*i + 1]*V_In[1] +
                T[3*i + 2]*V_In[2];
  } //   for(int i = 0; i < 3; i++) {

  return Prod;
} // Vector Tensor::operator*(const Vector & V_In) const {

Tensor Tensor::operator*(const double c) const {
  // We will store results in Prod
  Tensor Prod;

  // Scale T by c, store in Prod.
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      Prod(i,j) = T[3*i + j]*c;
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  return Prod;
} // Tensor Tensor::operator*(const double c) const {

Tensor Tensor::operator/(const double c) const {
  Tensor Quotient;

  // Check that quotient is non-zero
  if(c == 0)
    printf("You tried dividing a tensor by zero!!!\n");

  // Divide the components of T by c, store the results in Quotient.
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
        Quotient(i,j) = T[3*i + j]/c;
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  return Quotient;
} //Tensor Tensor::operator/(const double c) const {



////////////////////////////////////////////////////////////////////////////////
// Compound Arithmetic operations

Tensor & Tensor::operator+=(const Tensor & T_In) {
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i+j] = T[3*i+ j] + T_In(i,j);
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  // Return this tensor (element wise summation is done)
  return *this;
} // Tensor & Tensor::operator+=(const Tensor & T_In) {

Tensor & Tensor::operator+=(const double T_In[9]) {
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i+j] = T[3*i+ j] + T_In[3*i+j];
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  // Return this tensor (element wise summation is done)
  return *this;
} // Tensor & Tensor::operator+=(const Tensor T_In) {

Tensor & Tensor::operator*=(const double c) {
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++)  {
      T[3*i + j] = T[3*i + j]*c;
    } // for(int j = 0; j < 3; j++)  {
  } // for(int i = 0; i < 3; i++) {

  // Return this tensor (all elements have been scaled by c)
  return *this;
} // Tensor & Tensor::operator*=(const double c) {

Tensor & Tensor::operator*=(const Tensor & T_In) {
  /* When we multiply two tensors together, we do so element by element. The
  issue with this is that once an element has been updated, there will be no
  way to access the origional value. We will need that origional value to
  calculate other components of the product. Therefore, we need a new tensor to
  temporarly hold onto the product. */
  Tensor Prod;

  /* Calcualate Tensor Tensor product using nested for loops. Store each element
   Of the product in the Tensor 'Prod'.

   Let A and B be tensors. The (i,j) element of AB is the dot product of the
   ith row of A and the jth column of B. This nromally means that we'd need
   three nested loops to calculate a Tensor Tensor product; one for the rows,
   one for the columns, and one for the the dot product for each element.
   However, since the Tensors we are working with will always be of dimension
   3x3, the dot product has a predictable form. Therefore, rather than include
   a third nested loop (which would add more overhead), we write out the 3 terms
   of the dot product explicitly.

   It should be noted that this small change will significantly reduce overhead.
   The reason is that the 3rd nested loop would run 9 times (once for each
   element) this means that the 3 iterations of the 3rd nested loop would run
   9 times each, leading to a total of 27 equivalent iterations. By contrast,
   the j loop only runs 3 times leading to 9 equivalent iterations and the i
   loop only runs once for 3 equivelent iterations. Therefore, the innermost
   loop is by far the most expensive. Removing it should have an appreciable
   impact on performance. */
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      Prod(i,j) = T[3*i + 0]*T_In(0,j) +
                  T[3*i + 1]*T_In(1,j) +
                  T[3*i + 2]*T_In(2,j);
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  // Copy Prod to T (in this Tensor)
  for(int i = 0; i < 3; i++) {
    T[3*i + 0] = Prod(i,0);
    T[3*i + 1] = Prod(i,1);
    T[3*i + 2] = Prod(i,2);
  } //   for(int i = 0; i < 3; i++) {

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
// Tensor equality

Tensor & Tensor::operator=(const double T_In[9]) {
  for(int i = 0; i < 9; i++) {
    for(int j = 0; j < 9; j++) {
      T[i] = T_In[i];
    }
  } //   for(int i = 0; i < 9; i++) {

  // return this tensor (element wise equality is done)
  return *this;
} // Tensor & Tensor::operator=(const double T_In[9]) {

Tensor & Tensor::operator=(const Tensor & T_In) {
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i + j] = T_In(i,j);
    }
  }

  // return this tensor (element wise equality is done)
  return *this;
} // Tensor & Tensor::operator=(const Tensor & T_In) {



////////////////////////////////////////////////////////////////////////////////
// Tensor element access: ()

double & Tensor::operator()(const uByte row, const uByte col) {
  /* Check if desired rows or col is > 3. Note that this function only accepts
  unsigned integer inputs. Thus, there is no possibility of row or col being
  negative.  Therefore we only need to check that row and col are < 3 */
  if(row >= 3 || col >=3)
      printf("Index out of bounds\n");

  // Return the specified element (treat it like matirx notation)
  return T[3*row + col];
} // double & Tensor::operator()(const uByte row, const uByte col) {

double Tensor::operator()( const uByte row, const uByte col) const {
  if(row >= 3 || col >=3)
      printf("Index out of bounds\n");

  // Return the specified element (treat it like matirx notation)
  return T[3*row + col];
} // double Tensor::operator()(const uByte row, const uByte col) const {



////////////////////////////////////////////////////////////////////////////////
// Other Tensor methods

Tensor Tensor::Inverse(void) const{
  /* This method calculates and returns the inverse of this tensor. We do this
  using the fomrula for the inverse of a 3x3 matrix: Let
      | a b c |
  A = | d e f |
      | g h i |

  Then, by solving for A^(-1) in A*A^(-1) = I by working through some awful
  algebra we get,
                      | -hf + ei    ch - bi   -ce + bf |
  A^(-1) = (1/Det(A))*|  fg - di   -cg + ai    cd - af |
                      | -eg + dh    bg - ah   -bd + ae |
  Where Det(A) = aei - ahf - bdi + bgf + cdh - cge.
  */

  // First, let's define the return tensor
  Tensor T_Inv;

  // We will store the determinant in a double
  double Det_T = (T[0*3 + 0]*(T[1*3 + 1]*T[2*3 + 2]
                             -T[2*3 + 1]*T[1*3 + 2])
                 +T[0*3 + 1]*(T[2*3 + 0]*T[1*3 + 2]
                             -T[1*3 + 0]*T[2*3 + 2])
                 +T[0*3 + 2]*(T[1*3 + 0]*T[2*3 + 1]
                             -T[2*3 + 0]*T[1*3 + 1]));

  /* If the determinant is non-zero, then the matrix is singular and there is no
  Inverse. Thus, we should hault computation if the determinant is zero.
  */

  if(Det_T == 0) {
    printf("This tensor is singular. No Inverse exists!\n");
    return T_Inv;
  }

  T_Inv(0,0) = -T[2*3 + 1]*T[1*3 + 2]
               +T[1*3 + 1]*T[2*3 + 2];      // -hf + ei

  T_Inv(0,1) =  T[0*3 + 2]*T[2*3 + 1]
               -T[0*3 + 1]*T[2*3 + 2];      // ch - bi

  T_Inv(0,2) = -T[0*3 + 2]*T[1*3 + 1]
               +T[0*3 + 1]*T[1*3 + 2];      // -ce + bf

  T_Inv(1,0) =  T[1*3 + 2]*T[2*3 + 0]
               -T[1*3 + 0]*T[2*3 + 2];      // fg - di

  T_Inv(1,1) = -T[0*3 + 2]*T[2*3 + 0]
               +T[0*3 + 0]*T[2*3 + 2];      // -cg + ai

  T_Inv(1,2) =  T[0*3 + 2]*T[1*3 + 0]
               -T[0*0 + 0]*T[1*3 + 2];      // cd - af

  T_Inv(2,0) = -T[1*3 + 1]*T[2*3 + 0]
               +T[1*3 + 0]*T[2*3 + 1];      // -eg + dh

  T_Inv(2,1) =  T[0*3 + 1]*T[2*3 + 0]
               -T[0*3 + 0]*T[2*3 + 1];      // bg - ah

  T_Inv(2,2) = -T[0*3 + 1]*T[1*3 + 0]
               +T[0*3 + 0]*T[1*3 + 1];      // -bd + ae

  T_Inv = (1/Det_T)*T_Inv;

  return T_Inv;
} // Tensor Inverse(void) const{

double Tensor::Determinant(void) const {
  double Det_T = (T[0*3 + 0]*(T[1*3 + 1]*T[2*3 + 2]
                             -T[2*3 + 1]*T[1*3 + 2])
                 +T[0*3 + 1]*(T[2*3 + 0]*T[1*3 + 2]
                             -T[1*3 + 0]*T[2*3 + 2])
                 +T[0*3 + 2]*(T[1*3 + 0]*T[2*3 + 1]
                             -T[2*3 + 0]*T[1*3 + 1]));

  return Det_T;
} // double Tensor::Determinant(void) const {

Tensor Tensor::Transpose(void) const {
  Tensor T_Transpose;                         // Will store Transpose

  // Set components of the transpose using the rule A(i,j) = A^T(i,j) = A^(j,i).
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T_Transpose(i,j) = T[3*j + i];
    }
  }

  return T_Transpose;
} // Tensor Tensor::Transpose(void) const {

void Tensor::Print(void) const {
  for(int i = 0; i < 3; i++) {
    printf("| %8.4f %8.4f %8.4f |\n",T[i*3], T[i*3+1], T[i*3+2]);
  }
} // void Tensor::Print(void) const {



////////////////////////////////////////////////////////////////////////////////
// Friend functions

Tensor operator*(const double c, const Tensor & T_In) {
  // Used to store the resulting product.
  Tensor T_Scaled;

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T_Scaled(i,j) = T_In(i,j)*c;
    } // for(int i = 0; i < 3; i++) {
  } // for(int j = 0; j < 3; j++) {

  // Return the scaled tensor.
  return T_Scaled;
} // Tensor operator*(const double c, const Tensor & T_In) {

Tensor Inverse(const Tensor & T_In) {
  return T_In.Inverse();
} // Tensor Inverse(const Tensor & T_In) {

double Determinant(const Tensor & T_In) {
  return T_In.Determinant();
} // double Determinant(const Tensor & T_In) {

Tensor Transpose(const Tensor & T_In) {
  return T_In.Transpose();
} // Tensor Transpose(const Tensor & T_In) {

double Tensor_Dot_Product(const Tensor & T1, const Tensor & T2) {
  // This returns T1:T2, the tensor dot product of T1 and T2
  double dot_prod = 0;

  for(int i = 0; i < 3; i++) {
    dot_prod += T1(i,0)*T2(i,0) +
                T1(i,1)*T2(i,1) +
                T1(i,2)*T2(i,2);
  }

  return dot_prod;
} // double Tensor_Dot_Product(const Tensor & T1, const Tensor & T2) {

void Print(const Tensor & T_In) {
  T_In.Print();
} // void Print(const Tensor & T_In) {
  
#endif
