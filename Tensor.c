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
  //printf("Tensor default constructor\n");

  // Default constructor: Set 9 components to zero.
  T[0] = T[1] = T[2] =
  T[3] = T[4] = T[5] =
  T[6] = T[7] = T[8] = 0;
} // Tensor::Tensor(void) {

Tensor::Tensor(double t11, double t12, double t13,
               double t21, double t22, double t23,
               double t31, double t32, double t33) {
    //printf("Tensor component constrctor\n");

    // Set the 9 individual components of T using inputs.
    T[0] = t11; T[1] = t12; T[2] = t13; // Row 1
    T[3] = t21; T[4] = t22; T[5] = t23; // Row 2
    T[6] = t31; T[7] = t32; T[8] = t33; // Row 3
} // Tensor:Tensor(double t11,.... double t33) {

Tensor::Tensor(const Tensor & T_In) {
  //printf("Tensor Copy contructor\n");

  /* Here we create one tensor using another. This is done element-by-element.
  Normally, such an operation would require a double nested loop (one for the
  rows and one for the columns).  However, loops have overhead. To eliminate
  this overhead, I wrote what would have been the 9 iterations of this double
  loop as 9 statemenets.

  To make this a little  more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */
  T[3*0 + 0] = T_In(0,0);                        // i = 0, j = 0
  T[3*0 + 1] = T_In(0,1);                        // i = 0, j = 1
  T[3*0 + 2] = T_In(0,2);                        // i = 0, j = 2

  T[3*1 + 0] = T_In(1,0);                        // i = 1, j = 0
  T[3*1 + 1] = T_In(1,1);                        // i = 1, j = 1
  T[3*1 + 2] = T_In(1,2);                        // i = 1, j = 2

  T[3*2 + 0] = T_In(2,0);                        // i = 2, j = 0
  T[3*2 + 1] = T_In(2,1);                        // i = 2, j = 1
  T[3*2 + 2] = T_In(2,2);                        // i = 2, j = 2

  /* Old double nested loop
  for(int i = 1; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i + j] = T_In(i,j);
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 1; i < 3; i++) {
  */
} // Tensor::Tensor(const Tensor & T_In) {

Tensor::~Tensor(void) {
  //printf("Tensor destroyed\n");
}


////////////////////////////////////////////////////////////////////////////////
// Simple arithmetic operators

Tensor Tensor::operator+(const Tensor & T_In) const {
  // Declare some Tensor to store the sum
  Tensor Sum;

  /* Here we comput (*this) + T_In. This is done component-by-componnet.
  Normally, such calculations would require a double nested loop (one for rows)
  and one for columns). However, loops have overhead, and we want this code to
  run fast. Thus, I choose to write out the 9 loop iterations as statements
  rather than in a double loop.

  I included comments that specify which loop iteration a given statement would
  correspond to (with i as the row index, j as the column index  */
  Sum(0,0) = T[3*0 + 0] + T_In(0,0);             // i = 0, j = 0
  Sum(0,1) = T[3*0 + 1] + T_In(0,1);             // i = 0, j = 1
  Sum(0,2) = T[3*0 + 2] + T_In(0,2);             // i = 0, j = 2

  Sum(1,0) = T[3*1 + 0] + T_In(1,0);             // i = 1, j = 0
  Sum(1,1) = T[3*1 + 1] + T_In(1,1);             // i = 1, j = 1
  Sum(1,2) = T[3*1 + 2] + T_In(1,2);             // i = 1, j = 2

  Sum(2,0) = T[3*2 + 0] + T_In(2,0);             // i = 2, j = 0
  Sum(2,1) = T[3*2 + 1] + T_In(2,1);             // i = 2, j = 1
  Sum(2,2) = T[3*2 + 2] + T_In(2,2);             // i = 2, j = 2

  /* Old, double nested loop.
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++) {
      Sum(i,j) = T[3*i + j] + T_In(i,j);
    }
  } // for(int i = 0; i < 9; i++){
  */

  return Sum;
} // Tensor Tensor::operator+(const Tensor & T_In) const {

Tensor Tensor::operator-(const Tensor & T_In) const {
  // Declare some Tensor to store the difference
  Tensor Diff;

  /* Here we comput (*this) - T_In. This is done component-by-component
  Normally, such calculations would require a double nested loop (one for rows)
  and one for columns). However, loops have overhead, and we want this code to
  run fast. Thus, I choose to write out the 9 loop iterations as statements
  rather than in a double loop.

  I included comments that specify which loop iteration a given statement would
  correspond to (with i as the row index, j as the column index  */

  Diff(0,0) = T[3*0 + 0] - T_In(0,0);             // i = 0, j = 0
  Diff(0,1) = T[3*0 + 1] - T_In(0,1);             // i = 0, j = 1
  Diff(0,2) = T[3*0 + 2] - T_In(0,2);             // i = 0, j = 2

  Diff(1,0) = T[3*1 + 0] - T_In(1,0);             // i = 1, j = 0
  Diff(1,1) = T[3*1 + 1] - T_In(1,1);             // i = 1, j = 1
  Diff(1,2) = T[3*1 + 2] - T_In(1,2);             // i = 1, j = 2

  Diff(2,0) = T[3*2 + 0] - T_In(2,0);             // i = 2, j = 0
  Diff(2,1) = T[3*2 + 1] - T_In(2,1);             // i = 2, j = 1
  Diff(2,2) = T[3*2 + 2] - T_In(2,2);             // i = 2, j = 2

  /* Old double nested loop
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++) {
      Diff(i,j) = T[3*i + j] - T_In(i,j);
    }
  } // for(int i = 0; i < 9; i++){
  */

  return Diff;
} // Tensor Tensor::operator-(const Tensor & T_In) const {

Tensor Tensor::operator*(const Tensor & T_In) const{
  // Declare product tensor
  Tensor Prod;

  /* Calcualate Tensor Tensor product using nested for loops. Store each element
  Of the product in the Tensor 'Prod'.

  Let A and B be tensors. The (i,j) element of AB is the dot product of the
  ith row of A and the jth column of B. In other words:
                            AB(i,j) = Sum(p = 0 to p = 2) A(i,p)*B(p,j)
  Naturally, this would lead to three nested loops.

  We will call the loop that changes the rows the 'i loop', the loop that
  changes the columns the 'j loop' and the loop that performs the dot product
  the 'k loop'

  The order in which we perform these loops does not change the resulting tensor
  product. However, the order will completely change performance. I have stored
  my tensors in Row-Major order  ( A(i,j), A(i,j+1) occur in adjacnet memory
  locations).

  To optimize cache usage, we want the loops that change rows to be on the
  outside. Notice that for both AB and A, the i index changes rows. We want this
  to be the outer most loop. The p index changes columns on the B matrix, so we
  want this loop to come next. Finally, the j index never changes the rows, so
  this index should go on the inside.

  Thus, to improve performance, I choose to order my loops in this way.

  However, that's not the end of the story; loops cary overhead. We want
  exceptional performnace, this means minimizing unecessairry overhead.
  I found that 'unrolling' my loops - literally writing out the 27 loop
  interations in the order that they would occur in (so that we still optimize
  cache line usage) can improve perofmrnce by 20%+. This makes the code a little
  difficult to understand, however. To make the code more readible, I have
  included a comment on each statement that explains which loop iteration that
  statement would have corresponded to (which value of the i,j,p indicies) */

  /* i = 0 (1st row of C) */
  Prod[3*0 + 0] += T[3*0 + 0]*T_In[3*0 + 0];     // i = 0, j = 0, p = 0
  Prod[3*0 + 1] += T[3*0 + 0]*T_In[3*0 + 1];     // i = 0, j = 1, p = 0
  Prod[3*0 + 2] += T[3*0 + 0]*T_In[3*0 + 2];     // i = 0, j = 2, p = 0

  Prod[3*0 + 0] += T[3*0 + 1]*T_In[3*1 + 0];     // i = 0, j = 0, p = 1
  Prod[3*0 + 1] += T[3*0 + 1]*T_In[3*1 + 1];     // i = 0, j = 1, p = 1
  Prod[3*0 + 2] += T[3*0 + 1]*T_In[3*1 + 2];     // i = 0, j = 2, p = 1

  Prod[3*0 + 0] += T[3*0 + 2]*T_In[3*2 + 0];     // i = 0, j = 0, p = 2
  Prod[3*0 + 1] += T[3*0 + 2]*T_In[3*2 + 1];     // i = 0, j = 1, p = 2
  Prod[3*0 + 2] += T[3*0 + 2]*T_In[3*2 + 2];     // i = 0, j = 2, p = 2


  /* i = 1 (2nd row of C) */
  Prod[3*1 + 0] += T[3*1 + 0]*T_In[3*0 + 0];     // i = 1, j = 0, p = 0
  Prod[3*1 + 1] += T[3*1 + 0]*T_In[3*0 + 1];     // i = 1, j = 1, p = 0
  Prod[3*1 + 2] += T[3*1 + 0]*T_In[3*0 + 2];     // i = 1, j = 2, p = 0

  Prod[3*1 + 0] += T[3*1 + 1]*T_In[3*1 + 0];     // i = 1, j = 0, p = 1
  Prod[3*1 + 1] += T[3*1 + 1]*T_In[3*1 + 1];     // i = 1, j = 1, p = 1
  Prod[3*1 + 2] += T[3*1 + 1]*T_In[3*1 + 2];     // i = 1, j = 2, p = 1

  Prod[3*1 + 0] += T[3*1 + 2]*T_In[3*2 + 0];     // i = 1, j = 0, p = 2
  Prod[3*1 + 1] += T[3*1 + 2]*T_In[3*2 + 1];     // i = 1, j = 1, p = 2
  Prod[3*1 + 2] += T[3*1 + 2]*T_In[3*2 + 2];     // i = 1, j = 2, p = 2


  /* i = 2 (3rd row) */
  Prod[3*2 + 0] += T[3*2 + 0]*T_In[3*0 + 0];     // i = 2, j = 0, p = 0
  Prod[3*2 + 1] += T[3*2 + 0]*T_In[3*0 + 1];     // i = 2, j = 1, p = 0
  Prod[3*2 + 2] += T[3*2 + 0]*T_In[3*0 + 2];     // i = 2, j = 2, p = 0

  Prod[3*2 + 0] += T[3*2 + 1]*T_In[3*1 + 0];     // i = 2, j = 0, p = 1
  Prod[3*2 + 1] += T[3*2 + 1]*T_In[3*1 + 1];     // i = 2, j = 1, p = 1
  Prod[3*2 + 2] += T[3*2 + 1]*T_In[3*1 + 2];     // i = 2, j = 2, p = 1

  Prod[3*2 + 0] += T[3*2 + 2]*T_In[3*2 + 0];     // i = 2, j = 0, p = 2
  Prod[3*2 + 1] += T[3*2 + 2]*T_In[3*2 + 1];     // i = 2, j = 1, p = 2
  Prod[3*2 + 2] += T[3*2 + 2]*T_In[3*2 + 2];     // i = 2, j = 2, p = 2

  /* Old, loop based Tensor-Tensor product
  for(int i = 0; i < 3; i++) {    // Row loops
    for(int p = 0; p < 3; p++) {    // Dot prod loops
      Prod(i,0) += T[3*i + j]*T_In(j,0);
      Prod(i,1) += T[3*i + j]*T_In(j,1);
      Prod(i,2) += T[3*i + j]*T_In(j,2);
    } // for(int p = 0; p < 3; p++) {
  } // for(int i = 0; i < 3; i++) {
  */

  return Prod;
} // Tensor Tensor::operator*(const Tensor & T_In) const{

Vector Tensor::operator*(const Vector & V_In) const {
  // Declare product vector (matrix vector product is a vector)
  Vector Prod;

  /* Here I compute the product Prod = (*this) * V_In. Normally, this would
  require a loop to go through the components of Prod (and possible a second
  to compute what goes in each componnet). However, Loops have overhead. To
  eliminate this overhead, I have written what would have been the loop
  iterations as a series of statements.*/

  Prod[0] =   T[3*0 + 0]*V_In[0] +
              T[3*0 + 1]*V_In[1] +
              T[3*0 + 2]*V_In[2];

  Prod[1] =   T[3*1 + 0]*V_In[0] +
              T[3*1 + 1]*V_In[1] +
              T[3*1 + 2]*V_In[2];

  Prod[2] =   T[3*2 + 0]*V_In[0] +
              T[3*2 + 1]*V_In[1] +
              T[3*2 + 2]*V_In[2];

  /* Old loop
  for(int i = 0; i < 3; i++) {
    Prod[i] =   T[3*i + 0]*V_In[0] +
                T[3*i + 1]*V_In[1] +
                T[3*i + 2]*V_In[2];
  } //   for(int i = 0; i < 3; i++) {
  */

  return Prod;
} // Vector Tensor::operator*(const Vector & V_In) const {

Tensor Tensor::operator*(const double c) const {
  // We will store results in Prod
  Tensor Prod;

  /* Here we return the product (*this)*c. Normally, doing this would require
  cycling through the row's and col's of T in a double nested loop. However,
  loops have overhead. To eliminate this overhead, I wrote what would have been
  the 9 iterations of this double loop as 9 statemenets.

  To make this a little  more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */
  Prod(0,0) = T[3*0 + 0]*c;
  Prod(0,1) = T[3*0 + 1]*c;
  Prod(0,2) = T[3*0 + 2]*c;

  Prod(1,0) = T[3*1 + 0]*c;
  Prod(1,1) = T[3*1 + 1]*c;
  Prod(1,2) = T[3*1 + 2]*c;

  Prod(2,0) = T[3*2 + 0]*c;
  Prod(2,1) = T[3*2 + 1]*c;
  Prod(2,2) = T[3*2 + 2]*c;

  /* Old double nested loop
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      Prod(i,j) = T[3*i + j]*c;
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {
  */

  return Prod;
} // Tensor Tensor::operator*(const double c) const {

Tensor Tensor::operator/(const double c) const {
  Tensor Quotient;

  // Check that quotient is non-zero
  if(c == 0)
    printf("You tried dividing a tensor by zero!!!\n");

    /* Here we return the quotient (*this)/c. Normally, doing this would require
    cycling through the row's and col's of T in a double nested loop. However,
    loops have overhead. To eliminate this overhead, I wrote what would have been
    the 9 iterations of this double loop as 9 statemenets.

    To make this a little  more readible, I have included a comment with each
    statement that identifies which loop iteration that statement would have
    corresponded to (with i as the row index and j as the column index) */
    Quotient(0,0) = T[3*0 + 0]/c;
    Quotient(0,1) = T[3*0 + 1]/c;
    Quotient(0,2) = T[3*0 + 2]/c;

    Quotient(1,0) = T[3*1 + 0]/c;
    Quotient(1,1) = T[3*1 + 1]/c;
    Quotient(1,2) = T[3*1 + 2]/c;

    Quotient(2,0) = T[3*2 + 0]/c;
    Quotient(2,1) = T[3*2 + 1]/c;
    Quotient(2,2) = T[3*2 + 2]/c;

  /* Old double nested loop
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
        Quotient(i,j) = T[3*i + j]/c;
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {
  */

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

Tensor & Tensor::operator-=(const Tensor & T_In) {
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i+j] = T[3*i+ j] - T_In(i,j);
    } // for(int j = 0; j < 3; j++) {
  } // for(int i = 0; i < 3; i++) {

  // Return this tensor (element wise subtraction is done)
  return *this;
} // Tensor & Tensor::operator-=(const Tensor & T_In) {

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
  // Declare product tensor
  Tensor Prod;

  /* Calcualate Tensor Tensor product using nested for loops. Store each element
  Of the product in the Tensor 'Prod'.

  Let A and B be tensors. The (i,j) element of AB is the dot product of the
  ith row of A and the jth column of B. In other words:
                            AB(i,j) = Sum(p = 0 to p = 2) A(i,p)*B(p,j)
  Naturally, this would lead to three nested loops.

  We will call the loop that changes the rows the 'i loop', the loop that
  changes the columns the 'j loop' and the loop that performs the dot product
  the 'k loop'

  The order in which we perform these loops does not change the resulting tensor
  product. However, the order will completely change performance. I have stored
  my tensors in Row-Major order  ( A(i,j), A(i,j+1) occur in adjacnet memory
  locations).

  To optimize cache usage, we want the loops that change rows to be on the
  outside. Notice that for both AB and A, the i index changes rows. We want this
  to be the outer most loop. The p index changes columns on the B matrix, so we
  want this loop to come next. Finally, the j index never changes the rows, so
  this index should go on the inside.

  Thus, to improve performance, I choose to order my loops in this way.

  However, that's not the end of the story; loops cary overhead. We want
  exceptional performnace, this means minimizing unecessairry overhead.
  I found that 'unrolling' my loops - literally writing out the 27 loop
  interations in the order that they would occur in (so that we still optimize
  cache line usage) can improve perofmrnce by 20%+. This makes the code a little
  difficult to understand, however. To make the code more readible, I have
  included a comment on each statement that explains which loop iteration that
  statement would have corresponded to (which value of the i,j,p indicies) */

  /* i = 0 (1st row of C) */
  Prod[3*0 + 0] += T[3*0 + 0]*T_In[3*0 + 0];     // i = 0, j = 0, p = 0
  Prod[3*0 + 1] += T[3*0 + 0]*T_In[3*0 + 1];     // i = 0, j = 1, p = 0
  Prod[3*0 + 2] += T[3*0 + 0]*T_In[3*0 + 2];     // i = 0, j = 2, p = 0

  Prod[3*0 + 0] += T[3*0 + 1]*T_In[3*1 + 0];     // i = 0, j = 0, p = 1
  Prod[3*0 + 1] += T[3*0 + 1]*T_In[3*1 + 1];     // i = 0, j = 1, p = 1
  Prod[3*0 + 2] += T[3*0 + 1]*T_In[3*1 + 2];     // i = 0, j = 2, p = 1

  Prod[3*0 + 0] += T[3*0 + 2]*T_In[3*2 + 0];     // i = 0, j = 0, p = 2
  Prod[3*0 + 1] += T[3*0 + 2]*T_In[3*2 + 1];     // i = 0, j = 1, p = 2
  Prod[3*0 + 2] += T[3*0 + 2]*T_In[3*2 + 2];     // i = 0, j = 2, p = 2


  /* i = 1 (2nd row of C) */
  Prod[3*1 + 0] += T[3*1 + 0]*T_In[3*0 + 0];     // i = 1, j = 0, p = 0
  Prod[3*1 + 1] += T[3*1 + 0]*T_In[3*0 + 1];     // i = 1, j = 1, p = 0
  Prod[3*1 + 2] += T[3*1 + 0]*T_In[3*0 + 2];     // i = 1, j = 2, p = 0

  Prod[3*1 + 0] += T[3*1 + 1]*T_In[3*1 + 0];     // i = 1, j = 0, p = 1
  Prod[3*1 + 1] += T[3*1 + 1]*T_In[3*1 + 1];     // i = 1, j = 1, p = 1
  Prod[3*1 + 2] += T[3*1 + 1]*T_In[3*1 + 2];     // i = 1, j = 2, p = 1

  Prod[3*1 + 0] += T[3*1 + 2]*T_In[3*2 + 0];     // i = 1, j = 0, p = 2
  Prod[3*1 + 1] += T[3*1 + 2]*T_In[3*2 + 1];     // i = 1, j = 1, p = 2
  Prod[3*1 + 2] += T[3*1 + 2]*T_In[3*2 + 2];     // i = 1, j = 2, p = 2


  /* i = 2 (3rd row) */
  Prod[3*2 + 0] += T[3*2 + 0]*T_In[3*0 + 0];     // i = 2, j = 0, p = 0
  Prod[3*2 + 1] += T[3*2 + 0]*T_In[3*0 + 1];     // i = 2, j = 1, p = 0
  Prod[3*2 + 2] += T[3*2 + 0]*T_In[3*0 + 2];     // i = 2, j = 2, p = 0

  Prod[3*2 + 0] += T[3*2 + 1]*T_In[3*1 + 0];     // i = 2, j = 0, p = 1
  Prod[3*2 + 1] += T[3*2 + 1]*T_In[3*1 + 1];     // i = 2, j = 1, p = 1
  Prod[3*2 + 2] += T[3*2 + 1]*T_In[3*1 + 2];     // i = 2, j = 2, p = 1

  Prod[3*2 + 0] += T[3*2 + 2]*T_In[3*2 + 0];     // i = 2, j = 0, p = 2
  Prod[3*2 + 1] += T[3*2 + 2]*T_In[3*2 + 1];     // i = 2, j = 1, p = 2
  Prod[3*2 + 2] += T[3*2 + 2]*T_In[3*2 + 2];     // i = 2, j = 2, p = 2

  /*
  // Old triple nested loop.
  for(int i = 0; i < 3; i++) {    // Row loop
    for(int p = 0; p < 3; p++) {    // Dot prod loop
      for(int j = 0; j < 3; j++) {    // Col loop
        Prod(i,j) += T[3*i + p]*T_In(p,j);
      } // for(int j = 0; j < 3; j++) {
    } // for(int p = 0; p < 3; p++) {
  } // for(int i = 0; i < 3; i++) {
  */

  /* Copy Prod to T (in this Tensor). To do this normally would require two
  nested loops (one for the rows and one for the cols). However, as before, we
  want to minimize overhead. This means unrolling loops. Thus, rather than using
  a double loop to cycle through the 9 elements of our tensor, I have witten out
  the 9 iterations. To make the code a little more readible, I included the
  origional iteration numbers as comments. */
  T[3*0 + 0] = Prod[3*0 + 0];                    // i = 0, j = 0
  T[3*0 + 1] = Prod[3*0 + 1];                    // i = 0, j = 1
  T[3*0 + 2] = Prod[3*0 + 2];                    // i = 0, j = 1

  T[3*1 + 0] = Prod[3*1 + 0];                    // i = 1, j = 0
  T[3*1 + 1] = Prod[3*1 + 1];                    // i = 1, j = 1
  T[3*1 + 2] = Prod[3*1 + 2];                    // i = 1, j = 2

  T[3*2 + 0] = Prod[3*2 + 0];                    // i = 2, j = 0
  T[3*2 + 1] = Prod[3*2 + 1];                    // i = 2, j = 1
  T[3*2 + 2] = Prod[3*2 + 2];                    // i = 2, j = 2

  /* Old double loop
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i + j] = Prod(i,j);
    } // for(int j = 0; j < 3; j++) {
  } //   for(int i = 0; i < 3; i++) {
  */
  return *this;
} // Tensor & Tensor::operator*=(const Tensor & T_In) {

////////////////////////////////////////////////////////////////////////////////
// Tensor equality

Tensor & Tensor::operator=(const Tensor & T_In) {
  /* Here we set one tensor equal to another. This is done element-by-element.
  Normally, such an operation would require a double nested loop (one for the
  rows and one for the columns).  However, loops have overhead. To eliminate
  this overhead, I wrote what would have been the 9 iterations of this double
  loop as 9 statemenets.

  To make this a little  more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */
  T[3*0 + 0] = T_In[3*0 + 0];                    // i = 0, j = 0
  T[3*0 + 1] = T_In[3*0 + 1];                    // i = 0, j = 1
  T[3*0 + 2] = T_In[3*0 + 2];                    // i = 0, j = 1

  T[3*1 + 0] = T_In[3*1 + 0];                    // i = 1, j = 0
  T[3*1 + 1] = T_In[3*1 + 1];                    // i = 1, j = 1
  T[3*1 + 2] = T_In[3*1 + 2];                    // i = 1, j = 2

  T[3*2 + 0] = T_In[3*2 + 0];                    // i = 2, j = 0
  T[3*2 + 1] = T_In[3*2 + 1];                    // i = 2, j = 1
  T[3*2 + 2] = T_In[3*2 + 2];                    // i = 2, j = 2

  /* Old double nested loop
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      T[3*i + j] = T_In(i,j);
    }
  }
  */

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

double & Tensor::operator[](const uByte N) {
  // Return the specified element (treat it like matirx notation)
  return T[N];
} // double & Tensor::operator()(const uByte row, const uByte col) {

double Tensor::operator[]( const uByte N) const {
  // Return the specified element (treat it like matirx notation)
  return T[N];
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
  return T_In*c;
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
