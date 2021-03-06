#include "Tensor.h"
#include "Vector/Vector.h"
#include "Errors.h"
#include "Diagnostics/Operation_Count.h"
#include "Quick_Math.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* In this file, I define methods for tensor objects. It should be noted that
tensor data is stored in a 9 element array rather than a 3x3 matrix. This is
done because single dimension arrays tend to perform better and are, in my
personal opinion, easier to work with. As a result, successive elements in a
column are 3 elements away from each other. Successive elements in a row are one
element away from each other. Therefore, Ar[3*i+j] is the same as T(i,j). */

////////////////////////////////////////////////////////////////////////////////
// Constuctors and destructor

Tensor::Tensor(void) { }



Tensor::Tensor(double t11, double t12, double t13,
               double t21, double t22, double t23,
               double t31, double t32, double t33) {

  // Set the 9 individual components of S using inputs.
  Ar[0] = t11; Ar[1] = t12; Ar[2] = t13; // Row 1
  Ar[3] = t21; Ar[4] = t22; Ar[5] = t23; // Row 2
  Ar[6] = t31; Ar[7] = t32; Ar[8] = t33; // Row 3
} // Tensor:Tensor(double t11,.... double t33) {



Tensor::Tensor(const Tensor & Tensor_In) {
  /* Here we create one tensor using another. This is done element-by-element.
  Normally, such an operation would require a double nested loop (one for the
  rows and one for the columns).  However, loops have overhead. To eliminate
  this overhead, I wrote what would have been the 9 iterations of this double
  loop as 9 statemenets.

  To make this a little  more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */

  (*this).Ar[3*0 + 0] = Tensor_In[0*3 + 0];                    // i = 0, j = 0
  (*this).Ar[3*0 + 1] = Tensor_In[0*3 + 1];                    // i = 0, j = 1
  (*this).Ar[3*0 + 2] = Tensor_In[0*3 + 2];                    // i = 0, j = 2

  (*this).Ar[3*1 + 0] = Tensor_In[1*3 + 0];                    // i = 1, j = 0
  (*this).Ar[3*1 + 1] = Tensor_In[1*3 + 1];                    // i = 1, j = 1
  (*this).Ar[3*1 + 2] = Tensor_In[1*3 + 2];                    // i = 1, j = 2

  (*this).Ar[3*2 + 0] = Tensor_In[2*3 + 0];                    // i = 2, j = 0
  (*this).Ar[3*2 + 1] = Tensor_In[2*3 + 1];                    // i = 2, j = 1
  (*this).Ar[3*2 + 2] = Tensor_In[2*3 + 2];                    // i = 2, j = 2
} // Tensor::Tensor(const Tensor & Tensor_In) {



Tensor::~Tensor(void) { }





////////////////////////////////////////////////////////////////////////////////
// Tensor equality

Tensor & Tensor::operator=(const Tensor & Tensor_In) {
  /* Here we set one tensor equal to another. This is done element-by-element.
  Normally, such an operation would require a double nested loop (one for the
  rows and one for the columns).  However, loops have overhead. To eliminate
  this overhead, I wrote what would have been the 9 iterations of this double
  loop as 9 statemenets.

  To make this a little  easier to read, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */

  (*this).Ar[3*0 + 0] = Tensor_In.Ar[3*0 + 0];             // i = 0, j = 0
  (*this).Ar[3*0 + 1] = Tensor_In.Ar[3*0 + 1];             // i = 0, j = 1
  (*this).Ar[3*0 + 2] = Tensor_In.Ar[3*0 + 2];             // i = 0, j = 1

  (*this).Ar[3*1 + 0] = Tensor_In.Ar[3*1 + 0];             // i = 1, j = 0
  (*this).Ar[3*1 + 1] = Tensor_In.Ar[3*1 + 1];             // i = 1, j = 1
  (*this).Ar[3*1 + 2] = Tensor_In.Ar[3*1 + 2];             // i = 1, j = 2

  (*this).Ar[3*2 + 0] = Tensor_In.Ar[3*2 + 0];             // i = 2, j = 0
  (*this).Ar[3*2 + 1] = Tensor_In.Ar[3*2 + 1];             // i = 2, j = 1
  (*this).Ar[3*2 + 2] = Tensor_In.Ar[3*2 + 2];             // i = 2, j = 2

  // return this tensor (element wise equality is done)
  return *this;
} // Tensor & Tensor::operator=(const Tensor & Tensor_In) {





////////////////////////////////////////////////////////////////////////////////
// Simple arithmetic operators

Tensor Tensor::operator+(const Tensor & Tensor_In) const {
  // Declare some Tensor to store the sum
  Tensor Sum;

  /* Here we comput (*this) + Tensor_In. This is done component-by-componnet.
  Normally, such calculations would require a double nested loop (one for rows)
  and one for columns). However, loops have overhead, and we want this code to
  run fast. Thus, I choose to write out the 9 loop iterations as statements
  rather than in a double loop.

  I included comments that specify which loop iteration a given statement would
  correspond to (with i as the row index, j as the column index  */

  Sum.Ar[0*3 + 0] = (*this).Ar[3*0 + 0] + Tensor_In.Ar[0*3 + 0];     // i = 0, j = 0
  Sum.Ar[0*3 + 1] = (*this).Ar[3*0 + 1] + Tensor_In.Ar[0*3 + 1];     // i = 0, j = 1
  Sum.Ar[0*3 + 2] = (*this).Ar[3*0 + 2] + Tensor_In.Ar[0*3 + 2];     // i = 0, j = 2

  Sum.Ar[1*3 + 0] = (*this).Ar[3*1 + 0] + Tensor_In.Ar[1*3 + 0];     // i = 1, j = 0
  Sum.Ar[1*3 + 1] = (*this).Ar[3*1 + 1] + Tensor_In.Ar[1*3 + 1];     // i = 1, j = 1
  Sum.Ar[1*3 + 2] = (*this).Ar[3*1 + 2] + Tensor_In.Ar[1*3 + 2];     // i = 1, j = 2

  Sum.Ar[2*3 + 0] = (*this).Ar[3*2 + 0] + Tensor_In.Ar[2*3 + 0];     // i = 2, j = 0
  Sum.Ar[2*3 + 1] = (*this).Ar[3*2 + 1] + Tensor_In.Ar[2*3 + 1];     // i = 2, j = 1
  Sum.Ar[2*3 + 2] = (*this).Ar[3*2 + 2] + Tensor_In.Ar[2*3 + 2];     // i = 2, j = 2

  #ifdef OPERATION_COUNT
    // 9 additions in the calculations above.
    #pragma omp atomic update
    OP_Count::Addition += 9;
  #endif

  return Sum;
} // Tensor Tensor::operator+(const Tensor & Tensor_In) const {



Tensor Tensor::operator-(const Tensor & Tensor_In) const {
  // Declare some Tensor to store the difference
  Tensor Diff;

  /* Here we comput (*this) - Tensor_In. This is done component-by-component
  Normally, such calculations would require a double nested loop (one for rows)
  and one for columns). However, loops have overhead, and we want this code to
  run fast. Thus, I choose to write out the 9 loop iterations as statements
  rather than in a double loop.

  I included comments that specify which loop iteration a given statement would
  correspond to (with i as the row index, j as the column index  */

  Diff.Ar[0*3 + 0] = (*this).Ar[3*0 + 0] - Tensor_In.Ar[0*3 + 0];    // i = 0, j = 0
  Diff.Ar[0*3 + 1] = (*this).Ar[3*0 + 1] - Tensor_In.Ar[0*3 + 1];    // i = 0, j = 1
  Diff.Ar[0*3 + 2] = (*this).Ar[3*0 + 2] - Tensor_In.Ar[0*3 + 2];    // i = 0, j = 2

  Diff.Ar[1*3 + 0] = (*this).Ar[3*1 + 0] - Tensor_In.Ar[1*3 + 0];    // i = 1, j = 0
  Diff.Ar[1*3 + 1] = (*this).Ar[3*1 + 1] - Tensor_In.Ar[1*3 + 1];    // i = 1, j = 1
  Diff.Ar[1*3 + 2] = (*this).Ar[3*1 + 2] - Tensor_In.Ar[1*3 + 2];    // i = 1, j = 2

  Diff.Ar[2*3 + 0] = (*this).Ar[3*2 + 0] - Tensor_In.Ar[2*3 + 0];    // i = 2, j = 0
  Diff.Ar[2*3 + 1] = (*this).Ar[3*2 + 1] - Tensor_In.Ar[2*3 + 1];    // i = 2, j = 1
  Diff.Ar[2*3 + 2] = (*this).Ar[3*2 + 2] - Tensor_In.Ar[2*3 + 2];    // i = 2, j = 2

  #ifdef OPERATION_COUNT
    // 9 subtractions in the calculations above.
    #pragma omp atomic update
    OP_Count::Subtraction += 9;
  #endif

  return Diff;
} // Tensor Tensor::operator-(const Tensor & Tensor_In) const {



Tensor Tensor::operator*(const Tensor & Tensor_In) const{
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
  Prod.Ar[3*0 + 0] = (*this).Ar[3*0 + 0]*Tensor_In.Ar[3*0 + 0];      // i = 0, j = 0, p = 0
  Prod.Ar[3*0 + 1] = (*this).Ar[3*0 + 0]*Tensor_In.Ar[3*0 + 1];      // i = 0, j = 1, p = 0
  Prod.Ar[3*0 + 2] = (*this).Ar[3*0 + 0]*Tensor_In.Ar[3*0 + 2];      // i = 0, j = 2, p = 0

  Prod.Ar[3*0 + 0] += (*this).Ar[3*0 + 1]*Tensor_In.Ar[3*1 + 0];     // i = 0, j = 0, p = 1
  Prod.Ar[3*0 + 1] += (*this).Ar[3*0 + 1]*Tensor_In.Ar[3*1 + 1];     // i = 0, j = 1, p = 1
  Prod.Ar[3*0 + 2] += (*this).Ar[3*0 + 1]*Tensor_In.Ar[3*1 + 2];     // i = 0, j = 2, p = 1

  Prod.Ar[3*0 + 0] += (*this).Ar[3*0 + 2]*Tensor_In.Ar[3*2 + 0];     // i = 0, j = 0, p = 2
  Prod.Ar[3*0 + 1] += (*this).Ar[3*0 + 2]*Tensor_In.Ar[3*2 + 1];     // i = 0, j = 1, p = 2
  Prod.Ar[3*0 + 2] += (*this).Ar[3*0 + 2]*Tensor_In.Ar[3*2 + 2];     // i = 0, j = 2, p = 2


  /* i = 1 (2nd row of C) */
  Prod.Ar[3*1 + 0] = (*this).Ar[3*1 + 0]*Tensor_In.Ar[3*0 + 0];      // i = 1, j = 0, p = 0
  Prod.Ar[3*1 + 1] = (*this).Ar[3*1 + 0]*Tensor_In.Ar[3*0 + 1];      // i = 1, j = 1, p = 0
  Prod.Ar[3*1 + 2] = (*this).Ar[3*1 + 0]*Tensor_In.Ar[3*0 + 2];      // i = 1, j = 2, p = 0

  Prod.Ar[3*1 + 0] += (*this).Ar[3*1 + 1]*Tensor_In.Ar[3*1 + 0];     // i = 1, j = 0, p = 1
  Prod.Ar[3*1 + 1] += (*this).Ar[3*1 + 1]*Tensor_In.Ar[3*1 + 1];     // i = 1, j = 1, p = 1
  Prod.Ar[3*1 + 2] += (*this).Ar[3*1 + 1]*Tensor_In.Ar[3*1 + 2];     // i = 1, j = 2, p = 1

  Prod.Ar[3*1 + 0] += (*this).Ar[3*1 + 2]*Tensor_In.Ar[3*2 + 0];     // i = 1, j = 0, p = 2
  Prod.Ar[3*1 + 1] += (*this).Ar[3*1 + 2]*Tensor_In.Ar[3*2 + 1];     // i = 1, j = 1, p = 2
  Prod.Ar[3*1 + 2] += (*this).Ar[3*1 + 2]*Tensor_In.Ar[3*2 + 2];     // i = 1, j = 2, p = 2


  /* i = 2 (3rd row) */
  Prod.Ar[3*2 + 0] = (*this).Ar[3*2 + 0]*Tensor_In.Ar[3*0 + 0];      // i = 2, j = 0, p = 0
  Prod.Ar[3*2 + 1] = (*this).Ar[3*2 + 0]*Tensor_In.Ar[3*0 + 1];      // i = 2, j = 1, p = 0
  Prod.Ar[3*2 + 2] = (*this).Ar[3*2 + 0]*Tensor_In.Ar[3*0 + 2];      // i = 2, j = 2, p = 0

  Prod.Ar[3*2 + 0] += (*this).Ar[3*2 + 1]*Tensor_In.Ar[3*1 + 0];     // i = 2, j = 0, p = 1
  Prod.Ar[3*2 + 1] += (*this).Ar[3*2 + 1]*Tensor_In.Ar[3*1 + 1];     // i = 2, j = 1, p = 1
  Prod.Ar[3*2 + 2] += (*this).Ar[3*2 + 1]*Tensor_In.Ar[3*1 + 2];     // i = 2, j = 2, p = 1

  Prod.Ar[3*2 + 0] += (*this).Ar[3*2 + 2]*Tensor_In.Ar[3*2 + 0];     // i = 2, j = 0, p = 2
  Prod.Ar[3*2 + 1] += (*this).Ar[3*2 + 2]*Tensor_In.Ar[3*2 + 1];     // i = 2, j = 1, p = 2
  Prod.Ar[3*2 + 2] += (*this).Ar[3*2 + 2]*Tensor_In.Ar[3*2 + 2];     // i = 2, j = 2, p = 2

  #ifdef OPERATION_COUNT
    // 18 additions, 27 multiplications in the calculations above.
    #pragma omp atomic update
    OP_Count::Multiplication += 27;
    #pragma omp atomic update
    OP_Count::Addition += 18;
  #endif

  return Prod;
} // Tensor Tensor::operator*(const Tensor & Tensor_In) const{



Vector Tensor::operator*(const Vector & V_In) const {
  // Declare product vector (matrix vector product is a vector)
  Vector Prod;

  /* Here I compute the product Prod = (*this) * V_In. Normally, this would
  require a loop to go through the components of Prod (and possible a second
  to compute what goes in each componnet). However, Loops have overhead. To
  eliminate this overhead, I have written what would have been the loop
  iterations as a series of statements.*/

  Prod[0] =   (*this).Ar[3*0 + 0]*V_In[0] +
              (*this).Ar[3*0 + 1]*V_In[1] +
              (*this).Ar[3*0 + 2]*V_In[2];

  Prod[1] =   (*this).Ar[3*1 + 0]*V_In[0] +
              (*this).Ar[3*1 + 1]*V_In[1] +
              (*this).Ar[3*1 + 2]*V_In[2];

  Prod[2] =   (*this).Ar[3*2 + 0]*V_In[0] +
              (*this).Ar[3*2 + 1]*V_In[1] +
              (*this).Ar[3*2 + 2]*V_In[2];

  #ifdef OPERATION_COUNT
    // 9 multiplications,  6 additions in the calculations above.
    #pragma omp atomic update
    OP_Count::Multiplication += 9;
    #pragma omp atomic update
    OP_Count::Addition += 6;
  #endif

  return Prod;
} // Vector Tensor::operator*(const Vector & V_In) const {



Tensor Tensor::operator*(const double c) const {
  // We will store results in Prod
  Tensor Prod;

  /* Here we return the product (*this)*c. Normally, doing this would require
  cycling through the row's and col's of S in a double nested loop. However,
  loops have overhead. To eliminate this overhead, I wrote what would have been
  the 9 iterations of this double loop as 9 statemenets.

  To make this a little  more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */

  Prod.Ar[0*3 + 0] = (*this).Ar[3*0 + 0]*c;                // i = 0, j = 0
  Prod.Ar[0*3 + 1] = (*this).Ar[3*0 + 1]*c;                // i = 0, j = 1
  Prod.Ar[0*3 + 2] = (*this).Ar[3*0 + 2]*c;                // i = 0, j = 2

  Prod.Ar[1*3 + 0] = (*this).Ar[3*1 + 0]*c;                // i = 1, j = 0
  Prod.Ar[1*3 + 1] = (*this).Ar[3*1 + 1]*c;                // i = 1, j = 1
  Prod.Ar[1*3 + 2] = (*this).Ar[3*1 + 2]*c;                // i = 1, j = 2

  Prod.Ar[2*3 + 0] = (*this).Ar[3*2 + 0]*c;                // i = 2, j = 0
  Prod.Ar[2*3 + 1] = (*this).Ar[3*2 + 1]*c;                // i = 2, j = 1
  Prod.Ar[2*3 + 2] = (*this).Ar[3*2 + 2]*c;                // i = 2, j = 2

  #ifdef OPERATION_COUNT
    // 9 multiplications in the calculations above.
    #pragma omp atomic update
    OP_Count::Multiplication += 9;
  #endif

  return Prod;
} // Tensor Tensor::operator*(const double c) const {



Tensor Tensor::operator/(const double c) const {
  /* Check that quotient is non-zero
  Note: we keep this as an exception rather than an assertion because dividing
  by zero can happen in particular simulations (rather than just in buggy code). */
  if(c == 0) {
    throw Divide_By_Zero("Divide by Zero Exception: thrown by Tensor::operator/\n"
                         "You tried dividing a tensor by zero. Bad!\n");
  } // if(c == 0) {


  #ifdef OPERATION_COUNT
    /* 1 division in the calculation below (the multiplication is
    operator overloading and is counted elsewhere) */
    #pragma omp atomic update
    OP_Count::Division += 1;
  #endif

  /* To improve performance, we return (*this)*(1/c). One division (to
  calculate 1/c) and 9 multiplications (to to tensor-scalar multiplication)
  is faster than 9 divisions. */
  return (*this)*(1./c);
} //Tensor Tensor::operator/(const double c) const {





////////////////////////////////////////////////////////////////////////////////
// Compound Arithmetic operations

Tensor & Tensor::operator+=(const Tensor & Tensor_In) {
  /* Here we calculate the sum (*this) + Tensor_In. This is done by comuging the
  element-by-element sum of the two tensors. Normally, doing this would require
  a double loop (one for the roes, one for the cols). However, loops have
  overhead. To eliminate this overhead, I wrote what would have been
  the 9 iterations of this double loop as 9 statemenets.

  To make this a little  more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */

  (*this).Ar[3*0 + 0] += Tensor_In.Ar[3*0 + 0];       // i = 0, j = 0
  (*this).Ar[3*0 + 1] += Tensor_In.Ar[3*0 + 1];       // i = 0, j = 1
  (*this).Ar[3*0 + 2] += Tensor_In.Ar[3*0 + 2];       // i = 0, j = 2

  (*this).Ar[3*1 + 0] += Tensor_In.Ar[3*1 + 0];       // i = 1, j = 0
  (*this).Ar[3*1 + 1] += Tensor_In.Ar[3*1 + 1];       // i = 1, j = 1
  (*this).Ar[3*1 + 2] += Tensor_In.Ar[3*1 + 2];       // i = 1, j = 2

  (*this).Ar[3*2 + 0] += Tensor_In.Ar[3*2 + 0];       // i = 2, j = 0
  (*this).Ar[3*2 + 1] += Tensor_In.Ar[3*2 + 1];       // i = 2, j = 1
  (*this).Ar[3*2 + 2] += Tensor_In.Ar[3*2 + 2];       // i = 2, j = 2

  #ifdef OPERATION_COUNT
    // 9 additions in the calculations above.
    #pragma omp atomic update
    OP_Count::Addition += 9;
  #endif

  // Return this tensor (element wise summation is done)
  return *this;
} // Tensor & Tensor::operator+=(const Tensor & Tensor_In) {



Tensor & Tensor::operator-=(const Tensor & Tensor_In) {
  /* Here we calculate the difference (*this) - Tensor_In. This is done by computing the
  element-by-element difference of the two tensors. Normally, doing this would require
  a double loop (one for the roes, one for the cols). However, loops have
  overhead. To eliminate this overhead, I wrote what would have been
  the 9 iterations of this double loop as 9 statemenets.

  To make this a little  more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */

  (*this).Ar[3*0 + 0] -= Tensor_In.Ar[3*0 + 0];       // i = 0, j = 0
  (*this).Ar[3*0 + 1] -= Tensor_In.Ar[3*0 + 1];       // i = 0, j = 1
  (*this).Ar[3*0 + 2] -= Tensor_In.Ar[3*0 + 2];       // i = 0, j = 2

  (*this).Ar[3*1 + 0] -= Tensor_In.Ar[3*1 + 0];       // i = 1, j = 0
  (*this).Ar[3*1 + 1] -= Tensor_In.Ar[3*1 + 1];       // i = 1, j = 1
  (*this).Ar[3*1 + 2] -= Tensor_In.Ar[3*1 + 2];       // i = 1, j = 2

  (*this).Ar[3*2 + 0] -= Tensor_In.Ar[3*2 + 0];       // i = 2, j = 0
  (*this).Ar[3*2 + 1] -= Tensor_In.Ar[3*2 + 1];       // i = 2, j = 1
  (*this).Ar[3*2 + 2] -= Tensor_In.Ar[3*2 + 2];       // i = 2, j = 2

  #ifdef OPERATION_COUNT
    // 9 subtractions in the calculations above.
    #pragma omp atomic update
    OP_Count::Subtraction += 9;
  #endif

  // Return this tensor (element wise subtraction is done)
  return *this;
} // Tensor & Tensor::operator-=(const Tensor & Tensor_In) {



Tensor & Tensor::operator*=(const Tensor & Tensor_In) {
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
  Prod.Ar[3*0 + 0] = (*this).Ar[3*0 + 0]*Tensor_In.Ar[3*0 + 0];      // i = 0, j = 0, p = 0
  Prod.Ar[3*0 + 1] = (*this).Ar[3*0 + 0]*Tensor_In.Ar[3*0 + 1];      // i = 0, j = 1, p = 0
  Prod.Ar[3*0 + 2] = (*this).Ar[3*0 + 0]*Tensor_In.Ar[3*0 + 2];      // i = 0, j = 2, p = 0

  Prod.Ar[3*0 + 0] += (*this).Ar[3*0 + 1]*Tensor_In.Ar[3*1 + 0];     // i = 0, j = 0, p = 1
  Prod.Ar[3*0 + 1] += (*this).Ar[3*0 + 1]*Tensor_In.Ar[3*1 + 1];     // i = 0, j = 1, p = 1
  Prod.Ar[3*0 + 2] += (*this).Ar[3*0 + 1]*Tensor_In.Ar[3*1 + 2];     // i = 0, j = 2, p = 1

  Prod.Ar[3*0 + 0] += (*this).Ar[3*0 + 2]*Tensor_In.Ar[3*2 + 0];     // i = 0, j = 0, p = 2
  Prod.Ar[3*0 + 1] += (*this).Ar[3*0 + 2]*Tensor_In.Ar[3*2 + 1];     // i = 0, j = 1, p = 2
  Prod.Ar[3*0 + 2] += (*this).Ar[3*0 + 2]*Tensor_In.Ar[3*2 + 2];     // i = 0, j = 2, p = 2


  /* i = 1 (2nd row of C) */
  Prod.Ar[3*1 + 0] = (*this).Ar[3*1 + 0]*Tensor_In.Ar[3*0 + 0];      // i = 1, j = 0, p = 0
  Prod.Ar[3*1 + 1] = (*this).Ar[3*1 + 0]*Tensor_In.Ar[3*0 + 1];      // i = 1, j = 1, p = 0
  Prod.Ar[3*1 + 2] = (*this).Ar[3*1 + 0]*Tensor_In.Ar[3*0 + 2];      // i = 1, j = 2, p = 0

  Prod.Ar[3*1 + 0] += (*this).Ar[3*1 + 1]*Tensor_In.Ar[3*1 + 0];     // i = 1, j = 0, p = 1
  Prod.Ar[3*1 + 1] += (*this).Ar[3*1 + 1]*Tensor_In.Ar[3*1 + 1];     // i = 1, j = 1, p = 1
  Prod.Ar[3*1 + 2] += (*this).Ar[3*1 + 1]*Tensor_In.Ar[3*1 + 2];     // i = 1, j = 2, p = 1

  Prod.Ar[3*1 + 0] += (*this).Ar[3*1 + 2]*Tensor_In.Ar[3*2 + 0];     // i = 1, j = 0, p = 2
  Prod.Ar[3*1 + 1] += (*this).Ar[3*1 + 2]*Tensor_In.Ar[3*2 + 1];     // i = 1, j = 1, p = 2
  Prod.Ar[3*1 + 2] += (*this).Ar[3*1 + 2]*Tensor_In.Ar[3*2 + 2];     // i = 1, j = 2, p = 2


  /* i = 2 (3rd row) */
  Prod.Ar[3*2 + 0] = (*this).Ar[3*2 + 0]*Tensor_In.Ar[3*0 + 0];      // i = 2, j = 0, p = 0
  Prod.Ar[3*2 + 1] = (*this).Ar[3*2 + 0]*Tensor_In.Ar[3*0 + 1];      // i = 2, j = 1, p = 0
  Prod.Ar[3*2 + 2] = (*this).Ar[3*2 + 0]*Tensor_In.Ar[3*0 + 2];      // i = 2, j = 2, p = 0

  Prod.Ar[3*2 + 0] += (*this).Ar[3*2 + 1]*Tensor_In.Ar[3*1 + 0];     // i = 2, j = 0, p = 1
  Prod.Ar[3*2 + 1] += (*this).Ar[3*2 + 1]*Tensor_In.Ar[3*1 + 1];     // i = 2, j = 1, p = 1
  Prod.Ar[3*2 + 2] += (*this).Ar[3*2 + 1]*Tensor_In.Ar[3*1 + 2];     // i = 2, j = 2, p = 1

  Prod.Ar[3*2 + 0] += (*this).Ar[3*2 + 2]*Tensor_In.Ar[3*2 + 0];     // i = 2, j = 0, p = 2
  Prod.Ar[3*2 + 1] += (*this).Ar[3*2 + 2]*Tensor_In.Ar[3*2 + 1];     // i = 2, j = 1, p = 2
  Prod.Ar[3*2 + 2] += (*this).Ar[3*2 + 2]*Tensor_In.Ar[3*2 + 2];     // i = 2, j = 2, p = 2

  /* Copy Prod to this Tensor. To do this normally would require two
  nested loops (one for the rows and one for the cols). However, as before, we
  want to minimize overhead. This means unrolling loops. Thus, rather than using
  a double loop to cycle through the 9 elements of our tensor, I have witten out
  the 9 iterations. To make the code a little more readible, I included the
  origional iteration numbers as comments. */

  (*this).Ar[3*0 + 0] = Prod.Ar[3*0 + 0];                    // i = 0, j = 0
  (*this).Ar[3*0 + 1] = Prod.Ar[3*0 + 1];                    // i = 0, j = 1
  (*this).Ar[3*0 + 2] = Prod.Ar[3*0 + 2];                    // i = 0, j = 1

  (*this).Ar[3*1 + 0] = Prod.Ar[3*1 + 0];                    // i = 1, j = 0
  (*this).Ar[3*1 + 1] = Prod.Ar[3*1 + 1];                    // i = 1, j = 1
  (*this).Ar[3*1 + 2] = Prod.Ar[3*1 + 2];                    // i = 1, j = 2

  (*this).Ar[3*2 + 0] = Prod.Ar[3*2 + 0];                    // i = 2, j = 0
  (*this).Ar[3*2 + 1] = Prod.Ar[3*2 + 1];                    // i = 2, j = 1
  (*this).Ar[3*2 + 2] = Prod.Ar[3*2 + 2];                    // i = 2, j = 2

  #ifdef OPERATION_COUNT
    // 27 multiplications, 18 additions in the calculations above.
    #pragma omp atomic update
    OP_Count::Multiplication += 27;
    #pragma omp atomic update
    OP_Count::Addition += 18;
  #endif

  return *this;
} // Tensor & Tensor::operator*=(const Tensor & Tensor_In) {





////////////////////////////////////////////////////////////////////////////////
// Tensor element access: ()

double & Tensor::operator()(const unsigned row, const unsigned col) {
  /* Tensors have 3 rows and 3 columns.
  Thus, valid row and col indicies are 0, 1, and 2.
  Therefore, we check that row and col are both less than 3 */
  assert(row < 3);
  assert(col < 3);

  // Return the specified element (treat it like matirx notation)
  return (*this).Ar[3*row + col];
} // double & Tensor::operator()(const unsigned row, const unsigned col) {



double & Tensor::operator[](const unsigned index) {
  /* Tensors only have 9 components. Thus, the valid indicies are 0-8
  Thus, we check that index < 9.  */
  assert(index < 9);

  return (*this).Ar[index];
} // double & Tensor::operator[](const unsigned index) {



double Tensor::operator()( const unsigned row, const unsigned col) const {
  /* Tensors have 3 rows and 3 columns.
  Thus, valid row and col indicies are 0, 1, and 2.
  Therefore, we check that row and col are both less than 3 */
  assert(row < 3);
  assert(col < 3);

  // Return the specified element (treat it like matirx notation)
  return (*this).Ar[3*row + col];
} // double Tensor::operator()(const unsigned row, const unsigned col) const {



double Tensor::operator[](const unsigned index) const {
  /* Tensors only have 9 components. Thus, the valid indicies are 0-8
  Thus, we check that index < 9.  */
  assert(index < 9);

  return (*this).Ar[index];
} // double Tensor::operator[](const unsigned index) const {





////////////////////////////////////////////////////////////////////////////////
// Boolean operators

bool Tensor::operator==(const Tensor & T_In) const {
  /* First, check if any components are more than Tensor Epsilon apart. If so
  then return false */
  for(unsigned i = 0; i < 9; i++) {
    double d = (*this).Ar[i] - T_In.Ar[i];
    if(d < -Tensor_Epsilon || d > Tensor_Epsilon) { return false; }
  } // for(unsigned i = 0; i < 9; i++) {

  #ifdef OPERATION_COUNT
    // 9 subtractions in the loop above.
    #pragma omp atomic update
    OP_Count::Subtraction += 9;
  #endif

  // If all components are sufficiently close together then return true.
  return true;
} // bool Tensor::operator==(const Tensor & T_In) const {



bool Tensor::operator!=(const Tensor & T_In) const {
  // Return the negation of (*this) == T_In
  return !((*this) == T_In);
} // bool Tensor::operator!=(const Tensor & T_In) const {






////////////////////////////////////////////////////////////////////////////////
// Inverse methods

Tensor Tensor::Inverse(void) const{
  /* This method calculates and returns the inverse of this tensor. We do this
  using the fomrula for the inverse of a 3x3 matrix: Let
      | a b c |
  A = | d e f |
      | g h i |

  Then, by solving for A^(-1) in A*A^(-1) = I (which involves working through
  some awful algebra) we get,
                      | -hf + ei    ch - bi   -ce + bf |
  A^(-1) = (1/Det(A))*|  fg - di   -cg + ai    cd - af |
                      | -eg + dh    bg - ah   -bd + ae |
  Where Det(A) = aei - ahf - bdi + bgf + cdh - cge.
  */

  // First, let's define the return tensor
  Tensor Tensor_Inv;

  // We will store the determinant in a double
  double Det_S = (*this).Determinant();

  /* If the determinant is non-zero, then the matrix is singular and there is no
  Inverse. Thus, we should hault computation if the determinant is zero. */

  if(Det_S == 0) {
    char Error_Message_Buffer[500];
    sprintf(Error_Message_Buffer,
            "Singular_Matrix Exception: Thrown by Tensor::Inverse\n"
            "| %9.2e %9.2e %9.2e |\n"
            "| %9.2e %9.2e %9.2e |\n"
            "| %9.2e %9.2e %9.2e |\n"
            "This tensor is singular. No Inverse exists!\n",
            (*this)[0*3 + 0], (*this)[0*3 + 1], (*this)[0*3 + 2],
            (*this)[1*3 + 0], (*this)[1*3 + 1], (*this)[1*3 + 2],
            (*this)[2*3 + 0], (*this)[2*3 + 1], (*this)[2*3 + 2]);
    throw Singular_Matrix(Error_Message_Buffer);
  } // if(Det_S == 0) {

  Tensor_Inv.Ar[3*0 + 0] = -(*this).Ar[2*3 + 1]*(*this).Ar[1*3 + 2]
                           +(*this).Ar[1*3 + 1]*(*this).Ar[2*3 + 2]; // -hf + ei

  Tensor_Inv.Ar[3*0 + 1] =  (*this).Ar[0*3 + 2]*(*this).Ar[2*3 + 1]
                           -(*this).Ar[0*3 + 1]*(*this).Ar[2*3 + 2]; // ch - bi

  Tensor_Inv.Ar[3*0 + 2] = -(*this).Ar[0*3 + 2]*(*this).Ar[1*3 + 1]
                           +(*this).Ar[0*3 + 1]*(*this).Ar[1*3 + 2]; // -ce + bf

  Tensor_Inv.Ar[3*1 + 0] =  (*this).Ar[1*3 + 2]*(*this).Ar[2*3 + 0]
                           -(*this).Ar[1*3 + 0]*(*this).Ar[2*3 + 2]; // fg - di

  Tensor_Inv.Ar[3*1 + 1] = -(*this).Ar[0*3 + 2]*(*this).Ar[2*3 + 0]
                           +(*this).Ar[0*3 + 0]*(*this).Ar[2*3 + 2]; // -cg + ai

  Tensor_Inv.Ar[3*1 + 2] =  (*this).Ar[0*3 + 2]*(*this).Ar[1*3 + 0]
                           -(*this).Ar[0*0 + 0]*(*this).Ar[1*3 + 2]; // cd - af

  Tensor_Inv.Ar[3*2 + 0] = -(*this).Ar[1*3 + 1]*(*this).Ar[2*3 + 0]
                           +(*this).Ar[1*3 + 0]*(*this).Ar[2*3 + 1]; // -eg + dh

  Tensor_Inv.Ar[3*2 + 1] =  (*this).Ar[0*3 + 1]*(*this).Ar[2*3 + 0]
                           -(*this).Ar[0*3 + 0]*(*this).Ar[2*3 + 1]; // bg - ah

  Tensor_Inv.Ar[3*2 + 2] = -(*this).Ar[0*3 + 1]*(*this).Ar[1*3 + 0]
                           +(*this).Ar[0*3 + 0]*(*this).Ar[1*3 + 1]; // -bd + ae

  Tensor_Inv = (1./Det_S)*Tensor_Inv;

  #ifdef OPERATION_COUNT
    /* This one looks weirder than it is.
    Each component involves 1 subtraction and 2 multiplications (there's funny
    buisness with the - signs because of the way that determinants work, but
    there isn't really multiplication by -1 and addition in half of the
    above lines).

    Thus, in total, 18 multiplications, 9 subtractions, and 1 division. */
    #pragma omp atomic update
    OP_Count::Multiplication += 18;
    #pragma omp atomic update
    OP_Count::Subtraction += 9;
    #pragma omp atomic update
    OP_Count::Division += 1;
  #endif

  return Tensor_Inv;
} // Tensor Inverse(void) const{



Tensor Tensor::operator^(const int exp) {
  /* This method is defined to find the Inverse and/or tranpose of this
  tensor. We have defined a global variable, T, which is equal to 2. If
  this function is called with a 2, we return the transpose, if it is called
  with a -1, we return the inverse, and if it is called with a -2 then we
  return the inverse transpose. */

  switch (exp) {
    case -1:
      // If exp = -1 then return Inverse of T
      return (*this).Inverse();
    case T:
      // If exp = T (T = 2) then return transpose
      return (*this).Transpose();
    case -T:
      // If exp = -T (-T = -2) then return inverse transpose
      return (*this).Inverse().Transpose();
    default:
      /* If we end up in this case then exp was none of -1, T, and -T. Since
      exponentiation is only defined for these values, we throw an exception. */
      char Error_Message_Buffer[500];
      sprintf(Error_Message_Buffer,
              "Undefined Exponent Exception: Thrown by Tensor::operator^\n"
              "If S is a tensor, then S^exp is currently only defined if exp = -1 (for inverse), \n"
              "T (for transpose), or -T (for inverse transpose).\n"
              "You used exp = %d\n",
              exp);
      throw Undefined_Exponent(Error_Message_Buffer);
  } // switch (exp) {
} // Tensor & operator^(const char exp) {





////////////////////////////////////////////////////////////////////////////////
// Other Tensor methods

double Tensor::Determinant(void) const {
  #ifdef OPERATION_COUNT
    // 9 multiplications, 3 subtractions, 2 additions below
    #pragma omp atomic update
    OP_Count::Multiplication += 9;
    #pragma omp atomic update
    OP_Count::Addition += 2;
    #pragma omp atomic update
    OP_Count::Subtraction += 3;
  #endif

  return ((*this).Ar[0*3 + 0]*((*this).Ar[1*3 + 1]*(*this).Ar[2*3 + 2]
                              -(*this).Ar[2*3 + 1]*(*this).Ar[1*3 + 2])
         +(*this).Ar[0*3 + 1]*((*this).Ar[2*3 + 0]*(*this).Ar[1*3 + 2]
                              -(*this).Ar[1*3 + 0]*(*this).Ar[2*3 + 2])
         +(*this).Ar[0*3 + 2]*((*this).Ar[1*3 + 0]*(*this).Ar[2*3 + 1]
                              -(*this).Ar[2*3 + 0]*(*this).Ar[1*3 + 1]));
} // double Tensor::Determinant(void) const {



Tensor Tensor::Transpose(void) const {
  Tensor T_Transpose;                         // Will store Transpose

  /* Let A be some tensor and A^T be its transpose. The rule for a transpose is
  A(i,j) = A^T(j,i). Thus, to find the transpose we would normally cycle through
  the rows and columns of S using a double loop. However, loops have overhead
  and we want to avoid overhead if at all possible. Thus, we have instead
  written out 9 statements that would have made up the 9 iterations of this
  dobule loop.

  To make this a little more readible, I have included a comment with each
  statement that identifies which loop iteration that statement would have
  corresponded to (with i as the row index and j as the column index) */

  T_Transpose.Ar[3*0 + 0] = (*this).Ar[3*0 + 0];           // i = 0, j = 0
  T_Transpose.Ar[3*0 + 1] = (*this).Ar[3*1 + 0];           // i = 0, j = 1
  T_Transpose.Ar[3*0 + 2] = (*this).Ar[3*2 + 0];           // i = 0, j = 2

  T_Transpose.Ar[3*1 + 0] = (*this).Ar[3*0 + 1];           // i = 1, j = 0
  T_Transpose.Ar[3*1 + 1] = (*this).Ar[3*1 + 1];           // i = 1, j = 1
  T_Transpose.Ar[3*1 + 2] = (*this).Ar[3*2 + 1];           // i = 1, j = 2

  T_Transpose.Ar[3*2 + 0] = (*this).Ar[3*0 + 2];           // i = 2, j = 0
  T_Transpose.Ar[3*2 + 1] = (*this).Ar[3*1 + 2];           // i = 2, j = 1
  T_Transpose.Ar[3*2 + 2] = (*this).Ar[3*2 + 2];           // i = 2, j = 2

  return T_Transpose;
} // Tensor Tensor::Transpose(void) const {



const Vector Tensor::Eigenvalues(void) const {
  double p, p_inv, p1, p2, q, r, phi;
  Vector Eig_Values;
  Tensor B;
  Tensor I = {1,0,0,
              0,1,0,
              0,0,1};

  // First, calculate p1
  p1 = (*this).Ar[3*0 + 1]*(*this).Ar[3*0 + 1]
     + (*this).Ar[3*0 + 2]*(*this).Ar[3*0 + 2]
     + (*this).Ar[3*1 + 2]*(*this).Ar[3*1 + 2];

  #ifdef OPERATION_COUNT
    // 3 multiplications, 2 additions in the calculation above
    #pragma omp atomic update
    OP_Count::Addition += 2;
    #pragma omp atomic update
    OP_Count::Multiplication += 3;
  #endif

  /* If p1 == 0 then the Tensor is diagional . In this case, the eivenvalues are
  simply the diagional entris of the tensor */
  if(p1 == 0) {
    Eig_Values[0] = (*this).Ar[3*0 + 0];
    Eig_Values[1] = (*this).Ar[3*1 + 1];
    Eig_Values[2] = (*this).Ar[3*2 + 2];
  } // if(p1 == 0) {

  else {
    q = (1./3.)*((*this).Ar[3*0 + 0] + (*this).Ar[3*1 + 1] + (*this).Ar[3*2 + 2]);

    p2 = ((*this).Ar[3*0 + 0] - q)*((*this).Ar[3*0 + 0] - q) +
         ((*this).Ar[3*1 + 1] - q)*((*this).Ar[3*1 + 1] - q) +
         ((*this).Ar[3*2 + 2] - q)*((*this).Ar[3*2 + 2] - q) +
         2*p1;

    p = sqrt(p2/6.);
    p_inv = 1./p;
    B = (*this) - q*I;
    r = (.5)*(p_inv)*(p_inv)*(p_inv)*B.Determinant();

    #ifdef OPERATION_COUNT
      /* q:  2 additions, 1 multiplication (1/3 uses only constants/is computed at compile time)
      p2:    6 subtractions, 4 multiplications, 3 additions
      p:     1 sqrt, 1 division
      p_inv: 1 division
      B:     no computations, everything here is operator overloading (which is counted elsewhere)
      r:     4 multiplications */
      #pragma omp atomic update
      OP_Count::Multiplication += 9;
      #pragma omp atomic update
      OP_Count::Addition += 5;
      #pragma omp atomic update
      OP_Count::Subtraction += 6;
      #pragma omp atomic update
      OP_Count::Division += 2;
      #pragma omp atomic update
      OP_Count::Sqrt += 1;
    #endif

    /* In theory, r should be in (-1,1). However, because of floating point
    rounding errors, it is possible for r to be bigger than 1 or smaller
    than -1. However, we need r to be in -1,1 to use our acos function. Thus,
    if r is out of bounds then we move r back into bounds. */
    if(r >= 1) {
      /* This corresponds to phi = 0. However, if phi = 0, then cos(phi) = 1
      and cos(phi + 2*pi/3) = -1/2. */
      Eig_Values[0] = q + 2.*p;
      Eig_Values[1] = q - p;
      Eig_Values[2] = Eig_Values[1];

      #ifdef OPERATION_COUNT
        // 1 addition, 1 subtraction, and 1 multiplication in the three lines above
        #pragma omp atomic update
        OP_Count::Multiplication += 1;
        #pragma omp atomic update
        OP_Count::Addition += 1;
        #pragma omp atomic update
        OP_Count::Subtraction += 1;
      #endif
    } // if(r >= 1) {

    else if(r <= -1) {
      /* This corresponds to phi = PI/3. However, if phi = PI/3 then cos(phi) = 1/2
      and cos(phi + 2*Pi/3) = -1 */
      Eig_Values[0] = q + p;
      Eig_Values[1] = q - 2.*p;
      Eig_Values[2] = Eig_Values[0];

      #ifdef OPERATION_COUNT
        // The three lines above contain 1 addition, 1 subtraction, and 1 multiplication
        #pragma omp atomic update
        OP_Count::Multiplication += 1;
        #pragma omp atomic update
        OP_Count::Addition += 1;
        #pragma omp atomic update
        OP_Count::Subtraction += 1;
      #endif
    } // else if(r <= -1) {

    else {
      phi = (1./3.)*acos(r);
      Eig_Values[0] = q + 2.*p*cos(phi);
      Eig_Values[1] = q + 2.*p*cos(phi + 2.*PI/3.);
      Eig_Values[2] = 3.*q - Eig_Values[0] - Eig_Values[1];

      #ifdef OPERATION_COUNT
        /* phi:        1 multiplication, 1 acos (1/3 uses costants/is calculated at compile time)
        Eig_Values[0]: 1 addition, 2 multiplications, 1 cos
        Eig_Values[1]: 2 additions, 2 multiplications, 1 cos (2.*PI/3 uses all constants/is calculated at compile time)
        Eig_Values[2]: 1 multiplication, 2 subtractions. */
        #pragma omp atomic update
        OP_Count::Multiplication += 6;
        #pragma omp atomic update
        OP_Count::Addition += 3;
        #pragma omp atomic update
        OP_Count::Subtraction += 2;
        #pragma omp atomic update
        OP_Count::Cos += 2;
        #pragma omp atomic update
        OP_Count::Acos += 1;
      #endif
    } // else {
  } // else {

  return Eig_Values;
} // const Vector Tensor::Eigenvalues(void) {



void Tensor::Print(void) const {
  for(unsigned i = 0; i < 3; i++)
    printf("| %9.2e %9.2e %9.2e |\n",(*this).Ar[i*3], (*this).Ar[i*3 + 1], (*this).Ar[i*3 + 2]);
} // void Tensor::Print(void) const {


const double* Tensor::Get_Ar(void) const {
  /* This function returns the address of the tensor's internal array.

  This can be used to bypass the operator access methods and, thereby, improve
  runtime. However, it is extremely risky (no checks at all). Only use this if
  you know what you're doing. */
  return Ar;
} // const double* Tensor::Get_Ar(void) const {





////////////////////////////////////////////////////////////////////////////////
// Functions of a tensor

Tensor operator*(const double c, const Tensor & Tensor_In) { return Tensor_In*c; }

Tensor Inverse(const Tensor & Tensor_In) { return Tensor_In.Inverse(); }

double Determinant(const Tensor & Tensor_In) { return Tensor_In.Determinant(); }

Tensor Transpose(const Tensor & Tensor_In) { return Tensor_In.Transpose(); }

const Vector Eigenvalues(const Tensor & T_In) { return T_In.Eigenvalues(); }

void Print(const Tensor & Tensor_In) { Tensor_In.Print(); }

double Dot_Product(const Tensor & T1, const Tensor & T2) {
  /* This returns T1:T2, the tensor dot product of T1 and T2. This is given by
  T1(0,0)*T2(0,0) + T1(0,1)*T2(0,1) + .... T1(2,2)*T2(2,2) */

  #ifdef OPERATION_COUNT
    // 9 multiplications, 8 additions below.
    #pragma omp atomic update
    OP_Count::Multiplication += 9;
    #pragma omp atomic update
    OP_Count::Addition += 8;
  #endif

  return T1.Ar[0]*T2.Ar[0] +
         T1.Ar[1]*T2.Ar[1] +
         T1.Ar[2]*T2.Ar[2] +
         T1.Ar[3]*T2.Ar[3] +
         T1.Ar[4]*T2.Ar[4] +
         T1.Ar[5]*T2.Ar[5] +
         T1.Ar[6]*T2.Ar[6] +
         T1.Ar[7]*T2.Ar[7] +
         T1.Ar[8]*T2.Ar[8];
} // double Dot_Product(const Tensor & T1, const Tensor & T2) {
