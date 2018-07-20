#if !defined(TESTS_SOURCE)
#define TESTS_SOURCE

#include "Tests.h"
#include "Particle.h"
#include "Tensor.h"
#include "Vector.h"
#include "List.h"
#include "VTK_File.h"

void Vector_Tests(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Test vector methods

  printf("Vector method tests: \n\n");

  // Check that Component constructor works
  Vector V1;
  Vector V2(1,2,3);
  printf("V2(1,2,3)    : ");
  V2.Print();

  // Check that Vector-Vector Equality works
  V1 = V2;
  printf("V1 = V2      : ");
  V1.Print();

  // Check that Vector-Array equality works
  V2 = {4,5,6};
  printf("V2 = {4,5,6} : ");
  V2.Print();

  // Check Vector-Vector addition
  Vector V3 = V1 + V2;
  printf("V3 = V1 + V2 : ");
  V3.Print();

  // Check that Vector-Vector subtraction works
  V3 = V1 - V2;
  printf("V3 = V1 - V2 : ");
  V3.Print();

  // Check that scalar multiplication works
  V3 = V1*5;
  printf("V3 = V1*5    : ");
  V3.Print();

  // Check that the other scalar multiplication works
  V3 = 5*V1;
  printf("V3 = 5*V1    : ");
  V3.Print();

  // Check that scalar division works
  V3 = V1/((float)5);
  printf("V3 = V1/5.   : ");
  V3.Print();

  // Check that compound Vector-Vector addition works
  V1 += V2;
  printf("V1 += V2     : ");
  V1.Print();

  // Check that compound Vector-Array addition works
  V1 += {1,2,3};
  printf("V1 += {1,2,3}: ");
  V1.Print();

  // Check that compound scalar multiplication works
  V1 *= 5;
  printf("V1 *= 5      : ");
  V1.Print();

  // Check that () access works
  double v = V1(1);
  printf("v = V1(1)    : %f\n",v);

  // Test that magnitude method works
  v = Magnitude(V1);
  printf("v = |V1|     : %f\n",v);

  // Test dyadic product
  Tensor S = Dyadic_Product(V1,V2);
  printf(" S = V1 dyad V2\n");
  S.Print();

  // Test Vector dot product
  V1 = {1,2,3}; V2 = {1,2,3};
  double dot_prod = Vector_Dot_Product(V1,V2);
  printf("V1 = V2 = {1,2,3}. V1 dot V2 = %f\n", dot_prod);
} // void Vector_Tests(void) {

void Tensor_Tests(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Test Tensor methods

  printf("\nTensor Method tests \n\n");

  // Test component constructor
  Tensor S1;
  Tensor S2(1,2,3,4,5,6,7,8,9);
  printf("\nS2(1,2,...8,9): \n");
  S2.Print();

  // Test Tensor-Tensor equality
  S1 = S2;
  printf("\nS1 = S2\n");
  S2.Print();

  // Test Tensor-Array Equality
  S2 = {9,8,7,6,5,4,3,2,1};
  printf("\nS2 = {9,8,....2,1}:\n");
  S2.Print();

  // Test Tensor-Tensor Addition
  Tensor S3 = S1 + S2;
  printf("\nS3 = S1 + S2\n");
  S3.Print();

  // Test Tensor-Tensor multiplication
  S3 = S1 - S2;
  printf("\nS3 = S1 - S2\n");
  Print(S3);

  // Test Tensor-Tensor multiplication
  S3 = S1*S2;
  printf("\nS3 = S1*S2\n");
  S3.Print();

  // Test Tensor-Vector multiplication
  Vector V = {1,2,3};
  V = S1*V;
  printf("\nV = {1,2,3}\nV = S1*V\n");
  V.Print();

  // Test scalar multiplication
  S3 = S1*5;
  printf("\nS3 = S1*5\n");
  S3.Print();

  // Test other scalar multiplication
  S3 = 5*S1;
  printf("\nS3 = 5*S1\n");
  S3.Print();

  // Test Scalar division
  S3 = S1/5.;
  printf("\nS3 = S1/5.\n");
  S3.Print();

  // Test compound Tensor-Tensor addition
  S3 += S1;
  printf("\nS3 += S1\n");
  S3.Print();

  S3 -= S1;
  printf("\nS3 -= S1\n");
  Print(S3);

  // Test compound Tensor-Array addition
  S3 += {1,1,1,1,1,1,1,1,1};
  printf("\nS3 += {1,1... 1,1}\n");
  S3.Print();

  // Test compund Tensor multiplication
  S3 = S1;
  S3 *= S2;                           // together, these two lines are the same as S3 = (S1*S2)
  printf("\nS3 = S1; S3 *= S2\n");
  S3.Print();

  // Test () component access
  double t = S1(1,1);
  printf("\nt = S1(1,1)  : %f\n",t);

  // Test inverse method
  S2 = {1,4, 9, 29, 4, 67, 10, 4, 0};
  S3 = S2.Inverse();
  printf("\nS2 = {1,4,9, 29, 4, 67, 10, 4, 0} \nS3 = S2.Inverse()\n");
  S3.Print();

  // Test inverse operator
  S2 = {1,4, 9, 29, 4, 67, 10, 4, 0};
  S3 = (S2^-1);
  printf("S3 = S2^-1\n");
  S3.Print();

  // Test that Determinant method works
  double Det_S = S2.Determinant();
  printf("\nS2.Determinant = %f\n",Det_S);

  // Test Transpose method
  Tensor S3_T = S3.Transpose();
  printf("\nS3.Transpose()\n");
  S3_T.Print();

  // Test Transpose operator
  S3_T = (S3^T);
  printf("S3^T\n");
  S3_T.Print();

  // Test Inverse-Transpose operator
  S3_T = (S3^-T);
  printf("\nS3^-T\n");
  S3_T.Print();

  // Test Tensor Dot Product function
  S1 = {1,2,3,4,5,6,7,8,9};
  S2 = {1,2,3,4,5,6,7,8,9};
  double dot_prod = Tensor_Dot_Product(S1, S2);
  printf("\nS1 = S2 = {1,2,3,...9}. S1 : S2 = %f\n", dot_prod);

  // A Dyadic-Product identity (S2 and S3 should be equal at the end)
  Vector V1{1,2,3}, V2{92.392,-203.29, 5.2039};
  S1 = {1, -20, 39,
        6 ,2.293, -32.3020,
        .20392, .592, -.0001993};
  printf("\nV1 = {1,2,3}, V2 = {92.392,-203.29, 5.2039}\n");
  printf("S1 = \n");
  S1.Print();
  printf("S2 = Dyadic_Product(V1, S1^-1 * V2)\n");
  S2 = Dyadic_Product(V1, ((S1^(-1))*V2));
  S2.Print();

  printf("S3 = Dyadic_Product(V1, V2)*S1^-T\n");
  S3 = Dyadic_Product(V1, V2)*(S1^(-T));
  S3.Print();
} // void Tensor_Tests(void) {

void List_Tests(void) {
  // Create a new Node, test that the two ends point to NULL
  List<unsigned int> L1;
  printf("Newly created node. \n");
  L1.Print_Node_Info();

  L1.Add_Back(1);
  printf("Added '1' to end of list\n");
  L1.Print_Node_Info();

  int returned_value = L1.Remove_Back();
  printf("Removed rear, got %d\n",returned_value);
  L1.Print_Node_Info();

  L1.Add_Front(2);
  printf("Added '2' to front of list\n");
  L1.Print_Node_Info();

  returned_value = L1.Remove_Front();
  printf("Removed Front, got %d\n",returned_value);
  L1.Print_Node_Info();

  L1.Add_Back(2);
  L1.Add_Front(1);
  L1.Add_Back(3);
  L1.Add_Front(0);
  L1.Add_Back(4);
  printf("Added {0,1,2,3,4} in strange order \n");
  L1.Print_Node_Info();

} // void List_Tests(void) {

void Particle_Tests(void) {
} // void Particle_Tests(void) {

void Timing_Tests(void) {
  printf("\nTiming tests\n\n");

  // Set up timing variables. Note: All times will be reported in ms
  #define CLOCKS_PER_MS (CLOCKS_PER_SEC/1000.)
  long Ms_Elapsed;
  clock_t timer;
  const unsigned long Num_Tests = 10000;         // Number of tests (# of times that we cycle through the arrays)
  const unsigned long Num_El = 10000;            // Number of elements in Vector, Tensor arrays
  unsigned long i,k;                             // index variables

  // Dynamically allocate tensor, vector arrays
  double * C1 = new double[Num_El];
  Tensor * S1 = new Tensor[Num_El];
  Tensor * S2 = new Tensor[Num_El];
  Tensor * S3 = new Tensor[Num_El];
  Vector * V1 = new Vector[Num_El];
  Vector * V2 = new Vector[Num_El];

  /////////////////////////////////////////////////////////////////////////////////////////
  /* Tensor-Vector product timing tests */
  /*
  // First, populate the V1 and S1 elements
  for(i = 0; i < Num_El; i++) {
    V1[i] = {1,1,1};
    S1[i] = {1,0,0,
            0,1,0,
            0,0,1};
  }

  // Cycle through the Num_EL tensors Num_Tests times.
  timer = clock();
  for(k = 0; k < Num_Tests; k++) {
    for(i = 0; i < Num_El; i++) {
      V2[i] = S1[i]*V1[i];
    }
  }
  timer = clock() - timer;

  Ms_Elapsed = (int)((double)timer / (double)CLOCKS_PER_MS);
  printf("It took %ld ms to compute %3.0e Tensor-Vector products \n",Ms_Elapsed, (double)Num_Tests*Num_El);
  */
  /////////////////////////////////////////////////////////////////////////////////////////
  /* Tensor addition + Multiplication by a vector test*/
  /*
  // Populate the S1, S2, and V1 arrays
  for(i = 0; i < Num_El; i++) {
    V1[i] = {1,1,1};
    S1[i] = {1,2,3,
             4,5,6,
             7,8,9};
    S2[i] = {1,0,0,
             0,1,0,
             0,0,1};
  }

  // Two random scalars to force the compiler to perform tensor-scalar or vector-scalar
  //multiplication
  double d1 = rand();
  double d2 = rand();

  // Cycle through the Num_El tensors Num_Tests times.
  timer = clock();
  for(k = 0; k < Num_Tests; k++) {
    for(i = 0; i < Num_El; i++) {
      V2[i] = (d1*d2)*((S1[i] + S2[i])*V1[i]);
    }
  }
  timer = clock() - timer;

  Ms_Elapsed = (int)((double)timer / (double)CLOCKS_PER_MS);
  printf("It took %ld ms to compute %3.0e (T+T)*Vs \n",Ms_Elapsed, (double)Num_Tests*Num_El);
  */

  ///////////////////////////////////////////////////////////////////////////////////////
  // Maximum Eigenvalues test

  for(i = 0; i < Num_El; i++) {
    S1[i] = {     1,     7,   2.3,
                  7, 2.329,  -4.2,
                2.3,  -4.2, 1.392};
  } // for(i = 0; i < Num_El; i++) {

  timer = clock();

  for(k = 0; k < Num_Tests; k++) {
    for(i = 0; i < Num_El; i++) {
      C1[i] = Max_Component(Eigenvalues(S1[i],'F'));
    }
  }
  timer = clock() - timer;
  Ms_Elapsed = (int)((double)timer / (double)CLOCKS_PER_MS);
  printf("It took %ld ms to compute %3.0e max eigenvalues \n",Ms_Elapsed, (double)Num_Tests*Num_El);


  // Free the dynamic tensor, vector arrays.
  delete [] C1;
  delete [] S1;
  delete [] S2;
  delete [] S3;
  delete [] V1;
  delete [] V2;
} // void Timing_Tests(void) {

#endif
