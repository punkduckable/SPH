#include "Timing_Test.h"

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
      C1[i] = Max_Component(Eigenvalues(S1[i]));
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
