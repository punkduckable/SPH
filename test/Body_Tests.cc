#include "Tensor/Tensor.h"
#include "Vector/Vector.h"
#include "Errors.h"
#include <stdlib.h>
#include <time.h>

static void Calculate_Internal_Or_Viscosity_Force(Vector & F,
                                                  const double V_j,
                                                  const Tensor & T1,
                                                  const Tensor & T2,
                                                  const Vector & GradW_j) {
  /* This function computes F += V_j*((T1 + T2)*GradW_j) without any
  operator overloading. The goal is to eliminate any use of temporary objects
  and (hopefully) improve runtime. */

  /* Note: Tensors are stored in ROW MAJOR ordering. Thus, we want to change
  rows as infrequently. */
  for(unsigned i = 0; i < 3; i++) {
    F[i] += + V_j*( (T1[i*3 + 0] + T2[i*3 + 0])*GradW_j[0] +
                    (T1[i*3 + 1] + T2[i*3 + 1])*GradW_j[1] +
                    (T1[i*3 + 2] + T2[i*3 + 2])*GradW_j[2] );
  } // for(unsigned i = 0; i < 3; i++) {
} // static void Calculate_Internal_Or_Viscosity_Force(Tensor & F,...


TEST_CASE("Update_x_Optimizations","[Body]") {
  /* Seed the random number generator with the current time */
  srand(time(0));

  /* Repeatedly fill some tensor with random values, compute their product and
  then check that the specialized function gives the same result. */

  Vector F1, F2, GradW_j;
  Tensor T1, T2;
  double V_j;
  const unsigned n_trials = 10000;
  unsigned n_times_not_equal = 0;

  for(unsigned i = 0; i < n_trials; i++) {
    // Set up F1, F2, T1, T2, GradW_j, and V_j.
    V_j = ((double)rand()/((double)RAND_MAX));
    for(unsigned i = 0; i < 3; i++) {
      F1[i]      = ((double)rand()/((double)RAND_MAX));
      GradW_j[i] = ((double)rand()/((double)RAND_MAX));

      for(unsigned j = 0; j < 3; j++) {
        T1[i*3 + j] = ((double)rand()/((double)RAND_MAX));
        T2[i*3 + j] = ((double)rand()/((double)RAND_MAX));
      } // for(unsigned j = 0; j < 3; j++) {
    } // for(unsigned i = 0; i < 3; i++) {
    F2 = F1;

    // Calculate the product using operator overloading and the specialized function
    F1 += V_j*((T1 + T2)*GradW_j);
    Calculate_Internal_Or_Viscosity_Force(F2, V_j, T1, T2, GradW_j);

    // Check that the two approaches give the same tensor.
    if(F1 != F2) { n_times_not_equal++; }
  } // for(unsigned i = 0; i < n_trials; i++) {

  REQUIRE( n_times_not_equal == 0 );
} // TEST_CASE("read_line_after","[IO]") {
