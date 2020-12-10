#include "Tensor/Tensor.h"
#include "Vector/Vector.h"
#include "Errors.h"
#include <stdlib.h>
#include <time.h>

const double CLOCKS_PER_MS = CLOCKS_PER_SEC/1000.;

static void Calculate_Force(Vector & F,
                            const double V_j,
                            const Tensor & T1,
                            const Tensor & T2,
                            const Vector & Grad_Wj) {
  /* This function is used to calculate Force_Internal and Force_Viscosity due
  to one of a particle's neighbors. These quantities are calculated by the
  following expression:
      Force_Internal  += V_j*((P_i + P_j)*Grad_W[j])
      Force_Viscosity += V_j*((Visc + Particles[Neighbor_ID].Visc)*Grad_W[j])
  Thus, this function computes the following:
      F += V_j*((T1 + T2)*GradW_j) without any
  The goal is to eliminate any use of temporary objects and (hopefully) improve
  runtime.

  update_x is the only thing that should call this function. */

  /* Note: Tensors are stored in ROW MAJOR ordering. Thus, we want to change
  rows as infrequently. */
  const double * T1_Ar = T1.Get_Ar();
  const double * T2_Ar = T2.Get_Ar();
  const double * Grad_Wj_Ar = Grad_Wj.Get_Ar();


  F[0] += + V_j*( (T1_Ar[0*3 + 0] + T2_Ar[0*3 + 0])*Grad_Wj_Ar[0] +
                  (T1_Ar[0*3 + 1] + T2_Ar[0*3 + 1])*Grad_Wj_Ar[1] +
                  (T1_Ar[0*3 + 2] + T2_Ar[0*3 + 2])*Grad_Wj_Ar[2] );

  F[1] += + V_j*( (T1_Ar[1*3 + 0] + T2_Ar[1*3 + 0])*Grad_Wj_Ar[0] +
                  (T1_Ar[1*3 + 1] + T2_Ar[1*3 + 1])*Grad_Wj_Ar[1] +
                  (T1_Ar[1*3 + 2] + T2_Ar[1*3 + 2])*Grad_Wj_Ar[2] );

  F[2] += + V_j*( (T1_Ar[2*3 + 0] + T2_Ar[2*3 + 0])*Grad_Wj_Ar[0] +
                  (T1_Ar[2*3 + 1] + T2_Ar[2*3 + 1])*Grad_Wj_Ar[1] +
                  (T1_Ar[2*3 + 2] + T2_Ar[2*3 + 2])*Grad_Wj_Ar[2] );
} // static void Calculate_Force(Tensor & F,...



static double Calculate_Delta(const Tensor & F,
                              const Vector & R_j,
                              const Vector & rj,
                              const double Mag_rj) {
  /* This function computes delta_ij and delta_ji. By definition,
        delta_ji = Dot_Product(F_j*R[j], rj)/(Mag_rj) - Mag_rj
        delta_ij = Dot_Product(F_i*R[j], rj)/(Mag_rj) - Mag_rj
  Thus, this function returns the following quantity:
        Dot_Product(F*R_j, rj)/(Mag_rj) - Mag_rj
  The goal is to eliminate any use of temporary objects.

  update_x is the ONLY thing that should call this function. */

  const double * F_Ar  = F.Get_Ar();
  const double * Rj_Ar = R_j.Get_Ar();
  const double * rj_Ar = rj.Get_Ar();


  double FRj_0 = F_Ar[0*3 + 0]*Rj_Ar[0] + F_Ar[0*3 + 1]*Rj_Ar[1] + F_Ar[0*3 + 2]*Rj_Ar[2];
  double FRj_1 = F_Ar[1*3 + 0]*Rj_Ar[0] + F_Ar[1*3 + 1]*Rj_Ar[1] + F_Ar[1*3 + 2]*Rj_Ar[2];
  double FRj_2 = F_Ar[2*3 + 0]*Rj_Ar[0] + F_Ar[2*3 + 1]*Rj_Ar[1] + F_Ar[2*3 + 2]*Rj_Ar[2];

  return (FRj_0*rj_Ar[0] + FRj_1*rj_Ar[1] + FRj_2*rj_Ar[2])/(Mag_rj) - Mag_rj;
} // static double Calculate_Delta(const Tensor & F,





TEST_CASE("Calculate_Force_Correctness","[Body]") {
  /* Seed the random number generator with the current time */
  srand(time(0));

  /* Repeatedly fill some tensor, vectors with random values, compute their
  product and then check that Calculate_Force gives the same result. */
  Vector F1, F2, GradW_j;
  Tensor T1, T2;
  double V_j;
  const unsigned n_trials = 10000;
  unsigned n_times_not_equal = 0;

  for(unsigned n = 0; n < n_trials; n++) {
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
    Calculate_Force(F2, V_j, T1, T2, GradW_j);

    // Check that the two approaches give the same tensor.
    if(F1 != F2) { n_times_not_equal++; }
  } // for(unsigned n = 0; n < n_trials; n++) {

  REQUIRE( n_times_not_equal == 0 );
} // TEST_CASE("Calculate_Force_Correctness","[Body]") {



TEST_CASE("Calculate_Force_Speed", "[Body]") {
  /* Seed the random number generator with the current time */
  srand(time(0));

  /* Populate 100 vectors/tensors with random values (for R_i, Rj and rj)
  We can then use these to calculate delta using the random quantities. We
  check that Calculate_Delta gives the same result. */
  Vector F[100], Grad_W[100];
  Tensor T[100], S[100];
  double V[100];

  for(unsigned n = 0; n < 100; n++) {
    // Set up F1, F2, T1, T2, GradW_j, and V_j.
    V[n] = ((double)rand()/((double)RAND_MAX));
    for(unsigned i = 0; i < 3; i++) {
      F[n][i]      = ((double)rand()/((double)RAND_MAX));
      Grad_W[n][i]  = ((double)rand()/((double)RAND_MAX));

      for(unsigned j = 0; j < 3; j++) {
        T[n][i*3 + j] = ((double)rand()/((double)RAND_MAX));
        S[n][i*3 + j] = ((double)rand()/((double)RAND_MAX));
      } // for(unsigned j = 0; j < 3; j++) {
    } // for(unsigned i = 0; i < 3; i++) {
  } // for(unsigned n = 0; n < 100; n++) {


  // Compute Force using the old method a bunch of times and time it.
  const unsigned n_trials = 1000000;

  clock_t timer = clock();
  double sum = 0;            // I store the sum of the 0,0 component of the F vector in here so that the compiler won't optimize those calculations out.
  for(unsigned n = 0; n < n_trials; n++) {
    unsigned j = rand()%100;
    F[j] += V[j]*(T[j] + S[j])*Grad_W[j];
    sum += F[j][0];
  } // for(unsigned n = 0; n < n_trials; n++) {
  clock_t time_old = clock() - timer;

  // Now calculate it using the new method a bunch of times and time it.
  timer = clock();
  for(unsigned n = 0; n < n_trials; n++) {
    unsigned j = rand()%100;
    Calculate_Force(F[j], V[j], T[j], S[j], Grad_W[j]);
    sum += F[j][0];
  } // for(unsigned n = 0; n < n_trials; n++) {
  clock_t time_new = clock() - timer;

  unsigned long MS_Old = (unsigned long)((double)time_old / (double)CLOCKS_PER_MS);
  unsigned long MS_New = (unsigned long)((double)time_new / (double)CLOCKS_PER_MS);
  printf("Old Runtime: %lu ms\n", MS_Old);
  printf("New Runtime: %lu ms\n", MS_New);
  printf("sum = %lf\n", sum);
} // TEST_CASE("Calculate_Delta_Speed", "[Body]") {



TEST_CASE("Calculate_Delta_Correctness", "[Body]") {
  /* Seed the random number generator with the current time */
  srand(time(0));

  /* Repeatedly fill some tensors and vectors with random values, compute delta
  from them, and then check that Calculate_Delta gives the same result. */
  Tensor F_i{};
  Vector Rj{}, rj{};
  double Mag_rj;
  const unsigned n_trials = 10000;
  unsigned n_times_not_equal = 0;

  for(unsigned n = 0; n < n_trials; n++) {
    // Populate Rj, rj, and F with random values.
    for(unsigned i = 0; i < 3; i++) {
      Rj[i] = ((double)rand()/((double)RAND_MAX));
      rj[i] = ((double)rand()/((double)RAND_MAX));
      for(unsigned j = 0; j < 3; j++) {
        F_i[i*3 + j] = ((double)rand()/((double)RAND_MAX));
      } // for(unsigned j = 0; j < 3; j++) {
    } // for(unsigned i = 0; i < 3; i++) {
    Mag_rj = Magnitude(rj);

    /* Use these random values to calculate delta_ij. Compare this to the value
    given by the Calculate_Delta function. */
    double delta_1 = Dot_Product(F_i*Rj, rj)/(Mag_rj) - Mag_rj;
    double delta_2 = Calculate_Delta(F_i, Rj, rj, Mag_rj);

    if(delta_1 != delta_2) { n_times_not_equal++; }
  } // for(unsigned n = 0; n < n_trials; n++) {

  REQUIRE( n_times_not_equal == 0 );
} // TEST_CASE("Calculate_Delta", "[Body]") {


TEST_CASE("Calculate_Delta_Speed", "[Body]") {
  /* Seed the random number generator with the current time */
  srand(time(0));

  /* Populate 100 vectors/tensors with random values (for R_i, Rj and rj)
  We can then use these to calculate delta using the random quantities. We
  check that Calculate_Delta gives the same result. */
  Tensor F[100];
  Vector R[100], r[100];
  double Mag_r[100];

  for(unsigned n = 0; n < 100; n++) {
    // Populate Rj, rj, and F with random values.
    for(unsigned i = 0; i < 3; i++) {
      R[n][i] = ((double)rand()/((double)RAND_MAX));
      r[n][i] = ((double)rand()/((double)RAND_MAX));
      for(unsigned j = 0; j < 3; j++) {
        F[n][i*3 + j] = ((double)rand()/((double)RAND_MAX));
      } // for(unsigned j = 0; j < 3; j++) {
    } // for(unsigned i = 0; i < 3; i++) {
    Mag_r[n] = Magnitude(r[n]);
  } // for(unsigned n = 0; n < 100; n++) {


  // Compute delta_ji using the old method a bunch of times and time it.
  const unsigned n_trials = 1000000;

  clock_t timer = clock();
  double sum = 0;            // I store the sum of the delta's in here so that the compiler won't optimize those calculations out.
  for(unsigned n = 0; n < n_trials; n++) {
    unsigned j = rand()%100;
    sum += Dot_Product(F[j]*R[j], r[j])/(Mag_r[j]) - Mag_r[j];
  } // for(unsigned n = 0; n < n_trials; n++) {
  clock_t time_old = clock() - timer;

  // Now calculate it using the new method a bunch of times and time it.
  timer = clock();
  for(unsigned n = 0; n < n_trials; n++) {
    unsigned j = rand()%100;
    sum += Calculate_Delta(F[j], R[j], r[j], Mag_r[j]);
  } // for(unsigned n = 0; n < n_trials; n++) {
  clock_t time_new = clock() - timer;

  unsigned long MS_Old = (unsigned long)((double)time_old / (double)CLOCKS_PER_MS);
  unsigned long MS_New = (unsigned long)((double)time_new / (double)CLOCKS_PER_MS);
  printf("Old Runtime: %lu ms\n", MS_Old);
  printf("New Runtime: %lu ms\n", MS_New);
  printf("sum = %lf\n", sum);
} // TEST_CASE("Calculate_Delta_Speed", "[Body]") {
