#if !defined(PARTICLE_FILE_SOURCE)
#define PARTICLE_FILE_SOURCE

#include "SPH_Diagnostics.h"
#include "Particle.h"

void Particle_Debugger::Get_File_Name(string & Str) {
  char Buf[6];
  sprintf(Buf,"%05d",File_Number);
  File_Number++;

  Str += "_variables_";
  Str += Buf;
  Str += ".txt";
} // void Particle_Debugger::Get_File_Name(string & Str) {

void Particle_Debugger::Export_Pariticle_Forces(const unsigned int Num_Particles, const Particle * Particles) {
  string File_Name = "Particle";
  Get_File_Name(File_Name);

  string File_Path = "./Force_Files/";
  File_Path += File_Name;

  FILE * File = fopen(File_Path.c_str(), "w");

  // Print header.
  fprintf(File," Particle Pos  |");
  fprintf(File,"        Internal Force        |");
  //fprintf(File,"        Viscous Force         |");
  fprintf(File,"        Contact Force         |");
  fprintf(File,"        Hourglass Force       |");
  fprintf(File,"\n");

  // Cycle through particles, print spacial positions of each particle
  for(unsigned int i = 0; i < Num_Particles; i++) {
    fprintf(File,"%4.1f,%4.1f,%4.1f | ",Particles[i].X[0], Particles[i].X[1], Particles[i].X[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",Particles[i].Force_Int[0], Particles[i].Force_Int[1], Particles[i].Force_Int[2]);
    //fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",Particles[i].Force_Visc[0], Particles[i].Force_Visc[1], Particles[i].Force_Visc[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",Particles[i].Force_Contact[0], Particles[i].Force_Contact[1], Particles[i].Force_Contact[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e>\n",Particles[i].Force_HG[0], Particles[i].Force_HG[1], Particles[i].Force_HG[2]);
  } // for(unsigned int i = 0; i < Num_Particles; i++) {

  fclose(File);
} // void Particle_Debugger::Export_Pariticle_Forces(const unsigned int Num_Particles, const Particle * Particles) {

void OP_Count::Reset_Counts(void) {
  /* This function, as the name would suggest, is designed to reset the
  operation counts. This is done by literally setting each of the operation
  count varialbes to 0. */

  // Tensors
  T_Default_Constructor = 0;
  T_Component_Constructor = 0;
  T_Copy_Constructor = 0;
  T_Equality = 0;
  T_T_Addition = 0;
  T_T_Subtraction = 0;
  T_T_Multiplication = 0;
  T_V_Multiplication = 0;
  T_S_Multiplication = 0;
  T_S_Division = 0;
  Compound_T_T_Addition = 0;
  Compound_T_T_Subtraction = 0;
  Compound_T_T_Multiplication = 0;
  T_Inverse = 0;
  T_Determinant = 0;
  T_Transpose = 0;
  T_Dot_Product = 0;

  // Vectors
  V_Default_Constructor = 0;
  V_Component_Constructor = 0;
  V_Copy_Constructor = 0;
  V_Equality = 0;
  V_V_Addition = 0;
  V_V_Subtraction  = 0;
  V_S_Multiplication = 0;
  V_S_Division = 0;
  Compound_V_V_Addition = 0;
  Compound_V_V_Subtraction = 0;
  Compound_V_S_Multiplication = 0;
  V_Magnitude = 0;
  V_Dot_Product = 0;

  // Other
  Dyadic_Product = 0;
} // void OP_Count::Reset_Counts(void) {

void OP_Count::Print_Counts(void) {
  printf("\nOperation count:\n\n");

  // Tensors
  printf("Tensor Default Constructor =           %u\n",T_Default_Constructor);
  printf("Tensor Component Constructor =         %u\n",T_Component_Constructor);
  printf("Tensor Copy Constructor =              %u\n",T_Copy_Constructor);
  printf("Tensor Equality =                      %u\n",T_Equality);
  printf("Tensor-Tensor Addition =               %u\n",T_T_Addition);
  printf("Tensor_Tensor_Subtraction =            %u\n",T_T_Subtraction);
  printf("Tensor-Tensor Multiplication =         %u\n",T_T_Multiplication);
  printf("Tensor-Vector Multiplication =         %u\n",T_V_Multiplication);
  printf("Tensor-Scalar Multiplication =         %u\n",T_S_Multiplication);
  printf("Tensor-Scalar Division =               %u\n",T_S_Division);
  printf("Compound Tensor-Tensor Addition =      %u\n",Compound_T_T_Addition);
  printf("Compound Tensor-Tensor Subtraction =   %u\n",Compound_T_T_Subtraction);
  printf("Compound Tensor-Tensor Multiplication =%u\n",Compound_T_T_Multiplication);
  printf("Tensor Inverse =                       %u\n",T_Inverse);
  printf("Tensor Determinant =                   %u\n",T_Determinant);
  printf("Tensor Transpose =                     %u\n",T_Transpose);
  printf("Tensor Dot Product =                   %u\n",T_Dot_Product);

  // Vectors
  printf("Vector Default Constructor =           %u\n",V_Default_Constructor);
  printf("Vector Component Constructor =         %u\n",V_Component_Constructor);
  printf("Vector Copy Constructor =              %u\n",V_Copy_Constructor);
  printf("Vector Equality =                      %u\n",V_Equality);
  printf("Vector-Vector Addition =               %u\n",V_V_Addition);
  printf("Vector-Vector Subtraction =            %u\n",V_V_Subtraction );
  printf("Vector-Scalar Multiplication =         %u\n",V_S_Multiplication);
  printf("Vector-Scalar Division =               %u\n",V_S_Division);
  printf("Compound Vector-Vector Addition =      %u\n",Compound_V_V_Addition);
  printf("Compound Vector-Vector Subtraction =   %u\n",Compound_V_V_Subtraction);
  printf("Compound Vector-Scalar Multiplication =%u\n",Compound_V_S_Multiplication);
  printf("Vector Magnitude =                     %u\n",V_Magnitude);
  printf("Vector Dot Product =                   %u\n",V_Dot_Product);

  // Other
  printf("Dyadic Product =                       %u\n",Dyadic_Product);
} // void OP_Count::Print_Counts(void) {

#endif
