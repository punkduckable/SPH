#include "Operation_Count.h"
#include <stdio.h>

/* Initialize OP_Count variables. */
namespace OP_Count {
  unsigned long Addition       = 0;
  unsigned long Subtraction    = 0;
  unsigned long Multiplication = 0;
  unsigned long Division       = 0;
  unsigned long Modulus        = 0;
  unsigned long Sqrt           = 0;
  unsigned long Exp            = 0;
  unsigned long Log            = 0;
  unsigned long Cos            = 0;
  unsigned long Acos           = 0;
} // namespace OP_Count {



void OP_Count::Print(void) {
  #ifndef OPERATION_COUNT
    printf("OPERATION_COUNT not defined. Nothing to report\n");
    return;
  #endif

  printf("\nOperation count:\n\n");
  printf("Addition:           %lu\n", OP_Count::Addition);
  printf("Subtraction:        %lu\n", OP_Count::Subtraction);
  printf("Multiplication:     %lu\n", OP_Count::Multiplication);
  printf("Division:           %lu\n", OP_Count::Division);
  printf("Modulus:            %lu\n", OP_Count::Modulus);
  printf("Sqrt:               %lu\n", OP_Count::Sqrt);
  printf("Exp:                %lu\n", OP_Count::Exp);
  printf("Log:                %lu\n", OP_Count::Log);
  printf("Cos:                %lu\n", OP_Count::Cos);
  printf("Acos:               %lu\n", OP_Count::Acos);
  printf("\n");
} // void OP_Count::Pint(void) {
