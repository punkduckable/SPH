//#define NDEBUG
#include "Simulation/Simulation.h"
#include "Errors.h"

int main(void) {
  try { Simulation::Run(); }

  catch(const Exception & Error_In) {
    printf("%s\n", Error_In.what());
  } // catch(Exception & Error_In) {

  return 0;
} // catch.hpp
