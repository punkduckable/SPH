#include "Tensor/Tensor.h"
#include "Errors.h"

TEST_CASE("Tensor_Base", "[Tensor][Constructor]") {
  Tensor T1{0,1,2,
            3,4,5,
            6,7,8};

  // Check that all components are what they should be
  for(unsigned i = 0; i < 9; i++) {
    REQUIRE( T1[i] == i );
  } // for(unsigned i = 0; i < 9; i++) {

  // Now check that () and [] give the same results
  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      REQUIRE( T1[3*i + j] == T1(i,j) );
    } // for(unsigned j = 0; j < 3; j++) {
  } // for(unsigned i = 0; i < 3; i++) {


  // Now check that = works
  Tensor T2 = T1;
  for(unsigned i = 0; i < 9; i++) {
    REQUIRE( T2[i] == T1[i] );
  } // for(unsigned i = 0; i < 9; i++) {


  // Now check that == and != work
  REQUIRE( T1 == T2 );
  REQUIRE_FALSE( T1 != T2 );

  T1[5] = 293;
  REQUIRE_FALSE( T1 == T2 );
  REQUIRE( T1 != T2 );
} // TEST_CASE("Tensor_Base", "[Tensor][Constructor]") {



TEST_CASE("Tensor_Algebra", "[Tensor][Addition][Subtraction][Multiplication][Scaling]") {

} // TEST_CASE("Tensor_Algebra", "[Tensor][Addition][Subtraction][Multiplication][Scaling]") {



TEST_CASE("Tensor_Other","[Tensor][Dot_Product][Inverse][Transpose][Determinant][Eigenvalues]") {

} // TEST_CASE("Tensor_Other","[Tensor][Dot_Product][Inverse][Transpose][Determinant][Eigenvalues]") {
