#include "Vector/Vector.h"
#include "Errors.h"
#include "Catch.hpp"

TEST_CASE("Vector_Base","[Vector][Constructor][Access_Operator][Boolean_Operators]") {
  // First, check that the component constructor and access operators work.
  Vector V1{1,2,3};
  REQUIRE(V1[0] == 1);
  REQUIRE(V1[1] == 2);
  REQUIRE(V1[2] == 3);

  REQUIRE(V1[0] == V1(0));
  REQUIRE(V1[1] == V1(1));
  REQUIRE(V1[2] == V1(2));

  // Now check that the == and != operators work
  Vector V2{1,2,3};
  REQUIRE(V1 == V2);
  REQUIRE_FALSE(V1 != V2);

  V2[0] = 0;
  REQUIRE_FALSE(V1 == V2);
  REQUIRE(V1 != V2);

  V2[0] = 1;
  V2[1] = 5;
  REQUIRE_FALSE(V1 == V2);
  REQUIRE(V1 != V2);

  V2[1] = 2;
  V2[2] = 10;
  REQUIRE_FALSE(V1 == V2);
  REQUIRE(V1 != V2);
} // TEST_CASE("Vector_Base","[Vector][Constructor][Access_Operator][Boolean_Operators]") {

TEST_CASE("Vector_Algebra","[Vector][Addition][Subtraction][Scaling]") {
  Vector V1{0,1,2};
  Vector V2{1,1,1};

  Vector V3 = V1 + V2;
  REQUIRE( V3 == Vector{1,2,3});
} // TEST_CASE("Vector_Algebra","[Vector][Addition][Subtraction][Scaling]") {

/*
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
*/
