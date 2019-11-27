#include "Vector/Vector.h"
#include "Errors.h"

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
  double c;

  // Addition
  Vector V3 = V1 + V2;
  REQUIRE( V3 == Vector{1,2,3});

  // Subtraction
  V3 = V1 - V2;
  REQUIRE( V3 == Vector{-1, 0, 1});

  // Scalar Multiplication
  c = 17;
  V3 = V1*c;
  REQUIRE( V3 == Vector{0, 1*c, 2*c} );

  c = 82.2;
  V3 = V1*c;
  REQUIRE( V3 == Vector{0, c*1, c*2} );

  // Scalar Division
  c = 2.39;
  V3 = V1/c;
  REQUIRE( V3 == Vector{0, 1./c, 2./c} );

  // Compound Vector Addition
  V1 = {11930, 2932, 293920};
  V2 = {229, 293, 2039};
  V2 += V1;
  REQUIRE( V2 == Vector{11930+229, 2932+293, 293920+2039} );

  // Compound Vector Subtraction
  V1 = {15829, 29382, 22};
  V2 = {67839, 28, 6902};
  V2 -= V1;
  REQUIRE( V2 == Vector{67839-15829, 28-29382, 6902-22} );

  // Compound Scalar Multiplication
  V3 = {1,2,3};
  c = 128.29;
  V3 *= c;
  REQUIRE( V3 == Vector{1*c, 2*c, 3*c} );
} // TEST_CASE("Vector_Algebra","[Vector][Addition][Subtraction][Scaling]") {



TEST_CASE("Vector_Other","[Vector][Magnitude][Max_Component][Dot_Product]") {
  Vector V1{2., 4., 4.};
  Vector V2{29., 192., 4290.};

  // Mangnitude
  double c = V1.Magnitude();
  REQUIRE(c == 6.);

  c = Magnitude(V1);
  REQUIRE(c == 6.);

  // Maximum component
  c = V2.Max_Component();
  REQUIRE(c == 4290.);

  c = Max_Component(V2);
  REQUIRE(c == 4290.);

  // Dot Product
  c = Dot_Product(V1, V2);
  REQUIRE(c == (2.*29. + 4.*192. + 4.*4290.));
} // TEST_CASE("Vector_Other","[Vector][Magnitude][Max_Component][Dot_Product]") {
