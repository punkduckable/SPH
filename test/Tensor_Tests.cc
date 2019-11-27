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
  Tensor T1{0,1,2,
            3,4,5,
            6,7,8};


  Tensor T2{2.9, -3, 69.2,
            0.5, .3, -0.9,
            2.4, 92, -2.3};

  //////////////////////////////////////////////////////////////////////////////
  // Addition
  Tensor T3 = T1 + T2;
  REQUIRE( T3 == Tensor{T1[0]+T2[0], T1[1]+T2[1], T1[2]+T2[2],
                        T1[3]+T2[3], T1[4]+T2[4], T1[5]+T2[5],
                        T1[6]+T2[6], T1[7]+T2[7], T1[8]+T2[8]} );

  // Compound Addition
  T3 = T1;
  T3 += T2;
  REQUIRE( T3 == Tensor{T1[0]+T2[0], T1[1]+T2[1], T1[2]+T2[2],
                        T1[3]+T2[3], T1[4]+T2[4], T1[5]+T2[5],
                        T1[6]+T2[6], T1[7]+T2[7], T1[8]+T2[8]} );


  //////////////////////////////////////////////////////////////////////////////
  // Subtraction
  T3 = T1 - T2;
  REQUIRE( T3 == Tensor{T1[0]-T2[0], T1[1]-T2[1], T1[2]-T2[2],
                        T1[3]-T2[3], T1[4]-T2[4], T1[5]-T2[5],
                        T1[6]-T2[6], T1[7]-T2[7], T1[8]-T2[8]} );

  T3 = T2;
  T3 -= T1;
  REQUIRE( T3 == Tensor{T2[0]-T1[0], T2[1]-T1[1], T2[2]-T1[2],
                        T2[3]-T1[3], T2[4]-T1[4], T2[5]-T1[5],
                        T2[6]-T1[6], T2[7]-T1[7], T2[8]-T1[8]} );


  //////////////////////////////////////////////////////////////////////////////
  // Scalar Multiplication
  double c = 29.3;
  T3 = T2*c;
  REQUIRE( T3 == Tensor{c*T2[0], c*T2[1], c*T2[2],
                        c*T2[3], c*T2[4], c*T2[5],
                        c*T2[6], c*T2[7], c*T2[8]} );

  T3 = c*T2;
  REQUIRE( T3 == Tensor{c*T2[0], c*T2[1], c*T2[2],
                        c*T2[3], c*T2[4], c*T2[5],
                        c*T2[6], c*T2[7], c*T2[8]} );


  //////////////////////////////////////////////////////////////////////////////
  // Scalar Division
  c = 692.3;
  T3 = T1/c;
  REQUIRE( T3 == Tensor{T1[0]*(1./c), T1[1]*(1./c), T1[2]*(1./c),
                        T1[3]*(1./c), T1[4]*(1./c), T1[5]*(1./c),
                        T1[6]*(1./c), T1[7]*(1./c), T1[8]*(1./c)} );


  //////////////////////////////////////////////////////////////////////////////
  // Vector-Tensor Multiplication
  Vector V1{1,2,3};
  Vector V2 = T1*V1;
  REQUIRE( V2 == Vector{T1[0]*V1[0] + T1[1]*V1[1] + T1[2]*V1[2],
                        T1[3]*V1[0] + T1[4]*V1[1] + T1[5]*V1[2],
                        T1[6]*V1[0] + T1[7]*V1[1] + T1[8]*V1[2]} );


  //////////////////////////////////////////////////////////////////////////////
  // Tensor-Tensor Multiplication
  T3 = T1*T2;
  double p00, p01, p02, p10, p11, p12, p20, p21, p22;
  p00 = T1[0]*T2[0] + T1[1]*T2[3] + T1[2]*T2[6];
  p01 = T1[0]*T2[1] + T1[1]*T2[4] + T1[2]*T2[7];
  p02 = T1[0]*T2[2] + T1[1]*T2[5] + T1[2]*T2[8];

  p10 = T1[3]*T2[0] + T1[4]*T2[3] + T1[5]*T2[6];
  p11 = T1[3]*T2[1] + T1[4]*T2[4] + T1[5]*T2[7];
  p12 = T1[3]*T2[2] + T1[4]*T2[5] + T1[5]*T2[8];

  p20 = T1[6]*T2[0] + T1[7]*T2[3] + T1[8]*T2[6];
  p21 = T1[6]*T2[1] + T1[7]*T2[4] + T1[8]*T2[7];
  p22 = T1[6]*T2[2] + T1[7]*T2[5] + T1[8]*T2[8];

  REQUIRE( T3 == Tensor{p00, p01, p02,
                        p10, p11, p12,
                        p20, p21, p22} );

  // Compound Tensor-Tensor Multiplication
  T3 = T1;
  T3 *= T2;
  REQUIRE( T3 == Tensor{p00, p01, p02,
                        p10, p11, p12,
                        p20, p21, p22} );
} // TEST_CASE("Tensor_Algebra", "[Tensor][Addition][Subtraction][Multiplication][Scaling]") {



TEST_CASE("Tensor_Other","[Tensor][Dot_Product][Inverse][Transpose][Determinant]") {
  Tensor T1{1,2,3,
            4,5,6,
            7,8,9};

  Tensor T2{59.202, 214.29, -3.582,
            -29.29, 293.39, -.0293,
            0.0293, -394.2, 1784.9};

  //////////////////////////////////////////////////////////////////////////////
  // Dot Product
  double c = Dot_Product(T1, T2);
  REQUIRE( c == T1[0]*T2[0] +
                T1[1]*T2[1] +
                T1[2]*T2[2] +
                T1[3]*T2[3] +
                T1[4]*T2[4] +
                T1[5]*T2[5] +
                T1[6]*T2[6] +
                T1[7]*T2[7] +
                T1[8]*T2[8] );


  //////////////////////////////////////////////////////////////////////////////
  // Inverse (basic)
  T1 = Tensor{1,0,0,
              0,1,0,
              0,0,1};
  T2 = T1.Inverse();
  REQUIRE( T2 == Tensor{1,0,0,
                        0,1,0,
                        0,0,1} );

  T2 = Inverse(T1);
  REQUIRE( T2 == Tensor{1,0,0,
                        0,1,0,
                        0,0,1} );

  // Inverse (singular)
  T1 = {1,2,3,
        4,5,6,
        7,8,9};
  REQUIRE_THROWS(T1.Inverse());
  REQUIRE_THROWS(Inverse(T1));

  // Inverse (complex)
  T1 = {3, 0, 2,
        2, 0,-2,
        0, 1, 1};
  REQUIRE(T1.Inverse() == Tensor{.2, .2, 0,
                                -.2, .3, 1,
                                 .2,-.3, 0});
  REQUIRE(Inverse(T1) == Tensor{.2, .2, 0,
                                -.2, .3, 1,
                                 .2,-.3, 0});


  //////////////////////////////////////////////////////////////////////////////
  // Transpose
  T1 = Tensor{1,2,3,
              4,5,6,
              7,8,9};
  T2 = T1.Transpose();
  REQUIRE( T2 == Tensor{1,4,7,
                        2,5,8,
                        3,6,9} );

  T2 = Transpose(T1);
  REQUIRE( T2 == Tensor{1,4,7,
                        2,5,8,
                        3,6,9} );


  //////////////////////////////////////////////////////////////////////////////
  // Determinant
  T1 = Tensor{1,2,3,
              4,5,6,
              7,8,9};
  c = T1.Determinant();
  REQUIRE( c == T1[0]*(T1[4]*T1[8] - T1[7]*T1[5])
               -T1[1]*(T1[3]*T1[8] - T1[6]*T1[5])
               +T1[2]*(T1[3]*T1[7] - T1[6]*T1[4]) );

  c = Determinant(T1);
  REQUIRE( c == T1[0]*(T1[4]*T1[8] - T1[7]*T1[5])
               -T1[1]*(T1[3]*T1[8] - T1[6]*T1[5])
               +T1[2]*(T1[3]*T1[7] - T1[6]*T1[4]) );
} // TEST_CASE("Tensor_Other","[Tensor][Dot_Product][Inverse][Transpose][Determinant]") {
