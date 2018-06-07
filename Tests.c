#if !defined(_TESTS_SOURCE)
#define _TESTS_SOURCE

void Vector_Tests(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Test vector methods

  printf("Vector method tests: \n\n");

  // Check that default constructor initializes all components to zero
  Vector V1;
  printf("V1           : ");
  V1.Print();

  // Check that Component constructor works
  Vector V2(1,2,3);
  printf("V2(1,2,3)    : ");
  V2.Print();

  // Check that Vector-Vector Equality works
  V1 = V2;
  printf("V1 = V2      : ");
  V1.Print();

  // Check that Vector-Array equality works
  V2 = {1,2,3};
  printf("V2 = {1,2,3} : ");
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

  // check that [] access works
  v = V1[2];
  printf("v = V1[2]    : %f\n",v);

  // Test that magnitude function works
  v = V2.Magnitude();
  printf("v = V2.Mag.(): %f\n",v);
} // void Vector_Tests(void) {

void Tensor_Tests(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Test Tensor methods

  printf("\nTensor Method tests \n\n");

  // Test default tensor Constructor
  Tensor T1;
  printf("T1: \n");
  T1.Print();

  // Test component constructor
  Tensor T2(1,2,3,4,5,6,7,8,9);
  printf("T2(1,2,...8,9): \n");
  T2.Print();

  // Test Tensor-Tensor equality
  T1 = T2;
  printf("T1 = T2\n");
  T2.Print();

  // Test Tensor-Array Equality
  T2 = {9,8,7,6,5,4,3,2,1};
  printf("T2 = {9,8,....2,1}:\n");
  T2.Print();

  // Test Tensor-Tensor Addition
  Tensor T3 = T1 + T2;
  printf("T3 = T1 + T2\n");
  T3.Print();

  // Test Tensor-Tensor multiplication
  T3 = T1*T2;
  printf("T3 = T1*T2\n");
  T3.Print();

  // Test Tensor-Vector multiplication
  Vector V = {1,2,3};
  V = T1*V;
  printf("V = {1,2,3}\n V = T1*V\n");
  V.Print();

  // Test scalar multiplication
  T3 = T1*5;
  printf("T3 = T1*5\n");
  T3.Print();

  // Test other scalar multiplication
  T3 = 5*T1;
  printf("T3 = 5*T1\n");
  T3.Print();

  // Test Scalar division
  T3 = T1/5.;
  printf("T3 = T1/5.\n");
  T3.Print();

  // Test compound Tensor-Tensor addition
  T3 += T1;
  printf("T3 += T1\n");
  T3.Print();

  // Test compound Tensor-Array addition
  T3 += {1,1,1,1,1,1,1,1,1};
  printf("T3 += {1,1... 1,1}\n");
  T3.Print();

  // Test compound Scalar multiplication
  T3 *= 5;
  printf("T3 *= 5\n");
  T3.Print();

  // Test () component access
  double t = T1(1,1);
  printf("t = T(1,1)   : %f",t);

  // Test inverse method
  T2 = {1,4, 9, 29, 4, 67, 10, 4, 0};
  T3 = T2.Inverse();
  printf("T1 = {1,4,9, 29, 4, 67, 10, 4, 0} \n T3 = T2.Inverse()\n");
  T3.Print();
} // void Tensor_Tests(void) {

#endif
