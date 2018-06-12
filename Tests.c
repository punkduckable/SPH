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

  // Test that magnitude method works
  v = V2.Magnitude();
  printf("v = V2.Mag.(): %f\n",v);

  // Test dyadic product
  Tensor T = Dyadic_Product(V1,V2);
  printf(" T = V1 dyad V2\n");
  T.Print();

  // Test Vector dot product
  V1 = {1,2,3}; V2 = {1,2,3};
  double dot_prod = Vector_Dot_Product(V1,V2);
  printf("V1 = V2 = {1,2,3}. V1 dot V2 = %f\n", dot_prod);
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
  printf("V = {1,2,3}\nV = T1*V\n");
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

  // Test compund Tensor multiplication
  T3 = T1;
  T3 *= T2;                           // together, these two lines are the same as T3 = (T1*T2)
  printf("T3 = T1; T3 *= T2\n");
  T3.Print();

  // Test () component access
  double t = T1(1,1);
  printf("t = T(1,1)   : %f\n",t);

  // Test inverse method
  T2 = {1,4, 9, 29, 4, 67, 10, 4, 0};
  T3 = T2.Inverse();
  printf("T1 = {1,4,9, 29, 4, 67, 10, 4, 0} \nT3 = T2.Inverse()\n");
  T3.Print();

  // Test that Determinant method works
  double Det_T = T2.Determinant();
  printf("T2.Determinant = %f\n",Det_T);

  // Test Transpose method
  Tensor T3_T = T3.Transpose();
  printf("T3^T\n");
  T3_T.Print();

  // Test Tensor Dot Product function
  T1 = {1,2,3,4,5,6,7,8,9};
  T2 = {1,2,3,4,5,6,7,8,9};
  double dot_prod = Tensor_Dot_Product(T1, T2);
  printf("T1 = T2 = {1,2,3,...9}. T1 : T2 = %f\n", dot_prod);
} // void Tensor_Tests(void) {

void List_Tests(void) {
  // Create a new Node, test that the two ends point to NULL
  List L1;
  printf("Newly created node. \n");
  L1.Cycle_Nodes();

  L1.Add_End(1);
  printf("Added '1' to end of list\n");
  L1.Cycle_Nodes();

  int returned_value = L1.Remove_End();
  printf("Removed rear, got %d\n",returned_value);
  L1.Cycle_Nodes();

  L1.Add_Front(2);
  printf("Added '2' to front of list\n");
  L1.Cycle_Nodes();

  returned_value = L1.Remove_Front();
  printf("Removed Front, got %d\n",returned_value);
  L1.Cycle_Nodes();

  L1.Add_End(2);
  L1.Add_Front(1);
  L1.Add_End(3);
  L1.Add_Front(0);
  L1.Add_End(4);
  printf("Added {0,1,2,3,4} in strange order \n");
  L1.Cycle_Nodes();
}
#endif
