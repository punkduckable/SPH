#if !defined(LIST_HEADER)
#define LIST_HEADER

template <typename Type>
class List {
  private:
    struct Node {
      Type Value;
      Node *Next_Node;
      Node *Previous_Node;
    };

    Node* Beginning;                             // Start of the list
    Node* End;                                   // End of the list
    unsigned int Num_Nodes;                      // Number of nodes in the list

  public:
    // Constructors, destructor
    List();                                      // Default constructor
    List(const List & L_In);                     // Copy constructor (to prevent inadvertent shallow copying)
    ~List();                                     // Destructor

    // Operator overloading
    List & operator=(const List & L_In);         // Explicit = operator overload (to prevent inadvertent shallow copying)
    Type & operator[](const unsigned int i);     // Acces a particular element of the list.
    const Type operator[](const unsigned int i) const;


    // Method to add or remove items from the list
    void Add_Back(const Type ID_In);             // Add an item to the start of the list
    void Add_Front(const Type ID_In);            // Add an item to the end of the list
    Type Remove_Back(void);                      // Remove an item from the end of the list
    Type Remove_Front(void);                     // Add an item to the start of the list

    // Methods to get list info
    unsigned int Node_Count(void) const;         // Returns number of nodes
    void Print_Node_Info(void) const;            // For testing
};

////////////////////////////////////////////////////////////////////////////////
// Constructors, destructor
template <typename Type>
List<Type>::List(void) {
  Num_Nodes = 0;
  Beginning = NULL;
  End = NULL;
} // List<Type>::List() {



template <typename Type>
List<Type>::List(const List & L_In) {
  /* Since list type objects dynamically allocate memory, we need to perform a
  deep copy. To do this, we cycle throught the nodes of L_In. For each node,
  we acquire its value and then add a new node with the same value to the
  end of our list. */
  Node* Node_Ptr = L_In.Beginning;
  Type value;

  Beginning = Node_Ptr;
  End = Beginning = NULL;
  Num_Nodes = 0;

  while(Node_Ptr != NULL) {
    value = Node_Ptr->Value;                         // Get latest node's value
    Add_Back(value);                                        // Append a new node onto our list
    Num_Nodes++;                                           // Increment number of nodes
    Node_Ptr = Node_Ptr->Next_Node;                        // Point to next node in the list
  } //   while(Node_Ptr != NULL) {
} // List<Type>::List(const List & L_In) {



template <typename Type>
List<Type>::~List(void) {
  /* To delete this node, we need to free all memory that has bee nallocated by
  our list. To do this, we need to delete our nodes one by one. */

  Node * temp;
  while(Beginning != NULL) {
    temp = Beginning;                  // Store current Beginning in temp (so we can free it!)
    Beginning = Beginning->Next_Node;  // Move Beginning forward one node
    delete temp;                       // free the old Beginning node

    Num_Nodes--;
    //printf("%d nodes remaining\n",Num_Nodes);
  } // while(Beginning != NULL) {
} // List<Type>::~List(void) {

////////////////////////////////////////////////////////////////////////////////
// Operator overloading

// List equality: Deep copy
template <typename Type>
List<Type> & List<Type>::operator=(const List<Type> & L_In) {
  /* Since list type objects dynamically allocate memory, we need to perform a
  deep copy. To do this, we cycle throught the nodes of L_In. For each node,
  we acquire its value and then add a new node with the same value to the
  end of our list. */
  Node* Node_Ptr = L_In.Beginning;
  Type value;

  Beginning = End = NULL;
  Num_Nodes = 0;

  while(Node_Ptr != NULL) {
    value = Node_Ptr->Value;                         // Get latest node's value
    Add_Back(value);                                        // Append a new node onto our list
    Num_Nodes++;                                           // Increment number of nodes
    Node_Ptr = Node_Ptr->Next_Node;                        // Point to next node in the list
  } //   while(Node_Ptr != NULL) {

  return *this;
} // List & List<Type>::operator=(const List & L_In) {



template <typename Type>
Type & List<Type>::operator[](const unsigned int i) {
  // Check that the requested index is in the bounds of our list
  if(i >= Num_Nodes) {
    printf("Error! requested list index is out of bounds! Returning value of last node.\n");
    return End->Value;
  } // if(i >= Num_Nodes) {

  // Traverse the list until we get to the requested node.
  Node * Current_Node = Beginning;

  for(unsigned int j = 0; j < i; j++)
    Current_Node = Current_Node->Next_Node;

  // Now return the value of the requested node.
  return Current_Node->Value;
} // Type & List<Type>::operator[](const unsigned int i) {



template <typename Type>
const Type List<Type>::operator[](const unsigned int i) const {
  // Check that the requested index is in the bounds of our list
  if(i >= Num_Nodes) {
    printf("Error! requested list index is out of bounds! Returning value of last node.\n");
    return End->Value;
  } // if(i >= Num_Nodes) {

  // Traverse the list until we get to the requested node.
  Node * Current_Node = Beginning;

  for(unsigned int j = 0; j < i; j++)
    Current_Node = Current_Node->Next_Node;

  // Now return the value of the requested node.
  return Current_Node->Value;
} // const Type List<Type>::operator[](const unsigned int i) const

////////////////////////////////////////////////////////////////////////////////
// Method to add or remove items from the list

template <typename Type>
void List<Type>::Add_Back(const Type ID_In) {
  /* This method is used to add a node onto the end of our list. We do this by
  dynamically allocating a node, adding it into our list (old End, if there is
  one, points to the new node, make End point to new node). */
  Node * add = new Node;                       // Dynamically allocate the new node
  add->Value = ID_In;                    // New Node's ID is the input ID
  add->Next_Node = NULL;                       // New node is End, so there is no next node

  /* Check if there are any nodes in our list.

  If the list is empty (no nodes), then the added element will be both the
  Beginning and the End. However, this means that there is no previous node.

  If our list already has elements, then add is the new End. The old End
  is the Previous_Node for add.*/
  if(Num_Nodes == 0) {
    Beginning = add;
    add->Previous_Node = NULL;                 // Add is new Beginning (and End), so there is no previous node.
  }
  else {
    add->Previous_Node = End;                  // add is new End. Previous node previous node for add (new end)
    End->Next_Node = add;                      // old end's new node is add (new end)
  }

  End = add;                                   // add is new End
  Num_Nodes++;                                 // We added a new node, incremenet number of nodes
} // void List<Type>::Add_Back(const Type ID_In) {



template <typename Type>
void List<Type>::Add_Front(const Type ID_In) {
  /* this method adds a node onto the Beginning of our list. We do this by
  dynamically allocating a node and adding it to the Beginning of our list (this
  new node will become the new Beginning/will point to the old Beginning) */

  Node *add = new Node;                        // Dynamically allocate new node
  add->Value = ID_In;                    // The new node's value is the input ID
  add->Previous_Node = NULL;                   // add will be new Beginning, so there is no node before add.

  /* Check if there are any nodes in our list:

  If the list is empty, then add is both the Beginning and the End. This means that
  there is no next node to point to.

  If the list is not empty, then add is the new Beginning. The old Beginning is the
  Next_Node for add */
  if(Num_Nodes == 0) {
    End = add;
    add->Next_Node = NULL;
  }
  else {
    add->Next_Node = Beginning;                  // Old Beginning is next_node for add (new Beginning).
    Beginning->Previous_Node = add;              // add (new Beginning) is previous node of old Beginning.
  }

  Beginning = add;                               // add is new Beginning
  Num_Nodes++;                                   // Increment number of nodes
} // void List<Type>::Add_Front(const Type ID_In) {



template <typename Type>
Type List<Type>::Remove_Back(void) {
  /* This method removes the last node on our list (if there is one). This is
  done by having End point to the Old End's Previous node and then freeing the
  old End. */

  // Check if there is a node to remove.
  if(End == NULL)
    return -1;

  Type ID_Out = End->Value;

  Node * temp = End;                  // temp node points to old End
  End = temp->Previous_Node;          // Move End backward (to old End's previous node)
  delete temp;                         // free memory allocated by old End

  Num_Nodes--;                         // decrement number of nodes
  //printf("Node removed from end. %d nodes remaining\n",Num_Nodes);      // For testing

  /* check if there are any nodes left in the list. If there are none, then
  Beginning should point to NULL */
  if(Num_Nodes == 0)
    Beginning = NULL;

  return ID_Out;
} // Type List<Type>::Remove_Back(void) {



template <typename Type>
Type List<Type>::Remove_Front(void) {
  /* This method removes the first node on our list (if there is one). This is
  done by having Beginning point to the Old Beginning's Next_Node, and then freeing the
  old Beginning. */

  // Check if there is a noe to remove
  if(Beginning == NULL)
    return -1;

  Type ID_Out = Beginning->Value;

  Node * temp = Beginning;                       // temp points to old Beginning.
  Beginning = temp->Next_Node;                   // Move Beginning forward
  delete temp;                                   // free memory allocated by old Beginning

  Num_Nodes--;                                   // decrement number of nodes
  //printf("Node remove from Beginning. %d nodes remaining\n",Num_Nodes);       // For testing

  /* Check if there are any nodes left in the list. If there are none, then End
  should point to NULL. */
  if(Num_Nodes == 0)
    End = NULL;

  return ID_Out;
} // Type List<Type>::Remove_Front(void) {

////////////////////////////////////////////////////////////////////////////////
// Methods to get list info

template <typename Type>
unsigned int List<Type>::Node_Count(void) const {
  return Num_Nodes;
} // void List<Type>::Node_Count(void) const {

template <typename Type>
void List<Type>::Print_Node_Info(void) const {
  Node * Current_Node = Beginning;
  int Node_Num = 0;

  // Print info on list
  printf("Num Nodes: %d\n",Num_Nodes);
  printf("Beginning: %p\n",Beginning);
  printf("End:       %p\n\n",End);

  while(Current_Node != NULL) {
    // Print info on current node
    printf("Node #%d\n",Node_Num);
    printf("Value:     %d\n",Current_Node->Value);
    printf("Current    %p\n",Current_Node);
    printf("Previous:  %p\n",Current_Node->Previous_Node);
    printf("Next:      %p\n\n",Current_Node->Next_Node);

    // Go to next node, cycle node count
    Current_Node = Current_Node->Next_Node;
    Node_Num++;
  } //   while(Current_Node != NULL) {
} // void List<Type>::Print_Node_Info(void) const {

#endif
