//************************************************************************
class MapSearchNode{
public:
  unsigned int x;	 // the (x,y) positions of the node
  unsigned int y;	
  
  MapSearchNode() { x = y = 0; }
  MapSearchNode( unsigned int px, unsigned int py ) { x=px; y=py; }
  
  float GoalDistanceEstimate( MapSearchNode &nodeGoal){
    //return (pow(((float)x - (float)nodeGoal.x),2)+
    //pow(((float)y - (float)nodeGoal.y),2));
    return MAP_HEIGHT-1-y;
  }

  bool IsGoal(MapSearchNode &nodeGoal){
    return ((x == nodeGoal.x) && (y == nodeGoal.y));
  }

  void GetSuccessors( AStarSearch<MapSearchNode> *astarsearch, MapSearchNode *parent_node );

  float GetCost(MapSearchNode &successor){
    return (float) GetMap( x, y );
  }

  bool IsSameState(MapSearchNode &rhs){
    return ((x == rhs.x) && (y == rhs.y));
  }

  void PrintNodeInfo(){printf("Node position : (%d,%d)\n", x,y);}
};

// This generates the successors to the given Node. 
// It uses a helper function called AddSuccessor to give the successors 
// to the AStar class. The A* specific initialisation is done for each node 
// internally, so here you just set the state information that
// is specific to the application
void MapSearchNode::GetSuccessors( AStarSearch<MapSearchNode> *astarsearch, MapSearchNode *parent_node ){
  int parent_x = -1; 
  int parent_y = -1; 
  if( parent_node ){
    parent_x = parent_node->x;
    parent_y = parent_node->y;
  }
  MapSearchNode NewNode;
  // push each possible move except allowing the search to go backwards
  if( (GetMap( x-1, y ) < WALL) && !((parent_x == x-1) && (parent_y == y))) {
    NewNode = MapSearchNode( x-1, y );
    astarsearch->AddSuccessor( NewNode );
  }	
  if( (GetMap( x, y-1 ) < WALL) && !((parent_x == x) && (parent_y == y-1))){
    NewNode = MapSearchNode( x, y-1 );
    astarsearch->AddSuccessor( NewNode );
  }	
  if( (GetMap( x+1, y ) < WALL) && !((parent_x == x+1) && (parent_y == y))){
    NewNode = MapSearchNode( x+1, y );
    astarsearch->AddSuccessor( NewNode );
  }	
  if( (GetMap( x, y+1 ) < WALL) && !((parent_x == x) && (parent_y == y+1))){
    NewNode = MapSearchNode( x, y+1 );
    astarsearch->AddSuccessor( NewNode );
  }	
}

