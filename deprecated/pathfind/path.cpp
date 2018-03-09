#include "mathfun.h"
#define X 8
#define Y 8
#define Z 8
#define WALL 1
#include <vector>
#include <algorithm>
#include <list>
#include <cstdio>
#include <fstream>
//global:
float grid[X][Y][Z];

class node{
public:
  int x,y,z, px, py, pz;
  float f,g,h;
  bool is_parent(int i, int j, int k){return (i==px && j==py && k==pz);}
  node(){x=y=z=0; px=py=pz=-1; f=g=h=0;}
  node(int i, int j, int k){x=i; y=j; z=k; f=g=h=0; px=py=pz=-1;}
  friend ostream& operator << (ostream& out, const node& v){
    out<<v.x<<","<<v.y<<","<<v.z;
    return out;
  }
  bool is_same(node* n){return x==n->x && y==n->y && z==n->z;}
  node * parent;
  void set_parent(node * n){parent=n; px=n->x; py=n->y; pz=n->z;}
  float cost(node * n){ //cost to get to n from this
    //return grid[x][y][z];
    return grid[n->x][n->y][n->z]-grid[x][y][z];
  }
};

float dist(node* n1){
  return n1->z;
}

class Heap_f {
public:
  bool operator() (const node* x, const node* y) const{
    return x->f > y->f;
  }
};

void pushnode(vector<node*> & succs, int x, int y, int z){
  node * newnode = new node(x,y,z);
  succs.push_back(newnode);
}

void getsuccs(vector<node*> & succs, node * n){
  int x=n->x;
  int xl=(x-1+X)%X;
  int xr=(x+1)%X;
  int y=n->y;
  int yl=(y-1+Y)%Y;
  int yr=(y+1)%Y;
  int z=n->z;
  int zl=z-1;
  int zr=z+1;

  //plane above
  if (zl >= 0){
    if (grid[x][y][zl]<WALL && !n->is_parent(x,y,zl)) 
      pushnode(succs, x,y,zl);
    if (grid[x][yr][zl]<WALL && !n->is_parent(x,yr,zl)) 
      pushnode(succs, x,yr,zl);
    if (grid[x][yl][zl]<WALL && !n->is_parent(x,yl,zl)) 
      pushnode(succs, x,yl,zl);
    if (grid[xr][y][zl]<WALL && !n->is_parent(xr,y,zl)) 
      pushnode(succs, xr,y,zl);
    if (grid[xr][yr][zl]<WALL && !n->is_parent(xr,yr,zl)) 
      pushnode(succs, xr,yr,zl);
    if (grid[xr][yl][zl]<WALL && !n->is_parent(xr,yl,zl)) 
      pushnode(succs, xr,yl,zl);
    if (grid[xl][y][zl]<WALL && !n->is_parent(xl,y,zl)) 
      pushnode(succs, xl,y,zl);
    if (grid[xl][yr][zl]<WALL && !n->is_parent(xl,yr,zl)) 
      pushnode(succs, xl,yr,zl);
    if (grid[xl][yl][zl]<WALL && !n->is_parent(xl,yl,zl)) 
      pushnode(succs, xl,yl,zl);
  }
  
  //same plane
  if (grid[x][yr][z]<WALL && !n->is_parent(x,yr,z)) 
    pushnode(succs, x,yr,z);
  if (grid[x][yl][z]<WALL && !n->is_parent(x,yl,z)) 
    pushnode(succs, x,yl,z);
  if (grid[xr][y][z]<WALL && !n->is_parent(xr,y,z)) 
    pushnode(succs, xr,y,z);
  if (grid[xr][yr][z]<WALL && !n->is_parent(xr,yr,z)) 
    pushnode(succs, xr,yr,z);
  if (grid[xr][yl][z]<WALL && !n->is_parent(xr,yl,z)) 
    pushnode(succs, xr,yl,z);
  if (grid[xl][y][z]<WALL && !n->is_parent(xl,y,z)) 
    pushnode(succs, xl,y,z);
  if (grid[xl][yr][z]<WALL && !n->is_parent(xl,yr,z)) 
    pushnode(succs, xl,yr,z);
  if (grid[xl][yl][z]<WALL && !n->is_parent(xl,yl,z)) 
    pushnode(succs, xl,yl,z);
    
  //plane below
  if (zr < Y){
    if (grid[x][y][zr]<WALL && !n->is_parent(x,y,zr)) 
      pushnode(succs, x,y,zr);
    if (grid[x][yr][zr]<WALL && !n->is_parent(x,yr,zr)) 
      pushnode(succs, x,yr,zr);
    if (grid[x][yl][zr]<WALL && !n->is_parent(x,yl,zr)) 
      pushnode(succs, x,yl,zr);
    if (grid[xr][y][zr]<WALL && !n->is_parent(xr,y,zr)) 
      pushnode(succs, xr,y,zr);
    if (grid[xr][yr][zr]<WALL && !n->is_parent(xr,yr,zr)) 
      pushnode(succs, xr,yr,zr);
    if (grid[xr][yl][zr]<WALL && !n->is_parent(xr,yl,zr)) 
      pushnode(succs, xr,yl,zr);
    if (grid[xl][y][zr]<WALL && !n->is_parent(xl,y,zr)) 
      pushnode(succs, xl,y,zr);
    if (grid[xl][yr][zr]<WALL && !n->is_parent(xl,yr,zr)) 
      pushnode(succs, xl,yr,zr);
    if (grid[xl][yl][zr]<WALL && !n->is_parent(xl,yl,zr)) 
      pushnode(succs, xl,yl,zr);
  }
}

int main(int argc, char* argv[]){
//    int seed=0;
//    if (argc>1)
//      seed=atoi(argv[argc-1]);
//    else
//      seed=TimeSeed();
//    cout<<"Seed: "<<seed<<endl;
//    srand(seed);
  ifstream fin("ugrid.out", ios::in);
  int x,y,z;
  fin>>x>>y>>z;
  for (z=0; z<Z; z++)
    for (y=0; y<Y; y++)
      for (x=0; x<X; x++){
	fin>>grid[x][y][z];
      }
  
  vector<node*> open,closed, succs, trash;
  node *start=NULL;
  start=new node(4,4,6);

  list<node *> path;

  start->g=grid[4][4][6];
  start->h=dist(start);
  start->f=start->g+start->h;
  open.push_back(start);
  push_heap(open.begin(), open.end(), Heap_f());
  
  //*********
  int steps=0;
  int maxsteps=X*Y*Z;
  
  while(!open.empty() && steps<maxsteps){
    steps++;
    node *current = open.front();
    //cout<<*current<<endl; 
    pop_heap(open.begin(), open.end(), Heap_f());
    open.pop_back();

    if (current->z==0){
      node* pathmem=current;
      while (pathmem!=start){
	path.push_front(pathmem);
    	pathmem=pathmem->parent;
      }
      path.push_front(start);
      break;
    }
    else{
      succs.clear();
      getsuccs(succs, current);
      vector<node*>::iterator vni,o_result, c_result;
      
      for (vni=succs.begin(); vni!=succs.end(); vni++){
	float newg=current->g+current->cost(*vni);
	for (o_result=open.begin(); o_result!=open.end(); o_result++)
	  if ((*o_result)->is_same(*vni)) break;
	if (o_result !=open.end())
	  if ((*o_result)->g <=newg){
	    delete *vni;
	    continue;
	  }
	
       	for (c_result=closed.begin(); c_result!=closed.end(); c_result++)
	  if ((*c_result)->is_same(*vni)) break;
	if (c_result !=closed.end())
	  if ((*c_result)->g <= newg){
	    delete *vni;
	    continue;
	  }
	
	(*vni)->g=newg;
	(*vni)->h=dist(*vni);
	(*vni)->f=(*vni)->g+(*vni)->h;
	(*vni)->set_parent(current);
	
	// Remove successor from closed if it was on it
	// but put it on the trash list, so it isn't lost
	if(c_result != closed.end()){
	  trash.push_back(*c_result);
	  closed.erase(c_result);
	}
	// Update old version of this node
	if(o_result != open.end()){
	  delete *o_result;
	  open.erase(o_result);
	  sort_heap(open.begin(), open.end(), Heap_f());
	}
	open.push_back(*vni); //heap now unsorted
	push_heap(open.begin(), open.end(), Heap_f());
      }
      
      // push n onto Closed, as we have expanded it now
      closed.push_back(current);
    }
  }
  if (!path.empty()){
    //cout<<"Path found!!\n";
    for (list<node*>::iterator lni=path.begin(); lni!=path.end(); lni++)
      cout<<grid[(*lni)->x][(*lni)->y][(*lni)->z]<<"\t"<<(*lni)->g<<"\t"<<**lni<<endl;
  }
  else{
    cout<<"No path found. ";
    if (steps >=maxsteps) cout<<"Max steps exceeded. ";
    cout<<endl;
  }
  return 0;
}
