// === Code Name: CGP/FEM  Project                                                            
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                          

// === Class Name: Mesh
// === mesh.h
// === Description:
// The class constructs the 2D unstructured mesh based on P1_P1 elements.

#ifndef MESH_H
#define MESH_H

#include <math.h>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std ;

struct Coordinate{
	
	public :
	
	int node_number ;
	double x_coordinate ;
	double y_coordinate ;

};

struct Boundary{
	
	public :

	int node_number ;
	double x_coordinate ;
	double y_coordinate ;

};

struct Element{

	public :
	
	int numeral ; // element number 
	
	// global coordinate
	double x1_global , y1_global , x2_global , y2_global , x3_global , y3_global ; // Placed on P1 and P2 // in corners
	double x4_global , y4_global , x5_global , y5_global , x6_global , y6_global ; // Placed just on P2

	// global node numbers
	int node_1 , node_2 , node_3 , node_4 , node_5 , node_6 ;
	
	// The numbers should be determined by reading the mesh.
	
	// Global Number of Neighborhood   
	
	int neighborhood_1 ; // The Neighborhood_1 adapted with adjust_1_2
	int neighborhood_2 ; // The Neighborhood_2 adapted with adjust_2_3
	int neighborhood_3 ; // The Neighborhood_3 adapted with adjust_3_1

	int number_of_neighborhoods ; // It can be equal to 1,2,3
	
};

class Mesh{

	friend class Computation ;
	friend class Green ;

	public :
		
		int polynomial_degree ;

		int pressure_node_number ;
		int velocity_node_number ;

		int number_of_elements ;

		Element *domain ;
		Coordinate *velocity ;
		Coordinate *unstructured_mesh ;

		Boundary *boundary_nodes ;
		int boundary_nodes_counter ;

		bool mesh_type ; // It is true if we use GMSH , and It is false if we do not use GMSH!
		string Name ; // The name of a mesh file, such as "MeshP1.txt"
		
		int N_center ; // the number of circles inside the domain

	public :

		Mesh (string Mesh , int Polynomial_Degree = 1 , int n_center = 0) ;

		void generateMesh() ;
		void generateElements() ;
		void numberNodes() ;

		bool checkExistence(int flag  , int counter , int S[]) ;

		double getXCoordinateVelocity(int i) ;
		double getYCoordinateVelocity(int i) ; 

		int findNodeNumber(double x , double y) ; 
		
		void resetBoundaryNodes() ;

		void(Mesh :: *generate_mesh_using_GMSH[2])() ;
		void generateMeshUsingGMSHP1() ;
		void generateMeshUsingGMSHP2() ;

		void setNeighborhoodVertexNumber() ;
		void(Mesh :: *set_neighborhood_vertex_number[2])() ;
		void setNeighborhoodVertexNumberP1() ;
		void setNeighborhoodVertexNumberP2() ;

		void setBoundaryNodes() ;
		void(Mesh :: *set_boundary_nodes[4])() ;
		void setBoundaryNodesP1() ;
		void setBoundaryNodesP2() ;

		void printVelocityGrid() ;
		void printPressureGrid() ;
		void printElementInformation() ;
		void printElementsUsingTecPlot() ;
		void printVelocityCoordinates() ;
		void printBoundaryNodes() ;

		void Destroy() ;

		~Mesh() ; 
};

#endif


