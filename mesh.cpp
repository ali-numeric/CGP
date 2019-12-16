// === Code Name: CGP/FEM  Project                                                           
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                         

// === Class Name: Mesh														  
// === mesh.cpp

#include "mesh.h"

bool compareBoundary(const Boundary & a , const Boundary & b) ;

// -------------------- Constructor --------------------

Mesh :: Mesh(string Mesh , int Polynomial_Degree , int n_center){

	N_center = n_center ;

	mesh_type = true ;
	Name = Mesh ;

	polynomial_degree = Polynomial_Degree ;

	if (polynomial_degree == 1){ 
		
		ifstream fin(Name.c_str()) ;

			fin >> pressure_node_number ; 
			fin >> number_of_elements ;
			velocity_node_number = pressure_node_number ; // It's for P1-P1
			fin.close() ;
	}

	if (polynomial_degree == 2){ 
		
		ifstream fin(Name.c_str()) ;

			fin >> pressure_node_number ; 
			fin >> number_of_elements ;
			velocity_node_number = pressure_node_number ; // It's for P2-P2
			fin.close() ;
	}
}

// -------------------- Generate Mesh --------------------

void Mesh :: generateMesh(){
	
	if (mesh_type == true){
	
		generate_mesh_using_GMSH[0] = &Mesh :: generateMeshUsingGMSHP1 ;
		generate_mesh_using_GMSH[1] = &Mesh :: generateMeshUsingGMSHP2 ;
		
		(this->*generate_mesh_using_GMSH[polynomial_degree - 1])() ;

	}
}

// -------------------- Generate Mesh Using GMSH P1 --------------------

void Mesh :: generateMeshUsingGMSHP1(){

	ifstream fin (Name.c_str());
	
	fin >> pressure_node_number ; // It's for P1-P1 
	fin >> number_of_elements ;
	
	int node_number = pressure_node_number ; // It's for P1-P1

	unstructured_mesh = new Coordinate[node_number] ; 

	double z_coordinate = 0.0 ;

	for (int i = 0 ; i < node_number ; i ++){
		
		unstructured_mesh[i].node_number = i + 1 ;
		fin >> unstructured_mesh[i].x_coordinate ;  
		fin >> unstructured_mesh[i].y_coordinate ;
		fin >> z_coordinate ;			
	}
	
	domain = new Element [number_of_elements] ;
	
	int extra_index = 3 ;
	int index1 , index2 , index3 ;

	for (int i = 0 ; i < number_of_elements ; i ++){
	
		domain[i].numeral = i + 1 ;

		fin >> extra_index ;

//=================
//=================
// For Flow Inside a Circle:
// index1 + 1 --> index1
// index2 + 1 --> index2
// index3 + 1 --> index3
// Note: In this case, the FEM_Information does not match with the figure! But it is not an issue!
//=================
//=================
		
		fin >> index1 ;
		domain[i].node_1 = index1 + 1 - N_center ;
		fin >> index2 ;
		domain[i].node_2 = index2 + 1 - N_center ;
		fin >> index3 ;
		domain[i].node_3 = index3 + 1 - N_center ;
	}

	for (int i = 0 ; i < number_of_elements ; i ++ ){

		for (int j = 0 ; j < node_number ; j ++ ){

			if (domain[i].node_1 == unstructured_mesh[j].node_number){
			
				domain[i].x1_global = unstructured_mesh[j].x_coordinate ;
				domain[i].y1_global = unstructured_mesh[j].y_coordinate ;
			 }

			if (domain[i].node_2 == unstructured_mesh[j].node_number){
			
				domain[i].x2_global = unstructured_mesh[j].x_coordinate ;
				domain[i].y2_global = unstructured_mesh[j].y_coordinate ;
			 }

			if (domain[i].node_3 == unstructured_mesh[j].node_number){
			
				domain[i].x3_global = unstructured_mesh[j].x_coordinate ;
				domain[i].y3_global = unstructured_mesh[j].y_coordinate ;
			 }
		}
	}
	
	delete [] unstructured_mesh ;

	setNeighborhoodVertexNumber() ;
	setBoundaryNodes() ;
}

// -------------------- Generate Mesh Using GMSH P2 --------------------

void Mesh :: generateMeshUsingGMSHP2(){

	int corner ;
	int extra_information ;

	ifstream fin1 (Name.c_str());

	fin1 >> corner ; 
	fin1 >> extra_information ;
	
	Coordinate *p1 ; 
	p1 = new Coordinate[corner] ; 

	double z_coordinate = 0.0 ;

	for (int i = 0 ; i < corner ; i ++){
		
		p1[i].node_number = i + 1 ;
		fin1 >> p1[i].x_coordinate ;  
		fin1 >> p1[i].y_coordinate ;
		fin1 >> z_coordinate ;			
	}

	ifstream fin2 ("MeshP2.txt") ; // Modify in the future, right now we are just working with P1-P1
	
	fin2 >> pressure_node_number ; // It's for P1-P1  
	fin2 >> number_of_elements ;
	
	int node_number = pressure_node_number ; // It's for P1-P1

	unstructured_mesh = new Coordinate[node_number] ; 

	for (int i = 0 ; i < node_number ; i ++){
		
		unstructured_mesh[i].node_number = i + 1 ;
		
		fin2 >> unstructured_mesh[i].x_coordinate ;  
		fin2 >> unstructured_mesh[i].y_coordinate ;
		fin2 >> z_coordinate ;	
	}
	
	domain = new Element [number_of_elements] ;
	
	int extra_index = 3 ;
	int index1 , index2 , index3 , index4 , index5 , index6 ;
	
	int counter1 = 1 ;
	int counter2 = 1 ;
	int search = 1 ;

	for (int i = 0 ; i < number_of_elements ; i ++){
	
		domain[i].numeral = i + 1 ;

		fin2 >> extra_index ;

		fin2 >> index1 ;
		fin2 >> index2 ;
		fin2 >> index3 ;
		fin2 >> index4 ;
		fin2 >> index5 ;
		fin2 >> index6 ;
		
		int *Index ;
		Index = new int[6] ;
		
		Index[0] = index1 ;
		Index[1] = index2 ;
		Index[2] = index3 ;
		Index[3] = index4 ;
		Index[4] = index5 ;
		Index[5] = index6 ;
		
		for (int k = 0 ; k < 6 ; k ++){
			
			search = 1 ;

			for (int j = 0 ; j < corner ; j ++){

				if ((unstructured_mesh[Index[k]].x_coordinate == p1[j].x_coordinate) && (unstructured_mesh[Index[k]].y_coordinate == p1[j].y_coordinate)){
				
					if(counter1 == 1){
							domain[i].node_1 = Index[k] + 1 ;
							domain[i].x1_global = unstructured_mesh[Index[k]].x_coordinate ;
							domain[i].y1_global = unstructured_mesh[Index[k]].y_coordinate ;
							search = 2 ; 
						}

					if(counter1 == 2){
							domain[i].node_2 = Index[k] + 1 ; 
							domain[i].x2_global = unstructured_mesh[Index[k]].x_coordinate ;
							domain[i].y2_global = unstructured_mesh[Index[k]].y_coordinate ;
							search = 2 ; 
						}

					if(counter1 == 3){
							domain[i].node_3 = Index[k] + 1 ;
							domain[i].x3_global = unstructured_mesh[Index[k]].x_coordinate ;
							domain[i].y3_global = unstructured_mesh[Index[k]].y_coordinate ;
							search = 2 ; 
						}

					counter1 ++ ;
				}
			}

					if(counter2 == 1 && search == 1){
							domain[i].node_4 = Index[k] + 1 ;
							domain[i].x4_global = unstructured_mesh[Index[k]].x_coordinate ;
							domain[i].y4_global = unstructured_mesh[Index[k]].y_coordinate ;
					
						}

					if(counter2 == 2 && search == 1){
							domain[i].node_5 = Index[k] + 1 ;
							domain[i].x5_global = unstructured_mesh[Index[k]].x_coordinate ;
							domain[i].y5_global = unstructured_mesh[Index[k]].y_coordinate ;
					
						}

					if(counter2 == 3 && search == 1){
							domain[i].node_6 = Index[k] + 1 ;
							domain[i].x6_global = unstructured_mesh[Index[k]].x_coordinate ;
							domain[i].y6_global = unstructured_mesh[Index[k]].y_coordinate ;
						
						}
					
					if (search == 1) {counter2 ++ ;}
								
		}

		counter1 = 1 ;
		counter2 = 1 ;
		search = 1 ;
	}

	delete [] p1 ;
	delete [] unstructured_mesh ;

	setNeighborhoodVertexNumber() ;
	setBoundaryNodes() ;
}

// -------------------- Set Neighborhood Vertex Number --------------------

void Mesh :: setNeighborhoodVertexNumber(){

	set_neighborhood_vertex_number[0] = &Mesh ::  setNeighborhoodVertexNumberP1 ;
	set_neighborhood_vertex_number[1] = &Mesh ::  setNeighborhoodVertexNumberP2 ;

	(this->*set_neighborhood_vertex_number[polynomial_degree - 1])() ;
}

// -------------------- Set Neighborhood Vertex Number P1 --------------------

void Mesh :: setNeighborhoodVertexNumberP1(){

	for(int i = 0 ; i < number_of_elements ; i ++){
	
		domain[i].neighborhood_1 = 0 ;
		domain[i].neighborhood_2 = 0 ;
		domain[i].neighborhood_3 = 0 ;
	}

	int node1 , node2 , node3 ;
	bool flag1 , flag2 , flag3 ;
	
	for (int i = 0 ; i < number_of_elements ; i ++){
		
		flag1 = flag2 = flag3 = false ;
		node1 = domain[i].node_1 ; 
		node2 = domain[i].node_2 ; 
		node3 = domain[i].node_3 ; 
	
		for (int j = 0 ; j < number_of_elements ; j ++ ){
			
			if (j == i){continue ;}

			if (((node1 == domain[j].node_1) || (node1 == domain[j].node_2) || (node1 == domain[j].node_3)) && ((node2 == domain[j].node_1) || (node2 == domain[j].node_2) || (node2 == domain[j].node_3))){
				
					domain[i].neighborhood_1 = domain[j].numeral ;
					flag1 = true ;
				}

			if (((node1 == domain[j].node_1) || (node1 == domain[j].node_2) || (node1 == domain[j].node_3)) && ((node3 == domain[j].node_1) || (node3 == domain[j].node_2) || (node3 == domain[j].node_3))){
				
					domain[i].neighborhood_3 = domain[j].numeral ;
					flag2 = true ;
				}

			if (((node2 == domain[j].node_1) || (node2 == domain[j].node_2) || (node2 == domain[j].node_3)) && ((node3 == domain[j].node_1) || (node3 == domain[j].node_2) || (node3 == domain[j].node_3))){
				
					domain[i].neighborhood_2 = domain[j].numeral ;
					flag3 = true ;
				}
			}
	
			domain[i].number_of_neighborhoods = 1 ;
			if((flag1 && flag2) || (flag1 && flag3) || (flag2 && flag3)){domain[i].number_of_neighborhoods = 2 ;}
			if(flag1 && flag2 && flag3){domain[i].number_of_neighborhoods = 3 ;}		
	}
}

// -------------------- Set Neighborhood Vertex Number P2 --------------------

void Mesh :: setNeighborhoodVertexNumberP2(){

	for(int i = 0 ; i < number_of_elements ; i ++){

		domain[i].neighborhood_1 = 0 ;
		domain[i].neighborhood_2 = 0 ;
		domain[i].neighborhood_3 = 0 ;	
	}

	int node1 , node2 , node3 ;
	bool flag1 , flag2 , flag3 ;
	
	for (int i = 0 ; i < number_of_elements ; i ++){
		
		flag1 = flag2 = flag3 = false ;
		node1 = domain[i].node_1 ; 
		node2 = domain[i].node_2 ; 
		node3 = domain[i].node_3 ;
	
		for (int j = 0 ; j < number_of_elements ; j ++ ){
			
			if (j == i){continue ;}

			if (((node1 == domain[j].node_1) || (node1 == domain[j].node_2) || (node1 == domain[j].node_3)) && ((node2 == domain[j].node_1) || (node2 == domain[j].node_2) || (node2 == domain[j].node_3))){
				
					domain[i].neighborhood_1 = domain[j].numeral ;
					flag1 = true ;
				}

			if (((node1 == domain[j].node_1) || (node1 == domain[j].node_2) || (node1 == domain[j].node_3)) && ((node3 == domain[j].node_1) || (node3 == domain[j].node_2) || (node3 == domain[j].node_3))){
				
					domain[i].neighborhood_3 = domain[j].numeral ;
					flag2 = true ;
				}

			if (((node2 == domain[j].node_1) || (node2 == domain[j].node_2) || (node2 == domain[j].node_3)) && ((node3 == domain[j].node_1) || (node3 == domain[j].node_2) || (node3 == domain[j].node_3))){
				
					domain[i].neighborhood_2 = domain[j].numeral ;
					flag3 = true ;
				}
			}

			domain[i].number_of_neighborhoods = 1 ;
			if((flag1 && flag2) || (flag1 && flag3) || (flag2 && flag3)){domain[i].number_of_neighborhoods = 2 ;}
			if(flag1 && flag2 && flag3){domain[i].number_of_neighborhoods = 3 ;}
		}
}

// -------------------- Set Boundary Nodes --------------------

void Mesh :: setBoundaryNodes(){

	set_boundary_nodes[0] = &Mesh :: setBoundaryNodesP1 ;
	set_boundary_nodes[1] = &Mesh :: setBoundaryNodesP2 ;

	(this->*set_boundary_nodes[polynomial_degree - 1])() ;
	
	sort (boundary_nodes , boundary_nodes + boundary_nodes_counter , compareBoundary) ;
	resetBoundaryNodes() ;
}

// -------------------- Set Boundary Nodes P1 --------------------

void Mesh :: setBoundaryNodesP1(){

	boundary_nodes_counter = 0 ;

	for (int i = 0 ; i < number_of_elements ; i ++){
	
		if (domain[i].number_of_neighborhoods == 1){ boundary_nodes_counter = boundary_nodes_counter + 3 ; continue ;}
		if (domain[i].number_of_neighborhoods == 2){ boundary_nodes_counter = boundary_nodes_counter + 2 ; continue ;}
	}

	boundary_nodes = new Boundary[boundary_nodes_counter] ; // double numbering

	int counter = 0 ;

	for (int i = 0 ; i < number_of_elements ; i ++){
	
		if (domain[i].number_of_neighborhoods == 3){ continue ;}
			
		if (domain[i].number_of_neighborhoods == 2){
			
			if (domain[i].neighborhood_1 == 0){
				
					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;
					
					counter ++ ;
				}

			if (domain[i].neighborhood_2 == 0){
					
					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;
					
					counter ++ ;
				}

			if (domain[i].neighborhood_3 == 0){

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;
					
					counter ++ ;
				}
			}

		if (domain[i].number_of_neighborhoods == 1){
					
					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;
					
					counter ++ ;
			
					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;
					
					counter ++ ;
		}
	}
}

// -------------------- Set Boundary Nodes P2 --------------------

void Mesh :: setBoundaryNodesP2(){

	boundary_nodes_counter = 0 ;

	for (int i = 0 ; i < number_of_elements ; i ++){
	
		if (domain[i].number_of_neighborhoods == 1){ boundary_nodes_counter = boundary_nodes_counter + 5 ; continue ;}
		if (domain[i].number_of_neighborhoods == 2){ boundary_nodes_counter = boundary_nodes_counter + 3 ; continue ;}
	}

	boundary_nodes = new Boundary[boundary_nodes_counter] ;

	int counter = 0 ;

	for (int i = 0 ; i < number_of_elements; i ++){
	
		if (domain[i].number_of_neighborhoods == 3){ continue ;}
		
		if (domain[i].number_of_neighborhoods == 2){
			
			if (domain[i].neighborhood_1 == 0){
					
					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_6 ;
					boundary_nodes[counter].x_coordinate = domain[i].x6_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y6_global ;
					
					counter ++ ;
				}

			if (domain[i].neighborhood_2 == 0){
					
					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_4 ;
					boundary_nodes[counter].x_coordinate = domain[i].x4_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y4_global ;

					counter ++ ;
				}

			if (domain[i].neighborhood_3 == 0){

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_5 ;
					boundary_nodes[counter].x_coordinate = domain[i].x5_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y5_global ;
					
					counter ++ ;
				}
			}

		if (domain[i].number_of_neighborhoods == 1 && domain[i].neighborhood_1 != 0){
					
					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;
					
					counter ++ ;
			
					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_4 ;
					boundary_nodes[counter].x_coordinate = domain[i].x4_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y4_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_5 ;
					boundary_nodes[counter].x_coordinate = domain[i].x5_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y5_global ;
					
					counter ++ ;
			}

		if (domain[i].number_of_neighborhoods == 1 && domain[i].neighborhood_2 != 0){
					
					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;
					
					counter ++ ;
			
					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_5 ;
					boundary_nodes[counter].x_coordinate = domain[i].x5_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y5_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_6 ;
					boundary_nodes[counter].x_coordinate = domain[i].x6_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y6_global ;
					
					counter ++ ;
			}
				
		if (domain[i].number_of_neighborhoods == 1 && domain[i].neighborhood_3 != 0){
					
					boundary_nodes[counter].node_number = domain[i].node_1 ;
					boundary_nodes[counter].x_coordinate = domain[i].x1_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y1_global ;
					
					counter ++ ;
			
					boundary_nodes[counter].node_number = domain[i].node_2 ;
					boundary_nodes[counter].x_coordinate = domain[i].x2_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y2_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_3 ;
					boundary_nodes[counter].x_coordinate = domain[i].x3_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y3_global ;
					
					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_4 ;
					boundary_nodes[counter].x_coordinate = domain[i].x4_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y4_global ;

					counter ++ ;

					boundary_nodes[counter].node_number = domain[i].node_6 ;
					boundary_nodes[counter].x_coordinate = domain[i].x6_global ;
					boundary_nodes[counter].y_coordinate = domain[i].y6_global ;
					
					counter ++ ;
			}
	}
}

// -------------------- Reset Boundary Nodes --------------------

void Mesh :: resetBoundaryNodes(){

	int first = 0 ;
	int second = 0 ;
	int new_counter = 1 ;

	for (int i = 0 ; i < boundary_nodes_counter - 1 ; i ++){
		
		first = boundary_nodes[i].node_number ;
		second = boundary_nodes[i + 1].node_number ;
		
		if (second > first){ new_counter ++ ;}
	}

	Boundary *Reset ;
	Reset = new Boundary [new_counter] ;
	
	int counter = 0 ;
	Reset[counter].node_number = boundary_nodes[0].node_number ;
	Reset[counter].x_coordinate = boundary_nodes[0].x_coordinate ;
	Reset[counter].y_coordinate = boundary_nodes[0].y_coordinate ;
	
	counter ++ ;

	for (int i = 0 ; i < boundary_nodes_counter - 1 ; i ++){
		
		first = boundary_nodes[i].node_number ;
		second = boundary_nodes[i + 1].node_number ;
		
		if (second > first){
			
			Reset[counter].node_number = boundary_nodes[i + 1].node_number ;
			Reset[counter].x_coordinate = boundary_nodes[i + 1].x_coordinate ;
			Reset[counter].y_coordinate = boundary_nodes[i + 1].y_coordinate ;
			
			counter ++ ;
		}
	}

	delete [] boundary_nodes ;
	
	boundary_nodes_counter = new_counter ;
	boundary_nodes = new Boundary [boundary_nodes_counter] ;

	for (int i = 0 ; i < boundary_nodes_counter ; i ++){
		
		boundary_nodes[i].node_number = Reset[i].node_number ;
		boundary_nodes[i].x_coordinate = Reset[i].x_coordinate ;
		boundary_nodes[i].y_coordinate = Reset[i].y_coordinate ;
	}

	delete [] Reset ;
}

// -------------------- Print Velocity Grid --------------------

void Mesh :: printVelocityGrid(){

	int *storage ;
	storage = new int[velocity_node_number] ;
	for (int i = 0 ; i < velocity_node_number ; i ++ ){

		storage[i] = 0 ;
	}

	int flag = 0 ;
	int counter = 0 ;
	//cout << domain[0].node_1 << "  Annie" <<endl ;
	int i = 0 ;
	
	for (int j = 0 ; j < number_of_elements ; j ++){
	
		flag = domain[j].node_1 ;
		if (checkExistence( flag , counter + 1 , storage)){

			storage[counter] = flag ;
			counter ++ ;
			velocity[i].node_number = domain[j].node_1 ;
			velocity[i].x_coordinate = domain[j].x1_global ;
			velocity[i].y_coordinate = domain[j].y1_global ;
			
			i ++ ;
		}
		
		flag = domain[j].node_2 ;
		if (checkExistence(flag  , counter + 1 , storage)){
			
			storage[counter] = flag ;
			counter ++ ;
			velocity[i].node_number = domain[j].node_2 ;
			velocity[i].x_coordinate = domain[j].x2_global ;
			velocity[i].y_coordinate = domain[j].y2_global ;
			i ++ ;
		}

		flag = domain[j].node_3 ;
		if (checkExistence(flag  , counter + 1 , storage)){

			storage[counter] = flag ;
			counter ++ ;
			velocity[i].node_number = domain[j].node_3 ;
			velocity[i].x_coordinate = domain[j].x3_global ;
			velocity[i].y_coordinate = domain[j].y3_global ;
			i ++ ;
		}

		flag = domain[j].node_4 ;
		if (checkExistence(flag  , counter + 1 , storage)){

			storage[counter] = flag ;
			counter ++ ;
			velocity[i].node_number = domain[j].node_4 ;
			velocity[i].x_coordinate = domain[j].x4_global ;
			velocity[i].y_coordinate = domain[j].y4_global ;
			i ++ ;
		}

		flag = domain[j].node_5 ;
		if (checkExistence(flag  , counter + 1 , storage)){

			storage[counter] = flag ;
			counter ++ ;
			velocity[i].node_number = domain[j].node_5 ;
			velocity[i].x_coordinate = domain[j].x5_global ;
			velocity[i].y_coordinate = domain[j].y5_global ;
			i ++ ;
		}
		
		flag = domain[j].node_6 ;
		if (checkExistence(flag  , counter + 1 , storage)){

			storage[counter] = flag ;
			counter ++ ;
			velocity[i].node_number = domain[j].node_6 ;
			velocity[i].x_coordinate = domain[j].x6_global ;
			velocity[i].y_coordinate = domain[j].y6_global ;
			i ++ ;
		}
	}

	delete [] storage ;
	return ;
}

// -------------------- Check Existence--------------------

bool Mesh :: checkExistence(int flag  , int counter , int S[]){

	for (int i = 0 ; i < counter  ; i ++){
		
		if(S[i] == flag){return false ;}
	}
	return true ;
}

// -------------------- Print Boundary Nodes --------------------

void Mesh :: printBoundaryNodes(){

	ofstream fout("BoundaryNodes.txt") ;
	
	for (int i = 0 ; i < boundary_nodes_counter ; i ++){
		
		fout << "Node Number# " << boundary_nodes[i].node_number << "  " << "( " << boundary_nodes[i].x_coordinate << " , " << boundary_nodes[i].y_coordinate << " )" << endl ;
	}
}

// -------------------- Print Element Information --------------------

void Mesh :: printElementInformation(){

	ofstream fout1 ("ElementInformation.txt") ;
	
	for( int i = 0 ; i < number_of_elements ; i ++){
			
		fout1 << " Element Number = " << domain[i].numeral << endl ;
		fout1 << endl ;
		fout1 << " Node#1 = " << domain[i].node_1 << "  "<< "( " << domain[i].x1_global <<" , "<< domain[i].y1_global <<" )" << endl ;
		fout1 << " Node#2 = " << domain[i].node_2 << "  "<< "( " << domain[i].x2_global <<" , "<< domain[i].y2_global <<" )" << endl ;
		fout1 << " Node#3 = " << domain[i].node_3 << "  "<< "( " << domain[i].x3_global <<" , "<< domain[i].y3_global <<" )" << endl ;
	//	fout1 << " Node#4 = " << domain[i].node_4 << "  "<< "( " << domain[i].x4_global <<" , "<< domain[i].y4_global <<" )" << endl ;
	//	fout1 << " Node#5 = " << domain[i].node_5 << "  "<< "( " << domain[i].x5_global <<" , "<< domain[i].y5_global <<" )" << endl ;
	//	fout1 << " Node#6 = " << domain[i].node_6 << "  "<< "( " << domain[i].x6_global <<" , "<< domain[i].y6_global <<" )" << endl ;
		fout1 << endl ;
		fout1 << " ---------------------------------------------------------- " << endl ;
		fout1 << endl ;
	}
}

// -------------------- Print Elements Using TecPlot --------------------

void Mesh :: printElementsUsingTecPlot(){
	
	// Note: the function is designed for P1-P1

	ofstream fout("Figure.dat") ;

	fout << " VARIABLES = X , Y" << endl ;
	fout << "ZONE N = " << pressure_node_number << ",E=" << number_of_elements << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 

	double stuff , x_mesh , y_mesh ;
	int rouf , index1 , index2 , index3 ;
	
	ifstream fin (Name.c_str());

	fin >> stuff ; 
	fin >> stuff ;
	
	for (int i = 0 ; i < pressure_node_number ; i ++){
		
		fin >> x_mesh ;  
		fin >> y_mesh ;
		fin >> stuff ;

		fout << x_mesh << " " << y_mesh << endl ;
	}
	
	for (int i = 0 ; i < number_of_elements ; i ++){

		fin >> rouf ;
		fin >> index1 ;
		fin >> index2 ;
		fin >> index3 ;

//================
// General

	fout << index1 + 1 - N_center << " " << index2 + 1 - N_center << " " << index3 + 1 - N_center << endl ;

//================

//***************//

//================
// Circles

	// fout << index1 << " " << index2 << " " << index3 << endl ;

//================

	}

	fin.close() ;
}

// -------------------- Return X Coordinate Velocity --------------------

double Mesh :: getXCoordinateVelocity(int i){
	
	for (int j = 0 ; j < velocity_node_number ; j ++){

		if (velocity[j].node_number == i){return velocity[j].x_coordinate ;}
	}
}

// -------------------- Return Y Coordinate Velocity --------------------

double Mesh :: getYCoordinateVelocity(int i){
		
	for (int j = 0 ; j < velocity_node_number ; j ++){

		if (velocity[j].node_number == i){return velocity[j].y_coordinate ;}
	}
}

// -------------------- Find Node Number--------------------

int Mesh :: findNodeNumber(double x , double y){

	for (int i = 0 ; i < velocity_node_number ; i ++){
		
		if (velocity[i].x_coordinate == x && velocity[i].y_coordinate == y){return velocity[i].node_number ; }
	}
}

// -------------------- Print Velocity Nodes Coordinate --------------------

void Mesh :: printVelocityCoordinates(){

	for (int i = 0 ; i < velocity_node_number ; i ++){
		
		cout << "Number  " << i + 1 << "  " << getXCoordinateVelocity(i + 1) <<"  " << getYCoordinateVelocity(i + 1) << endl ;
		//cout << velocity[i].node_number << "  " << velocity[i].x_coordinate << "  " << velocity[i].y_coordinate << endl ; 
		//put real node number to get real x coordinate and y coordinate, there is no necessiate to input i-1 or something else.
	}
}

// -------------------- Destroy --------------------

void Mesh :: Destroy(){
	
	delete [] domain ;
}

// -------------------- Distructor --------------------

Mesh :: ~Mesh(){

	//cout << "The Mesh Class is destroyed." << endl ;
}

// ===============================
// === Out of Class Functions  ===
// ===============================

// -------------------- Make a Comparison (Mapp Version) --------------------

bool compareBoundary(const Boundary & a , const Boundary & b){

	return a.node_number < b.node_number ;	
}