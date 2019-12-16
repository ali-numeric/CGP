// === Code Name: CGP/FEM  Project                                                            
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                          

// === Class Name: Shape														   
// === shape.cpp

#include "shape.h"

// -------------------- Constructor --------------------

Shape :: Shape(Element e){

		//Taylor and Hood (P2_P1)
		this->pressure_node_number = 3 ; //P1 element
		this->velocity_node_number = 3 ; //P1 element

		critical_value = pow(10.0 , -12.0) ;

		number_guass_point  = 7 ;

		x1_global = e.x1_global ; 
		y1_global = e.y1_global ;
		x2_global = e.x2_global ;
		y2_global = e.y2_global ;
		x3_global = e.x3_global ;
		y3_global = e.y3_global ;

		x4_global = e.x4_global ; 
		y4_global = e.y4_global ;
		x5_global = e.x5_global ;
		y5_global = e.y5_global ;
		x6_global = e.x6_global ;
		y6_global = e.y6_global ;

		node_1 = e.node_1 ;
		node_2 = e.node_2 ;
		node_3 = e.node_3 ;
		node_4 = e.node_4 ;
		node_5 = e.node_5 ;
		node_6 = e.node_6 ;
	

		guass_points = new double *[number_guass_point] ;
			for(int j = 0 ; j < number_guass_point ; j++ ){
				guass_points[j] = new double[3] ;
			}
		
		pressure_shape_function = new double *[number_guass_point] ;
			for(int j = 0 ; j < number_guass_point ; j++ ){
				pressure_shape_function[j] = new double[3] ;
			}
		
		pressure_shape_function_s_derivatives = new double *[number_guass_point] ;
			for(int j = 0 ; j < number_guass_point ; j++ ){
				pressure_shape_function_s_derivatives[j] = new double[3] ;
			}

		pressure_shape_function_r_derivatives = new double *[number_guass_point] ;
			for(int j = 0 ; j < number_guass_point ; j++ ){
				pressure_shape_function_r_derivatives[j] = new double[3] ;
			} 
		
		fi_1 = new double [ number_guass_point] ;
		fi_2 = new double [ number_guass_point] ;
		fi_3 = new double [ number_guass_point] ;

		fi_1_r = new double [ number_guass_point] ;
		fi_2_r = new double [ number_guass_point] ;
		fi_3_r = new double [ number_guass_point] ;

		fi_1_s = new double [ number_guass_point] ;
		fi_2_s = new double [ number_guass_point] ;
		fi_3_s = new double [ number_guass_point] ;

		fi_1_x = new double [ number_guass_point] ;
		fi_2_x = new double [ number_guass_point] ;
		fi_3_x = new double [ number_guass_point] ; 

		fi_1_y = new double [ number_guass_point] ;
		fi_2_y = new double [ number_guass_point] ;
		fi_3_y = new double [ number_guass_point] ;
		
		weight = new double [number_guass_point] ;

		fi_i = new Function[pressure_node_number] ;
		fi_i_x = new Function[pressure_node_number] ;
		fi_i_y = new Function[pressure_node_number] ;
		
		velocity_mass_matrix = new double *[velocity_node_number] ;
			for(int j = 0 ; j < velocity_node_number ; j++ ){
				velocity_mass_matrix [j] = new double[velocity_node_number] ;
			}
		
		pressure_laplace_operator_matrix = new double *[pressure_node_number] ;
			for(int j = 0 ; j < pressure_node_number ; j++ ){
				pressure_laplace_operator_matrix [j] = new double[pressure_node_number] ;
			}
		
		velocity_laplace_operator_matrix = new double *[velocity_node_number] ;
			for(int j = 0 ; j < velocity_node_number ; j++ ){
				velocity_laplace_operator_matrix [j] = new double[velocity_node_number] ;
			}

		H = new double *[velocity_node_number] ;  // For P2_P1, check the size of the matrix
			for(int j = 0 ; j < velocity_node_number ; j++ ){
				H[j] = new double[pressure_node_number] ;
			}

		J = new double *[velocity_node_number] ;  // For P2_P1, check the size of the matrix
			for(int j = 0 ; j < velocity_node_number ; j++ ){
				J[j] = new double[pressure_node_number] ;
			}

		mapp = new int[3] ;
}

// -------------------- Set Integration Weight --------------------

void Shape :: setWeight(){

	w1=0.225000000000000;
	w2=0.132394152788506;
	w3=0.125939180544827;


	weight[0] = w1 ; 
	weight [1] = weight [2] = weight [3] = w2 ;
	weight [4] = weight [5] = weight [6] = w3 ;

}

// -------------------- Set Guass Point --------------------

void Shape :: setGuassPoint(){

	a1=0.333333333333333;
	a2=0.059715871789770;
	a3=0.797426985353087;
	b2=(1-a2)/2;
	b3=(1-a3)/2;

	// guass_point [i][2] = 1 - r - s = 1 - eta - zeta = 1 - L2 - L3
	// guass_point [i][0] = r = zeta = L3
	// guass_point [i][1] = s = eta  = L2

	for (int i = 0 ; i < number_guass_point ; i ++ ){
		if (i == 0){
			guass_points [i][0] = a1 ;
			guass_points [i][1] = a1 ;
			guass_points [i][2] = 1.0 - guass_points [i][0] - guass_points [i][1] ;  
			}
		if (i == 1){
			guass_points [i][0] = a2 ; 
			guass_points [i][1] = b2 ;
			guass_points [i][2] = 1.0 - guass_points [i][0] - guass_points [i][1] ;
	
			}
		if (i == 2){
			guass_points [i][0] = b2 ; 
			guass_points [i][1] = a2 ;
			guass_points [i][2] = 1.0 - guass_points [i][0] - guass_points [i][1] ;
	
			}
		if (i == 3){
			guass_points [i][0] = b2 ; 
			guass_points [i][1] = b2 ;
			guass_points [i][2] = 1.0 - guass_points [i][0] - guass_points [i][1] ;
	
			}
		if (i == 4){
			guass_points [i][0] = a3 ; 
			guass_points [i][1] = b3 ;
			guass_points [i][2] = 1.0 - guass_points [i][0] - guass_points [i][1] ;

			}
		if (i == 5){
			guass_points [i][0] = b3 ; 
			guass_points [i][1] = a3 ;
			guass_points [i][2] = 1.0 - guass_points [i][0] - guass_points [i][1] ;
	
			}
		if (i == 6){
			guass_points [i][0] = b3 ; 
			guass_points [i][1] = b3 ;
			guass_points [i][2] = 1.0 - guass_points [i][0] - guass_points [i][1] ;
	
			}
		continue ;
	}
}

// ----- Set Shape Functions and Their Local Derivatives ------

void Shape :: setShapeFunctions(){

	setWeight() ;
	setGuassPoint() ; 
	
	for (int i = 0 ; i < number_guass_point ; i ++) { 

		pressure_shape_function[i][0] = guass_points[i][2];
		pressure_shape_function[i][1] = guass_points[i][0];
		pressure_shape_function[i][2] = guass_points[i][1];

		pressure_shape_function_r_derivatives[i][0] = -1.0 ;
		pressure_shape_function_r_derivatives[i][1] = 1.0 ;
		pressure_shape_function_r_derivatives[i][2] = 0.0 ;

		pressure_shape_function_s_derivatives[i][0] = -1.0 ;
		pressure_shape_function_s_derivatives[i][1] = 0.0 ;
		pressure_shape_function_s_derivatives[i][2] = 1.0 ;
	}

		for (int i = 0 ; i < number_guass_point ; i ++){

			fi_1[i] = pressure_shape_function[i][0] ; 
			fi_2[i] = pressure_shape_function[i][1] ; 
			fi_3[i] = pressure_shape_function[i][2] ;

			fi_1_r[i] = pressure_shape_function_r_derivatives[i][0] ;
			fi_2_r[i] = pressure_shape_function_r_derivatives[i][1] ;
			fi_3_r[i] = pressure_shape_function_r_derivatives[i][2] ;

			fi_1_s[i] = pressure_shape_function_s_derivatives[i][0] ;
			fi_2_s[i] = pressure_shape_function_s_derivatives[i][1] ; 
			fi_3_s[i] = pressure_shape_function_s_derivatives[i][2] ;
		}

		for (int i = 0 ; i < pressure_node_number ; i ++){

				fi_i[i].f() ;
				fi_i_x[i].f() ;
				fi_i_y[i].f() ;
			}
}

// -------------------- Set Local Coordinate --------------------

void Shape :: setLocalCoordinate() {

	x1_local = x1_global ; 
	y1_local = y1_global ;

	mapp[0] = node_1 ;

	if (((x2_global - x1_local)*(y3_global - y1_local) - (x3_global - x1_local)*(y2_global - y1_local)) >= 0 ){
		
		
		x2_local =  x2_global ;
		y2_local =  y2_global ; 
		x3_local =  x3_global ;
		y3_local =	y3_global ;

		mapp[1] = node_2 ;
		mapp[2] = node_3 ;
	}
	
	if (((x2_global - x1_local)*(y3_global - y1_local) - (x3_global - x1_local)*(y2_global - y1_local)) < 0 ){
	
		x2_local =  x3_global ;
		y2_local =  y3_global ; 
		x3_local =  x2_global ;
		y3_local =	y2_global ;

		mapp[1] = node_3 ;
		mapp[2] = node_2 ;

	}
}

// -------------------- Plot Line Between Point4 and Point5  --------------------

double Shape :: plot_4_5( double test ){

	return ((y4_global - y5_global)/(x4_global - x5_global))*(test - x4_global) + y4_global ; 
} 

// -------------------- Plot Line Between Point4 and Point6  --------------------

double Shape :: plot_4_6( double test ){

	return ((y4_global - y6_global)/(x4_global - x6_global))*(test - x4_global) + y4_global ; 
}

// -------------------- Plot Line Between Point5 and Point6 --------------------

double Shape :: plot_5_6( double test ){

	return ((y6_global - y5_global)/(x6_global - x5_global))*(test - x6_global) + y6_global ; 
}

// -------------------- Set Local Coeficient --------------------

void Shape ::setCoeficient() {

		setLocalCoordinate() ;

		alfa1 = x2_local*y3_local - x3_local*y2_local ;
		alfa2 = x3_local*y1_local - x1_local*y3_local ;
		alfa3 = x1_local*y2_local - x2_local*y1_local ;

		beta1 = y2_local - y3_local ;
		beta2 = y3_local - y1_local ;
		beta3 = y1_local - y2_local ;

		gama1 = - (x2_local - x3_local) ;
		gama2 = - (x3_local - x1_local) ; 
		gama3 = - (x1_local - x2_local) ;
	
}

// -------------------- Set Jacubian --------------------

void Shape :: setJacubian(){

	setCoeficient() ;
	jacobian = (beta2*gama3 - gama2*beta3) ;

}

// -------------------- Set Global Derivative--------------------

void Shape :: setGlobalDerivative(){
			
			setJacubian() ;
			setShapeFunctions() ;

			for (int i = 0 ; i < number_guass_point ; i ++){

				fi_1_x[i] = (1.0/jacobian)*((beta2 * fi_1_r[i]) + (beta3 * fi_1_s[i])) ;
				fi_2_x[i] = (1.0/jacobian)*((beta2 * fi_2_r[i]) + (beta3 * fi_2_s[i])) ; 
				fi_3_x[i] = (1.0/jacobian)*((beta2 * fi_3_r[i]) + (beta3 * fi_3_s[i])) ;

				fi_1_y[i] = (1.0/jacobian)*((gama2 * fi_1_r[i]) + (gama3 * fi_1_s[i])) ;
				fi_2_y[i] = (1.0/jacobian)*((gama2 * fi_2_r[i]) + (gama3 * fi_2_s[i])) ; 
				fi_3_y[i] = (1.0/jacobian)*((gama2 * fi_3_r[i]) + (gama3 * fi_3_s[i])) ;
	
				}

			for (int i = 0 ; i < this->number_guass_point ; i ++){

				fi_i[0].function[i] = this-> fi_1[i] ;
				fi_i[1].function[i] = this-> fi_2[i] ;
				fi_i[2].function[i] = this-> fi_3[i] ;

				fi_i_x[0].function[i] = this-> fi_1_x [i] ;
				fi_i_x[1].function[i] = this-> fi_2_x [i] ;
				fi_i_x[2].function[i] = this-> fi_3_x [i] ;

				fi_i_y[0].function[i] = this-> fi_1_y [i] ;
				fi_i_y[1].function[i] = this-> fi_2_y [i] ;
				fi_i_y[2].function[i] = this-> fi_3_y [i] ;
			
			}
}

// -------------------- Set Velocity Mass Matrix --------------------

void Shape :: setVelocityMassMatrix(){

	// Without considering the density

	for( int i = 0 ; i < velocity_node_number ; i ++){
		for( int j = 0 ; j < velocity_node_number ; j ++){

				velocity_mass_matrix[i][j] = integrate (fi_i[i] , fi_i[j]) ;
		}
	}
}

// -------------------- Mapping --------------------

int Shape :: mapping(int global_number){

	for(int i = 0 ; i < 3 ; i ++){

		if(mapp[i] == global_number){return i ;}
	
	}
	return 0 ;
}

// -------------------- Set Velocity Laplace Operator Matrix --------------------

void Shape :: setVelocityLaplaceOperatorMatrix(){

	// Without considering the viscosity

		for( int i = 0 ; i < velocity_node_number ; i ++){
			for( int j = 0 ; j < velocity_node_number ; j ++){
				
				velocity_laplace_operator_matrix[i][j] = -1.0*(integrate (fi_i_x[i] , fi_i_x[j]) + integrate (fi_i_y[i] , fi_i_y[j])) ;
		}
	}
}

// -------------------- Set Pressure Laplace Operator Matrix --------------------

void Shape :: setPressureLaplaceOperatorMatrix(){

		for( int i = 0 ; i < pressure_node_number ; i ++){
			for( int j = 0 ; j < pressure_node_number ; j ++){

				pressure_laplace_operator_matrix [i][j] = -1.0 * (integrate (fi_i_x[i] , fi_i_x[j]) + integrate (fi_i_y[i] , fi_i_y[j]));
		}
	}
}

// -------------------- Get NonLinear U --------------------

double Shape :: getNonLinearU(int k , int i , int j ){

	return integrate (fi_i[i] , fi_i_x[j] , fi_i[k]) ;
}

// -------------------- Get NonLinear V --------------------

double Shape :: getNonLinearV(int k , int i , int j ){
	
	return integrate (fi_i[i] , fi_i_y[j] , fi_i[k]) ; 
}

// -------------------- Set H Matrix --------------------

void Shape :: setH(){

		for( int j = 0 ; j < velocity_node_number ; j ++){
			for( int i = 0 ; i < pressure_node_number ; i ++){

				H[i][j] = integrate (fi_i_x[j] , fi_i[i]) ;
		}
	}
}

// -------------------- Set J Matrix --------------------

void Shape :: setJ(){

		for( int j = 0 ; j < velocity_node_number ; j ++){
			for( int i = 0 ; i < pressure_node_number ; i ++){

				J[i][j] = integrate (fi_i_y[j] , fi_i[i]) ;
		}
	}
}

// -------------------- Integration 3 arguments --------------------

double Shape :: integrate( Function f1 , Function f2 , Function f3){

	double integration = 0.0 ;
	double sigma = 0.0 ;

	for( int i = 0 ; i < number_guass_point ; i ++){

		sigma += (weight[i] * f1.function[i] * f2.function[i] * f3.function[i]) ;
			
		}
	
	integration = (0.5)* jacobian * sigma ;

	return integration ;
}

// -------------------- Integration 2 arguments --------------------

double Shape :: integrate( Function f1 , Function f2){

	double integration = 0.0 ;
	double sigma = 0.0 ;

	for( int i = 0 ; i < number_guass_point ; i ++){

			sigma += (weight[i] * f1.function[i] * f2.function[i]) ;
		}
	
	integration = (0.5)* jacobian * sigma ;

	return integration ;
}

// -------------------- Integration 1 argument --------------------

double Shape :: integrate( Function f1 ){

	double integration = 0.0 ;
	double sigma = 0.0 ;

	for( int i = 0 ; i < number_guass_point ; i ++){

			sigma += (weight[i] * f1.function[i]) ;
		}
	
	integration = (0.5)* jacobian * sigma ;
	return integration ;
}

// -------------------- Set Operators --------------------

void Shape :: setOperators(){
		
	setGlobalDerivative() ;

	setVelocityMassMatrix() ;
	setVelocityLaplaceOperatorMatrix() ;
	setPressureLaplaceOperatorMatrix() ;
	setH() ;
	setJ() ;
}

// -------------------- get H value-------------------

double Shape :: getH (int i , int j){
	
	return H[i][j] ;
}

// -------------------- get J value-------------------

double Shape :: getJ (int i , int j){
	
	return J[i][j] ;
}

// -------------------- get Velocity Laplace value -------------------

double Shape :: getVelocityLaplace(int i , int j){
	
	return velocity_laplace_operator_matrix[i][j] ;
}

// -------------------- get Pressure Laplace value -------------------

double Shape :: getPressureLaplace(int i , int j){

	return pressure_laplace_operator_matrix[i][j] ;
}

// -------------------- get Velocity Mass value -------------------

double Shape :: getVelocityMass(int i , int j){

	return velocity_mass_matrix[i][j] ;
}

// -------------------- Destruction --------------------

void Shape :: Destroy(){

	for(int j = 0 ; j < number_guass_point ; ++ j)
			
		delete [] guass_points[j] ;
		delete [] guass_points ;

		delete [] weight ;

	for (int j = 0 ; j < pressure_node_number ; j ++){
		
		delete [] velocity_mass_matrix[j] ;
		delete [] velocity_laplace_operator_matrix[j] ;
		delete [] pressure_laplace_operator_matrix[j] ;

		delete [] H[j] ;
		delete [] J[j] ;

	 }

		delete [] velocity_mass_matrix ;
		delete [] velocity_laplace_operator_matrix ;
		delete [] pressure_laplace_operator_matrix ;
		
		delete [] H ;
		delete [] J ;

		delete [] mapp ;

	for (int j = 0 ; j < number_guass_point ; j ++ ){
		
			delete [] pressure_shape_function[j] ;
			delete [] pressure_shape_function_s_derivatives[j] ;
			delete [] pressure_shape_function_r_derivatives[j] ;
		} 
		
	delete [] pressure_shape_function ;
	delete [] pressure_shape_function_s_derivatives ;
	delete [] pressure_shape_function_r_derivatives ;

	delete [] fi_1 , fi_2 , fi_3 ;
	delete [] fi_1_r , fi_2_r , fi_3_r ;
	delete [] fi_1_s , fi_2_s , fi_3_s ; 
	delete [] fi_1_x , fi_2_x , fi_3_x ; 
	delete [] fi_1_y , fi_2_y , fi_3_y ; 

	delete [] fi_i ;
	delete [] fi_i_x ;
	delete [] fi_i_y ;
}

// -------------------- Destruction --------------------

Shape :: ~Shape(){
//
}