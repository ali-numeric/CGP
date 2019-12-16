// === Code Name: CGP/FEM  Project                                                            
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                          

// === Flow Past A Cylinder
// === Class Name: Cylinder
// === cylinder.h
// === Description:
// The class developes various operators for a master P1_P1 element.  

#ifndef SHAPE_H
#define SHAPE_H

#include "mesh.h"
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std ;

struct Function{

	public :

		double *function ;
		void f(){ function = new double [7] ;}
};

class Shape{

	private :

		int velocity_node_number , pressure_node_number ; 
		double x1_local , y1_local , x2_local , y2_local , x3_local , y3_local ;
		double x4_local , y4_local , x5_local , y5_local , x6_local , y6_local ;
		double w1 , w2 , w3 , a1 , a2 , a3 , b2 , b3 ;
		double alfa1 , alfa2 , alfa3 ;
		double beta1 , beta2 , beta3 ;
		double gama1 , gama2 , gama3 ;
		double x1_global , y1_global , x2_global , y2_global , x3_global , y3_global ;
		double x4_global , y4_global , x5_global , y5_global , x6_global , y6_global ;
		int node_1 , node_2 , node_3 , node_4 , node_5 , node_6 ; //global node numberes ;
		double jacobian ;

		int number_guass_point ;

		double *weight ;

		double ** guass_points ;
		double ** pressure_shape_function ;
		double ** pressure_shape_function_s_derivatives ;
		double ** pressure_shape_function_r_derivatives ;
		
		double *fi_1 , *fi_2 , *fi_3  ;
		double *fi_1_r , *fi_2_r , *fi_3_r ; 
		double *fi_1_s , *fi_2_s , *fi_3_s ;
		double *fi_1_x , *fi_2_x , *fi_3_x ;
		double *fi_1_y , *fi_2_y , *fi_3_y ;

		Function *fi_i ;
		Function *fi_i_x ;
		Function *fi_i_y ;

		double **velocity_mass_matrix ;
		double **velocity_laplace_operator_matrix ;
		double **pressure_laplace_operator_matrix ;

		double **H ;
		double **J ;
		
		double *force ;

		double critical_value ; // indicating real zero

		int *mapp ;

		double plot_4_5( double test ) ; 
		double plot_4_6( double test ) ;
		double plot_5_6( double test ) ;

	public :

		Shape( Element e ) ;
		void setWeight() ;
		void setGuassPoint() ;
		void setJacubian() ;
		void setCoeficient() ;
		void setShapeFunctions() ;
		void setLocalCoordinate() ;
		void setGlobalDerivative() ;
		void setVelocityMassMatrix() ;
		void setPressureMassMatrix() ;
		
		void setVelocityLaplaceOperatorMatrix() ;
		void setPressureLaplaceOperatorMatrix() ;
	
		void setH() ;
		void setJ() ;
		
		double integrate( Function f1 , Function f2 , Function f3) ;
		double integrate( Function f1 , Function f2) ;
		double integrate( Function f1 ) ;
		
		void setOperators() ;

		int mapping(int global_number) ;

		double getVelocityLaplace(int i , int j) ;
		double getPressureLaplace(int i , int j) ;		
		double getVelocityMass(int i , int j) ;
		double getH (int i , int j) ;
		double getJ (int i , int j) ;
		double getNonLinearU(int k , int i , int j ) ;
		double getNonLinearV(int k , int i , int j ) ;

		void Destroy() ;

		~Shape() ;
};

#endif
