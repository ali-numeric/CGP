// === Code Name: CGP/FEM  Project                                                            
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                          

// === Flow Past A Cylinder
// === Class Name: Cylinder
// === cylinder.h

#ifndef CYLINDER_H
#define CYLINDER_H

#include "computation.h"
#include <sstream>
#include <iostream>
using namespace std ;

class Cylinder : public Computation{

	private :

		double pi ;
		double dynamic_viscosiy ;
		double density ;

		double U_stream , H , L_Box , X_Center , Y_Center , D , X_S ;

	public :

		Cylinder(double Viscosity , double Density , double Time , double Time_Step , string S , double Alfa = 1.0 , double penalty_term = 0.0 , double Beta_q = 3.0/2.0 , int level = 0 , int n_center = 0);

		void initializeVelocity() ;
		void initializePressure() ;
		void initializeForce() ;
		
		void constructNuemannDirichletBounds() ;
		void constructNuemannDirichletBoundsVkMinus1() ;
		void constructNuemannDirichletBoundsVkMinus2() ;
		void constructNuemannDirichletBoundsVkMinus3() ;
		
		~Cylinder() ;
};

#endif