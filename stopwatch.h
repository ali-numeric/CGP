// === Code Name: CGP/FEM  Project                                                           
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                          

// === Flow Past A Cylinder
// === Class Name: Cylinder
// === cylinder.h
// === Description:
// The class computes the time processing.								   

#ifndef StopWatch_H
#define StopWatch_H

#include <time.h>

class StopWatch{
	
	public :
		
		clock_t start ;
		clock_t stop ;
		clock_t full ;

		float final ;

		StopWatch(int random) ;
		~StopWatch();
		void startTime() ;
		void stopTime() ;
		void resetTime() ;
		int computeCPUTime() ;
		float finalTime() ;
};

#endif 
