#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#define id(width,row,col) (width*col+row)

extern "C" struct size_m {
	int m;
	int n;
};


/*
	NEARZERO is an interpretation of "zero"
*/
extern "C" const double NEARZERO;

/*
	TOLERANCE is the epsilon for the convergence error
*/

//extern const double TOLERANCE = 1.0e-10;
extern "C" const double TOLERANCE;


#endif /*PARAMETERS_H_*/

