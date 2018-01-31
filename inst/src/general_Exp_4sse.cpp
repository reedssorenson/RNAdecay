#include<TMB.hpp>
#include<math.h>
//C++ file defining the sum of squared errors.
//I will need one function for each of the two general models.

template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(t);
	DATA_VECTOR(m);
	PARAMETER(a1);
	PARAMETER(a2);
	PARAMETER(a3);
	PARAMETER(a4);

	int n = t.size();
	int set = n/4;
	Type obj = 0.0;
	for (int i=0; i<set; i++){
	   obj += pow(m[i] - exp(-a1*t[i]),2) + pow(m[set+i] - exp(-a2*t[i]),2) + pow(m[2*set+i] - exp(-a3*t[i]),2) + pow(m[3*set+i] - exp(-a4*t[i]),2);
	}

	return obj;
}
