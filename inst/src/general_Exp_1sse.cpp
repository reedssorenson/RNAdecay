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

	int set = t.size();
	Type obj = 0.0;
	for (int i=0; i<set; i++){
	   obj += pow(m[i] - exp(-a1*t[i]),2);
	}

	return obj;
}
