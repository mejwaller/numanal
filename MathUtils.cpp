#include "MathUtils.hpp"

long double factorial(unsigned long int num)
{    
    //numbers > 170 cause overflow...
    if(num == 0)
	{
        return 1;
    }

    return (long double)num*(long double)factorial(num-1);
}
