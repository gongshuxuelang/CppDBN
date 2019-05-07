#include "head.h"

int power(int base, int exponent)
{
	int result = 1;	
	if(exponent == 0)
	{
		return result;
	}
	for (int i = 0; i < exponent; ++i)
	{
		result *= base;
	}
 
	return result;
}