#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
	double *tTotal = (double *) malloc(2 * sizeof(double));

	SistLinear_t *SL;
	SL = lerSistLinear();
	//eliminacaoGauss(SL,SL->b,tTotal);
	//gaussJacobi(SL,SL->b,tTotal);
	gaussSeidel(SL,SL->b,tTotal);
	//prnSistLinear(SL);
}
