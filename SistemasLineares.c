#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "SistemasLineares.h"
#include <float.h>

#define EPS 1.0e-4


SistLinear_t* trocaLinhas(SistLinear_t *SL, int linha_pivo, int coluna);
int encontraMax (SistLinear_t *SL, int n,int linha, int coluna);
real_t normaMAX (real_t *next, real_t *prev, unsigned int n);
/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{

}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{

  tTotal[0] = timestamp();

  unsigned int n = SL->n;
  int linha_pivo, coluna, i ;
  real_t  soma, solucao[3];
  double m;


  for (coluna = 0; coluna < n-1; coluna++){
    linha_pivo = encontraMax(SL,n,coluna,coluna);
    
    if (linha_pivo != coluna){ //quando linha e coluna são iguais é pq é a diag principal
    	trocaLinhas(SL, linha_pivo, coluna);
    }

    for (i= coluna+1; i<n; i++){
    	m = SL->A[i][coluna] / SL->A[coluna][coluna];
     	SL->A[i][coluna] = 0.0;
      	for (int j = coluna+1; j<n; j++){
  			SL->A[i][j] = SL->A[i][j]-  m* SL->A[coluna][j];	
      	}
  		SL->b[i]-= SL->b[coluna]*m;		
    }

  }
  //retrosubtituição
  soma = 0;
  n = SL->n ;
  n--;
  x[n] = SL->b[n]/SL->A[n][n];
  for(int i = n - 1; i>= 0; i--){
  	soma = SL->b[i];
  	for (int j = i+1; j<= n;j++){
  		soma -= SL->A[i][j]*x[j];
  	}
  	x[i] = soma / SL->A[i][i];
  	//printf("%f\n",solucao[i] );
  }
  tTotal[1] = timestamp();

  printf("===>Eliminação de Gauss: %f ms \n",tTotal[1]-tTotal[0]);
  printf("-->X: ");
  prnVetor(x,SL->n);
  printf("-->Norma L2 do resíduo: \n");
  printf("\n");

}

int encontraMax(SistLinear_t *SL, int n,int linha, int coluna){
  real_t max;
  max = SL->A[linha][coluna];
  int linha_pivo = 0 ;


  for(int i = linha; i < n-1; ++i) {
    if (fabs(SL->A[i][coluna]) >= fabs(max)){
      max = SL->A[i][coluna];
      linha_pivo = i;
    }
  }
}//OK

SistLinear_t* trocaLinhas(SistLinear_t *SL, int linha_pivo, int coluna){
	double aux, tmp;

	for(int i =0; i< SL->n; i++){
        aux = SL->A[linha_pivo][i];
        SL->A[linha_pivo][i] = SL->A[coluna][i];
        SL->A[coluna][i] = aux;
    }

    tmp = SL->b[linha_pivo];
    SL->b[linha_pivo]= SL->b[coluna];
    SL->b[coluna] = tmp;

    //prnSistLinear(SL);
    //prnVetor(SL->b,SL->n);

    return SL;
}//OK

/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal)
{
  tTotal[0] = timestamp();

	int n = SL->n;
	int iteracoes =1;

	real_t *vetor_mult = (real_t *) malloc(n * sizeof(real_t));

	for(int i=0;i<n;i++){
		vetor_mult[i] =0;
	}

	while(iteracoes < MAXIT){
    for (int i = 0; i < SL->n; i++){
			real_t soma = 0;
			for(int j =0; j< SL->n ; j++){
				if (j!= i){
          soma = soma + SL->A[i][j]* -1;
					soma = soma * vetor_mult[j]; 
				}
      soma = soma + SL->b[i]; 
			x[i] = soma/SL->A[i][i];
			}
		}

    real_t erro = normaMAX(x,vetor_mult,n);
    if(erro > SL->erro){
      for(int k=0;k<SL->n;k++){
        vetor_mult[k]=x[k];
      }

    }else{
      tTotal[0] = timestamp();
      printf("===>Jacobi: %f ms --> %d Iterações \n",tTotal[1]-tTotal[0], iteracoes);
      printf("-->X: ");
      prnVetor(x,SL->n);
      printf("-->Norma L2 do resíduo: \n");
      printf("\n");
      return 0;
    }
    iteracoes++;
	}
  tTotal[1] = timestamp();
  printf("===>Jacobi: %f ms --> %d Iterações \n",tTotal[1]-tTotal[0], iteracoes);
  printf("-->X: ");
  prnVetor(x,SL->n);
  printf("-->Norma L2 do resíduo: \n");
  printf("\n");
  
}

real_t normaMAX (real_t *next, real_t *prev, unsigned int n)
{
  real_t maior = 0.0, sub;
  unsigned int i;
  sub = fabs(next[0]-prev[0]);
  real_t m = next[0];
  for (int j = 1; j< n-1;j++){
    if(next[j] > m){
      m = next[j];
    }
  }
  for (i = 1; i < n; ++i)
  {
    sub = fabs (next[i] - prev[i]);
    if (sub > maior)
      maior = sub/m;
  }
  return maior;
}; //TERMINADO

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal)
{
	tTotal[0] = timestamp();
  int n = SL->n;  //dimensão da matriz
  int iteracoes =1;

  real_t *tmp = (real_t *) malloc(n * sizeof(real_t));

  for(int i=0;i<n;i++){
    tmp[i] =0;
  }

  while(iteracoes < MAXIT){
    for (int i = 0; i < SL->n; i++){
      real_t soma = 0;
      for(int j =0; j< SL->n ; j++){
        if (j!= i){
          soma = soma + SL->A[i][j]* -1;
          soma = soma * tmp[j]; 
        }
      soma = soma + SL->b[i]; 
      x[i] = soma/SL->A[i][i];
      tmp[i] = x[i];
      }
    }

    real_t erro = normaMAX(x,tmp,n);
    if(erro > SL->erro){
      for(int k=0;k<SL->n;k++){
        tmp[k]=x[k];
      }

    }else{
      tTotal[1] = timestamp();
      printf("===>Gauss Seidel: %f ms --> %d Iterações \n",tTotal[1]-tTotal[0], iteracoes);
      printf("-->X: ");
      prnVetor(x,SL->n);
      printf("-->Norma L2 do resíduo: \n");
      printf("\n");
      return 0;
    }
    iteracoes++;
  }
  tTotal[1] = timestamp();
  printf("===>Jacobi: %f ms --> %d Iterações \n",tTotal[1]-tTotal[0], iteracoes);
  printf("-->X: ");
  prnVetor(x,SL->n);
  printf("-->Norma L2 do resíduo: \n");
  printf("\n");
}





/*!
  \brief Método de Refinamento
  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal)
{


}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (unsigned int n)
{


  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  if ( SL ) 
  {
    SL->A = (real_t**) malloc(n * sizeof(real_t*)); //aloca coluna
    for(int i=0; i<n; i++){
      SL->A[i] = (real_t *) malloc(n * sizeof(real_t));
    }
    SL->b = (real_t *) malloc(n * sizeof(real_t));     //VETOR DE TERMOS INDEPENDENTES
    if (!(SL->A) || !(SL->b) )
      liberaSistLinear (SL);
  }


  return (SL);

} //OK

/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL)
{
  free(SL->A);
  free(SL->b);
  free(SL);

} //OK

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear ()
{
  int n;
  float erro;
  
  scanf("%d", &n);
  //printf("DIMENSÃO: %d\n",n );
  scanf("%f", &erro);
 // printf("O ERRO: %f\n", erro);

  SistLinear_t *SL = alocaSistLinear (n);
  SL->n = n;
  SL->erro = erro;

  for (int lin = 0; lin<n; lin++){
  	//printf("primeiro for\n");
    for (int col=0; col<n ; col++){
    	//rintf("segundo for\n");
      	scanf("%f", &SL->A[lin][col]);
      	//printf("%*.1f",5,SL->A[lin][col]);

    }
    //printf("\n");
  }
  for (int j = 0; j <n; ++j ){
    scanf("%f", &SL->b[j]);
  }

  
  if (SL){
    //prnSistLinear(SL);
    //prnVetor(SL->b,n);
    return SL;
    //printf("Terminou leitura\n");
  }else

  printf("NULL");
    return NULL;
} //OK

// Exibe SL na saída padrão
void prnSistLinear (SistLinear_t *SL)
{

  printf("Imprimindo o sistema linear:\n");
  for (int lin = 0; lin < SL->n; lin++){
    for (int col=0; col < SL->n ; col++){
      printf("%f ", SL->A[lin][col]);
    }
    printf("\n");
  }

} //OK

// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n)
{
  //printf("Imprimindo o vetor:\n");

  for (int i = 0; i < n; i++){
    printf("%f ", v[i] );
  }
  printf("\n");

} //OK

