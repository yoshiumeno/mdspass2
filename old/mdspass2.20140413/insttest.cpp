/*
 Test code for solving eigenvalue problems using LAPACK.
 To compile:
  g++ insttest.cpp -llapack -lblas
*/

/*
double precision :: AR(NM1,NM1),AI(NM1,NM1),EVR(NM1,MM),EVI(NM1,MM) &
  ,E(MM),FV1(NM1),FV2(NM1),FV3(NM1),FV4(NM1),FV5(NM1) &
  ,FV6(NM1),FV7(NM1),FV8(NM1),FM1(2,NM1),IV1(MM)

 call DGEEV('N','V',nn,ar,nn,e,ei,evi,nn,evr,nn,work,nn*4,INFO)
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

extern "C"
  long int dgeev_(char *jobvl, char *jobvr, long int *n, double *a, 
	     long int *lda, double *wr, double *wi, double *vl, 
	     long int *ldvl, double *vr, long int *ldvr, double *work, 
	     long int *lwork, long int *info);

int main(void)
{
  long int nn, nn4, info;
  nn = 3; nn4 = nn*nn*nn;
  //  double ar[3][3], e[3], ei[3], evi[3][3], evr[3][3], work[12];
  //  double *e, *ei, *evi, *evr, *work;
  
  //  double *ar = (double *)calloc(sizeof(double), nn * nn); 
  double *ar = (double *)calloc(nn*nn, sizeof(double)); 
  double *aro = (double *)calloc(nn*nn, sizeof(double)); 
  double *e  = (double *)calloc(nn, sizeof(double)); 
  double *ei = (double *)calloc(nn, sizeof(double)); 
  double *evi = (double *)calloc(nn*nn, sizeof(double)); 
  double *evr = (double *)calloc(nn*nn, sizeof(double)); 
  double *work = (double *)calloc(nn*4, sizeof(double)); 

  ar[0]=0.0; ar[3]=1.0; ar[6]=1.0;
  ar[1]=1.0; ar[4]=0.0; ar[7]=1.0;
  ar[2]=1.0; ar[5]=1.0; ar[8]=0.0;
  for (int i=0; i<9; i++) { aro[i]=ar[i]; }

  dgeev_("N","V",&nn,ar,&nn,e,ei,evi,&nn,evr,&nn,work,&nn4,&info);

  printf("eigenvalues\n");
  printf("%20.15f %20.15f %20.15f\n",e[0],e[1],e[2]);

  printf("eigenvector1\n");
  printf("%20.15f\n",evr[0]);
  printf("%20.15f\n",evr[1]);
  printf("%20.15f\n",evr[2]);
  printf("eigenvector2\n");
  printf("%20.15f\n",evr[3]);
  printf("%20.15f\n",evr[4]);
  printf("%20.15f\n",evr[5]);
  printf("eigenvector3\n");
  printf("%20.15f\n",evr[6]);
  printf("%20.15f\n",evr[7]);
  printf("%20.15f\n",evr[8]);

}
