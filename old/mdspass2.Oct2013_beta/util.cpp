#include <stdio.h>
#include <math.h>
#include "myheader.h"

void outproduct (double a[3], double b[3], double x[3]);

void inverse(double h[3][3], double hi[3][3])
{
  double deth
    = h[0][0]*h[1][1]*h[2][2] + h[0][1]*h[1][2]*h[2][0] 
    + h[0][2]*h[1][0]*h[2][1] - h[0][0]*h[1][2]*h[2][1] 
    - h[0][1]*h[1][0]*h[2][2] - h[0][2]*h[1][1]*h[2][0];
    
  hi[0][0]=(h[1][1]*h[2][2]-h[1][2]*h[2][1])/deth;
  hi[0][1]=(h[0][2]*h[2][1]-h[0][1]*h[2][2])/deth;
  hi[0][2]=(h[0][1]*h[1][2]-h[0][2]*h[1][1])/deth;
  hi[1][0]=(h[1][2]*h[2][0]-h[1][0]*h[2][2])/deth;
  hi[1][1]=(h[0][0]*h[2][2]-h[0][2]*h[2][0])/deth;
  hi[1][2]=(h[0][2]*h[1][0]-h[0][0]*h[1][2])/deth;
  hi[2][0]=(h[1][0]*h[2][1]-h[1][1]*h[2][0])/deth;
  hi[2][1]=(h[0][1]*h[2][0]-h[0][0]*h[2][1])/deth;
  hi[2][2]=(h[0][0]*h[1][1]-h[0][1]*h[1][0])/deth;
}

void resetmat(double a[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      a[i][j]=0; } }
}
void resetmat6(double a[6][6])
{
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      a[i][j]=0; } }
}

void transpose(double a[3][3])
{
  double x;
  for (int i=0; i<3; i++) {
    for (int j=i+1; j<3; j++) {
      x=a[i][j]; a[i][j]=a[j][i]; a[j][i]=x; } }
}
void transpose(double a[3][3], double b[3][3])
{
  double x;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b[i][j]=a[j][i]; } }
}

//c = product(a,b)
void matmul(double a[3][3], double b[3][3], double c[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      c[i][j]=0.0;
      for (int k=0; k<3; k++) {
	c[i][j] += a[i][k]*b[k][j];
      } } }
}
void matadd(double a[3][3], double b[3][3], double c[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      c[i][j]=a[i][j]+b[i][j];
    } }
}
void matsub(double a[3][3], double b[3][3], double c[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      c[i][j]=a[i][j]-b[i][j];
    } }
}
void matcpy(double a[3][3], double b[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b[i][j] = a[i][j];
    } }
}

void recipro(double h[3][3], double r[3][3])
{
  double a[3], b[3], c[3], x[3];
  for (int i=0; i<3; i++) {
    a[i] = h[i][0];
    b[i] = h[i][1];
    c[i] = h[i][2];
  }
  outproduct(b,c,x);
  for (int i=0; i<3; i++) {
    r[i][0] = x[i];
  }
  outproduct(c,a,x);
  for (int i=0; i<3; i++) {
    r[i][1] = x[i];
  }
  outproduct(a,b,x);
  for (int i=0; i<3; i++) {
    r[i][2] = x[i];
  }
}

void outproduct (double a[3], double b[3], double x[3])
{
  x[0] = a[1]*b[2]-a[2]*b[1];
  x[1] = a[2]*b[0]-a[0]*b[2];
  x[2] = a[0]*b[1]-a[1]*b[0];
}


double veclength(double x[3])
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}
