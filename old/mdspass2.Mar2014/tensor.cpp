#include <string.h>
#include <iostream>
#include <GL/glui.h>
#include <fstream>
#include <math.h>
//#include <unistd.h>
#if defined __linux__ || defined __APPLE__
#include <unistd.h>
#else
#include <direct.h>
#endif
#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

void tensor_rotation(double s[3][3],  double x1[3],  double x2[3],  double x3[3],
		     double sp[3][3], double xp1[3], double xp2[3], double xp3[3]);
void resetvec(double a[3]);
void resetmat(double a[3][3]);

int main(int argc, char* argv[])
{
  double s[3][3], sp[3][3];
  double x1[3],  x2[3],  x3[3];
  double xp1[3], xp2[3], xp3[3];
  resetmat(s); resetmat(sp);
  resetvec(x1);  resetvec(x2);  resetvec(x3);
  resetvec(xp1); resetvec(xp2); resetvec(xp3);
  s[0][0] =  100;
  s[1][1] = -100;
  x1[0]=1;
  x2[1]=1;
  x3[2]=1;
  xp1[0]=1;
  xp2[1]=1;
  xp3[2]=1;
  for (double theta=0; theta <=360; theta += 5) {
    double rad = theta/180*M_PI;
    xp1[0] =  cos(rad); xp1[1] =  sin(rad);
    xp2[0] = -sin(rad); xp2[1] =  cos(rad);
    tensor_rotation(s,x1,x2,x3,sp,xp1,xp2,xp3);
    printf("%8.2f  ",theta);
    printf("%8.2f %8.2f %8.2f  ",sp[0][0],sp[0][1],sp[0][2]);
    printf("%8.2f %8.2f %8.2f  ",sp[1][0],sp[1][1],sp[1][2]);
    printf("%8.2f %8.2f %8.2f\n",sp[2][0],sp[2][1],sp[2][2]);
  }

}
