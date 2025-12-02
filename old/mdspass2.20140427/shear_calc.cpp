#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "myheader.h"
#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

extern int itolfor; extern float tolfor;

void stretch_celladjust(float c1x, float c1y, float c1z,
			float c2x, float c2y, float c2z,
			float c3x, float c3y, float c3z);
void potential();
void md(); void bookkeep();
void normal_strain_adjust(float sigx0, float sigy0, float sigz0, float range);

void shear_calc()
{

  //cell3x += 1.0;
  printf("cell3x = %f\n",cell3x);
  stretch_celladjust(cell1x,cell1y,cell1z,
		     cell2x,cell2y,cell2z,
		     cell3x,cell3y,cell3z);
  //normal_strain_adjust(0,0,0,10);
  if (mdmotion==0) { mdmotion = 1; }
  else { mdmotion = 0; }
}

void normal_strain_adjust(float sigx0, float sigy0, float sigz0, float range)
{
  
  ensemble = 2; 
  mdmotion = 1;
  //  itolfor = 1; tolfor = 0.001;
  float sigx,sigy,sigz;

  int icnt = 0;
  bookkeep();
  while (mdmotion == 1) {
    md(); icnt++; if (icnt>100) { mdmotion = 0; }
    sigx=cell.sgmmat[0][0]*1e-6; sigy=cell.sgmmat[1][1]*1e-6; sigz=cell.sgmmat[2][2]*1e-6; 
    if ((fabs(sigx-sigx0)<range)&&(fabs(sigy-sigy0)<range)&&(fabs(sigz-sigz0)<range)) {
      mdmotion = 0;
    }
  }
  printf("%d\n",icnt);
  printf("Normal stress (xx,yy,zz): %f %f %f\n",
	 cell.sgmmat[0][0]*1e-6,cell.sgmmat[1][1]*1e-6,cell.sgmmat[2][2]*1e-6);
  printf("Shear  stress (xy,yz,zx): %f %f %f\n",
	 cell.sgmmat[0][1]*1e-6,cell.sgmmat[1][2]*1e-6,cell.sgmmat[2][0]*1e-6);
  printf("Energy: %f\n",epotatom);

  mdmotion = 0;
}
