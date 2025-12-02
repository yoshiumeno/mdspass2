#include <iostream>
#include "myheader.h"

void stretch(double x, double y, double z)
{
  cell.hmat[0][0]=cell.hmat[0][0]*x;
  cell.hmat[1][1]=cell.hmat[1][1]*y;
  cell.hmat[2][2]=cell.hmat[2][2]*z;
  for (int i=1; i<=atom.natom; i++)
    {
      atom.rx[i]=atom.rx[i]*x;
      atom.ry[i]=atom.ry[i]*y;
      atom.rz[i]=atom.rz[i]*z;
    }
}

void stretch_celladjust(float x, float y, float z)
{
  double strx, stry, strz;
  strx = x*ang/cell.hmat[0][0];
  stry = y*ang/cell.hmat[1][1];
  strz = z*ang/cell.hmat[2][2];
  stretch(strx,stry,strz);
}
