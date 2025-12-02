#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "myheader.h"

void bond_set();
void bookkeep_alloc(int size);
void outproduct (double a[3], double b[3], double x[3]);
double veclength(double x[3]);

void bookkeep()
{
 START:
  if (book.alloc) {
    bookkeep_alloc(book.nbook_new);
  }

  int ix, iy, iz, nrpx, nrpy, nrpz;
  double r2;
  if (rcut>0) { book.frc = ((double)rcut_f+frcmar)*1e-10; }
  //if (rcut>0) { book.frc = (rcut+frcmar)*1e-10; }
  book.frc2 = book.frc * book.frc;
  frc_f = book.frc*1e10;
  if (book.algo == 3) { return; }
  nrpx = 0; nrpy = 0; nrpz = 0;

  double vol = cell.hmat[0][0]*cell.hmat[1][1]*cell.hmat[2][2] + cell.hmat[0][1]*cell.hmat[1][2]*cell.hmat[2][0] 
    + cell.hmat[0][2]*cell.hmat[1][0]*cell.hmat[2][1] - cell.hmat[0][0]*cell.hmat[1][2]*cell.hmat[2][1] 
    - cell.hmat[0][1]*cell.hmat[1][0]*cell.hmat[2][2] - cell.hmat[0][2]*cell.hmat[1][1]*cell.hmat[2][0];

  double area;

  if (cell.pbcx)
  { 
	double b[3] = {cell.hmat[1][0], cell.hmat[1][1], cell.hmat[1][2]};
	double c[3] = {cell.hmat[2][0], cell.hmat[2][1], cell.hmat[2][2]};
	double a[3]; outproduct(b, c, a);
	area = veclength(a);
	nrpx=(int)(book.frc*2/(vol/area)+1);
  }

  if (cell.pbcy)
  { 
	double c[3] = {cell.hmat[2][0], cell.hmat[2][1], cell.hmat[2][2]};
	double a[3] = {cell.hmat[0][0], cell.hmat[0][1], cell.hmat[0][2]};
	double b[3]; outproduct(c, a, b);
	area = veclength(b);
	nrpy=(int)(book.frc*2/(vol/area)+1);
  }

  if (cell.pbcz) 
  { 
	double a[3] = {cell.hmat[0][0], cell.hmat[0][1], cell.hmat[0][2]};
	double b[3] = {cell.hmat[1][0], cell.hmat[1][1], cell.hmat[1][2]};
	double c[3]; outproduct(a, b, c);
	area = veclength(c);
	nrpz=(int)(book.frc*2/(vol/area)+1); 
  }

  //std::cout<<"Bookkeep: Range= "<<nrpx<<" x "<<nrpy<<" x "<<nrpz<<std::endl;
  for (int i=1; i<=atom.natom; i++) {
    book.alistnum[i] = 0;
  for (int j=1; j<=atom.natom; j++) {
    if ((nrpx<=1)&&(nrpy<=1)&&(nrpz<=1)) {
      if (atom.Dist2Closest(i,j,ix,iy,iz)<book.frc2) {
	if (atom.Dist2Closest(i,j,ix,iy,iz)>1.0e-30) {
	  book.alistnum[i]++;
	  book.alist[i][book.alistnum[i]][0] = j;
	  book.alist[i][book.alistnum[i]][1] = ix;
	  book.alist[i][book.alistnum[i]][2] = iy;
	  book.alist[i][book.alistnum[i]][3] = iz;
	  //	  std::cout<<"hoge"<<std::endl;
	}
      }
    } else {
      for (int ixx=-nrpx; ixx<=nrpx; ixx++) {
      for (int iyy=-nrpy; iyy<=nrpy; iyy++) {
      for (int izz=-nrpz; izz<=nrpz; izz++) {
	if (atom.Dist2(i,j,ixx,iyy,izz)<book.frc2) {
	  if (atom.Dist2(i,j,ixx,iyy,izz)>1.0e-30) {
	    book.alistnum[i]++;
	    if (book.alistnum[i]>=book.nbook) { goto OUT; }
	    book.alist[i][book.alistnum[i]][0] = j;
	    book.alist[i][book.alistnum[i]][1] = ixx;
	    book.alist[i][book.alistnum[i]][2] = iyy;
	    book.alist[i][book.alistnum[i]][3] = izz;
	  }
	}
      }
      }
      }
    }
  OUT:
    //if (book.alistnum[i]>=NBOOK) {
    if (book.alistnum[i]>=book.nbook) {
      //std::cout<<"Bookkeep Error: NBOOK too small"<<std::endl;
      book.nbook_new = book.nbook*1.5;
      book.alloc = true;
      std::cout<<"Bookkeep: Reallocate arrays (from ";
      std::cout<<book.nbook<<" to "<<book.nbook_new<<")"<<std::endl;
      goto START;
    }
  }
  }
  if (book.alloc) goto START;
  bond_set();
}

void bookkeep_alloc(int size)
{
  //alist[NMAX+1][NBOOK+1][4];
  if (book.alist) {
    for (int i=0; i<=book.natom; i++) {
      for (int j=0; j<=book.nbook; j++) {
	delete[] book.alist[i][j]; }
      delete[] book.alist[i]; }
    delete[] book.alist; book.alist = NULL;
  }
  book.natom = atom.natom;
  book.alist = new int**[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) {
    book.alist[i] = new int*[size+1];
    for (int j=0; j<=size; j++) {
      book.alist[i][j] = new int[4];
    }
  }
  book.nbook = size;
  book.alloc = false;
}

/*    
int search_neighbor(i, j, &ix, &iy, &iz)
{
  double sdx, sdy, sdz;
  sdx = atom.rx[i] - atom.rx[j];
  sdy = atom.ry[i] - atom.ry[j];
  sdz = atom.rz[i] - atom.rz[j];
}
*/
