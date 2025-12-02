#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "myheader.h"

void inverse(double mat[][3], double imat[][3]);
int count_num_species(int nesp[50], char **aesp);

void writeconfig(const char* fname)
{
  printf ("Write config file = %s\n", fname);
  /*
  std::ofstream foutconfig( "config.end", std::ios::out );
  foutconfig << "ATOM_TYPE" << std::endl;
  foutconfig << "POTENTIAL" << std::endl;
  foutconfig << atom.natom << std::endl;
  */
  FILE *fp;
  double xx, yy, zz;
  double hinmat[3][3];
  char species[3];
  //  fp = fopen("config.end", "w");
  fp = fopen(fname, "w");
  strcpy(species,atom.asp[1]);
  if (atom.natom > 1) {
    for (int i=2; i<=atom.natom; i++) {
      if(strcmp(atom.asp[i], species)==1) {strcpy(species,"M");} } }
  fprintf(fp," %s\n",species);
  fprintf(fp," %s %s\n",atom.potential_func,atom.potential_arg);
  fprintf(fp," %d\n",atom.natom);
  fprintf(fp," %e\n",cell.alat);
  xx=cell.hmat[0][0]/cell.alat/1.0e-10; yy=cell.hmat[1][0]/cell.alat/1.0e-10; zz=cell.hmat[2][0]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][1]/cell.alat/1.0e-10; yy=cell.hmat[1][1]/cell.alat/1.0e-10; zz=cell.hmat[2][1]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][2]/cell.alat/1.0e-10; yy=cell.hmat[1][2]/cell.alat/1.0e-10; zz=cell.hmat[2][2]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  fprintf(fp," %s\n","fra");
  inverse(cell.hmat, hinmat);
  for (int i=1; i<=atom.natom; i++) {
    xx = hinmat[0][0]*atom.rx[i] + hinmat[0][1]*atom.ry[i] + hinmat[0][2]*atom.rz[i];
    yy = hinmat[1][0]*atom.rx[i] + hinmat[1][1]*atom.ry[i] + hinmat[1][2]*atom.rz[i];
    zz = hinmat[2][0]*atom.rx[i] + hinmat[2][1]*atom.ry[i] + hinmat[2][2]*atom.rz[i];
    fprintf(fp," %20.15e %20.15e %20.15e %s\n",xx,yy,zz,atom.asp[i]);
  }

  fclose(fp);
}

void writeconfig_abs(const char* fname)
{
  printf ("Write config file = %s\n", fname);
  /*
  std::ofstream foutconfig( "config.end", std::ios::out );
  foutconfig << "ATOM_TYPE" << std::endl;
  foutconfig << "POTENTIAL" << std::endl;
  foutconfig << atom.natom << std::endl;
  */
  FILE *fp;
  double xx, yy, zz;
  double hinmat[3][3];
  char species[3];
  //  fp = fopen("config.end", "w");
  fp = fopen(fname, "w");
  strcpy(species,atom.asp[1]);
  if (atom.natom > 1) {
    for (int i=2; i<=atom.natom; i++) {
      if(strcmp(atom.asp[i], species)==1) {strcpy(species,"M");} } }
  fprintf(fp," %s\n",species);
  fprintf(fp," %s %s\n",atom.potential_func,atom.potential_arg);
  fprintf(fp," %d\n",atom.natom);
  fprintf(fp," %e\n",cell.alat);
  xx=cell.hmat[0][0]/cell.alat/1.0e-10; yy=cell.hmat[1][0]/cell.alat/1.0e-10; zz=cell.hmat[2][0]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][1]/cell.alat/1.0e-10; yy=cell.hmat[1][1]/cell.alat/1.0e-10; zz=cell.hmat[2][1]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][2]/cell.alat/1.0e-10; yy=cell.hmat[1][2]/cell.alat/1.0e-10; zz=cell.hmat[2][2]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  fprintf(fp," %s\n","abs");
  inverse(cell.hmat, hinmat);
  for (int i=1; i<=atom.natom; i++) {
    fprintf(fp," %20.15e %20.15e %20.15e %s\n",atom.rx[i]/ang,atom.ry[i]/ang,atom.rz[i]/ang,atom.asp[i]);
  }

  fclose(fp);
}

void writeposcar(const char* fname)
{
  printf ("Write POSCAR file = %s\n", fname);
  FILE *fp;
  double xx, yy, zz;
  double hinmat[3][3];
  char species[3];
  int nesp[50];
  char **aesp;
  aesp = new char*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) { aesp[i] = new char[3]; }

  fp = fopen(fname, "w");
  int nsp = count_num_species(nesp,aesp);
  for (int i=1; i<=nsp; i++) {
    fprintf(fp, " %s", aesp[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "1.000\n");
  for (int j=0; j<3; j++) {
    double xx = cell.hmat[0][j]/ang;
    double yy = cell.hmat[1][j]/ang;
    double zz = cell.hmat[2][j]/ang;
    fprintf(fp, " %20.15e  %20.15e  %20.15e \n",xx,yy,zz);
  }
  for (int i=1; i<=nsp; i++) {
    fprintf(fp, " %d", nesp[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, " direct\n");
  inverse(cell.hmat, hinmat);
  for (int j=1; j<=nsp; j++) {
    for (int i=1; i<=atom.natom; i++) {
      if (strcmp(atom.asp[i],aesp[j])==0) {
	double xx = hinmat[0][0]*atom.rx[i]+hinmat[0][1]*atom.ry[i]+hinmat[0][2]*atom.rz[i];
	double yy = hinmat[1][0]*atom.rx[i]+hinmat[1][1]*atom.ry[i]+hinmat[1][2]*atom.rz[i];
	double zz = hinmat[2][0]*atom.rx[i]+hinmat[2][1]*atom.ry[i]+hinmat[2][2]*atom.rz[i];
	fprintf(fp, " %20.15e  %20.15e  %20.15e \n",xx,yy,zz);
      }
    }
  }
  fclose(fp);
  for (int i=0; i<=atom.natom; i++) { delete[] aesp[i]; }
  delete[] aesp;
}

int count_num_species(int nesp[50], char **aesp)
{
  int isnewspecies, atom_number;
  int nsp = 1;
  //  int nesp[50];

  strcpy(aesp[1],atom.asp[1]);
  for (int i=2; i<=atom.natom; i++) {
    // is this new species of atom?
    isnewspecies=1;
    for (int j=1; j<=i-1; j++) {
      if (strcmp(atom.asp[j],atom.asp[i])==0) {
	isnewspecies=0;
	break;
      }
    }
    if (isnewspecies==1) {
      nsp++;
      strcpy(aesp[nsp],atom.asp[i]);
    }
  }
        
  for (int i=1; i<=nsp; i++) {
    nesp[i]=0;
    for (int j=1; j<=atom.natom; j++) {
      if (strcmp(atom.asp[j],aesp[i])==0) { nesp[i]++; }
    }
  }
  return nsp;
}
