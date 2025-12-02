#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"
#include <GL/glui.h>
#include <complex>

using namespace std;

extern GLuint objects;
//extern GLfloat *color[NOBJECTS];
extern GLfloat **color, **color0;
extern int *iatom, *repidx;
extern GLfloat red[], yellow[];
extern float ex;
extern char config_atom[];

void configmake_fcc();
void configmake_bcc();
void configmake_dia();
void cellmultiply(int ncell);
void configmake_cnt();
int ncount_cnt();
int lcm(int n1, int n2);
void set_atom_color();
void set_atom_weight();
void potential();
void bookkeep();
void pot_initialize_all();
void cnt_wall_discard();
int atom_number(char* at);
void deallocate_arrays();
void allocate_arrays();

void createconfig()
{
  deallocate_arrays();

  /*
  // Deallocate arrays
  if (color) { for (int i=1; i<=atom.natom*3; i++) { if (color[i]) free(color[i]); } free(color); }
  if (color0) { for (int i=1; i<=atom.natom*3; i++) { if (color0[i]) free(color0[i]); } free(color0); }
  if (iatom) { free(iatom); }
  if (repidx) { free(repidx); }
  if (atom.asp) { for (int i=1; i<=atom.natom; i++) { if (atom.asp[i]) free(atom.asp[i]); } free(atom.asp);  }
  if (atom.anum) {delete[] atom.anum; atom.anum = NULL; }
  if (atom.wm) { delete[] atom.wm; atom.wm = NULL;}
  if (atom.rx) {free(atom.rx);}; if (atom.ry) {free(atom.ry);}; if (atom.rz) {free(atom.rz);}
  if (atom.rx_float) {free(atom.rx_float);}; if (atom.ry_float) {free(atom.ry_float);}; if (atom.rz_float) {free(atom.rz_float);}
  if (atom.rx_org) {free(atom.rx_org);}; if (atom.ry_org) {free(atom.ry_org);}; if (atom.rz_org) {free(atom.rz_org);}
  if (atom.rx_p) {free(atom.rx_p);}; if (atom.ry_p) {free(atom.ry_p);}; if (atom.rz_p) {free(atom.rz_p);}
  if (atom.vx) {free(atom.vx);}; if (atom.vy) {free(atom.vy);}; if (atom.vz) {free(atom.vz);}
  if (atom.fx) {free(atom.fx);}; if (atom.fy) {free(atom.fy);}; if (atom.fz) {free(atom.fz);}
  if (atom.ax) { delete[] atom.ax; atom.ax = NULL; }
  if (atom.ay) { delete[] atom.ay; atom.ay = NULL; }
  if (atom.az) { delete[] atom.az; atom.az = NULL; }
  if (atom.bx) { delete[] atom.bx; atom.bx = NULL; }
  if (atom.by) { delete[] atom.by; atom.by = NULL; }
  if (atom.bz) { delete[] atom.bz; atom.bz = NULL; }
  if (atom.cx) { delete[] atom.cx; atom.cx = NULL; }
  if (atom.cy) { delete[] atom.cy; atom.cy = NULL; }
  if (atom.cz) { delete[] atom.cz; atom.cz = NULL; }
  if (atom.fx_l) { delete[] atom.fx_l; atom.fx_l = NULL; }
  if (atom.fy_l) { delete[] atom.fy_l; atom.fy_l = NULL; }
  if (atom.fz_l) { delete[] atom.fz_l; atom.fz_l = NULL; }
  if (atom.fx_float) {free(atom.fx_float);}; if (atom.fy_float) {free(atom.fy_float);}; if (atom.fz_float) {free(atom.fz_float);}
  if (atom.mfx) {free(atom.mfx);}; if (atom.mfy) {free(atom.mfy);}; if (atom.mfz) {free(atom.mfz);}
  if (atom.epot) {free(atom.epot);}
  if (atom.epot_float) {free(atom.epot_float);}
  if (atom.repatom) {free(atom.repatom);}
  if (atom.elem_id) {free(atom.elem_id);}
  if (atom.evecx) { for (int i=1; i<=atom.natom; i++) { free(atom.evecx[i]); } free(atom.evecx);}
  if (atom.evecy) { for (int i=1; i<=atom.natom; i++) { free(atom.evecy[i]); } free(atom.evecy);}
  if (atom.evecz) { for (int i=1; i<=atom.natom; i++) { free(atom.evecz[i]); } free(atom.evecz);}
  if (atom.satom) {
    for (int i=1; i<=atom.natom; i++) {
      for (int j=0; j<3; j++) {	free(atom.satom[i][j]); }  free(atom.satom[i]); } 
    free(atom.satom); }
  if (atom.nneighbor) { delete[] atom.nneighbor; atom.nneighbor = NULL;
    for (int i=0; i<=atom.natom; i++) { delete[] atom.neighbor[i]; }
    delete[] atom.neighbor; atom.neighbor = NULL; }
  delete[] book.alistnum;
 // Deallocate arrays end
 */

  //if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); printf("glDeleteLists\n"); }
  if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); }

  if (config_type==0) {
    atom.natom = 4 * irepx * irepy * irepz;
  } else if (config_type==1) {
    atom.natom = 2 * irepx * irepy * irepz;
  } else if (config_type==2) {
    atom.natom = 8 * irepx * irepy * irepz;
  } else if (config_type==3) {
    atom.natom = ncount_cnt();
  }

  allocate_arrays();

  /*
  // Allocate arrays
  //  color = (GLfloat **)malloc(sizeof(GLfloat*)*(atom.natom*3+1));
  //  for (int i=1; i<=atom.natom*3; i++) { color[i] = (GLfloat *)malloc(sizeof(GLfloat)*4); }
  color = (GLfloat **)malloc(sizeof(GLfloat*)*(atom.natom*3+1));
  color0 = (GLfloat **)malloc(sizeof(GLfloat*)*(atom.natom*3+1));
  for (int i=1; i<=atom.natom*3; i++) { color[i] = (GLfloat *)malloc(sizeof(GLfloat)*4); }
  for (int i=1; i<=atom.natom*3; i++) { color0[i] = (GLfloat *)malloc(sizeof(GLfloat)*4); }
  iatom  = (int *)calloc(atom.natom*2+1, sizeof(int));
  repidx = (int *)calloc(atom.natom*2+1, sizeof(int));
  atom.asp = (char **)malloc(sizeof(char*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.asp[i] = (char *)malloc(sizeof(char)*3); }
  atom.anum = new int[atom.natom+1];
  atom.wm = new double[atom.natom+1];
  atom.rx = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ry = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rz = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ax = new double[atom.natom+1];
  atom.ay = new double[atom.natom+1];
  atom.az = new double[atom.natom+1];
  atom.bx = new double[atom.natom+1];
  atom.by = new double[atom.natom+1];
  atom.bz = new double[atom.natom+1];
  atom.cx = new double[atom.natom+1];
  atom.cy = new double[atom.natom+1];
  atom.cz = new double[atom.natom+1];
  atom.rx_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.ry_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.rz_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.rx_org = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ry_org = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rz_org = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rx_p = (double *)calloc(atom.natom+1, sizeof(double));
  atom.ry_p = (double *)calloc(atom.natom+1, sizeof(double));
  atom.rz_p = (double *)calloc(atom.natom+1, sizeof(double));
  atom.vx = (double *)calloc(atom.natom+1, sizeof(double));
  atom.vy = (double *)calloc(atom.natom+1, sizeof(double));
  atom.vz = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fx = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fy = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fz = (double *)calloc(atom.natom+1, sizeof(double));
  atom.fx_l = new double[atom.natom+1];
  atom.fy_l = new double[atom.natom+1];
  atom.fz_l = new double[atom.natom+1];
  atom.fx_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.fy_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.fz_float = (FLOAT *)calloc(atom.natom+1, sizeof(FLOAT));
  atom.mfx = (bool *)calloc(atom.natom+1, sizeof(bool));
  atom.mfy = (bool *)calloc(atom.natom+1, sizeof(bool));
  atom.mfz = (bool *)calloc(atom.natom+1, sizeof(bool));
  atom.epot = (double *)calloc(atom.natom+1, sizeof(double));
  atom.epot_float = (FLOAT2 *)calloc(atom.natom+1, sizeof(FLOAT2));
  atom.repatom = (int *)calloc(atom.natom+1, sizeof(int));
  atom.elem_id = (int *)calloc(atom.natom+1, sizeof(int));
  atom.evecx = (double **)malloc(sizeof(double*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.evecx[i] = (double *)malloc(sizeof(double)*20); }
  atom.evecy = (double **)malloc(sizeof(double*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.evecy[i] = (double *)malloc(sizeof(double)*20); }
  atom.evecz = (double **)malloc(sizeof(double*)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) { atom.evecz[i] = (double *)malloc(sizeof(double)*20); }
  atom.satom = (double ***)malloc(sizeof(double**)*(atom.natom+1));
  for (int i=1; i<=atom.natom; i++) {
    atom.satom[i] = (double **)malloc(sizeof(double*)*3);
    for (int j=0; j<3; j++) { atom.satom[i][j] = (double *)malloc(sizeof(double)*3); } }
  atom.nneighbor = new int[atom.natom+1]; atom.neighbor = new int*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) { atom.neighbor[i] = new int[MAXNEIGHBOR]; }
  book.alistnum = new int[atom.natom+1]; book.alloc = true;
  
  // Allocate arrays end
  */

  objects = glGenLists(atom.natom*3); //printf("objects = %d\n",objects);
  for (int i=1; i<=atom.natom*3; i++) { memcpy(color[i],yellow,sizeof(GLfloat)*4); }

  printf("======= CREATION OF NEW CONFIG ======\n");
  if (config_type==0) {
    printf("FCC %d x %d x %d x, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_fcc();
  } else if (config_type==1) {
    printf("BCC %d x %d x %d x, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_bcc();
  } else if (config_type==2) {
    printf("Diamond %d x %d x %d x, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_dia();
  } else if (config_type==3) {
    printf("Nanotube (%d, %d) x %d, N = %d\n",icntm,icntn,irepz,atom.natom);
    configmake_cnt(); cnt_wall_discard();
  }

  /*
  if (config_type ==3) {
    for (int i=1; i<=atom.natom; i++) { strcpy (atom.asp[i], "C"); }
  } else {
    for (int i=1; i<=atom.natom; i++) { strcpy (atom.asp[i], "Cu"); }
  }
  */

  if (strcmp(config_atom,"M")==0) {
    printf("Multi atom system\n");
  } else{
    printf("Atom type: %s\n",config_atom);
    for (int i=1; i<=atom.natom; i++) {
      strcpy(atom.asp[i], config_atom);
    }
  }
  printf("Cell size (ang): %f x %f x %f\n",
	 cell.hmat[0][0]/ang,cell.hmat[1][1]/ang,cell.hmat[2][2]/ang);
  printf("=====================================\n");
  cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;

  // Atom color
  set_atom_color();
  // Atom weight
  set_atom_weight();
  // Atom number
  for (int i=1; i<=atom.natom; i++) {
    atom.anum[i] = atom_number(atom.asp[i]);
  }
  for (int i=1; i<=atom.natom; i++) {
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
    atom.vx[i] = 0.0; atom.vy[i] = 0.0; atom.vz[i] = 0.0;
    atom.ax[i] = 0.0; atom.ay[i] = 0.0; atom.az[i] = 0.0;
    atom.bx[i] = 0.0; atom.by[i] = 0.0; atom.bz[i] = 0.0;
    atom.cx[i] = 0.0; atom.cy[i] = 0.0; atom.cz[i] = 0.0;
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.hmat_org[i][j]=cell.hmat[i][j];
      cell.hvmat[i][j] = 0; cell.hamat[i][j] = 0; cell.hbmat[i][j] = 0;
      cell.hcmat[i][j] = 0; cell.sgmmat_set[i][j] = 0;
    }
  }
  istep = 0;
  ex = 0.0;
  pot_initialize_all();
  bookkeep();
  potential();
}

void configmake_fcc()
{
  double ao=alat*1.0e-10;
  atom.rx[1]=   0.0;atom.ry[1]=   0.0;atom.rz[1]=   0.0;
  atom.rx[2]=   0.0;atom.ry[2]=ao/2.0;atom.rz[2]=ao/2.0;
  atom.rx[3]=ao/2.0;atom.ry[3]=   0.0;atom.rz[3]=ao/2.0;
  atom.rx[4]=ao/2.0;atom.ry[4]=ao/2.0;atom.rz[4]=   0.0;
  int ncell=4;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) { cell.hmat[i][j]=0.0; }
    cell.hmat[i][i]=ao;
  }
  cellmultiply(ncell);
}
void configmake_bcc()
{
  double ao=alat*1.0e-10;
  atom.rx[1]=   0.0;atom.ry[1]=   0.0;atom.rz[1]=   0.0;
  atom.rx[2]=ao/2.0;atom.ry[2]=ao/2.0;atom.rz[2]=ao/2.0;
  int ncell=2;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) { cell.hmat[i][j]=0.0; }
    cell.hmat[i][i]=ao;
  }
  cellmultiply(ncell);
}
void configmake_dia()
{
  double ao=alat*1.0e-10;
  atom.rx[1]=   0.0;atom.ry[1]=   0.0;atom.rz[1]=   0.0;
  atom.rx[2]=   0.0;atom.ry[2]=ao/2.0;atom.rz[2]=ao/2.0;
  atom.rx[3]=ao/2.0;atom.ry[3]=   0.0;atom.rz[3]=ao/2.0;
  atom.rx[4]=ao/2.0;atom.ry[4]=ao/2.0;atom.rz[4]=   0.0;
  for (int i=5; i<=8; i++) {
    atom.rx[i]=atom.rx[i-4]+ao/4.0;
    atom.ry[i]=atom.ry[i-4]+ao/4.0;
    atom.rz[i]=atom.rz[i-4]+ao/4.0;
  }
  int ncell=8;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) { cell.hmat[i][j]=0.0; }
    cell.hmat[i][i]=ao;
  }
  cellmultiply(ncell);
}

void cellmultiply(int ncell)
{
  int iii=ncell;
  int jjj=iii;
  if (irepx > 1) {
    for (int k=2; k<=irepx; k++) {
      for (int j=1; j<=jjj; j++) {
	iii++;
	atom.rx[iii]=atom.rx[iii-jjj]+cell.hmat[0][0];
	atom.ry[iii]=atom.ry[iii-jjj]+cell.hmat[1][0];
	atom.rz[iii]=atom.rz[iii-jjj]+cell.hmat[2][0];
      }
    }
    jjj=iii;
  }
  if (irepy > 1) {
    for (int k=2; k<=irepy; k++) {
      for (int j=1; j<=jjj; j++) {
	iii++;
	atom.rx[iii]=atom.rx[iii-jjj]+cell.hmat[0][1];
	atom.ry[iii]=atom.ry[iii-jjj]+cell.hmat[1][1];
	atom.rz[iii]=atom.rz[iii-jjj]+cell.hmat[2][1];
      }
    }
    jjj=iii;
  }
  if (irepz > 1) {
    for (int k=2; k<=irepz; k++) {
      for (int j=1; j<=jjj; j++) {
	iii++;
	atom.rx[iii]=atom.rx[iii-jjj]+cell.hmat[0][2];
	atom.ry[iii]=atom.ry[iii-jjj]+cell.hmat[1][2];
	atom.rz[iii]=atom.rz[iii-jjj]+cell.hmat[2][2];
      }
    }
  }
  for (int i=0; i<3; i++) {
    cell.hmat[i][0]=cell.hmat[i][0]*(double)irepx;
    cell.hmat[i][1]=cell.hmat[i][1]*(double)irepy;
    cell.hmat[i][2]=cell.hmat[i][2]*(double)irepz;
  }
}

void configmake_cnt()
{
  int na=4;
  double x[5],y[5],z[5];
  double a1[4],a2[4],a3[4];
  complex<double> cr, crcnj, ctube, ctubecnj, ctrans;

  int natom = 1;
  int iatom = 0;
  double ao = 1.42e-10;
  int nwall = 0;

  int iwall = 1;
  int n1 = icntm;
  int n2 = icntn;
  int nlayer = irepz;

  double pi = 4.0*atan(1.0);

  if (n1+n2 < 2) { return; }

  a1[1]=ao*sqrt(3.0);
  a1[2]=0.0;
  a1[3]=0.0;
  a2[1]=0.0;
  a2[2]=ao*3.0;
  a2[3]=0.0;
  a3[1]=0.0;
  a3[2]=0.0;
  a3[3]=6.696/1.42*ao;
  
  x[1]=0.0;
  y[1]=0.0;
  z[1]=0.0;
  x[2]=0.0;
  y[2]=ao;
  z[2]=0.0;
  x[3]=ao*sqrt(3.0)/2.0;
  y[3]=-ao/2.0;
  z[3]=0.0;
  x[4]=ao*sqrt(3.0)/2.0;
  y[4]=ao*3.0/2.0;
  z[4]=0.0;

  double at11=a1[1];
  double at12=0.0;
  double at21=a1[1]/2.0;
  double at22=-ao*3.0/2.0;

  double atube1=at11*(double)n1+at21*(double)n2;
  double atube2=at12*(double)n1+at22*(double)n2;
  double rtubex=sqrt(atube1*atube1+atube2*atube2);
  
  double tubecos=atube1/rtubex;
  double tubesin=atube2/rtubex;
  ctube=complex<double>(tubecos,tubesin);
  ctubecnj=complex<double>(tubecos,-tubesin);
  for (int i=1; i<=na; i++) {
    cr=complex<double>(x[i],y[i])*ctubecnj;
    x[i]=cr.real();
    y[i]=cr.imag();
  }

  int nn1=2*n1+n2;
  int nn2=2*n2+n1;
  //  call lcm(nn1,nn2,nc)
  int nc = lcm(nn1,nn2);
  
  double as1=at11*(double)(nc/nn1)-at21*(double)(nc/nn2);
  double as2=at12*(double)(nc/nn1)-at22*(double)(nc/nn2);
  double rtubey=sqrt(as1*as1+as2*as2);
  //  printf("%e %e\n",as1,as2);

  ctrans=complex<double>(a1[1],a1[2]);
  ctrans=ctrans*ctubecnj;
  a1[1]=ctrans.real();
  a1[2]=ctrans.imag();
  ctrans=complex<double>(a2[1],a2[2]);
  ctrans=ctrans/ctube;
  a2[1]=ctrans.real();
  a2[2]=ctrans.imag();

  for (int i=1; i<=na; i++) {
    for (int ix=-100; ix<=100; ix++) {
      for (int iy=-100; iy<=100; iy++) {
	double xx=x[i]+(double)ix*a1[1]+(double)iy*a2[1];
	double yy=y[i]+(double)ix*a1[2]+(double)iy*a2[2];
	if ((xx >= -1.0e-12)&&(xx < rtubex-1.0e-12)) {
	  if ((yy >= -1.0e-12)&&(yy < rtubey-1.0e-12)) {
	    iatom++;
	    if (iatom > atom.natom) { printf("Warning in configmake_cnt\n"); }
	    atom.rx[iatom]=xx;
	    atom.ry[iatom]=yy;
	    atom.rz[iatom]=0.0;
	  }
	}
      }
    }
  }

  for (int i=natom; i<=iatom; i++) {
    atom.rz[i]=atom.ry[i];
    double xx=atom.rx[i];
    atom.rx[i]=rtubex/(2.0*pi)*cos(2.0*pi*xx/rtubex+2.0*pi/96.0);
    atom.ry[i]=rtubex/(2.0*pi)*sin(2.0*pi*xx/rtubex+2.0*pi/96.0);
  }
  natom=iatom+1;

  int katom=iatom;
  for (int i=1; i<=nlayer-1; i++) {
    for (int j=1; j<=katom; j++) {
      iatom++;
      if (iatom > atom.natom) { printf("Warning in configmake_cnt (2)\n"); }
      atom.rx[iatom]=atom.rx[j];
      atom.ry[iatom]=atom.ry[j];
      atom.rz[iatom]=atom.rz[j]+rtubey*(double)i;
    }
  }

  // Rotation
  double rad = sqrt(atom.rx[1]*atom.rx[1]+atom.ry[1]*atom.ry[1]);
  for (int i=1; i<=iatom; i++) {
    double theta = 180.0/pi*atan(atom.ry[i]/atom.rx[i]);
    if ((atom.ry[i] > 0.0)&&(theta < 0.0)) { theta = theta + 180.0; }
    if ((atom.ry[i] < 0.0)&&(theta > 0.0)) { theta = theta + 180.0; }
    if ((atom.ry[i] < 0.0)&&(theta < 0.0)) { theta = theta + 360.0; }
    if ((atom.ry[i] == 0.0)&&(atom.rx[i]>0.0)) { theta = 0; }
    if ((atom.ry[i] == 0.0)&&(atom.rx[i]<0.0)) { theta = 180; }
    if ((atom.rx[i] == 0.0)&&(atom.ry[i]>0.0)) { theta = 90; }
    if ((atom.rx[i] == 0.0)&&(atom.ry[i]<0.0)) { theta = 270; }
    atom.rx[i]=rad*cos((theta+rotz)/180.0*pi);
    atom.ry[i]=rad*sin((theta+rotz)/180.0*pi);
  }
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { cell.hmat[i][j]=0.0; } }
  cell.hmat[0][0]=cscnt*1e-10;
  cell.hmat[1][1]=cscnt*1e-10;
  cell.hmat[2][2]=rtubey*(double)nlayer;
  // Shift
  for (int i=1; i<=iatom; i++) {
    atom.rx[i]=atom.rx[i]+cell.hmat[0][0]/2.0;
    atom.ry[i]=atom.ry[i]+cell.hmat[1][1]/2.0;
    atom.rz[i]=atom.rz[i]+cell.hmat[2][2]*shiftz/100.0;
    if (atom.rz[i] > cell.hmat[2][2]) { atom.rz[i] = atom.rz[i] - cell.hmat[2][2]; }
  }
  printf("CNT diameter (ang) = %f\n",rad*2/ang);
}

int lcm(int n1, int n2)
{
  int m1, m2, mm1, nc;
  int dif = n1-n2;
  if (dif < 0) {
    m1=n1;
    m2=n2;
  } else {
    m1=n2;
    m2=n1;
  }
  mm1=m1;
  nc=1;

  int ifg=1;
  while(ifg == 1)
    {
      for (int i=2; i<=mm1; i++) {
	if ((m1%i == 0)&&(m2%i == 0)) {
	  m1=m1/i;
	  m2=m2/i;
	  nc=nc*i;
	  ifg = 1;
	  break;
	} else {
	  ifg = 0;
	}
      }
    }
  nc=nc*m1*m2;
  return nc;
}

int ncount_cnt()
{
  int na=4;
  double x[5],y[5],z[5];
  double a1[4],a2[4],a3[4];
  complex<double> cr, crcnj, ctube, ctubecnj, ctrans;

  int natom = 1;
  int iatom = 0;
  double ao = 1.42e-10;
  int nwall = 0;

  int iwall = 1;
  int n1 = icntm;
  int n2 = icntn;
  int nlayer = irepz;

  double pi = 4.0*atan(1.0);

  if (n1+n2 < 2) { return 1; }

  a1[1]=ao*sqrt(3.0);
  a1[2]=0.0;
  a1[3]=0.0;
  a2[1]=0.0;
  a2[2]=ao*3.0;
  a2[3]=0.0;
  a3[1]=0.0;
  a3[2]=0.0;
  a3[3]=6.696/1.42*ao;
  
  x[1]=0.0;
  y[1]=0.0;
  z[1]=0.0;
  x[2]=0.0;
  y[2]=ao;
  z[2]=0.0;
  x[3]=ao*sqrt(3.0)/2.0;
  y[3]=-ao/2.0;
  z[3]=0.0;
  x[4]=ao*sqrt(3.0)/2.0;
  y[4]=ao*3.0/2.0;
  z[4]=0.0;

  double at11=a1[1];
  double at12=0.0;
  double at21=a1[1]/2.0;
  double at22=-ao*3.0/2.0;

  double atube1=at11*(double)n1+at21*(double)n2;
  double atube2=at12*(double)n1+at22*(double)n2;
  double rtubex=sqrt(atube1*atube1+atube2*atube2);
  
  double tubecos=atube1/rtubex;
  double tubesin=atube2/rtubex;
  ctube=complex<double>(tubecos,tubesin);
  ctubecnj=complex<double>(tubecos,-tubesin);
  for (int i=1; i<=na; i++) {
    cr=complex<double>(x[i],y[i])*ctubecnj;
    x[i]=cr.real();
    y[i]=cr.imag();
  }

  int nn1=2*n1+n2;
  int nn2=2*n2+n1;
  int nc = lcm(nn1,nn2);
  
  double as1=at11*(double)(nc/nn1)-at21*(double)(nc/nn2);
  double as2=at12*(double)(nc/nn1)-at22*(double)(nc/nn2);
  double rtubey=sqrt(as1*as1+as2*as2);

  ctrans=complex<double>(a1[1],a1[2]);
  ctrans=ctrans*ctubecnj;
  a1[1]=ctrans.real();
  a1[2]=ctrans.imag();
  ctrans=complex<double>(a2[1],a2[2]);
  ctrans=ctrans/ctube;
  a2[1]=ctrans.real();
  a2[2]=ctrans.imag();

  for (int i=1; i<=na; i++) {
    for (int ix=-100; ix<=100; ix++) {
      for (int iy=-100; iy<=100; iy++) {
	double xx=x[i]+(double)ix*a1[1]+(double)iy*a2[1];
	double yy=y[i]+(double)ix*a1[2]+(double)iy*a2[2];
	if ((xx >= -1.0e-12)&&(xx < rtubex-1.0e-12)) {
	  if ((yy >= -1.0e-12)&&(yy < rtubey-1.0e-12)) {
	    iatom++;
	  }
	}
      }
    }
  }
  natom=iatom+1;

  int katom=iatom;
  for (int i=1; i<=nlayer-1; i++) {
    for (int j=1; j<=katom; j++) {
      iatom++;
    }
  }
  return iatom;
}

