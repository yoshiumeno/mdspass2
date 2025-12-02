#include <string.h>
#include <iostream>
#include <GL/glui.h>
#include<fstream>
#include<math.h>
#include "myheader.h"

extern float obj_pos[3];
extern float rotate[16];
extern int wireframe, mdspeed, radius, show_axis, show_only_elem, show_cell;
extern int segments, ibase, icnt, draw_replica, draw_bond, draw_bond_pbc, draw_force;
extern int draw_load, draw_aux;
extern int *iatom, *repidx, ievec, ievec_num, show_cnt_wall, show_cnt_wall_num, show_cnt_ring;
extern int show_cnt_wallv;
extern float scl, evec_len, vscl, vscl_force, cnt_ring_radius;
extern double mat[3][3];
extern GLuint objects;
extern int edit_elem_mode, select_atom[10], select_atom_repidx[10];
extern GLfloat white[], red[], yellow[], gray[], blue[], blue2[], green[], green2[];
extern GLfloat **color, **color0;
extern GLUI *glui;
extern int createconfig_mode, iatom_pick;

void md();
void bookkeep();
void stretch(double x, double y, double z);
void writedata();
void drawCell(); void drawYZPlane(float x);
void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3]);
void putSphere( float x, float y, float z, float rad, int seg, int solid);
//void putArrow(float xx0, float yy0, float zz0, float vx, float vy, float vz, float vscale);
void glDrawArrowd(float x0, float y0, float z0, float x1, float y1, float z1);
void glDrawArrowd2(float x, float y, float z, float vx, float vy, float vz, float u);
void glDrawPiped(float x0, float y0, float z0, float x1, float y1, float z1, float r);
void glDrawAxisd(float length);
void out_of_cell();
void change_atom_color();

void putSphere( float x, float y, float z, float rad, int seg, int solid)
{
  glTranslated(x*scl, y*scl, z*scl);
  if (solid==0) {
    glutWireSphere( rad*scl, seg, seg );
  } else {
    glutSolidSphere( rad*scl, seg, seg );
  }
  glTranslated(-x*scl, -y*scl, -z*scl);
}

void drawCell()
{
  float xx1, yy1, zz1, xx2, yy2, zz2;
  glDisable( GL_LIGHTING );  
  glColor3d(0.0, 0.0, 0.0);
  glBegin(GL_LINES);
  xx1=0.0; yy1=0.0; zz1=0.0;
  xx2=cell.hmat[0][0]*1e9; yy2=cell.hmat[1][0]*1e9; zz2=cell.hmat[2][0]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx1=0.0; yy1=0.0; zz1=0.0;
  xx2=cell.hmat[0][1]*1e9; yy2=cell.hmat[1][1]*1e9; zz2=cell.hmat[2][1]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx1=0.0; yy1=0.0; zz1=0.0;
  xx2=cell.hmat[0][2]*1e9; yy2=cell.hmat[1][2]*1e9; zz2=cell.hmat[2][2]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);

  xx1=cell.hmat[0][0]*1e9; yy1=cell.hmat[1][0]*1e9; zz1=cell.hmat[2][0]*1e9; 
  xx2=xx1+cell.hmat[0][1]*1e9; yy2=yy1+cell.hmat[1][1]*1e9; zz2=zz1+cell.hmat[2][1]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx2=xx1+cell.hmat[0][2]*1e9; yy2=yy1+cell.hmat[1][2]*1e9; zz2=zz1+cell.hmat[2][2]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx1=xx2; yy1=yy2; zz1=zz2;
  xx2=xx1+cell.hmat[0][1]*1e9; yy2=yy1+cell.hmat[1][1]*1e9; zz2=zz1+cell.hmat[2][1]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx1=xx2; yy1=yy2; zz1=zz2;
  xx2=xx1-cell.hmat[0][2]*1e9; yy2=yy1-cell.hmat[1][2]*1e9; zz2=zz1-cell.hmat[2][2]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);

  xx1=cell.hmat[0][1]*1e9; yy1=cell.hmat[1][1]*1e9; zz1=cell.hmat[2][1]*1e9; 
  xx2=xx1+cell.hmat[0][0]*1e9; yy2=yy1+cell.hmat[1][0]*1e9; zz2=zz1+cell.hmat[2][0]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx2=xx1+cell.hmat[0][2]*1e9; yy2=yy1+cell.hmat[1][2]*1e9; zz2=zz1+cell.hmat[2][2]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx1=xx2; yy1=yy2; zz1=zz2;
  xx2=xx1+cell.hmat[0][0]*1e9; yy2=yy1+cell.hmat[1][0]*1e9; zz2=zz1+cell.hmat[2][0]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);

  xx1=cell.hmat[0][2]*1e9; yy1=cell.hmat[1][2]*1e9; zz1=cell.hmat[2][2]*1e9; 
  xx2=xx1+cell.hmat[0][0]*1e9; yy2=yy1+cell.hmat[1][0]*1e9; zz2=zz1+cell.hmat[2][0]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  xx2=xx1+cell.hmat[0][1]*1e9; yy2=yy1+cell.hmat[1][1]*1e9; zz2=zz1+cell.hmat[2][1]*1e9; 
  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  glEnd();
  glEnable( GL_LIGHTING );  
}

void drawYZPlane(float x)
{
  float xx1, yy1, zz1, xx2, yy2, zz2;
  glDisable( GL_LIGHTING );  
  glColor3d(1.0, 0.0, 1.0);
  glBegin(GL_LINES);
  for (int i=0; i<=10; i++) {
    xx1=cell.hmat[0][0]/2 + x*ang; yy1=0.0; zz1=cell.hmat[2][2]/10*(float)i;
    xx2=cell.hmat[0][0]/2 + x*ang; yy2=cell.hmat[1][1]; zz2=cell.hmat[2][2]/10*(float)i;
    xx1 *= 1e9; yy1 *= 1e9; zz1 *= 1e9; xx2 *= 1e9; yy2 *= 1e9; zz2 *= 1e9; 
    glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  }
  for (int i=0; i<=10; i++) {
    xx1=cell.hmat[0][0]/2 + x*ang; yy1=cell.hmat[1][1]/10*(float)i; zz1=0.0;
    xx2=cell.hmat[0][0]/2 + x*ang; yy2=cell.hmat[1][1]/10*(float)i; zz2=cell.hmat[2][2];
    xx1 *= 1e9; yy1 *= 1e9; zz1 *= 1e9; xx2 *= 1e9; yy2 *= 1e9; zz2 *= 1e9; 
    glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
  }

  glEnd();
  glEnable( GL_LIGHTING );  
}
  
void putLine(float x1, float y1, float z1, float x2, float y2, float z2)
{
  double vec1[3], vec2[3];
  glBegin(GL_LINES);
  vec1[0]=x1*scl; vec1[1]=y1*scl; vec1[2]=z1*scl;
  vec2[0]=x2*scl; vec2[1]=y2*scl; vec2[2]=z2*scl;
  glVertex3dv(vec1);  glVertex3dv(vec2);
  glEnd();
}

void myGlutDisplay( void )
{
  float xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, xx4, yy4, zz4;
  glClearColor( .9f, .9f, .9f, 1.0f );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  //  glMatrixMode( GL_PROJECTION );
  //  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] ); 
  glMultMatrixf( rotate );
  int solid = 1;
  if (wireframe == 1) { solid = 0; } else { solid = 1; }
  /*
  if (mdmotion == 1) {
    for (int i=1; i<=mdspeed; i++) { // MD: numerical integral of EoM
      if (istep % book.nbk ==0) {
	if (incell) { out_of_cell(); }
	bookkeep();
      }
      //if (istep % book.nbk ==0) { bookkeep(); }
      md(); if (mdmotion == 0) { break; }
      writedata();
      //printf("Step= %d, Epot= %f, Ekin= %e, Temp= %f Lx= %f Fmax= %e\n",
      //istep,atom.epotsum/ev/atom.natom,atom.ekinsum/ev/atom.natom,
      //tempc,cell.hmat[1][1]/ang,atom.Fmax()/ev*ang);
      istep++;
    }
    //    writedata();
  }
  */
  
  float rad = (float)radius/100.0;
    
  // draw cell
  if (show_cell) {
    drawCell();
  }
  // draw axis
  if (show_axis) {
    glDisable( GL_LIGHTING );  
    glDrawAxisd(0.5);
  }
  if (ievec) {
    if (ievec_num<1) { ievec_num=1; }
    if (ievec_num>MAXMODE) { ievec_num=MAXMODE; }
    //glDisable( GL_LIGHTING );  glColor3d(1.0, 0.0, 1.0);
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, green); 
    for (int i=1; i<=atom.natom; i++) {
      //glBegin(GL_LINES);
      //      float vlen = 10.0;
      float xx1=atom.rx[i]*1e9; float yy1=atom.ry[i]*1e9; float zz1=atom.rz[i]*1e9;
      float xx2=atom.evecx[i][ievec_num-1]*evec_len;
      float yy2=atom.evecy[i][ievec_num-1]*evec_len;
      float zz2=atom.evecz[i][ievec_num-1]*evec_len;
      //glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
      //glDrawArrowd(xx1*scl,yy1*scl,zz1*scl,xx2*scl,yy2*scl,zz2*scl);
      glDrawArrowd2(xx1,yy1,zz1,xx2,yy2,zz2,scl);
      //glEnd();
    }
    //glEnable( GL_LIGHTING );
 }

  // draw elements
  if (atom.QC) {
    glDisable( GL_LIGHTING );  
    glColor3d(0.0, 0.5, 0.5);
    glBegin(GL_LINES);
    for (int i=1; i<=atom.nelem; i++) {
      mk_helmat(i, atom.rx, atom.ry, atom.rz, mat);
      int ix0=0; int iy0=0; int iz0=0;
      if (atom.elem_v_rep[i][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
      if ((atom.elem_v_rep[i][1] % 4 == 2)||(atom.elem_v_rep[i][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
      if (atom.elem_v_rep[i][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
      xx1=atom.rx[atom.elem_v[i][1]]*1e9; yy1=atom.ry[atom.elem_v[i][1]]*1e9; zz1=atom.rz[atom.elem_v[i][1]]*1e9; 
      xx1=xx1+(cell.hmat[0][0]*ix0+cell.hmat[0][1]*iy0+cell.hmat[0][2]*iz0)*1e9;
      yy1=yy1+(cell.hmat[1][0]*ix0+cell.hmat[1][1]*iy0+cell.hmat[1][2]*iz0)*1e9;
      zz1=zz1+(cell.hmat[2][0]*ix0+cell.hmat[2][1]*iy0+cell.hmat[2][2]*iz0)*1e9;
      xx2=xx1+mat[0][0]*1e9; yy2=yy1+mat[1][0]*1e9; zz2=zz1+mat[2][0]*1e9;
      xx3=xx1+mat[0][1]*1e9; yy3=yy1+mat[1][1]*1e9; zz3=zz1+mat[2][1]*1e9;
      xx4=xx1+mat[0][2]*1e9; yy4=yy1+mat[1][2]*1e9; zz4=zz1+mat[2][2]*1e9;

      glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
      glVertex3d(xx2*scl,yy2*scl,zz2*scl); glVertex3d(xx3*scl,yy3*scl,zz3*scl);
      glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx3*scl,yy3*scl,zz3*scl);
      
      glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx4*scl,yy4*scl,zz4*scl);
      glVertex3d(xx2*scl,yy2*scl,zz2*scl); glVertex3d(xx4*scl,yy4*scl,zz4*scl);
      glVertex3d(xx3*scl,yy3*scl,zz3*scl); glVertex3d(xx4*scl,yy4*scl,zz4*scl);

    }
    glEnd();
    glEnable( GL_LIGHTING );  
  }

  // draw walls (for CNT)
  if (show_cnt_wall) {
    glEnable( GL_LIGHTING );  
    glColor3d(0.0, 0.5, 0.5);
    float cube3 = (float)cell.hmat[2][2]*1e9;
    for (int i=1; i<=atom.nwall; i++) {
      int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
      if ((i0!=0)&&(i1!=0)&&(i2!=0)) {
	glBegin(GL_TRIANGLES);
	xx1=atom.rx[abs(i0)]*1e9; yy1=atom.ry[abs(i0)]*1e9; zz1=atom.rz[abs(i0)]*1e9;
	xx2=atom.rx[abs(i1)]*1e9; yy2=atom.ry[abs(i1)]*1e9; zz2=atom.rz[abs(i1)]*1e9;
	xx3=atom.rx[abs(i2)]*1e9; yy3=atom.ry[abs(i2)]*1e9; zz3=atom.rz[abs(i2)]*1e9;
	if (zz2>zz1+cube3/2) {zz2=zz2-cube3;} if (zz2<zz1-cube3/2) {zz2=zz2+cube3;}
	if (zz3>zz1+cube3/2) {zz3=zz3-cube3;} if (zz3<zz1-cube3/2) {zz3=zz3+cube3;}
	if (cube3 < 3.2) {
	  if (i1<0) { if (zz2<cube3/2) { zz2 += cube3; } else { zz2 -= cube3; } }
	  if (i2<0) { if (zz3<cube3/2) { zz3 += cube3; } else { zz3 -= cube3; } }
	}
	glVertex3d(xx1*scl,yy1*scl,zz1*scl);
	glVertex3d(xx2*scl,yy2*scl,zz2*scl);
	glVertex3d(xx3*scl,yy3*scl,zz3*scl);
	glEnd();
      } }
    glDisable( GL_LIGHTING );  
    glColor3d(0.5, 0.5, 0.5);
    for (int i=1; i<=atom.nwall; i++) {
      int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
      if ((i0!=0)&&(i1!=0)&&(i2!=0)) {
	glBegin(GL_LINES);
	xx1=atom.rx[abs(i0)]*1e9; yy1=atom.ry[abs(i0)]*1e9; zz1=atom.rz[abs(i0)]*1e9;
	xx2=atom.rx[abs(i1)]*1e9; yy2=atom.ry[abs(i1)]*1e9; zz2=atom.rz[abs(i1)]*1e9;
	xx3=atom.rx[abs(i2)]*1e9; yy3=atom.ry[abs(i2)]*1e9; zz3=atom.rz[abs(i2)]*1e9;
	if (zz2>zz1+cube3/2) {zz2=zz2-cube3;} if (zz2<zz1-cube3/2) {zz2=zz2+cube3;}
	if (zz3>zz1+cube3/2) {zz3=zz3-cube3;} if (zz3<zz1-cube3/2) {zz3=zz3+cube3;}
	if (cube3 < 3.2) {
	  if (i1<0) { if (zz2<cube3/2) { zz2 += cube3; } else { zz2 -= cube3; } }
	  if (i2<0) { if (zz3<cube3/2) { zz3 += cube3; } else { zz3 -= cube3; } }
	}
	glVertex3d(xx1*scl,yy1*scl,zz1*scl);glVertex3d(xx2*scl,yy2*scl,zz2*scl);
	glVertex3d(xx2*scl,yy2*scl,zz2*scl);glVertex3d(xx3*scl,yy3*scl,zz3*scl);
	glVertex3d(xx3*scl,yy3*scl,zz3*scl);glVertex3d(xx1*scl,yy1*scl,zz1*scl);
	glEnd();
	if (show_cnt_wallv) {
	  float xx0, yy0, zz0, vx, vy, vz, vscale=vscl;
	  xx0 = (xx1+xx2+xx3)/3; yy0 = (yy1+yy2+yy3)/3; zz0 = (zz1+zz2+zz3)/3;
	  vx = atom.wall_nvec[i][0]; vy = atom.wall_nvec[i][1]; vz = atom.wall_nvec[i][2]; 
	  vx = xx0 + vx*vscale; vy = yy0 + vy*vscale; vz = zz0 + vz*vscale; 
	  glDrawArrowd(xx0*scl,yy0*scl,zz0*scl,vx*scl,vy*scl,vz*scl); }
      } }
  }
  if ((show_cnt_wall_num != 0)&&(show_cnt_wall_num<=atom.nwall)) {
    glEnable( GL_LIGHTING );  
    glColor3d(1.0, 1.0, 1.0);
    float cube3 = (float)cell.hmat[2][2]*1e9;
    int i = show_cnt_wall_num;
    int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
    if ((i0!=0)&&(i1!=0)&&(i2!=0)) {
      glBegin(GL_TRIANGLES);
      xx1=atom.rx[abs(i0)]*1e9; yy1=atom.ry[abs(i0)]*1e9; zz1=atom.rz[abs(i0)]*1e9;
      xx2=atom.rx[abs(i1)]*1e9; yy2=atom.ry[abs(i1)]*1e9; zz2=atom.rz[abs(i1)]*1e9;
      xx3=atom.rx[abs(i2)]*1e9; yy3=atom.ry[abs(i2)]*1e9; zz3=atom.rz[abs(i2)]*1e9;
      if (zz2>zz1+cube3/2) {zz2=zz2-cube3;} if (zz2<zz1-cube3/2) {zz2=zz2+cube3;}
      if (zz3>zz1+cube3/2) {zz3=zz3-cube3;} if (zz3<zz1-cube3/2) {zz3=zz3+cube3;}
      if (cube3 < 3.2) {
	if (i1<0) { if (zz2<cube3/2) { zz2 += cube3; } else { zz2 -= cube3; } }
	if (i2<0) { if (zz3<cube3/2) { zz3 += cube3; } else { zz3 -= cube3; } }
      }
      glVertex3d(xx1*scl,yy1*scl,zz1*scl);
      glVertex3d(xx2*scl,yy2*scl,zz2*scl);
      glVertex3d(xx3*scl,yy3*scl,zz3*scl);
      glEnd();
      float xx0, yy0, zz0, vx, vy, vz, vscale=vscl;
      xx0 = (xx1+xx2+xx3)/3; yy0 = (yy1+yy2+yy3)/3; zz0 = (zz1+zz2+zz3)/3;
      vx = atom.wall_nvec[i][0]; vy = atom.wall_nvec[i][1]; vz = atom.wall_nvec[i][2]; 
      vx = xx0 + vx*vscale; vy = yy0 + vy*vscale; vz = zz0 + vz*vscale; 
      glDrawArrowd(xx0*scl,yy0*scl,zz0*scl,vx*scl,vy*scl,vz*scl);
    }
    glDisable( GL_LIGHTING );  
    glColor3d(0.5, 0.5, 0.5);
    //int i = show_cnt_wall_num;
    //int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
    if ((i0!=0)&&(i1!=0)&&(i2!=0)) {
      glBegin(GL_LINES);
      xx1=atom.rx[abs(i0)]*1e9; yy1=atom.ry[abs(i0)]*1e9; zz1=atom.rz[abs(i0)]*1e9;
      xx2=atom.rx[abs(i1)]*1e9; yy2=atom.ry[abs(i1)]*1e9; zz2=atom.rz[abs(i1)]*1e9;
      xx3=atom.rx[abs(i2)]*1e9; yy3=atom.ry[abs(i2)]*1e9; zz3=atom.rz[abs(i2)]*1e9;
      if (zz2>zz1+cube3/2) {zz2=zz2-cube3;} if (zz2<zz1-cube3/2) {zz2=zz2+cube3;}
      if (zz3>zz1+cube3/2) {zz3=zz3-cube3;} if (zz3<zz1-cube3/2) {zz3=zz3+cube3;}
      if (cube3 < 3.2) {
	if (i1<0) { if (zz2<cube3/2) { zz2 += cube3; } else { zz2 -= cube3; } }
	if (i2<0) { if (zz3<cube3/2) { zz3 += cube3; } else { zz3 -= cube3; } }
      }
      glVertex3d(xx1*scl,yy1*scl,zz1*scl);glVertex3d(xx2*scl,yy2*scl,zz2*scl);
      glVertex3d(xx2*scl,yy2*scl,zz2*scl);glVertex3d(xx3*scl,yy3*scl,zz3*scl);
      glVertex3d(xx3*scl,yy3*scl,zz3*scl);glVertex3d(xx1*scl,yy1*scl,zz1*scl);
      glEnd();
    }
  }
  // draw auxiliary constraint
  if (draw_aux) {
    if (yzplane_punch) {
      drawYZPlane(yzplane_punch_d/2); drawYZPlane(-yzplane_punch_d/2);
    }
  }
  // draw ring (for CNT)
  if (show_cnt_ring) {
    GLUquadricObj *ring = gluNewQuadric();
    gluQuadricDrawStyle(ring, GLU_LINE);
    float cx = (float)cell.hmat[0][0];
    float cy = (float)cell.hmat[1][1];
    float cz = (float)cell.hmat[2][2];
    float rad = (float)cnt_ring_radius*ang;
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, red); 
    glPushMatrix();
    glTranslated( cx/2*1e9*scl, cy/2*1e9*scl, 0);
    //glTranslated( 0,0,0);
    gluCylinder(ring, rad*1e9*scl, rad*1e9*scl, cz*1e9*scl, 50, 20);
    gluDeleteQuadric(ring);
    glPopMatrix();
  }
  if (draw_bond) {
    glDisable( GL_LIGHTING );  
    glColor3d(0.5, 0.5, 0.5);
    float cube1 = (float)cell.hmat[0][0]*1e9;
    float cube2 = (float)cell.hmat[1][1]*1e9;
    float cube3 = (float)cell.hmat[2][2]*1e9;
    float cube1x=(float)cell.hmat[0][0]*1e9;
    float cube1y=(float)cell.hmat[1][0]*1e9;
    float cube1z=(float)cell.hmat[2][0]*1e9;
    float cube2x=(float)cell.hmat[0][1]*1e9;
    float cube2y=(float)cell.hmat[1][1]*1e9;
    float cube2z=(float)cell.hmat[2][1]*1e9;
    float cube3x=(float)cell.hmat[0][2]*1e9;
    float cube3y=(float)cell.hmat[1][2]*1e9;
    float cube3z=(float)cell.hmat[2][2]*1e9;
    if (book.algo == 3) {
      //for (int i=1; i<=bre.np; i++) {
      for (int i=1; i<=atom.natom; i++) {
	if (atom.nneighbor[i] > MAXNEIGHBOR) { continue; }
	for (int j=0; j<atom.nneighbor[i]; j++) {
	  int jn = atom.neighbor[i][j];
	  if (jn > atom.natom) { continue; }
	  xx1=atom.rx[i]*1e9; yy1=atom.ry[i]*1e9; zz1=atom.rz[i]*1e9;
	  xx2=atom.rx[jn]*1e9; yy2=atom.ry[jn]*1e9; zz2=atom.rz[jn]*1e9;
	  int drawflg = 1;
	  if (cell.pbcx) {
	    if (xx2>xx1+cube1/1.5) {
	      xx2=xx2-cube1; if (draw_bond_pbc==0) { drawflg = 0; } }
	    if (xx2<xx1-cube1/1.5) {
	      xx2=xx2+cube1; if (draw_bond_pbc==0) { drawflg = 0; } } }
	  if (cell.pbcy) {
	    if (yy2>yy1+cube2/1.5) {
	      yy2=yy2-cube2; if (draw_bond_pbc==0) { drawflg = 0; } }
	    if (yy2<yy1-cube2/1.5) {
	      yy2=yy2+cube2; if (draw_bond_pbc==0) { drawflg = 0; } } }
	  if (cell.pbcz) {
	    if (zz2>zz1+cube3/1.5) {
	      zz2=zz2-cube3; if (draw_bond_pbc==0) { drawflg = 0; } }
	    if (zz2<zz1-cube3/1.5) {
	      zz2=zz2+cube3; if (draw_bond_pbc==0) { drawflg = 0; } } }
	  if (drawflg ==1) {
	    glBegin(GL_LINES);
	    glVertex3d(xx1*scl,yy1*scl,zz1*scl);glVertex3d(xx2*scl,yy2*scl,zz2*scl);
	    glEnd(); }
	}
      }
    } else if (book.algo = 1) {
      //std::cout<<"[myglutdisplay.cpp]"<<std::endl;
      //std::cout<<"Setting bond using BK algo #1"<<std::endl;
      //std::cout<<"Not supported yet."<<std::endl;
      for (int i=1; i<=atom.natom; i++) {
	if (atom.nneighbor[i] > MAXNEIGHBOR) { continue; }
	for (int j=0; j<atom.nneighbor[i]; j++) {
	  int jn = atom.neighbor[i][j];
	  xx1=atom.rx[i]*1e9; yy1=atom.ry[i]*1e9; zz1=atom.rz[i]*1e9;
	  xx2=atom.rx[jn]*1e9; yy2=atom.ry[jn]*1e9; zz2=atom.rz[jn]*1e9;
	  xx3=atom.rx[jn]-atom.rx[i]; yy3=atom.ry[jn]-atom.ry[i]; zz3=atom.rz[jn]-atom.rz[i];
	  xx4=cell.hinmat[0][0]*xx3 + cell.hinmat[0][1]*yy3 + cell.hinmat[0][2]*zz3;
	  yy4=cell.hinmat[1][0]*xx3 + cell.hinmat[1][1]*yy3 + cell.hinmat[1][2]*zz3;
	  zz4=cell.hinmat[2][0]*xx3 + cell.hinmat[2][1]*yy3 + cell.hinmat[2][2]*zz3;
	  int drawflg = 1; int chk=0;
	  if (cell.pbcx) {
	    //if (xx2>xx1+cube1/2) {xx2=xx2-cube1;chk=1;}
	    //if (xx2<xx1-cube1/2) {xx2=xx2+cube1;chk=1;}
	    if (xx4> 0.5) {xx2-=cube1x;yy2-=cube1y;zz2-=cube1z;chk=1;}
	    if (xx4<-0.5) {xx2+=cube1x;yy2+=cube1y;zz2+=cube1z;chk=1;}
	  } else {
	    //if (xx2>xx1+cube1/2) {drawflg=0;}
	    //if (xx2<xx1-cube1/2) {drawflg=0;} }
	    if (xx4> 0.5) {drawflg=0;}
	    if (xx4<-0.5) {drawflg=0;} }
	  if (cell.pbcy) {
	    //if (yy2>yy1+cube2/2) {yy2=yy2-cube2;chk=1;}
	    //if (yy2<yy1-cube2/2) {yy2=yy2+cube2;chk=1;}
	    if (yy4> 0.5) {xx2-=cube2x;yy2-=cube2y;zz2-=cube2z;chk=1;}
	    if (yy4<-0.5) {xx2+=cube2x;yy2+=cube2y;zz2+=cube2z;chk=1;}
	  } else {
	    //if (yy2>yy1+cube2/2) {drawflg=0;}
	    //if (yy2<yy1-cube2/2) {drawflg=0;} }
	    if (yy4> 0.5) {drawflg=0;}
	    if (yy4<-0.5) {drawflg=0;} }
	  if (cell.pbcz) {
	    //if (zz2>zz1+cube3/2) {zz2=zz2-cube3;chk=1;}
	    //if (zz2<zz1-cube3/2) {zz2=zz2+cube3;chk=1;}
	    if (zz4> 0.5) {xx2-=cube3x;yy2-=cube3y;zz2-=cube3z;chk=1;}
	    if (zz4<-0.5) {xx2+=cube3x;yy2+=cube3y;zz2+=cube3z;chk=1;}
	  } else {
	    //if (zz2>zz1+cube3/2) {drawflg=0;}
	    //if (zz2<zz1-cube3/2) {drawflg=0;} }
	    if (zz4> 0.5) {drawflg=0;}
	    if (zz4<-0.5) {drawflg=0;} }
	  if ((chk==1)&&(draw_bond_pbc==0)) { drawflg = 0; }
	  if (drawflg ==1) {
	    glBegin(GL_LINES);
	    glVertex3d(xx1*scl,yy1*scl,zz1*scl);glVertex3d(xx2*scl,yy2*scl,zz2*scl);
	    glEnd(); }
	}
      }
    } 
  } // end of draw_bond

  glEnable( GL_LIGHTING );  

  // draw force
  if (draw_force) {
    glEnable( GL_LIGHTING );  
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, green2); 
    for (int i=1; i<=atom.natom; i++) {
      float xx0, yy0, zz0, vx, vy, vz, vscale;
      xx0=atom.rx[i]*1e9; yy0=atom.ry[i]*1e9; zz0=atom.rz[i]*1e9;
      vx =atom.fx[i];     vy =atom.fy[i];     vz =atom.fz[i];
      double vlen=sqrt(vx*vx+vy*vy+vz*vz);
      vscale=0.5/(0.1*ev/ang)*vscl_force; //0.1eV/A -> length 0.5
      //if (vlen>0.1*ev/ang) { vscale=0.5/vlen*vscl_force; } //Force over 0.1eV/A -> length 0.5
      //vscale=0.5/atom.Fmax(); //Fmax -> length 0.5
      //vx = xx0 + vx*vscale; vy = yy0 + vy*vscale; vz = zz0 + vz*vscale; 
      vx = vx*vscale; vy = vy*vscale; vz = vz*vscale; 
      //glDrawArrowd(xx0*scl,yy0*scl,zz0*scl,vx*scl,vy*scl,vz*scl);
      glDrawArrowd2(xx0,yy0,zz0,vx,vy,vz,scl);
    }
  } // end of draw force
  if (draw_load) {
    glEnable( GL_LIGHTING );  
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, blue2); 
    for (int i=1; i<=atom.natom; i++) {
      float xx0, yy0, zz0, vx, vy, vz, vscale;
      xx0=atom.rx[i]*1e9; yy0=atom.ry[i]*1e9; zz0=atom.rz[i]*1e9;
      vx =atom.fx_l[i];   vy =atom.fy_l[i];   vz =atom.fz_l[i];
      double vlen=sqrt(vx*vx+vy*vy+vz*vz);
      vscale=0.5/(0.1*ev/ang)*vscl_force; //0.1eV/A -> length 0.5
      //if (vlen>0.1*ev/ang) { vscale=0.5/vlen*vscl_force; } //Force over 0.1eV/A -> length 0.5
      //vscale=0.5/atom.Fmax(); //Fmax -> length 0.5
      //vx = xx0 + vx*vscale; vy = yy0 + vy*vscale; vz = zz0 + vz*vscale; 
      vx = vx*vscale; vy = vy*vscale; vz = vz*vscale; 
      //glDrawArrowd(xx0*scl,yy0*scl,zz0*scl,vx*scl,vy*scl,vz*scl);
      glDrawArrowd2(xx0,yy0,zz0,vx,vy,vz,scl);
    }
  } // end of draw force
  glPopMatrix();

  // draw atoms
  //  float rad = (float)radius/100.0;
  if (mdmotion) {
    change_atom_color();
  }
  if (createconfig_mode>0) {
    if ((iatom_pick>0)&&(iatom_pick<=atom.natom)) {
      for (int i=1; i<atom.natom; i++) {
	if (i==iatom_pick) {
	  memcpy(color[i],green,sizeof(GLfloat)*4);
	} else {
	  memcpy(color[i],color0[i],sizeof(GLfloat)*4);
	}
      }
    }
  }
  for (int i=1; i<=atom.natom; i++) {
    glNewList(objects+i, GL_COMPILE);
    glPushMatrix();
    glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] ); 
    glMultMatrixf( rotate );
    float xx=atom.rx[i]*1e9;
    float yy=atom.ry[i]*1e9;
    float zz=atom.rz[i]*1e9;
    if (atom.QC) {
      int iflg = 0;
      if (edit_elem_mode>0) {
	for (int ii=0; ii<4; ii++) { if ((select_atom[ii]==i)&&(select_atom_repidx[ii]==0)) { iflg=1; } }
      }
      if (show_only_elem) {
	if (atom.repatom[i]==0) {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, white);
	  putSphere(xx, yy, zz, rad, segments, solid);
	}
      } else if (atom.repatom[i]==1) {
	if (iflg==0) {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, red); 
	} else {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color[i]); 
	}
	putSphere(xx, yy, zz, rad, segments, solid);
      } else {
	if (iflg==0) {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, white);
	} else {
	  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color[i]); 
	}
	putSphere(xx, yy, zz, rad, segments, solid);
      }
    } else { //if not atom.QC
      // glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, red);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color[i]);
      putSphere(xx, yy, zz, rad, segments, solid);
    }
    glPopMatrix();
    glEndList();
    glCallList(objects+i);
  } // end of loop for real atoms
  ibase=objects+atom.natom;
  if (draw_replica>0) { // replica atoms
    icnt=0; int ix0, iy0, iz0;
    for (int rbit=1; rbit<=7; rbit++) {
      if (rbit>=4) { ix0 = 1; } else { ix0 = 0; }
      if ((rbit % 4 == 2)||(rbit % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
      if (rbit % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
      if ((draw_replica==2)&&(iy0==1)) { continue; }
      for (int i=1; i<=atom.natom; i++) {
	float xx=(atom.rx[i]+cell.hmat[0][0]*ix0)*1e9;
	float yy=(atom.ry[i]+cell.hmat[1][1]*iy0)*1e9;
	float zz=(atom.rz[i]+cell.hmat[2][2]*iz0)*1e9;
	if ((xx/1e9/cell.hmat[0][0]<1.1)&&(yy/1e9/cell.hmat[1][1]<1.1)&&(zz/1e9/cell.hmat[2][2]<1.1)&&(icnt<atom.natom*2)) {
	  icnt++;
	  glNewList(ibase+icnt, GL_COMPILE);
	  glPushMatrix();
	  glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] ); 
	  glMultMatrixf( rotate );
	  //glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color[ibase+icnt]); //!!!!
	  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color[atom.natom+icnt]); //!!!!
	  iatom[icnt]=i; repidx[icnt]=rbit;
	  putSphere(xx, yy, zz, rad, segments, 0);
	  glPopMatrix(); // fixed 2013.10.30 (previously out of if clause)
	}
	glEndList();
	//	glCallList(ibase+i);
	glCallList(ibase+icnt);
      }
    }
    // Draw lines between selected atoms
    glPushMatrix();
    glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] ); 
    glMultMatrixf( rotate );
    glDisable( GL_LIGHTING ); glColor3d(0.0, 0.0, 0.0); glBegin(GL_LINES);
    for (int ii=0; ii<3; ii++) {
      for (int jj=ii+1; jj<4; jj++) {
	int si=select_atom[ii]; int sj=select_atom[jj];
	if ((si!=0)&&(sj!=0)) {
	  int ix0, iy0, iz0, rbit;
	  rbit = select_atom_repidx[ii];
	  if (rbit>=4) { ix0 = 1; } else { ix0 = 0; }
	  if ((rbit % 4 == 2)||(rbit % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
	  if (rbit % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
	  float xx1=(atom.rx[si]+cell.hmat[0][0]*ix0)*1e9;
	  float yy1=(atom.ry[si]+cell.hmat[1][1]*iy0)*1e9;
	  float zz1=(atom.rz[si]+cell.hmat[2][2]*iz0)*1e9;
	  rbit = select_atom_repidx[jj];
	  if (rbit>=4) { ix0 = 1; } else { ix0 = 0; }
	  if ((rbit % 4 == 2)||(rbit % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
	  if (rbit % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
	  float xx2=(atom.rx[sj]+cell.hmat[0][0]*ix0)*1e9;
	  float yy2=(atom.ry[sj]+cell.hmat[1][1]*iy0)*1e9;
	  float zz2=(atom.rz[sj]+cell.hmat[2][2]*iz0)*1e9;
	  glVertex3d(xx1*scl,yy1*scl,zz1*scl); glVertex3d(xx2*scl,yy2*scl,zz2*scl);
	}
      }
    }
    glEnd(); glEnable( GL_LIGHTING );  
    glPopMatrix();
  } // end of draw_replica

  glLoadIdentity();
  gluLookAt(0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

  glutSwapBuffers(); 
  //glui->sync_live(); 
  GLUI_Master.sync_live_all();

}

void glDrawArrowd(float x0, float y0, float z0,
		  float x1, float y1, float z1)
{
  GLUquadricObj *arrows[2];
  float x2, y2, z2, len, angle;
  x2 = x1-x0; y2 = y1-y0; z2 = z1-z0;
  len = sqrt(x2*x2 + y2*y2 + z2*z2);
  if(len != 0.0){
    angle = acos(z2*len/(sqrt(x2*x2+y2*y2+z2*z2)*len))/M_PI*180.0;
    glPushMatrix();
    glTranslated( x0, y0, z0);
    glRotated( angle, -y2*len, x2*len, 0.0);
    arrows[0] = gluNewQuadric();
    gluQuadricDrawStyle(arrows[0], GLU_FILL);
    //gluCylinder(arrows[0], len/80, len/80, len*0.9, 8, 8);
    gluCylinder(arrows[0], 1e-2, 1e-2, len*0.9, 8, 8);
    glPushMatrix();
    glTranslated( 0.0, 0.0, len*0.9);
    arrows[1] = gluNewQuadric();
    gluQuadricDrawStyle(arrows[1], GLU_FILL);
    gluCylinder(arrows[1], len/30, 0.0f, len/10, 8, 8);
    gluDeleteQuadric(arrows[0]); gluDeleteQuadric(arrows[1]);
    glPopMatrix(); //glPopMatrix();
  }
}
void glDrawArrowd2(float x, float y, float z,
		   float vx, float vy, float vz, float u)
{
  GLUquadricObj *arrows[2];
  float len = sqrt(vx*vx + vy*vy + vz*vz)*u;
  if(len != 0.0){
    float angle = acos(vz*u/len)/M_PI*180.0;
    glPushMatrix();
    glTranslated( x*u, y*u, z*u);
    glRotated( angle, -vy, vx, 0.0);
    arrows[0] = gluNewQuadric();
    gluQuadricDrawStyle(arrows[0], GLU_FILL);
    //gluCylinder(arrows[0], len/80, len/80, len*0.9, 8, 8);
    gluCylinder(arrows[0], 1e-2, 1e-2, len*0.9, 8, 8);
    glPushMatrix();
    glTranslated( 0.0, 0.0, len*0.9);
    arrows[1] = gluNewQuadric();
    gluQuadricDrawStyle(arrows[1], GLU_FILL);
    //gluCylinder(arrows[1], len/30, 0.0f, len/10, 8, 8);
    gluCylinder(arrows[1], 3e-2, 0.0f, len/10, 8, 8);
    gluDeleteQuadric(arrows[0]); gluDeleteQuadric(arrows[1]);
    glPopMatrix(); //glPopMatrix();
  }
}

// Draw Pipe
void glDrawPiped(float x0, float y0, float z0,
 float x1, float y1, float z1,
 float r)
{
	GLUquadricObj *arrow;
	float x2, y2, z2, len, angle;

	x2 = x1-x0; y2 = y1-y0; z2 = z1-z0;
	len = sqrt(x2*x2 + y2*y2 + z2*z2);
	if(len != 0.0){
		angle = acos(z2*len/(sqrt(x2*x2+y2*y2+z2*z2)*len))/M_PI*180.0;

		glPushMatrix();
			glTranslated( x0, y0, z0);
			glRotated( angle, -y2*len, x2*len, 0.0f);
			arrow = gluNewQuadric();
			gluQuadricDrawStyle(arrow, GLU_FILL);
			gluCylinder(arrow, r, r, len, 8, 8);
			gluDeleteQuadric(arrow);
		glPopMatrix();
	}
}

void glDrawAxisd(float length) // Draw Coordinate Axis
{
	GLUquadricObj *arrows[3];

	// Draw X-axis
	glColor3ub(255, 0, 0);
	glBegin(GL_LINES);
		glVertex3d(-length*scl*0.5, -0.1*scl, -0.1*scl);
		glVertex3d( length*scl, -0.1*scl, -0.1*scl);
	glEnd();
	glPushMatrix();
		arrows[0] = gluNewQuadric();
		gluQuadricDrawStyle(arrows[0], GLU_FILL);
		glTranslated(length*scl, -0.1*scl, -0.1*scl);
		glRotated(90.0f, 0,1,0);
		gluCylinder(arrows[0], length/30*scl, 0.0f, length/5*scl, 8, 8);
		gluDeleteQuadric(arrows[0]);
	glPopMatrix();

	// Draw Y-axis
	glColor3ub(  0,255, 0);
	glBegin(GL_LINES);
		glVertex3d( -0.1*scl,-length*scl*0.5, -0.1*scl);
		glVertex3d( -0.1*scl, length*scl, -0.1*scl);
	glEnd();
	glPushMatrix();
		arrows[1] = gluNewQuadric();
		gluQuadricDrawStyle(arrows[1], GLU_FILL);
		glTranslated(-0.1*scl, length*scl, -0.1*scl);
		glRotated(-90.0f, 1,0,0);
		gluCylinder(arrows[1], length/30*scl, 0.0f, length/5*scl, 8, 8);
		gluDeleteQuadric(arrows[1]);
	glPopMatrix();

	// Draw Z-axis
	glColor3ub(  0, 0,255);
	glBegin(GL_LINES);
		glVertex3d( -0.1*scl, -0.1*scl,-length*scl*0.5);
		glVertex3d( -0.1*scl, -0.1*scl, length*scl);
	glEnd();
	glPushMatrix();
		arrows[2] = gluNewQuadric();
		gluQuadricDrawStyle(arrows[2], GLU_FILL);
		glTranslated(-0.1*scl, -0.1*scl, length*scl);
		gluCylinder(arrows[2], length/30*scl, 0.0f, length/5*scl, 8, 8);
		gluDeleteQuadric(arrows[2]);
	glPopMatrix();

	glColor3ub(255, 255, 255);
	glDisable(GL_DEPTH_TEST);
	//glPrint3d(length*1.1*scl, 0.0, 0.0, (void *)font, "X");
	//glPrint3d(0.0, length*1.1*scl, 0.0, (void *)font, "Y");
	//glPrint3d(0.0, 0.0, length*1.1*scl, (void *)font, "Z");
	glEnable(GL_DEPTH_TEST);

}
