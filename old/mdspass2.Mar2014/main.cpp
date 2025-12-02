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

#define OPEN_FILE_ID 100
#define OPEN_FILE_QCELEMENT_ID 103
#define FILE_CLOSE_ID 101
#define QCELEMENT_CLOSE_ID 104
#define RESET_ID 102
#define WRITECONFIG_ID 1061
#define CREATECONFIG_ID 1062
#define CREATECONFIG_CLOSE_ID 1063
#define CREATECONFIG_DO_ID 1064
#define MULTIPLYCELL_DO_ID 1069
#define REMOVEATOM_DO_ID 1070
#define ADDATOM_DO_ID 1071
#define MOVEATOM_DO_ID 1072
#define CNTWALL_DO_ID 1065
#define CNTWALL_READ_ID 1066
#define CNTWALL_WRITE_ID 1067
#define CNTWALL_CLOSE_ID 1068
#define WRITEQCELEMENT_ID 116
#define INST_ID 105
//#define SETDEXDT_ID 107
//#define SETDEXDT_CLOSE_ID 108
#define MDSWITCH_ID 109
#define EDIT_ELEM_ID 110
#define EDIT_ELEM_CLOSE_ID 111
#define EDIT_ELEM_TAKE_ID 112
#define EDIT_ELEM_XZ_ID 113
#define EDIT_ELEM_XZ_CLOSE_ID 114
#define EDIT_ELEM_XZ_TAKE_ID 115
#define MEASURE_ID 301
#define MEASURE_CLOSE_ID 302
#define MEASURE_CALC2_ID 303
#define MEASURE_CALC3_ID 304
#define MEASURE_CALC4_ID 305
#define MEASURE_CLEAR_ID 306
#define SETPARAM_ID 200
#define SETPARAM_CLOSE_ID 201
#define SETPARAM_READ_ID 204
#define CELLSIZE_APPLY_ID 202
#define CAPTURE_ID 203
#define STRS_ID 205
#define STRS_CLOSE_ID 206
#define STRSCHK_ID 207
#define STRSSET_ID 211
#define EXTRA_ID 212
#define EXTRA_CLOSE_ID 213
#define PHONON_CALC_ID 214
#define NEB_CALC_ID 215
#define CNTSHELL_TEST_ID 216
#define CNTSHELL_CALC_ID 217
#define CALC_ID 208
#define POTFILE_ID 209
#define POTFILE_CLOSE_ID 210
#define CB_BOND 1001
#define CB_ENSEMBLE 1004
#define CB_CFG_FB 1005
#define CB_QCELM_FB 1006
#define CB_CNTWALL_FB 1007
#define CB_COLOR_MODE 1008
#define CB_POTFILE_FB 1009
#define CB_INST 1010
#define CB_PHONON 1011
#define CB_EDITCONFIG 1012
#define CB_NEB 1013
#define CB_NEB_SHOW 1014
#define CB_ROTATE 1015
#define CB_ROTATE_REV 1016
#define CB_RESET_VIEW 1017
#define CB_CELL_DIM 1018
#define CB_SLICE 1019
#define CB_FIRE 1020
#define CB_CG 1023
#define CB_REL 1024
#define CB_ELASC 1021
#define CB_FCHECK 1022
#define CB_MARKED_ATOM 1025
//using namespace std;
#include "myheader.h"
#include "myclass.h"
int istep = 0, istep0 = 0;
int mdmotion = 0;
int mdspeed = 10;
double dt = 1.0e-15;
int confwrint = 0, autocap = 0;
int itolfor = 1; float tolfor = 0.02;
int itolstep = 0; int tolstep = 10000;
int itolstress = 0; float tolstress = 10.0;
int irecipe = 0;
int inogui = 0;
//double rcut = 10.3176e-10; // cut off radius, in angstrom
//double rcut = 13.0e-10; // cut off radius, in angstrom
//double rcut2 = rcut * rcut;
double rcut, rcut2;
float temp_set = 100.0;
float tempc = 0.0, cellx,celly,cellz,f_max,epotatom,dtm;
float cell1x,cell1y,cell1z;
float cell2x,cell2y,cell2z;
float cell3x,cell3y,cell3z;
float slice_1min=0,slice_1max=1,slice_2min=0,slice_2max=1,slice_3min=0,slice_3max=1;
int ifix_atoms;
FILE *energyfile = fopen("energy.d","w");
FILE *stressfile = fopen("stress.d","w");
FILE *ssnorfile = fopen("ssnormal.d","w");
FILE *cellfile = fopen("cell.d","w");
FILE *ecnorfile = fopen("ecnormal.d","w");
FILE *ecallfile = fopen("ecall.d","w");
float strs_xx,strs_yy,strs_zz,strs_xy,strs_yz,strs_zx;
float strs_set_xx,strs_set_yy,strs_set_zz,strs_set_xy,strs_set_yz,strs_set_zx;
int cellfix_xx, cellfix_yy, cellfix_zz, cellfix_xy, cellfix_yz, cellfix_zx;
float rcut_f, frcmar = 0.2, frc_f;
float fp_alph_ini=0.1, fp_ffinc=1.1, fp_ffdec=0.5, fp_ffalph=0.99, fp_ffdtmax=10.0;
int fp_nfmin = 5;

Atom atom;
Cell cell;
Book book;
Tersoff tersoff;
GEAM geam;
Eammis eammis;
Adp adp;
Dipole dipole;
Cond cond;
//Bre *bre = new Bre;
Bre bre;
Fire fire;
Cg cg;

#ifdef CUDA
int *arrayDint;
int *arrayDrepatom;
double *arrayDrx, *arrayDry, *arrayDrz, *arrayDfx, *arrayDfy, *arrayDfz;
double *arrayDepot, *arrayDhmat, arrayHhmat[9];
int *arrayDalistnum, *arrayDalist;
int iarray[(NMAX+1)*(NBOOK+1)*4];
#endif

/*
void Atom::Distance(int i, int j)
{

}
*/

void out_of_cell();
void md(), md_set();
//void vscale();
void bookkeep();
//void e_force();
void readconfig(const char* fname); 
void readconfig(const char* fname, int mode); 
void readqcelement(const char* fname); 
void readqcelement(); 
void writedata();
void writedata_initialize();
void writeconfig(const char* fname); void writeconfig_abs(const char* fname);
void writeposcar(const char* fname);
void writeqcelement(const char* fname);
void readsetdat(const char* fname);
//void instability();
//void instability_dsy();
void instability_atomcell();
void instability_atom();
void instability_atom_noremove();
void instability_QC();
void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3]);
void myGlutReshape(int x, int y);
void idle( void );
void createconfig();
void cnt_wall_set();
void cnt_wall_read(const char* fname); void cnt_wall_write(const char* fname);
void cnt_wall_discard();
void myGlutDisplay(void);
void bond_set();
void stretch_celladjust(float x, float y, float z);
void stretch_celladjust(float c1x, float c1y, float c1z,
			float c2x, float c2y, float c2z,
			float c3x, float c3y, float c3z);
void slice(float r1min, float r1max, float r2min, float r2max, float r3min, float r3max);
void rotate_cell(float x, float y, float z, int reverse);
void relax_fire_reset();
void relax_cg_reset();
void get_first_arg(std::string &line, std::string &arg1);
int count_arg_number(std::string line);
#ifndef NOPNG
void capture();
#endif
void phonon_calc();
void cntshell_calc_test(double dmu, double dnu, int n);
void cntshell_calc_1(); void cntshell_calc_2(); void cntshell_calc_all();
void neb_calc(int neb_num, const char* fini, const char* fend);
void stresscheck(double eps_strschk);
void stress_set();
void force_check(double fcheck_disp);
void calc_elastic_const();
void potential_set (int control);
void change_atom_color();
double csp(int i, int mi);
void potential();
void loading();
void set_potfile();
void multiply_cell(int ix, int iy, int iz);
void remove_atom(int ia); void move_atom(int ia, float arx, float ary, float arz);
void move_atom(int ia, const char *aasp, float arx, float ary, float arz);
void set_atom_color(); void set_atom_color(int ia);
void add_atom(const char *aasp, float arx, float ary, float arz);
void calc_distance(int select_atom[], int select_atom_repidx[]);
void calc_angle(int select_atom[], int select_atom_repidx[]);
void calc_dihedral(int select_atom[], int select_atom_repidx[]);
void (*integral_a[NINTEGRAL])(), (*integral_b[NINTEGRAL])();
int integral_type;

float xy_aspect;
int   last_x, last_y;
float rotationX = 0.0, rotationY = 0.0;
float rotate[16] = {
  1,0,0,0,
  0,1,0,0,
  0,0,1,0,
  0,0,0,1
};
float obj_pos[] = {0.0, 0.0, 0.0};
float scl = 1.0;
float vscl = 0.2, vscl_force = 1.0;
int   main_window;


GLfloat red[] = {1.0, 0.0, 0.0, 1.0};
GLfloat yellow[] = {1.0, 1.0, 0.0, 1.0};
GLfloat white[] = {0.8, 0.8, 0.8, 1.0};
GLfloat gray[] = { 0.7, 0.7, 0.7 };
GLfloat blue[] = { 0.1, 0.1, 0.9 };
GLfloat blue2[] = { 0.0, 0.0, 1.0 };
GLfloat green[] = { 0.0, 1.0, 0.0 };
GLfloat green2[] = { 0.0, 0.7, 0.0 };
GLfloat black[] = {0.0, 0.0, 0.0, 1.0};

//GLfloat *color[NOBJECTS];
GLfloat **color, **color0;
//int iatom[NREPLICA]; int repidx[NREPLICA]; int ibase, icnt;
int *iatom; int *repidx; int ibase, icnt;


/** These are the live variables passed into GLUI ***/
int   ortho = 1;
int   wireframe = 0;
int   draw_bond = 0, draw_bond_pbc = 0, draw_force = 0, draw_load = 0;
int   ensemble = 0;
int   relax_algo = 0, relax_accel = 1, relax_accel_interval = 10;
float relax_accel_threshold = 0.2;
int   config_type = 0;
char  config_atom[3] = "Al";
int   irepx=1,irepy=1,irepz=1;
int   icntm=8,icntn=0;
float cscnt = 50.0, rotz = 0.0, shiftz = 0.0;
float alat = 4.0;
int   imerge = 0;
char  current_config_name[100], current_qcelement_name[100];
char  cwdname[80] = "aaa";
const char *potstring_list[] = {"Morse", "GEAM", "Tersoff", "Brenner", "EAM Mishin", "Dipole", "ADP", "AIREBO" };
const char *ensemble_list[] = {"NVE", "NVT", "Relaxation (atom)", "NPH", "NPT", "NPH + PRdamper", "NPT + PRdamper", "Full relaxation"};
const char *kpstring_list[] = {"fcc/dia", "bcc", "sc", "hcp", "chain"};
const char *color_list[] = {"Red", "Blue", "Black", "White"};
int   ipottype = 0;
int   notrans = 0;
int   incell = 0;
int   segments = 8;
int   radius = 8;
int   show_only_elem = 0;
int   show_axis = 0;
int   show_cell = 1;
int   show_cnt_wall = 0, show_cnt_wallv = 0, show_cnt_ring = 0;
int   show_cnt_wall_num = 1;
int   mode_cnt_corrugation = 0, read_cntwall = 0;
int   cnt_load_algo = 0;
float dexdt = 0.0, ex = 0.0;
float deydt = 0.0, ey = 0.0;
float dezdt = 0.0, ez = 0.0;
int repeat_lz; float repeat_lz_min, repeat_lz_max;
float cnt_pressure = 0.1, cnt_pressure_ftot = 0, cnt_pressure_gpa;
float outermost_radius_f, cylinder_side_area_f;
float cnt_ring_radius = 15.0, cnt_ring_fmax = 10.0, cnt_ring_sharpness = 2.0;
float cntshell_dmu = 0.0, cntshell_dnu = 0.0, cntshell_eps = 0.0001;
int cntshell_n = 2;
int   ievec = 0;
int   ievec_num = 1;
float evec_len = 10.0, eigval = 0;
double mat[3][3];
float org_x = 0.0, org_y = 0.0, org_z = 0.0;
int size_w, size_h;
int b_state = 1;
float xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, xx4, yy4, zz4;
GLuint objects;
int edit_elem_mode = 0, select_atom[10], select_atom_repidx[10];
int measure_mode = 0;
int draw_replica = 0;
int hoge = 0;
float bondlength = 1.52; // in ang
int capture_count = 0;
float eps_strschk = 0.001;
float fcheck_disp = 0.0001;
int color_mode = 0, color_mode_auto = 0;
float color_mode_vmin, color_mode_vmax;
float prmass_scale = 0;
bool iprdamper = false, ivscale = false, irelax = false;
float prdamper_val1 = 1, prdamper_val2 = 0, prlimit = 0;
int hessian_read = 0, hessian_write = 0, inst_mode;
int phonon_rep = 2, phonon_kp = 0, phonon_knum = 50;
int neb_num = 10, neb_node = 0, neb_ite = 100, neb_init_read = 0;
float neb_tol_fac = -11;
int iatom_pick = 0, createconfig_mode = 0;
float atomrx, atomry, atomrz;
int milX1=1, milX2=0, milX3=0;
int milY1=0, milY2=1, milY3=0;
int milZ1=0, milZ2=0, milZ3=1;
float cellrot_x=0,cellrot_y=0,cellrot_z=0;
int marked_atom = 0, marked_atom_color = 0;
int initial_velocity = 0;

// Using a std::string as a live variable is safe.
std::string text = "Test string..";
std::string potfile;

GLUI *glui, *filename_config_glui=0;
GLUI *filename_qcelement_glui=0;
GLUI *filename_cntwall_glui=0;
GLUI *filename_potfile_glui=0;
//GLUI *setdexdt_glui=0;
GLUI *edit_elem_glui=0, *measure_glui=0;
GLUI *setparam_glui=0, *strs_glui=0, *extra_glui=0;
GLUI *createconfig_glui=0;
GLUI_Checkbox   *checkbox;
//GLUI_Spinner    *spinner;
GLUI_Checkbox   *checkbox_pbcx;
GLUI_Checkbox   *checkbox_pbcy;
GLUI_Checkbox   *checkbox_pbcz;
GLUI_Checkbox   *checkbox_QC;
GLUI_Checkbox   *checkbox_Trans;
GLUI_Checkbox   *checkbox_show_only_elem;
//GLUI_Checkbox   *checkbox_axis;
GLUI_Checkbox   *checkbox_fatom;
GLUI_Checkbox   *checkbox_evec;
//GLUI_Checkbox   *checkbox_ortho;
GLUI_Spinner    *spinner_evec_num;
GLUI_Spinner    *spinner_evec_len;
GLUI_Spinner    *spinner_instcenter_num;
GLUI_RadioGroup *radio;
GLUI_RadioGroup *createconfig_radio;
GLUI_EditText   *edittext;
GLUI_CommandLine *filename_config, *filename_qcelement, *filename_cntwall;
//GLUI_CommandLine *setdexdt;
GLUI_Button *open_file_btn, *open_file_qcelement_btn, *cntwall_btn;
GLUI_Button *inst_btn;
GLUI_Button *reset_btn;
GLUI_Button *writeconfig_btn;
GLUI_Button *createconfig_btn;
GLUI_Button *writeqcelement_btn;
//GLUI_Button *setdexdt_btn;
GLUI_Button *mdswitch_btn;
GLUI_Button *edit_elem_btn;
GLUI_Button *edit_elem_xz_btn;
GLUI_Button *measure_btn;
GLUI_Button *setparam_btn, *strs_btn, *calc_btn, *potfile_btn, *extra_btn;
GLUI_FileBrowser *config_fb, *qcelm_fb, *cntwall_fb, *potfile_fb;
GLUI_EditText *status_lx, *status_ly, *status_lz, *status_dt;

/* GLUI control callback                                                 */
void control_cb( int control )
{
  //printf( "callback: %d\n", control );
  /*
  printf( "             checkbox: %d\n", checkbox->get_int_val() );
  printf( "              spinner: %d\n", spinner->get_int_val() );
  printf( "             spinner2: %d\n", spinner2->get_int_val() );
  printf( "          radio group: %d\n", radio->get_int_val() );
  printf( "             ensemble: %d\n", ensemble );
  printf( "                 text: %s\n", edittext->get_text() );
  */
  //  std::cout<<ensemble<<std::endl;
  //  if (ievec_num>20) { ievec_num = 20; }
  //printf("control_cb %d\n",control);
  glutPostRedisplay();
  if (control == 0) { myGlutReshape( size_w, size_h ); }
  if (show_cnt_wall_num > atom.nwall) { show_cnt_wall_num = atom.nwall; GLUI_Master.sync_live_all(); }
  if (control == CB_INST) {
    if (ievec_num < 1) { ievec_num = 1; }
    if (ievec_num > MAXMODE) { ievec_num = MAXMODE; }
    eigval = atom.eigval[ievec_num-1];
  }
  if (control == CB_PHONON) {
    //if (phonon_rep > 4) { phonon_rep = 4; }
    if (phonon_rep < 1) { phonon_rep = 1; }
  }
  if (control == CB_NEB) {
    if (neb_num < 3)  { neb_num = 3;  }
    if (neb_num > 99) { neb_num = 99; }
  }
  if (control == CB_NEB_SHOW) {
    if (neb_node >= neb_num) { neb_node = neb_num - 1; }
    if (neb_node < 0)       { neb_node = 0;       }
    char filepath[80] = "CONFIG.NEB."; char numc[10];
    if (neb_node == neb_num - 1) { strcpy(numc,"END");
    } else if (neb_node == 0) { strcpy(numc,"INI");
    } else { snprintf(numc, sizeof(numc), "%02d", neb_node); }
    strcat(filepath,numc);
    readconfig(filepath, 1);
  }
  if (control == CB_EDITCONFIG) {
    if (iatom_pick > atom.natom) { iatom_pick = atom.natom; }
    if (iatom_pick < 0 ) { iatom_pick = 0; }
    if ((iatom_pick >0)&&(iatom_pick <= atom.natom)) {
      strcpy(config_atom,atom.asp[iatom_pick]);
      atomrx=atom.rx[iatom_pick]/ang; atomry=atom.ry[iatom_pick]/ang; atomrz=atom.rz[iatom_pick]/ang;
    }
  }
  if (control == CB_BOND) { bond_set(); }
  if (control == CB_COLOR_MODE) { change_atom_color(); }
  if (control == CB_MARKED_ATOM) {
    if (marked_atom < 0) { marked_atom = 0; }
    else if (marked_atom > atom.natom) { marked_atom = atom.natom; }
    change_atom_color(); }
  if (control == CB_ENSEMBLE) { md_set(); }//if (ensemble == 1) { velocity_random(temp_set); } }
  /* Relocated to pointer_cb...
  if (control == CB_CFG_FB) {
    char fname[60] = "aaa"; std::string fnames; fnames = config_fb->get_file();
    strcpy (fname, fnames.c_str());
    readconfig( fname ); readqcelement( current_qcelement_name );
    glutPostRedisplay();
    writedata_initialize(); }
  if (control == CB_QCELM_FB) {
    std::string text = qcelm_fb->get_file();
    char fname[60] = "aaa";
    strcpy (fname, text.c_str());
    readqcelement( fname );  }
  if (control == CB_CNTWALL_FB) {
    std::string text = cntwall_fb->get_file();
    char fname[60] = "aaa";
    strcpy (fname, text.c_str());
    cnt_wall_read( fname ); }
  */
  if (control == CB_POTFILE_FB) {
    potfile = cwdname;
#if defined __linux__ || defined __APPLE__
    potfile += "/pot/";
#else
    potfile += "\\pot\\";
#endif
    potfile += potfile_fb->get_file();
    set_potfile();
    /*
    if (ipottype == 5) {
      printf("Dipole parameter file is set to %s\n",potfile.c_str());
      strcpy(dipole.fname,potfile.c_str());
      dipole.initialize = true;
    }
    */
  }
  if (control == CB_ROTATE) { rotate_cell(cellrot_x,cellrot_y,cellrot_z,0); }
  if (control == CB_ROTATE_REV) { rotate_cell(cellrot_x,cellrot_y,cellrot_z,1); }
  if (control == CB_RESET_VIEW) {
    org_x=0.0; org_y=0.0; org_z=0.0;
    rotationX=0.0; rotationY=0.0;
    for (int i=0; i<=15; i++) { rotate[i]=0; }
    rotate[0]=1;rotate[6]=-1;rotate[9]=1;rotate[15]=1; }
  if (control == CB_CELL_DIM) {
    stretch_celladjust(cell1x,cell1y,cell1z,cell2x,cell2y,cell2z,cell3x,cell3y,cell3z);
    glui->sync_live();
  }
  if (control == CB_SLICE) {
    slice(slice_1min,slice_1max,slice_2min,slice_2max,slice_3min,slice_3max);
    glui->sync_live();
  }
  if ((control == CB_FIRE)||(control == CB_REL)) {
    relax_fire_reset(); relax_cg_reset();
    glui->sync_live();
  }
  if (control == CB_ELASC) {
    calc_elastic_const();
  }
  if (control == CB_FCHECK) {
    force_check((double)fcheck_disp);
  }
}

void pointer_cb ( GLUI_Control* control )
{
  if (control->get_id() == OPEN_FILE_ID) {
    filename_config_glui = GLUI_Master.create_glui( "Enter filename:", 0, 600, 150 );
    filename_config = new GLUI_CommandLine( filename_config_glui, "Config file:", NULL, -1, pointer_cb );
    filename_config->set_w( 300 );
    //config_fb = new GLUI_FileBrowser(filename_config_glui, "", false, CB_CFG_FB, control_cb);
    config_fb = new GLUI_FileBrowser(filename_config_glui, "", false, CB_CFG_FB, pointer_cb);
    config_fb->set_w(300);
    GLUI_Panel *panel = new GLUI_Panel(filename_config_glui, "", GLUI_PANEL_NONE);
    new GLUI_Checkbox(panel, "Merge", &imerge, 1, control_cb );
    new GLUI_Column(panel, false );
    new GLUI_Button(panel, "Close", FILE_CLOSE_ID, pointer_cb);
    filename_config_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  // Relocated from control_cb...
  else if (control->get_id() == CB_CFG_FB) {
    char fname[60] = "aaa"; std::string fnames; fnames = config_fb->get_file();
    strcpy (fname, fnames.c_str());
    readconfig( fname ); readqcelement( current_qcelement_name );
    glutPostRedisplay();
    writedata_initialize(); 
    open_file_btn->enable(); control->glui->close(); // After reading, pop-up window closes
  }
  else if ( control->get_id() == CB_QCELM_FB) {
    std::string text = qcelm_fb->get_file();
    char fname[60] = "aaa";
    strcpy (fname, text.c_str());
    readqcelement( fname ); open_file_qcelement_btn->enable(); control->glui->close(); }
  else if ( control->get_id() == CB_CNTWALL_FB) {
    std::string text = cntwall_fb->get_file();
    char fname[60] = "aaa";
    strcpy (fname, text.c_str());
    cnt_wall_read( fname ); cntwall_btn->enable(); control->glui->close(); }
   //
  else if ( control->get_id() == FILE_CLOSE_ID ) {
    open_file_btn->enable();
    control->glui->close();
  }
  else if (control->get_id() == OPEN_FILE_QCELEMENT_ID) {
    filename_qcelement_glui = GLUI_Master.create_glui( "Enter filename:", 0, 600, 150 );
    filename_qcelement = new GLUI_CommandLine( filename_qcelement_glui, "File:", NULL, -1, pointer_cb );
    filename_qcelement->set_w( 300 );
    qcelm_fb = new GLUI_FileBrowser(filename_qcelement_glui, "", false, CB_QCELM_FB, pointer_cb);
    qcelm_fb->set_w(300);
    GLUI_Panel *panel_qcelement = new GLUI_Panel(filename_qcelement_glui, "", GLUI_PANEL_NONE);
    new GLUI_Button(panel_qcelement, "Close", QCELEMENT_CLOSE_ID, pointer_cb);
    filename_qcelement_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == QCELEMENT_CLOSE_ID ) {
    open_file_qcelement_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == RESET_ID ) {
    //    org_x=0.0; org_y=0.0; org_z=0.0;
    //    rotationX=0.0; rotationY=0.0;
    //    for (int i=0; i<=15; i++) { rotate[i]=0; }
    //    rotate[0]=1;rotate[5]=1;rotate[10]=1;rotate[15]=1;
    relax_fire_reset(); relax_cg_reset();
    for (int i=1; i<=atom.natom; i++) {
      atom.rx[i]=atom.rx_org[i];
      atom.ry[i]=atom.ry_org[i];
      atom.rz[i]=atom.rz_org[i];
      atom.vx[i]=0.0;atom.vy[i]=0.0;atom.vz[i]=0.0;
      atom.fx[i]=0.0;atom.fy[i]=0.0;atom.fz[i]=0.0;
      atom.ax[i]=0.0;atom.ay[i]=0.0;atom.az[i]=0.0;
      atom.bx[i]=0.0;atom.by[i]=0.0;atom.bz[i]=0.0;
      atom.cx[i]=0.0;atom.cy[i]=0.0;atom.cz[i]=0.0;
    }
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	cell.hmat[i][j]=cell.hmat_org[i][j];
      }
    }
    //cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
    cellx=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
    celly=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
    cellz=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
    cell1x=cell.hmat[0][0]/ang;cell1y=cell.hmat[1][0]/ang;cell1z=cell.hmat[2][0]/ang;
    cell2x=cell.hmat[0][1]/ang;cell2y=cell.hmat[1][1]/ang;cell2z=cell.hmat[2][1]/ang;
    cell3x=cell.hmat[0][2]/ang;cell3y=cell.hmat[1][2]/ang;cell3z=cell.hmat[2][2]/ang;
    cell.Reset(); cell.volume = cell.Getvolume();
    // give 'dummy' value for istep and then reset it (to sync display)
    istep = -1; glui->sync_live();
    istep =  0; glui->sync_live();
    ex = 0.0; dexdt = 0.0; ey = 0.0; deydt = 0.0; ez = 0.0; dezdt = 0.0;
    GLUI_Master.sync_live_all();
    bond_set();
    writedata_initialize();
  }
  else if ( control->get_id() == INST_ID ) {
    if (atom.QC) {
      instability_QC();
    } else {
      //if (inst_mode == 0) {
      //instability();
      //} else if (inst_mode == 1) {
      //instability_dsy();
      //} else if (inst_mode == 2) {
      if (inst_mode == 0) {
	instability_atomcell();
      } else if (inst_mode == 1) {
	instability_atom();
      } else if (inst_mode == 2) {
	instability_atom_noremove();
      }
    }
  }
  else if ( control->get_id() == WRITECONFIG_ID ) {
    chdir(cwdname);
    writeconfig("CONFIG.OUT"); writeconfig_abs("CONFIG.OUT.ABS");
    writeposcar("POSCAR");
  }
  else if ( control->get_id() == WRITEQCELEMENT_ID ) {
    chdir(cwdname);
    writeqcelement("QCELEMENT.OUT");
  }
  else if ( control->get_id() == MDSWITCH_ID ) {
    if (mdmotion == 0) { glutPostRedisplay(); mdmotion = 1; istep0 = istep; }
    else { glutPostRedisplay();mdmotion = 0; }
  }
  else if ( control->get_id() == STRSCHK_ID ) {
    stresscheck((double)eps_strschk);
  }
  else if ( control->get_id() == PHONON_CALC_ID ) {
    phonon_calc();
  }
  else if ( control->get_id() == NEB_CALC_ID ) {
    neb_calc(neb_num, "CONFIG.NEB.INI", "CONFIG.NEB.END");
  }
  else if ( control->get_id() == CNTSHELL_TEST_ID ) {
    cntshell_calc_test((double)cntshell_dmu, (double)cntshell_dnu, cntshell_n);
  }
  else if ( control->get_id() == CNTSHELL_CALC_ID ) {
    cntshell_calc_1();
    //cntshell_calc_2();
  }
  else if ( control->get_id() == POTFILE_ID ) {
    filename_potfile_glui = GLUI_Master.create_glui("Potential parameter filename", 0, 600, 150);
    GLUI_Listbox *potlist3 = new GLUI_Listbox(filename_potfile_glui, "Potential", &ipottype,0,potential_set);
    for (int i=0; i<MAXPOTTYPE; i++) { potlist3->add_item(i, potstring_list[i]); }
#ifdef __GNUC__
    chdir("pot");
#endif
#ifdef _WIN32
    //SetCurrentDirectory("pot");
	chdir("pot");
#endif
    potfile_fb = new GLUI_FileBrowser(filename_potfile_glui, "", false, CB_POTFILE_FB, control_cb);
    potfile_fb->set_w(300);
    GLUI_Panel *panel_potfile = new GLUI_Panel(filename_potfile_glui, "", GLUI_PANEL_NONE);
    new GLUI_Button(panel_potfile, "Close", POTFILE_CLOSE_ID, pointer_cb);
    filename_potfile_glui->set_main_gfx_window(main_window);
    control->disable();
  }
  else if ( control->get_id() == POTFILE_CLOSE_ID ) {
#ifdef __GNUC__
    chdir(cwdname);
#endif
#ifdef _WIN32
    //SetCurrentDirectory(cwdname);
	chdir(cwdname);
#endif
    potfile_btn->enable();
    control->glui->close();
  }
  else if ( control == filename_config ) {
    std::string text = filename_config->get_text();
    char fname[60] = "aaa";
    strcpy (fname, text.c_str());
    readconfig( fname ); readqcelement( current_qcelement_name );
    glutPostRedisplay();
    writedata_initialize();
  }
  else if ( control == filename_qcelement ) {
    std::string text = filename_qcelement->get_text();
    char fname[60] = "aaa";
    strcpy (fname, text.c_str());
    readqcelement( fname );
  }
  else if ( control == filename_cntwall ) {
    std::string text = filename_cntwall->get_text();
    char fname[60] = "aaa";
    strcpy (fname, text.c_str());
    cnt_wall_read( fname );
  }
  //else if ( control == setdexdt ) {
  //  std::string text = setdexdt->get_text();
  //  char line[60] = "aaa"; strcpy (line, text.c_str());
  //  printf("Set dex/dt = %e ang\n",atof(line));
  //  dexdt = atof(line);
  //}
  else if ( control->get_id() == MEASURE_ID ) {
    measure_glui = GLUI_Master.create_glui("Measure tool", 0, 700, 50);
    for (int i=0; i<10; i++) { select_atom[i]=0; select_atom_repidx[i]=0; }
    measure_mode=1;
    draw_replica=1;
    glutPostRedisplay();
    GLUI_Panel *panel_measure = new GLUI_Panel(measure_glui, "", GLUI_PANEL_NONE);
    new GLUI_EditText(panel_measure, "1:", &select_atom[0], 5, control_cb );
    new GLUI_EditText(panel_measure, "1-idx:", &select_atom_repidx[0], 5, control_cb );
    new GLUI_EditText(panel_measure, "2:", &select_atom[1], 5, control_cb );
    new GLUI_EditText(panel_measure, "2-idx:", &select_atom_repidx[1], 5, control_cb );
    new GLUI_EditText(panel_measure, "3:", &select_atom[2], 5, control_cb );
    new GLUI_EditText(panel_measure, "3-idx:", &select_atom_repidx[2], 5, control_cb );
    //new GLUI_EditText(panel_measure, "4:", &select_atom[3], 5, control_cb );
    //new GLUI_EditText(panel_measure, "4-idx:", &select_atom_repidx[3], 5, control_cb );
    new GLUI_Button(panel_measure, "Distance 1-2", MEASURE_CALC2_ID, pointer_cb);
    new GLUI_Button(panel_measure, "Angle 1-2-3", MEASURE_CALC3_ID, pointer_cb);
    //new GLUI_Button(panel_measure, "Dihedral 1-2-3-4", MEASURE_CALC4_ID, pointer_cb);
    new GLUI_Button(panel_measure, "Clear", MEASURE_CLEAR_ID, pointer_cb);
    new GLUI_Button(panel_measure, "Close", MEASURE_CLOSE_ID, pointer_cb);
    measure_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == MEASURE_CLEAR_ID ) {
    for (int i=0; i<3; i++) { select_atom[i] = 0; select_atom_repidx[i] = 0; }
    set_atom_color();
  }
  else if ( control->get_id() == MEASURE_CALC2_ID ) {
    calc_distance(select_atom, select_atom_repidx);
    //GLUI_Master.sync_live_all(); glutPostRedisplay();
  }
  else if ( control->get_id() == MEASURE_CALC3_ID ) {
    calc_angle(select_atom, select_atom_repidx);
    //GLUI_Master.sync_live_all(); glutPostRedisplay();
  }
  else if ( control->get_id() == MEASURE_CALC4_ID ) {
    calc_dihedral(select_atom, select_atom_repidx);
    //GLUI_Master.sync_live_all(); glutPostRedisplay();
  }
  else if ( control->get_id() == MEASURE_CLOSE_ID ) {
    measure_mode=0;
    draw_replica=0;
    set_atom_color();
    //for (int i=1; i<=atom.natom+icnt; i++) { memcpy(color[i],yellow,sizeof(GLfloat)*4); }
    glutPostRedisplay();
    measure_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == EDIT_ELEM_ID ) {
    edit_elem_glui = GLUI_Master.create_glui("Edit element", 0, 700, 50);
    for (int i=0; i<10; i++) { select_atom[i]=0; select_atom_repidx[i]=0; }
    edit_elem_mode=1;
    draw_replica=1;
    glutPostRedisplay();
    GLUI_Panel *panel_edit_elem = new GLUI_Panel(edit_elem_glui, "", GLUI_PANEL_NONE);
    //    new GLUI_EditText(panel_edit_elem, "1", GLUI_EDITTEXT_INT, NULL, 1, 
    new GLUI_EditText(panel_edit_elem, "1:", &select_atom[0], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "1-idx:", &select_atom_repidx[0], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "2:", &select_atom[1], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "2-idx:", &select_atom_repidx[1], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "3:", &select_atom[2], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "3-idx:", &select_atom_repidx[2], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "4:", &select_atom[3], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "4-idx:", &select_atom_repidx[3], 5, control_cb );
    new GLUI_Button(panel_edit_elem, "Take this", EDIT_ELEM_TAKE_ID, pointer_cb);
    new GLUI_Button(panel_edit_elem, "Close", EDIT_ELEM_CLOSE_ID, pointer_cb);
    edit_elem_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == EDIT_ELEM_TAKE_ID ) {
    if ((select_atom[0]!=0)&&(select_atom[1]!=0)&&(select_atom[2]!=0)&&(select_atom[3]!=0)) {
      printf("New element: %d %d %d %d  %d %d %d %d\n",
	     select_atom[0],select_atom[1],select_atom[2],select_atom[3],
	     select_atom_repidx[0],select_atom_repidx[1],select_atom_repidx[2],select_atom_repidx[3]);
      int i = atom.nelem+1;
      atom.elem_v     = (int **)realloc(atom.elem_v,     sizeof(int*)*(i+1));
      atom.elem_v_rep = (int **)realloc(atom.elem_v_rep, sizeof(int*)*(i+1));
      atom.elem_v[i]     = (int *)malloc(sizeof(int)*5);
      atom.elem_v_rep[i] = (int *)malloc(sizeof(int)*5);
      atom.elem_v[i][1]=select_atom[0]; atom.elem_v_rep[i][1]=select_atom_repidx[0];
      atom.elem_v[i][2]=select_atom[1]; atom.elem_v_rep[i][2]=select_atom_repidx[1];
      atom.elem_v[i][3]=select_atom[2]; atom.elem_v_rep[i][3]=select_atom_repidx[2];
      atom.elem_v[i][4]=select_atom[3]; atom.elem_v_rep[i][4]=select_atom_repidx[3];
      atom.nelem = i;
      readqcelement();
      for (int ii=0; ii<10; ii++) { select_atom[ii]=0; select_atom_repidx[ii]=0; }
      for (int ii=1; ii<=atom.natom+icnt; ii++) { memcpy(color[ii],yellow,sizeof(GLfloat)*4); }
      GLUI_Master.sync_live_all();
      glutPostRedisplay();
    }
  }
  else if ( control->get_id() == EDIT_ELEM_CLOSE_ID ) {
    edit_elem_mode=0;
    draw_replica=0;
    //for (int i=1; i<=atom.natom+icnt; i++) { memcpy(color[i],yellow,sizeof(GLfloat)*4); }
    glutPostRedisplay();
    edit_elem_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == EDIT_ELEM_XZ_ID ) {
    edit_elem_glui = GLUI_Master.create_glui("Edit element (xz-plane mode)", 0, 700, 50);
    for (int i=0; i<10; i++) { select_atom[i]=0; select_atom_repidx[i]=0; }
    edit_elem_mode=2;
    draw_replica=2;
    glutPostRedisplay();
    GLUI_Panel *panel_edit_elem = new GLUI_Panel(edit_elem_glui, "", GLUI_PANEL_NONE);
    new GLUI_EditText(panel_edit_elem, "1:", &select_atom[0], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "1-idx:", &select_atom_repidx[0], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "2:", &select_atom[1], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "2-idx:", &select_atom_repidx[1], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "3:", &select_atom[2], 5, control_cb );
    new GLUI_EditText(panel_edit_elem, "3-idx:", &select_atom_repidx[2], 5, control_cb );
    new GLUI_Button(panel_edit_elem, "Take this", EDIT_ELEM_XZ_TAKE_ID, pointer_cb);
    new GLUI_Button(panel_edit_elem, "Close", EDIT_ELEM_XZ_CLOSE_ID, pointer_cb);
    edit_elem_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == EDIT_ELEM_XZ_TAKE_ID ) {
    if ((select_atom[0]!=0)&&(select_atom[1]!=0)&&(select_atom[2]!=0)) {
      //      printf("New element: %d %d %d %d  %d %d %d %d\n",
      //	     select_atom[0],select_atom[1],select_atom[2],select_atom[3],
      //	     select_atom_repidx[0],select_atom_repidx[1],select_atom_repidx[2],select_atom_repidx[3]);
      int i = atom.nelem+3;
      atom.elem_v     = (int **)realloc(atom.elem_v,     sizeof(int*)*(i+1));
      atom.elem_v_rep = (int **)realloc(atom.elem_v_rep, sizeof(int*)*(i+1));
      for (int ii=atom.nelem+1; ii<=atom.nelem+3; ii++) {
	atom.elem_v[ii]     = (int *)malloc(sizeof(int)*5);
	atom.elem_v_rep[ii] = (int *)malloc(sizeof(int)*5);
      }
      i = atom.nelem+1;
      atom.elem_v[i][1]=select_atom[0]; atom.elem_v_rep[i][1]=select_atom_repidx[0];
      atom.elem_v[i][2]=select_atom[1]; atom.elem_v_rep[i][2]=select_atom_repidx[1];
      atom.elem_v[i][3]=select_atom[2]; atom.elem_v_rep[i][3]=select_atom_repidx[2];
      atom.elem_v[i][4]=select_atom[2]; atom.elem_v_rep[i][4]=select_atom_repidx[2]+2;
      i = atom.nelem+2;
      atom.elem_v[i][1]=select_atom[0]; atom.elem_v_rep[i][1]=select_atom_repidx[0];
      atom.elem_v[i][2]=select_atom[1]; atom.elem_v_rep[i][2]=select_atom_repidx[1];
      atom.elem_v[i][3]=select_atom[1]; atom.elem_v_rep[i][3]=select_atom_repidx[1]+2;
      atom.elem_v[i][4]=select_atom[2]; atom.elem_v_rep[i][4]=select_atom_repidx[2]+2;
      i = atom.nelem+3;
      atom.elem_v[i][1]=select_atom[0]; atom.elem_v_rep[i][1]=select_atom_repidx[0];
      atom.elem_v[i][2]=select_atom[0]; atom.elem_v_rep[i][2]=select_atom_repidx[0]+2;
      atom.elem_v[i][3]=select_atom[1]; atom.elem_v_rep[i][3]=select_atom_repidx[1]+2;
      atom.elem_v[i][4]=select_atom[2]; atom.elem_v_rep[i][4]=select_atom_repidx[2]+2;
      atom.nelem = i;
      readqcelement();
      for (int ii=0; ii<10; ii++) { select_atom[ii]=0; select_atom_repidx[ii]=0; }
      for (int ii=1; ii<=atom.natom+icnt; ii++) { memcpy(color[ii],yellow,sizeof(GLfloat)*4); }
      GLUI_Master.sync_live_all();
      glutPostRedisplay();
    }
  }
  else if ( control->get_id() == EDIT_ELEM_XZ_CLOSE_ID ) {
    edit_elem_mode=0;
    draw_replica=0;
    for (int i=1; i<=atom.natom+icnt; i++) { memcpy(color[i],yellow,sizeof(GLfloat)*4); }
    glutPostRedisplay();
    edit_elem_xz_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == SETPARAM_ID ) {
    setparam_glui = GLUI_Master.create_glui("Set parameters", 0, 700, 50);
    glutPostRedisplay();
    GLUI_Rollout *setparam_panel010 = new GLUI_Rollout(setparam_glui, "Setup for drawing", false);
    GLUI_Panel *setparam_panel011 = new GLUI_Panel(setparam_panel010, "", GLUI_PANEL_NONE);
    new GLUI_Checkbox(setparam_panel011, "Wireframe", &wireframe, 1, control_cb );
    new GLUI_Checkbox(setparam_panel011, "Ortho", &ortho, 0, control_cb);
    new GLUI_Column(setparam_panel011, false );
    new GLUI_Checkbox(setparam_panel011, "Show axis", &show_axis);
    new GLUI_Checkbox(setparam_panel011, "Show cell", &show_cell);
    new GLUI_Column(setparam_panel011, false );
    GLUI_Spinner *spinner2  = new GLUI_Spinner(setparam_panel011, "Radius:", &radius, 5, control_cb );
    spinner2->set_int_limits( 1, 20 );
    GLUI_Panel *setparam_panel012 = new GLUI_Panel(setparam_panel010, "", GLUI_PANEL_NONE);
    new GLUI_Checkbox(setparam_panel012, "Draw bond", &draw_bond, CB_BOND, control_cb);
    new GLUI_Checkbox(setparam_panel012, "Bond-PBC", &draw_bond_pbc, CB_BOND, control_cb);
    GLUI_EditText *spinner5 = new GLUI_EditText(setparam_panel012,"Bond length",GLUI_EDITTEXT_FLOAT,&bondlength,CB_BOND,control_cb);
    spinner5->set_float_limits( 0.0, 10.0 );
    new GLUI_Column(setparam_panel012, false );
    new GLUI_Checkbox(setparam_panel012, "Draw Force", &draw_force, 10, control_cb);
    new GLUI_Checkbox(setparam_panel012, "Draw Load", &draw_load, 10, control_cb);
    new GLUI_EditText(setparam_panel012, "F-arw length", &vscl_force);
    GLUI_Panel *setparam_panel0121 = new GLUI_Panel(setparam_panel010, "", GLUI_PANEL_NONE);
    GLUI_Spinner *spinner3 = new GLUI_Spinner(setparam_panel0121, "Redraw itvl:", &mdspeed);
    spinner3->set_int_limits( 1, 50 );
    GLUI_Rollout *setparam_panel014 = new GLUI_Rollout(setparam_glui, "Atom color", false);
    GLUI_RadioGroup *color_radio
      = new GLUI_RadioGroup(setparam_panel014,&color_mode,CB_COLOR_MODE,control_cb );
    new GLUI_RadioButton(color_radio, "Species");
    new GLUI_RadioButton(color_radio, "Energy");
    new GLUI_RadioButton(color_radio, "CSP");
    new GLUI_Checkbox(setparam_panel014, "Autorange", &color_mode_auto,CB_COLOR_MODE,control_cb);
    new GLUI_Column(setparam_panel014, false );
    new GLUI_EditText(setparam_panel014, "min", &color_mode_vmin,CB_COLOR_MODE,control_cb);
    new GLUI_EditText(setparam_panel014, "max", &color_mode_vmax,CB_COLOR_MODE,control_cb);
    new GLUI_Column(setparam_panel014, false );
    new GLUI_Spinner(setparam_panel014, "Marked atom", &marked_atom,CB_MARKED_ATOM,control_cb);
    GLUI_Listbox *colorlist = new GLUI_Listbox(setparam_panel014, "Color      ", &marked_atom_color,CB_MARKED_ATOM,control_cb);
    for (int i=0; i<4; i++) { colorlist->add_item(i, color_list[i]); }
    GLUI_Rollout *setparam_panel013 = new GLUI_Rollout(setparam_glui, "Cutoff/Bookkeep", false);
    new GLUI_EditText(setparam_panel013, "Cutoff (A)", &rcut_f);
    new GLUI_Column(setparam_panel013, false );
    new GLUI_EditText(setparam_panel013, "B-keep margin (A)", &frcmar);
    new GLUI_EditText(setparam_panel013, "B-keep cutoff (A)", &frc_f);
    new GLUI_Column(setparam_panel013, false );
    GLUI_Spinner *spinner4 = new GLUI_Spinner(setparam_panel013, "B-keep itvl:", &book.nbk);
    spinner4->set_int_limits( 5, 500 );
    //GLUI_Panel *setparam_panel02 = new GLUI_Panel(setparam_glui, "", GLUI_PANEL_NONE );
    //new GLUI_Column(setparam_panel02, false );
    GLUI_Panel *setparam_panel015 = new GLUI_Panel(setparam_glui, "", GLUI_PANEL_NONE);
    GLUI_Panel *setparam_panel020 = new GLUI_Panel(setparam_panel015, "Relax algo");
    GLUI_RadioGroup *relax_radio = new GLUI_RadioGroup(setparam_panel020,&relax_algo,CB_REL,control_cb );
    new GLUI_RadioButton(relax_radio, "GLOC" );
    new GLUI_RadioButton(relax_radio, "FIRE" );
    new GLUI_RadioButton(relax_radio, "CG" );
    new GLUI_Checkbox(setparam_panel020, "Accel(Gloc)", &relax_accel);
    new GLUI_EditText(setparam_panel020, "Interval", &relax_accel_interval);
    new GLUI_EditText(setparam_panel020, "Threshold", &relax_accel_threshold);
    new GLUI_Column(setparam_panel015, false );
    new GLUI_EditText(setparam_panel015, "CONF Wr itvl", &confwrint);
    new GLUI_Checkbox(setparam_panel015, "Auto capture", &autocap);
    new GLUI_Column(setparam_panel015, false );
    new GLUI_Checkbox(setparam_panel015, "Follow recipe", &irecipe);
    new GLUI_Checkbox(setparam_panel015, "Stop MD by Fmax", &itolfor);
    new GLUI_EditText(setparam_panel015, "Fmax tlrc (eV/A)", &tolfor);
    new GLUI_Checkbox(setparam_panel015, "Stop MD by step", &itolstep);
    new GLUI_EditText(setparam_panel015, "# of steps", &tolstep);
    new GLUI_Checkbox(setparam_panel015, "Stop MD by stress", &itolstress);
    new GLUI_EditText(setparam_panel015, "Stress tlrc (MPa)", &tolstress);
    GLUI_Rollout *setparam_panel0150 = new GLUI_Rollout(setparam_glui, "Fire param",false);
    new GLUI_EditText(setparam_panel0150, "alpha init", &fp_alph_ini,CB_FIRE,control_cb);
    new GLUI_EditText(setparam_panel0150, "min accel steps", &fp_nfmin,CB_FIRE,control_cb);
    new GLUI_Column(setparam_panel0150, false );
    new GLUI_EditText(setparam_panel0150, "dt inc (>1)", &fp_ffinc,CB_FIRE,control_cb);
    new GLUI_EditText(setparam_panel0150, "dt dec (<1)", &fp_ffdec,CB_FIRE,control_cb);
    new GLUI_Column(setparam_panel0150, false );
    new GLUI_EditText(setparam_panel0150, "al rate (<1)", &fp_ffalph,CB_FIRE,control_cb);
    new GLUI_EditText(setparam_panel0150, "dt max (fs)", &fp_ffdtmax,CB_FIRE,control_cb);
    GLUI_Rollout *setparam_panel021 = new GLUI_Rollout(setparam_glui, "Deformation settings", false);
    GLUI_Panel *setparam_panel0211 = new GLUI_Panel(setparam_panel021, "", GLUI_PANEL_NONE);
    GLUI_EditText *text_setdexdt = new GLUI_EditText(setparam_panel0211, "ex(/ps):", &dexdt);
    text_setdexdt->set_float_limits(-1.0,1.0);
    new GLUI_Column(setparam_panel0211, false );
    GLUI_EditText *text_setdeydt = new GLUI_EditText(setparam_panel0211, "ey(/ps):", &deydt);
    text_setdeydt->set_float_limits(-1.0,1.0);
    new GLUI_Column(setparam_panel0211, false );
    GLUI_EditText *text_setdezdt = new GLUI_EditText(setparam_panel0211, "ez(/ps):", &dezdt);
    text_setdezdt->set_float_limits(-1.0,1.0);
    GLUI_Panel *setparam_panel0212 = new GLUI_Panel(setparam_panel021, "", GLUI_PANEL_NONE);
    new GLUI_Checkbox(setparam_panel0212, "Repeat Lz", &repeat_lz);
    new GLUI_Column(setparam_panel0212, false );
    GLUI_EditText *text_lz_min = new GLUI_EditText(setparam_panel0212, "Lz(min):", &repeat_lz_min);
    new GLUI_Column(setparam_panel0212, false );
    GLUI_EditText *text_lz_max = new GLUI_EditText(setparam_panel0212, "Lz(max):", &repeat_lz_max);
    /*
    GLUI_Panel *setparam_panel0212 = new GLUI_Panel(setparam_panel021, "", GLUI_PANEL_NONE);
    GLUI_EditText *text_setcellx = new GLUI_EditText(setparam_panel0212, "Lx", &cellx);
    new GLUI_Column(setparam_panel0212, false );
    GLUI_EditText *text_setcelly = new GLUI_EditText(setparam_panel0212, "Ly", &celly);
    new GLUI_Column(setparam_panel0212, false );
    GLUI_EditText *text_setcellz = new GLUI_EditText(setparam_panel0212, "Lz", &cellz);
    new GLUI_Button(setparam_panel0212, "Apply", CELLSIZE_APPLY_ID, pointer_cb);
    */
    GLUI_Rollout *setparam_panel03 = new GLUI_Rollout(setparam_glui, "Special settings for CNT", false);
    new GLUI_Checkbox(setparam_panel03, "Corrugation mode", &mode_cnt_corrugation);
    new GLUI_Checkbox(setparam_panel03, "Show CNT wall", &show_cnt_wall);
    new GLUI_Checkbox(setparam_panel03, "Show CNT n-vec", &show_cnt_wallv);
    new GLUI_Checkbox(setparam_panel03, "Show CNT ring", &show_cnt_ring);
    new GLUI_EditText(setparam_panel03, "n-vec length", &vscl);
    GLUI_Panel *corr_radio_panel = new GLUI_Panel(setparam_panel03, "Loading type" );
    GLUI_RadioGroup *corr_radio = new GLUI_RadioGroup(corr_radio_panel,&cnt_load_algo,0,control_cb );
    new GLUI_RadioButton(corr_radio, "Wall" );
    new GLUI_RadioButton(corr_radio, "Ring" );
    new GLUI_Column(setparam_panel03, false );
    GLUI_Panel *corr_wall_panel = new GLUI_Panel(setparam_panel03, "Wall" );
    GLUI_Spinner *spinner6 = new GLUI_Spinner(corr_wall_panel, "CNT wall #", &show_cnt_wall_num,5,control_cb);
    spinner6->set_speed( 0.01 ); //spinner6->set_int_limits( 0, NMAX*3 );
    GLUI_Spinner *spinner7 = new GLUI_Spinner(corr_wall_panel, "Pressure", &cnt_pressure,5,control_cb);
    spinner7->set_speed( 0.01 );
    GLUI_Panel *corr_ring_panel = new GLUI_Panel(setparam_panel03, "Ring" );
    new GLUI_EditText(corr_ring_panel, "Radius(A)", &cnt_ring_radius);
    new GLUI_EditText(corr_ring_panel, "Fmax(eV/A)",   &cnt_ring_fmax);
    new GLUI_EditText(corr_ring_panel, "Sharpness", &cnt_ring_sharpness);
    //spinner7->set_float_limits( 0.0, 100.0);

    GLUI_Panel *setparam_panel10 = new GLUI_Panel(setparam_glui, "", GLUI_PANEL_NONE );
    new GLUI_Button(setparam_panel10, "Read SETDAT", SETPARAM_READ_ID, pointer_cb);
    new GLUI_Column(setparam_panel10, false );
    potfile_btn = new GLUI_Button(setparam_panel10, "Pot file", POTFILE_ID, pointer_cb );
    new GLUI_Column(setparam_panel10, false );
    new GLUI_Button(setparam_panel10, "Close", SETPARAM_CLOSE_ID, pointer_cb);
    setparam_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == CALC_ID ) {
    bookkeep(); potential(); loading();
    // For output in "Status" area
    //cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
    cellx=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
    celly=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
    cellz=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
    cell1x=cell.hmat[0][0]/ang;cell1y=cell.hmat[1][0]/ang;cell1z=cell.hmat[2][0]/ang;
    cell2x=cell.hmat[0][1]/ang;cell2y=cell.hmat[1][1]/ang;cell2z=cell.hmat[2][1]/ang;
    cell3x=cell.hmat[0][2]/ang;cell3y=cell.hmat[1][2]/ang;cell3z=cell.hmat[2][2]/ang;
    f_max=atom.Fmax()/ev*ang;
    epotatom=atom.epotsum/ev/atom.natom;
  }
  else if ( control->get_id() == STRS_ID ) {
    strs_glui = GLUI_Master.create_glui("Stress", 0, 700, 50);
    glutPostRedisplay();
    GLUI_Panel *strs_panel01 = new GLUI_Panel(strs_glui, "", GLUI_PANEL_NONE);
    new GLUI_StaticText( strs_panel01, "Stress in MPa" );
    GLUI_Panel *strs_panel011 = new GLUI_Panel(strs_panel01, "", GLUI_PANEL_NONE);
    GLUI_EditText *str1 = new GLUI_EditText(strs_panel011, "xx", &strs_xx);
    new GLUI_Column(strs_panel011, false ); str1->set_w(150);
    GLUI_EditText *str2 = new GLUI_EditText(strs_panel011, "yy", &strs_yy);
    new GLUI_Column(strs_panel011, false );str2->set_w(150);
    GLUI_EditText *str3 = new GLUI_EditText(strs_panel011, "zz", &strs_zz);
    str3->set_w(150);
    GLUI_Panel *strs_panel012 = new GLUI_Panel(strs_panel01, "", GLUI_PANEL_NONE);
    GLUI_EditText *str4 = new GLUI_EditText(strs_panel012, "xy", &strs_xy);
    new GLUI_Column(strs_panel012, false ); str4->set_w(150);
    GLUI_EditText *str5 = new GLUI_EditText(strs_panel012, "yz", &strs_yz);
    new GLUI_Column(strs_panel012, false ); str5->set_w(150);
    GLUI_EditText *str6 = new GLUI_EditText(strs_panel012, "zx", &strs_zx);
    str6->set_w(150);
    GLUI_Panel *strs_panel015 = new GLUI_Panel(strs_panel01, "", GLUI_PANEL_NONE);
    new GLUI_EditText(strs_panel015, "Stress check with de", &eps_strschk);
    new GLUI_Column(strs_panel015, false );
    new GLUI_Button(strs_panel015, "Check", STRSCHK_ID, pointer_cb);
    GLUI_Rollout *strs_panel017 = new GLUI_Rollout(strs_panel01, "PR Setting (Stress in MPa)", true);
    new GLUI_EditText(strs_panel017, "PR mass x 10^", &prmass_scale);

    GLUI_Panel *strs_panel018 = new GLUI_Panel(strs_panel017, "", GLUI_PANEL_NONE);
    GLUI_EditText *str1set = new GLUI_EditText(strs_panel018, "xx", &strs_set_xx, STRSSET_ID, pointer_cb);
    new GLUI_Column(strs_panel018, false ); str1set->set_w(150);
    GLUI_EditText *str2set = new GLUI_EditText(strs_panel018, "yy", &strs_set_yy, STRSSET_ID, pointer_cb);
    new GLUI_Column(strs_panel018, false ); str2set->set_w(150);
    GLUI_EditText *str3set = new GLUI_EditText(strs_panel018, "zz", &strs_set_zz, STRSSET_ID, pointer_cb);
    str3set->set_w(150);
    GLUI_Panel *strs_panel019 = new GLUI_Panel(strs_panel017, "", GLUI_PANEL_NONE);
    GLUI_EditText *str4set = new GLUI_EditText(strs_panel019, "xy", &strs_set_xy, STRSSET_ID, pointer_cb);
    new GLUI_Column(strs_panel019, false ); str4set->set_w(150);
    GLUI_EditText *str5set = new GLUI_EditText(strs_panel019, "yz", &strs_set_yz, STRSSET_ID, pointer_cb);
    new GLUI_Column(strs_panel019, false ); str5set->set_w(150);
    GLUI_EditText *str6set = new GLUI_EditText(strs_panel019, "zx", &strs_set_zx, STRSSET_ID, pointer_cb);
    str6set->set_w(150);
    GLUI_Panel *strs_panel020 = new GLUI_Panel(strs_panel017, "", GLUI_PANEL_NONE);
    new GLUI_EditText(strs_panel020, "Damper", &prdamper_val1);
    new GLUI_Column(strs_panel020, false );
    //new GLUI_EditText(strs_panel020, "Damper2", &prdamper_val2);
    //new GLUI_Column(strs_panel020, false );
    new GLUI_EditText(strs_panel020, "Hv Limit", &prlimit);
    GLUI_Rollout *strs_panel0170 = new GLUI_Rollout(strs_panel01, "Cell constraint", true);
    GLUI_Panel *strs_panel0180 = new GLUI_Panel(strs_panel0170, "", GLUI_PANEL_NONE);
    new GLUI_Checkbox(strs_panel0180, "xx", &cellfix_xx); new GLUI_Column(strs_panel0180, false);
    new GLUI_Checkbox(strs_panel0180, "yy", &cellfix_yy); new GLUI_Column(strs_panel0180, false);
    new GLUI_Checkbox(strs_panel0180, "zz", &cellfix_zz); new GLUI_Column(strs_panel0180, false);
    //GLUI_Panel *strs_panel0190 = new GLUI_Panel(strs_panel0170, "", GLUI_PANEL_NONE);
    new GLUI_Checkbox(strs_panel0180, "xy", &cellfix_xy); new GLUI_Column(strs_panel0180, false);
    new GLUI_Checkbox(strs_panel0180, "yz", &cellfix_yz); new GLUI_Column(strs_panel0180, false);
    new GLUI_Checkbox(strs_panel0180, "zx", &cellfix_zx);

    GLUI_Rollout *strs_panel016 = new GLUI_Rollout(strs_panel01, "CNT corrugation", true);
    GLUI_EditText *cnt_text1 = new GLUI_EditText(strs_panel016, "Load(nN)", &cnt_pressure_ftot);
    GLUI_EditText *cnt_text2 = new GLUI_EditText(strs_panel016, "Outermost CNT radius (A)", &outermost_radius_f);
    new GLUI_Column(strs_panel016, false );
    GLUI_EditText *cnt_text3 = new GLUI_EditText(strs_panel016, "Pressure (GPa)", &cnt_pressure_gpa);
    GLUI_EditText *cnt_text4 = new GLUI_EditText(strs_panel016, "Outermost CNT area (A^2)", &cylinder_side_area_f);
    cnt_text1->set_w(130); cnt_text2->set_w(210); cnt_text3->set_w(150); cnt_text4->set_w(220);
    GLUI_Panel *strs_panel021 = new GLUI_Panel(strs_panel01, "", GLUI_PANEL_NONE);
    new GLUI_Button(strs_panel021, "Calc elastic coeff.", CB_ELASC, control_cb);
    new GLUI_Button(strs_panel01, "Close", STRS_CLOSE_ID, pointer_cb);
    strs_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == STRS_CLOSE_ID ) {
    strs_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == EXTRA_ID ) {
    extra_glui = GLUI_Master.create_glui("Extra analysis", 0, 700, 50);
    glutPostRedisplay();
    GLUI_Panel *measure_panel01 = new GLUI_Panel(extra_glui, "", GLUI_PANEL_NONE);
    measure_btn = new GLUI_Button(measure_panel01, "Measure", MEASURE_ID, pointer_cb );
    GLUI_Panel *fcheck_panel01 = new GLUI_Panel(extra_glui, "Force check");
    new GLUI_EditText(fcheck_panel01, "disp", &fcheck_disp, CB_FCHECK, control_cb);
    new GLUI_Button(fcheck_panel01, "Check", CB_FCHECK, control_cb);
    GLUI_Panel *phonon_panel01 = new GLUI_Panel(extra_glui, "Phonon");
    new GLUI_EditText(phonon_panel01, "# of replica", &phonon_rep, CB_PHONON, control_cb);
    GLUI_Listbox *kplist = new GLUI_Listbox(phonon_panel01,"k path",&phonon_kp,CB_PHONON,control_cb);
    for (int i=0; i<MAXKP; i++) { kplist->add_item(i, kpstring_list[i]); }
    new GLUI_EditText(phonon_panel01, "# of k-points", &phonon_knum, CB_PHONON, control_cb);
    new GLUI_Button(phonon_panel01, "Calc", PHONON_CALC_ID, pointer_cb);
    GLUI_Panel *neb_panel01 = new GLUI_Panel(extra_glui, "NEB");
    new GLUI_Checkbox(neb_panel01, "Init conf from file", &neb_init_read);
    new GLUI_EditText(neb_panel01, "# of nodes", &neb_num,CB_NEB,control_cb);
    new GLUI_EditText(neb_panel01, "Max iter", &neb_ite,CB_NEB,control_cb);
    new GLUI_EditText(neb_panel01, "tol 10^", &neb_tol_fac,CB_NEB,control_cb);
    new GLUI_Button(neb_panel01, "Calc", NEB_CALC_ID, pointer_cb);
    new GLUI_Spinner(neb_panel01, "Node to show", &neb_node,CB_NEB_SHOW,control_cb);
    GLUI_Panel *cntshell_panel01 = new GLUI_Panel(extra_glui, "CNT SHELL ANALYSIS");
    new GLUI_EditText(cntshell_panel01, "n", &cntshell_n);
    new GLUI_EditText(cntshell_panel01, "d_mu", &cntshell_dmu);
    new GLUI_EditText(cntshell_panel01, "d_nu", &cntshell_dnu);
    new GLUI_Button(cntshell_panel01, "Dfm test", CNTSHELL_TEST_ID, pointer_cb);
    new GLUI_EditText(cntshell_panel01, "eps", &cntshell_eps);
    new GLUI_Button(cntshell_panel01, "Calc", CNTSHELL_CALC_ID, pointer_cb);
    new GLUI_Button(extra_glui, "Close", EXTRA_CLOSE_ID, pointer_cb);
    extra_glui->set_main_gfx_window( main_window );
    control->disable();
  } 
  else if ( control->get_id() == EXTRA_CLOSE_ID ) {
    extra_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == SETPARAM_READ_ID ) {
    readsetdat("SETDAT");
  }
  else if ( control->get_id() == SETPARAM_CLOSE_ID ) {
    setparam_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == CELLSIZE_APPLY_ID ) {
    stretch_celladjust(cellx,celly,cellz);
    glui->sync_live();
  }
  else if ((control == status_lx)||(control == status_ly)||(control == status_lz)) {
    stretch_celladjust(cellx,celly,cellz);
    glui->sync_live();
  }
  else if (control == status_dt) {
    dt = (double)dtm*1e-15;
  }
  else if ( control->get_id() == STRSSET_ID ) {
    stress_set();
  }
  else if (control->get_id() == CREATECONFIG_ID) {
    createconfig_glui = GLUI_Master.create_glui("Create config", 0, 700, 50);
    glutPostRedisplay();
    createconfig_mode = 1;
    GLUI_Panel *createconfig_panel01 = new GLUI_Panel( createconfig_glui, "", GLUI_PANEL_NONE);
    GLUI_EditText *createconfig_asp = new GLUI_EditText(createconfig_panel01,"Atom", config_atom);
    GLUI_Listbox *potlist2 = new GLUI_Listbox(createconfig_panel01, "Potential", &ipottype,0,potential_set);
    for (int i=0; i<MAXPOTTYPE; i++) { potlist2->add_item(i, potstring_list[i]); }
    createconfig_radio = new GLUI_RadioGroup( createconfig_panel01,&config_type,0,control_cb );
    new GLUI_RadioButton( createconfig_radio, "FCC" );
    new GLUI_RadioButton( createconfig_radio, "BCC" );
    new GLUI_RadioButton( createconfig_radio, "Diamond" );
    new GLUI_RadioButton( createconfig_radio, "Nanotube" );
    GLUI_Panel *createconfig_panel011 = new GLUI_Panel( createconfig_glui, "Miller index");
    new GLUI_EditText(createconfig_panel011, "1 ", &milX1);
    new GLUI_EditText(createconfig_panel011, "2 ", &milY1);
    new GLUI_EditText(createconfig_panel011, "3 ", &milZ1);
    new GLUI_Column(createconfig_panel011, false);
    new GLUI_EditText(createconfig_panel011, "", &milX2);
    new GLUI_EditText(createconfig_panel011, "", &milY2);
    new GLUI_EditText(createconfig_panel011, "", &milZ2);
    new GLUI_Column(createconfig_panel011, false);
    new GLUI_EditText(createconfig_panel011, "", &milX3);
    new GLUI_EditText(createconfig_panel011, "", &milY3);
    new GLUI_EditText(createconfig_panel011, "", &milZ3);
    GLUI_Panel *createconfig_panel012 = new GLUI_Panel( createconfig_glui, "Rotation");
    new GLUI_EditText(createconfig_panel012, "X ", &cellrot_x);
    new GLUI_Column(createconfig_panel012, false);
    new GLUI_EditText(createconfig_panel012, "Y ", &cellrot_y);
    new GLUI_Button(createconfig_panel012, "Rotate", CB_ROTATE, control_cb);
    new GLUI_Column(createconfig_panel012, false);
    new GLUI_EditText(createconfig_panel012, "Z ", &cellrot_z);
    new GLUI_Button(createconfig_panel012, "Reverse", CB_ROTATE_REV, control_cb);
    GLUI_Rollout *createconfig_panel013 = new GLUI_Rollout(createconfig_glui, "Cell dimension", true);
    new GLUI_EditText(createconfig_panel013, "1 ", &cell1x, CB_CELL_DIM, control_cb);
    new GLUI_EditText(createconfig_panel013, "2 ", &cell2x, CB_CELL_DIM, control_cb);
    new GLUI_EditText(createconfig_panel013, "3 ", &cell3x, CB_CELL_DIM, control_cb);
    new GLUI_Checkbox(createconfig_panel013, "Fix atoms", &ifix_atoms);
    new GLUI_Column(createconfig_panel013, false);
    new GLUI_EditText(createconfig_panel013, "", &cell1y, CB_CELL_DIM, control_cb);
    new GLUI_EditText(createconfig_panel013, "", &cell2y, CB_CELL_DIM, control_cb);
    new GLUI_EditText(createconfig_panel013, "", &cell3y, CB_CELL_DIM, control_cb);
    new GLUI_Column(createconfig_panel013, false);
    new GLUI_EditText(createconfig_panel013, "", &cell1z, CB_CELL_DIM, control_cb);
    new GLUI_EditText(createconfig_panel013, "", &cell2z, CB_CELL_DIM, control_cb);
    new GLUI_EditText(createconfig_panel013, "", &cell3z, CB_CELL_DIM, control_cb);
    GLUI_Rollout *createconfig_panel014 = new GLUI_Rollout(createconfig_glui, "Slicer", true);
    new GLUI_EditText(createconfig_panel014, "1 min", &slice_1min);
    new GLUI_EditText(createconfig_panel014, "1 max", &slice_1max);
    new GLUI_Column(createconfig_panel014, false);
    new GLUI_EditText(createconfig_panel014, "2 min", &slice_2min);
    new GLUI_EditText(createconfig_panel014, "2 max", &slice_2max);
    new GLUI_Column(createconfig_panel014, false);
    new GLUI_EditText(createconfig_panel014, "3 min", &slice_3min);
    new GLUI_EditText(createconfig_panel014, "3 max", &slice_3max);
    new GLUI_Button(createconfig_panel014, "Slice", CB_SLICE, control_cb);

    GLUI_Spinner *createconfig_spinner10
      = new GLUI_Spinner(createconfig_panel01, "NT Rotation (z) [deg]", &rotz, 0, control_cb );
    createconfig_spinner10->set_float_limits( 0.0, 360.0 );
    GLUI_Spinner *createconfig_spinner11
      = new GLUI_Spinner(createconfig_panel01, "NT Shift (z) [%]", &shiftz, 0, control_cb );
    createconfig_spinner11->set_float_limits( 0.0, 100.0 );
    new GLUI_Column(createconfig_panel01, false );
    GLUI_Spinner *createconfig_spinner1
      = new GLUI_Spinner(createconfig_panel01, "# of rep in x", &irepx, 0, control_cb );
    createconfig_spinner1->set_int_limits( 1, 20 );
    GLUI_Spinner *createconfig_spinner2
      = new GLUI_Spinner(createconfig_panel01, "# of rep in y", &irepy, 0, control_cb );
    createconfig_spinner2->set_int_limits( 1, 20 );
    GLUI_Spinner *createconfig_spinner3
      = new GLUI_Spinner(createconfig_panel01, "# of rep in z", &irepz, 0, control_cb );
    createconfig_spinner3->set_int_limits( 1, 20 );
    GLUI_Spinner *createconfig_spinner4
      = new GLUI_Spinner(createconfig_panel01, "m of (m,n)", &icntm, 0, control_cb );
    createconfig_spinner4->set_int_limits( 0, 200 );
    GLUI_Spinner *createconfig_spinner5
      = new GLUI_Spinner(createconfig_panel01, "n of (m, n)", &icntn, 0, control_cb );
    createconfig_spinner5->set_int_limits( 0, 200 );
    GLUI_Spinner *createconfig_spinner7
      = new GLUI_Spinner(createconfig_panel01, "NT cell size (x/y) [ang]", &cscnt, 0, control_cb );
    createconfig_spinner7->set_float_limits( 10.0f, 200.0 );
    GLUI_Spinner *createconfig_spinner6
      = new GLUI_Spinner(createconfig_panel01, "Lattice const. [ang]", &alat, 0, control_cb );
    createconfig_spinner6->set_float_limits( 0.0f, 100.0 );

    GLUI_Rollout *createconfig_panel02 = new GLUI_Rollout(createconfig_glui, "CNT compression (wall mode)", false );
    new GLUI_Button(createconfig_panel02, "Set CNT wall", CNTWALL_DO_ID, pointer_cb);
    new GLUI_Column(createconfig_panel02, false );
    cntwall_btn = new GLUI_Button(createconfig_panel02, "Read wall data", CNTWALL_READ_ID, pointer_cb);
    new GLUI_Column(createconfig_panel02, false );
    new GLUI_Button(createconfig_panel02, "Write wall data", CNTWALL_WRITE_ID, pointer_cb);
    //new GLUI_Checkbox(createconfig_panel02, "Show CNT wall", &show_cnt_wall);
    GLUI_Rollout *createconfig_panel03 = new GLUI_Rollout(createconfig_glui, "Edit config", false);
    //GLUI_Panel *createconfig_panel03 = new GLUI_Panel(createconfig_glui, "Edit config" );
    new GLUI_Spinner(createconfig_panel03, "atom #", &iatom_pick, CB_EDITCONFIG, control_cb );
    new GLUI_EditText(createconfig_panel03,"Atom", config_atom);
    new GLUI_EditText(createconfig_panel03,"x (A)", &atomrx);
    new GLUI_EditText(createconfig_panel03,"y (A)", &atomry);
    new GLUI_EditText(createconfig_panel03,"z (A)", &atomrz);
    new GLUI_Column(createconfig_panel03, false );
    new GLUI_Button(createconfig_panel03, "Remove atom", REMOVEATOM_DO_ID, pointer_cb);
    new GLUI_Button(createconfig_panel03, "Add atom", ADDATOM_DO_ID, pointer_cb);
    new GLUI_Button(createconfig_panel03, "Move/Change", MOVEATOM_DO_ID, pointer_cb);

    GLUI_Panel *createconfig_panel10 = new GLUI_Panel(createconfig_glui, "", GLUI_PANEL_NONE );
    new GLUI_Button(createconfig_panel10, "Create", CREATECONFIG_DO_ID, pointer_cb);
    new GLUI_Column(createconfig_panel10, false );
    new GLUI_Button(createconfig_panel10, "Multiply", MULTIPLYCELL_DO_ID, pointer_cb);
    new GLUI_Column(createconfig_panel10, false );
    new GLUI_Button(createconfig_panel10, "Close", CREATECONFIG_CLOSE_ID, pointer_cb);
    createconfig_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == CREATECONFIG_CLOSE_ID ) {
    createconfig_mode = 0;
    set_atom_color();
    createconfig_btn->enable();
    control->glui->close();
  }
  else if ( control->get_id() == REMOVEATOM_DO_ID ) {
    remove_atom(iatom_pick);
    glutPostRedisplay();
  }
  else if ( control->get_id() == ADDATOM_DO_ID ) {
    add_atom(config_atom, atomrx, atomry, atomrz);
    iatom_pick = atom.natom;
    atomrx=atom.rx[iatom_pick]/ang; atomry=atom.ry[iatom_pick]/ang; atomrz=atom.rz[iatom_pick]/ang;
    glutPostRedisplay();
  }
  else if ( control->get_id() == MOVEATOM_DO_ID ) {
    //move_atom(iatom_pick, atomrx, atomry, atomrz);
    move_atom(iatom_pick, config_atom, atomrx, atomry, atomrz);
    glutPostRedisplay();
  }
  else if ( control->get_id() == CREATECONFIG_DO_ID ) {
    createconfig();
    glutPostRedisplay();
  }
  else if ( control->get_id() == MULTIPLYCELL_DO_ID ) {
    multiply_cell(irepx,irepy,irepz);
    glutPostRedisplay();
  }
  else if ( control->get_id() == CNTWALL_DO_ID ) {
    cnt_wall_set();
    glutPostRedisplay();
  }
  else if ( control->get_id() == CNTWALL_READ_ID ) {
    filename_cntwall_glui = GLUI_Master.create_glui( "CNTWALL filename", 0, 600, 150 );
    filename_cntwall = new GLUI_CommandLine( filename_cntwall_glui, "File:", NULL, -1, pointer_cb );
    filename_cntwall->set_w( 300 );
    cntwall_fb = new GLUI_FileBrowser(filename_cntwall_glui, "", false, CB_CNTWALL_FB, pointer_cb);
    cntwall_fb->set_w(300);
    GLUI_Panel *panel_cntwall = new GLUI_Panel(filename_cntwall_glui, "", GLUI_PANEL_NONE);
    new GLUI_Button(panel_cntwall, "Close", CNTWALL_CLOSE_ID, pointer_cb);
    filename_cntwall_glui->set_main_gfx_window( main_window );
    control->disable();
  }
  else if ( control->get_id() == CNTWALL_CLOSE_ID ) {
    cntwall_btn->enable();
    control->glui->close();
  }

  else if ( control->get_id() == CNTWALL_WRITE_ID ) {
    cnt_wall_write("CNTWALL.SAVE");
    glutPostRedisplay();
  }
#ifndef NOPNG
  else if ( control->get_id() == CAPTURE_ID) {
    capture();
    glutPostRedisplay();
  }
#endif
}

/**************************************** myGlutKeyboard() **********/
void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    // A few keys here to test the sync_live capability.
  case 'e':
    // Cycle through ensemble types
    ++ensemble %= 2;
    GLUI_Master.sync_live_all();
    break;
  case 'w':
    // Toggle wireframe mode
    wireframe = !wireframe;
    GLUI_Master.sync_live_all();
    break;
  case 'r':
    org_x=0.0; org_y=0.0; org_z=0.0;
    rotationX=0.0; rotationY=0.0;
    for (int i=0; i<=15; i++) { rotate[i]=0; }
    rotate[0]=1;rotate[6]=-1;rotate[9]=1;rotate[15]=1;
    /*
    for (int i=1; i<=atom.natom; i++) {
      atom.rx[i]=atom.rx_org[i]; atom.ry[i]=atom.ry_org[i]; atom.rz[i]=atom.rz_org[i];
      atom.vx[i]=0.0;atom.vy[i]=0.0;atom.vz[i]=0.0;
      atom.fx[i]=0.0;atom.fy[i]=0.0;atom.fz[i]=0.0;  }
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	cell.hmat[i][j]=cell.hmat_org[i][j];  } }
    istep = 0;
    */
    glui->sync_live();

    writedata_initialize();

    break;
  case 27: 
  case 'q':
    exit(0);
    break;
  };
  glutPostRedisplay();
}

/*
void potential_set (int control)
{
  printf("Potential: No. %d, %s\n",ipottype,potstring_list[ipottype]);
  //strcpy(atom.potential_func, potstring_list[ipottype]);
  std::string line = std::string(potstring_list[ipottype]);
  std::string arg1, arg2; int narg = count_arg_number(line);
  if (narg==1) { get_first_arg(line,arg1);
    strcpy(atom.potential_func, arg1.c_str());
    strcpy(atom.potential_arg, "");
  } else { get_first_arg(line,arg1); get_first_arg(line,arg2);
    strcpy(atom.potential_func, arg1.c_str());
    strcpy(atom.potential_arg,  arg2.c_str()); }
  book.algo = 1;
  if (strcmp(atom.potential_func,"Tersoff")==0) { book.algo = 2; }
  if (strcmp(atom.potential_func,"Brenner")==0) { book.algo = 3; }
}
*/

/***************************************** myGlutMenu() ***********/

void myGlutMenu( int value )
{
  myGlutKeyboard( value, 0, 0);
}

void idle( void )
{
  if (mdmotion == 1) {
    for (int i=1; i<=mdspeed; i++) {
      if (istep % book.nbk ==0) {
	if (incell) { out_of_cell(); }
	bookkeep();
      }
      md(); if (mdmotion == 0) { break; }
      writedata();
      istep++;
    }
  }
  glutPostRedisplay();

  glui->sync_live();
}

void myGlutIdle( void )
{
  if ( glutGetWindow() != main_window ) 
    glutSetWindow(main_window);
  idle();
  /*
  if (mdmotion == 1) {
    for (int i=1; i<=mdspeed; i++) { // MD: numerical integral of EoM
      if (istep % book.nbk ==0) {
	if (incell) { out_of_cell(); }
	bookkeep();
      }
      md(); if (mdmotion == 0) { break; }
      writedata();
      istep++;
    }
  }
  glutPostRedisplay();
  glui->sync_live();
  */
}

/***************************************** myGlutMouse() **********/

void myGlutMouse(int button, int button_state, int x, int y )
{
  static GLuint selection[SELECTIONS];
  static GLint hits = 0;

  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (button_state == GLUT_DOWN) {
      GLuint *ptr; GLint vp[4];
      glSelectBuffer(SELECTIONS, selection);
      glRenderMode(GL_SELECT);
      glInitNames();
      glPushName(-1);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      glGetIntegerv(GL_VIEWPORT, vp);
      gluPickMatrix(x, vp[3] - y - 1, 1, 1, vp);
      if (ortho) {
	glOrtho(-size_w/200.0, size_w/200.0, -size_h/200.0, size_h/200.0, -100.0, 100.0);
      } else {
	gluPerspective(30.0, (double)vp[2] / (double)vp[3], 1.0, 100.0);
	glTranslated(0.0, 0.0, -1.0);
      }
      glMatrixMode(GL_MODELVIEW);
      //      for (int i = 0; i < atom.natom; i++) {
      for (int i = 1; i <= atom.natom; i++) {
	glLoadName(i);
	glCallList(objects + i);
      }
      if (draw_replica>0) {
	for (int i = 1; i <= icnt; i++) {
	  glLoadName(ibase+i);
	  glCallList(ibase + i);
	}
      }
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      hits = glRenderMode(GL_RENDER);
      //printf("hits = %d\n",hits);
      ptr = selection;
      //      for (int i = 0; i < hits; i++) {
      //      printf("selection is %d\n",*ptr);
      if (hits>0) {
	unsigned int j, n = ptr[0];
	double near = (double)ptr[1] / (double)0x7fffffff;
	double far  = (double)ptr[2] / (double)0x7fffffff;
	ptr += 3;
	for (j = 0; j < n; j++) {
	  if (*ptr<=atom.natom) {
	    printf("Atom = %d [%s] Pos (A): %f %f %f\n",
		   *ptr,atom.asp[*ptr],atom.rx[*ptr]/ang,atom.ry[*ptr]/ang,atom.rz[*ptr]/ang);
	    printf("               F (eV/A): %f %f %f (total)\n",
		   atom.fx[*ptr]/ev*ang,atom.fy[*ptr]/ev*ang,atom.fz[*ptr]/ev*ang);
	    printf("                       : %f %f %f (from loading)\n",
		   atom.fx_l[*ptr]/ev*ang,atom.fy_l[*ptr]/ev*ang,atom.fz_l[*ptr]/ev*ang);
	    printf("               Ene (eV): %f\n",
		   atom.epot[*ptr]/ev);
	    printf("               V  (m/s): %f %f %f\n",
		   atom.vx[*ptr],atom.vy[*ptr],atom.vz[*ptr]);
	    printf("           Stress (MPa): xx = %f yy = %f zz = %f\n",
		   atom.satom[*ptr][0][0]/1e6,atom.satom[*ptr][1][1]/1e6,
		   atom.satom[*ptr][2][2]/1e6);
	    printf("                         xy = %f yz = %f zx = %f\n",
		   atom.satom[*ptr][0][1]/1e6,atom.satom[*ptr][1][2]/1e6,
		   atom.satom[*ptr][2][0]/1e6);
	    if (createconfig_mode) {
	      iatom_pick = *ptr;
	      strcpy(config_atom,atom.asp[iatom_pick]);
	      atomrx=atom.rx[iatom_pick]/ang;
	      atomry=atom.ry[iatom_pick]/ang;
	      atomrz=atom.rz[iatom_pick]/ang;
	    }
	  } else {
	    int ia=iatom[*ptr-atom.natom-objects];
	    printf("Replica Atom = %d  Pos: %f %f %f  Rep-index: %d\n",
		   ia,atom.rx[ia]*1e9,atom.ry[ia]*1e9,atom.rz[ia]*1e9,
		   repidx[*ptr-objects]);
	  }
	    
	  if (edit_elem_mode>0) {
	    int iimax=4; if (edit_elem_mode==2) { iimax=3; }
	    int ifnd = 0;
	    for (int ii=0; ii<iimax; ii++) {
	      int ix = *ptr; int iy = 0;
	      if (ix>atom.natom) { ix=iatom[*ptr-atom.natom-objects]; iy = repidx[*ptr-atom.natom-objects]; }
	      if ((select_atom[ii] == ix)&&(select_atom_repidx[ii] == iy)) {
		if (*ptr<=atom.natom) {
		  memcpy(color[*ptr],yellow,sizeof(GLfloat)*4);
		} else {
		  memcpy(color[*ptr-objects],yellow,sizeof(GLfloat)*4);
		}
		select_atom[ii] = 0; select_atom_repidx[ii] = 0;
		ifnd = 1;
		break;
	      }
	    }
	    if (ifnd==0) {
	      for (int ii=0; ii<iimax; ii++) {
		if (select_atom[ii] == 0) {
		  if (*ptr<=atom.natom) {
		    memcpy(color[*ptr],green,sizeof(GLfloat)*4);
		  } else {
		    memcpy(color[*ptr-objects],green,sizeof(GLfloat)*4);
		  }
		  if (*ptr<=atom.natom) {
		    select_atom[ii]=*ptr; select_atom_repidx[ii]=0;
		  } else {
		    select_atom[ii]=iatom[*ptr-atom.natom-objects];
		    select_atom_repidx[ii]=repidx[*ptr-atom.natom-objects];
		  }
		  printf("Selected Atom %d = %d \n", ii, select_atom[ii]);
		  break;
		}
	      }
	    }
	  } else if (measure_mode>0) {
	    int iimax=3;
	    int ifnd = 0;
	    for (int ii=0; ii<iimax; ii++) {
	      int ix = *ptr; int iy = 0;
	      if (ix>atom.natom) { ix=iatom[*ptr-atom.natom-objects]; iy = repidx[*ptr-atom.natom-objects]; }
	      if ((select_atom[ii] == ix)&&(select_atom_repidx[ii] == iy)) {
		if (*ptr<=atom.natom) { set_atom_color(*ptr);
		} else { set_atom_color(*ptr-objects); }

		select_atom[ii] = 0; select_atom_repidx[ii] = 0;
		ifnd = 1; break;
	      }
	    }
	    if (ifnd==0) {
	      for (int ii=0; ii<iimax; ii++) {
		if (select_atom[ii] == 0) {
		  if (*ptr<=atom.natom) {
		    memcpy(color[*ptr],green,sizeof(GLfloat)*4);
		  } else {
		    memcpy(color[*ptr-objects],green,sizeof(GLfloat)*4);
		  }
		  if (*ptr<=atom.natom) {
		    select_atom[ii]=*ptr; select_atom_repidx[ii]=0;
		  } else {
		    select_atom[ii]=iatom[*ptr-atom.natom-objects];
		    select_atom_repidx[ii]=repidx[*ptr-atom.natom-objects];
		  }
		  printf("Selected Atom %d = %d \n", ii, select_atom[ii]);
		  break;
		}
	      }
	    }

	  } else { // if not edit_elem_mode nor measure_mode
	    memcpy(color[*ptr],blue,sizeof(GLfloat)*4);
	  } // end if edit_elem_mode & measure_mode
	}
	GLUI_Master.sync_live_all();
      }
    }
    else {
      GLuint *ptr = selection;
      int i;
      //      for (i = 0; i < hits; i++) {
      if (hits>0) {
	unsigned int j, n = ptr[0];
	ptr += 3;
	for (j = 0; j < n; j++) {
	  if ((edit_elem_mode>0)||(measure_mode>0)) {

	  } else { // return to original color
	    //memcpy(color[*ptr],yellow,sizeof(GLfloat)*4);
	    memcpy(color[*ptr],color0[*ptr],sizeof(GLfloat)*4);
	  }
	}
      }
    }
    glutPostRedisplay();
    break;

  case GLUT_MIDDLE_BUTTON:
    glutIdleFunc(idle);
    mdmotion = 1;
    break;
  case GLUT_RIGHT_BUTTON:
    glutIdleFunc(0);
    mdmotion = 0;
    break;
    /*  
  case GLUT_MIDDLE_BUTTON:
    if (button_state == GLUT_DOWN) {
      glutIdleFunc(0);
      glutPostRedisplay();
      mdmotion = 1;
    } else {
      //      glutPostRedisplay();
      glutIdleFunc(0);
      mdmotion = 0;
    }      
    break;
    */
  }
}

/***************************************** myGlutMotion() **********/
void myGlutMotion(int x, int y )
{
  //glutPostRedisplay(); 
}

/**************************************** myGlutReshape() *************/
void myGlutReshape( int x, int y )
{
  xy_aspect = (float)x / (float)y;
  glViewport( 0, 0, x, y);
  glMatrixMode(GL_PROJECTION);// <==
  glLoadIdentity();
  if (ortho) {
    GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
    GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
    GLfloat light0_position[] = {200.0f, 200.0f, 200-6.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    glOrtho(-x/200.0, x/200.0, -y/200.0, y/200.0, -100.0, 100.0);
  } else {
    GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
    GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
    GLfloat light0_position[] = {20.0f, 20.0f, 200-4.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    gluPerspective(30.0, (double)xy_aspect, 1.0, 100.0);
    glTranslated(0.0, 0.0, -1.0);
  }
  glMatrixMode(GL_MODELVIEW);// <==
  glutPostRedisplay();
  size_w = x; size_h = y;
}


/**************************************** main() ********************/
int main(int argc, char* argv[])
{
  /*   Initialize GLUT and create window  */
  // This must be done before READCONFIG, otherwise it crashes in some systems.
  glutInit(&argc, argv);
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition( 50, 50 );
  glutInitWindowSize( 600, 600 );
  readsetdat("SETDAT000");
  main_window = glutCreateWindow( "MD viewer" );
  //

#if defined __linux__ || defined __APPLE__
  if (argc>1) {
    printf("Change directory to %s\n",argv[1]); chdir(argv[1]); }
  // default initial configuration
  strcpy(atom.potential_func, "GEAM"); book.algo = 1;
  if (argc<=2) {
    readsetdat("SETDAT"); readconfig("CONFIG"); //readqcelement("QCELEMENT");
  } else if (argc==3) {
    readsetdat("SETDAT"); readconfig(argv[2]); //readqcelement("QCELEMENT");
  } else if (argc==4) {
    readsetdat(argv[3]); readconfig(argv[2]); //readqcelement("QCELEMENT");
  }
#else
  readsetdat("SETDAT"); readconfig("CONFIG");
#endif
  if ((mode_cnt_corrugation)&&(read_cntwall)&&(cnt_load_algo==0)) {
    cnt_wall_read("CNTWALL");
  }
  getcwd(cwdname,80);
  printf("Current Working Directory = %s\n",cwdname);
  writedata_initialize(); md_set();

  for (int i=1; i<=atom.natom; i++) {
    atom.vx[i] = 0.0e3;   atom.vy[i] = 0.0e3  ; atom.vz[i] = 0.0e3;
  }
  for (int i=1; i<=atom.natom; i++) {
    atom.rx_org[i]=atom.rx[i];atom.ry_org[i]=atom.ry[i];atom.rz_org[i]=atom.rz[i];
  }

  // if SETDAT000 has "nogui yes"
  if (inogui>0) {
    while(mdmotion==1) { md(); istep++; }
    exit(0);
  }

  //  std::ofstream fouttraj("trajectory.d");
  //  std::ofstream fout("energy.d");

  // Then GLUT functions contiune... 
  glutDisplayFunc( myGlutDisplay );
  glutReshapeFunc( myGlutReshape );  
  glutKeyboardFunc( myGlutKeyboard );
  glutMotionFunc( myGlutMotion );
  glutMouseFunc( myGlutMouse );

  /*       Set up OpenGL lights           */

  //GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
  //GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
  //GLfloat light0_position[] = {200.0f, 200.0f, -4.0f, 1.0f};
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  //glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  //glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

  /*          Enable z-buferring          */
  glEnable(GL_DEPTH_TEST);

  //objects = glGenLists(NOBJECTS); printf("objects = %d\n",objects);

  /*         Here's the GLUI code         */
  //  GLUI *glui = GLUI_Master.create_glui( "Control", 0, 600, 50 ); // name, flags, x, and y
  glui = GLUI_Master.create_glui( "Control", 0, 600, 50 ); // name, flags, x, and y

  new GLUI_StaticText( glui, "MDSPASS ver.2.0" );
  new GLUI_Separator( glui );

  GLUI_Panel *panel000 = new GLUI_Panel( glui, "", GLUI_PANEL_NONE );
  checkbox_pbcx  = new GLUI_Checkbox( panel000, "PBC x", &cell.pbcx);
  new GLUI_Column( panel000, false );
  checkbox_pbcy  = new GLUI_Checkbox( panel000, "PBC y", &cell.pbcy);
  new GLUI_Column( panel000, false );
  checkbox_pbcz  = new GLUI_Checkbox( panel000, "PBC z", &cell.pbcz);
  new GLUI_Column( panel000, false );
  GLUI_Listbox *potlist = new GLUI_Listbox( panel000, "Potential:", &ipottype,0,potential_set);
  for (int i=0; i<MAXPOTTYPE; i++) { potlist->add_item(i, potstring_list[i]); }
  GLUI_Panel *panel001 = new GLUI_Panel( glui, "", GLUI_PANEL_NONE );
  //setdexdt_btn = new GLUI_Button( panel001, "Set dex/dt", SETDEXDT_ID, pointer_cb );
  //new GLUI_Column( panel001, false );
  setparam_btn = new GLUI_Button( panel001, "Set param", SETPARAM_ID, pointer_cb );
  calc_btn = new GLUI_Button(panel001, "Calc", CALC_ID, pointer_cb);
  new GLUI_Column( panel001, false );
  mdswitch_btn = new GLUI_Button( panel001, "MD on/off", MDSWITCH_ID, pointer_cb );
  strs_btn = new GLUI_Button(panel001, "Stress", STRS_ID, pointer_cb);
  new GLUI_Column( panel001, false );
  GLUI_EditText *spinner_temp_set = new GLUI_EditText( panel001, "Temp set:", &temp_set );
  spinner_temp_set->set_float_limits( 0.0f, 1000.0 );
  spinner_temp_set->set_alignment(GLUI_ALIGN_RIGHT);
  extra_btn = new GLUI_Button(panel001, "Extra", EXTRA_ID, pointer_cb);

  GLUI_Rollout *panel003 = new GLUI_Rollout(glui, "QC settings", false);
  GLUI_Panel *panel0031 = new GLUI_Panel(panel003, "", GLUI_PANEL_NONE );
  checkbox_QC  = new GLUI_Checkbox( panel0031, "QC", &atom.QC);
  new GLUI_Column(panel0031, false );
  checkbox_show_only_elem  = new GLUI_Checkbox( panel0031, "Show only QC element", &show_only_elem);
  GLUI_Panel *panel0032 = new GLUI_Panel(panel003, "", GLUI_PANEL_NONE );
  edit_elem_btn = new GLUI_Button( panel0032, "Edit elem", EDIT_ELEM_ID, pointer_cb );
  new GLUI_Column( panel0032, false );
  edit_elem_xz_btn = new GLUI_Button( panel0032, "Edit elem xz", EDIT_ELEM_XZ_ID, pointer_cb );
  GLUI_Panel *panel0033 = new GLUI_Panel(panel003, "", GLUI_PANEL_NONE );
  open_file_qcelement_btn = new GLUI_Button(panel0033, "Read qcelement", OPEN_FILE_QCELEMENT_ID, pointer_cb );
  new GLUI_Column( panel0033, false );
  writeqcelement_btn = new GLUI_Button( panel0033, "Write qcelement", WRITEQCELEMENT_ID, pointer_cb );

  //  edittext = new GLUI_EditText( glui, "Text:", text, 3, control_cb );
  GLUI_Panel *panel005 = new GLUI_Panel( glui, "", GLUI_PANEL_NONE );
  GLUI_Listbox *enslist
    = new GLUI_Listbox( panel005, "Algorithm:", &ensemble,CB_ENSEMBLE, control_cb);
  for (int i=0; i<MAXENSTYPE; i++) { enslist->add_item(i, ensemble_list[i]); }
  /*
  GLUI_Panel *ensemble_panel = new GLUI_Panel( panel005, "Ensemble Type" );
  radio = new GLUI_RadioGroup( ensemble_panel,&ensemble,CB_ENSEMBLE,control_cb );
  new GLUI_RadioButton( radio, "NVE" );
  new GLUI_RadioButton( radio, "NVT" );
  new GLUI_RadioButton( radio, "Relaxation" );
  new GLUI_RadioButton( radio, "NPH" );
  new GLUI_RadioButton( radio, "NPV" );
  new GLUI_RadioButton( radio, "NPH + PRdamper" );
  new GLUI_RadioButton( radio, "NPV + PRdamper" );
  new GLUI_RadioButton( radio, "Full relaxation" );
  */
  //  new GLUI_RadioButton( radio, "diamond" );

  new GLUI_Column( panel005, false );
  //  new GLUI_Column( panel005, false );
  checkbox_Trans  = new GLUI_Checkbox( panel005, "No translation", &notrans);
  new GLUI_Checkbox(panel005, "Keep in cell", &incell);

  GLUI_Rollout *panel00 = new GLUI_Rollout(glui, "Instability analysis", false);
  //  GLUI_Panel *panel00 = new GLUI_Panel( glui, "", GLUI_PANEL_NONE );
  checkbox_evec = new GLUI_Checkbox( panel00, "Evector", &ievec);
  //  new GLUI_Column( panel00, false);
  spinner_instcenter_num  = new GLUI_Spinner( panel00, "Center atom (inst)", &atom.instcenter, 2, control_cb);
  //spinner_instcenter_num->set_int_limits( 1, NMAX );
  spinner_instcenter_num->set_speed( 0.01 );
  GLUI_RadioGroup *radio_inst =  new GLUI_RadioGroup(panel00, &inst_mode, 0, control_cb);
  //new GLUI_RadioButton(radio_inst, "instability (DGEEV)");
  //new GLUI_RadioButton(radio_inst, "instability_dsy");
  new GLUI_RadioButton(radio_inst, "instability_atomcell");
  new GLUI_RadioButton(radio_inst, "instability_atom");
  new GLUI_RadioButton(radio_inst, "instability_atom_noremove");
  new GLUI_Checkbox(panel00, "Read Hessian", &hessian_read);
  new GLUI_Checkbox(panel00, "Write Hessian", &hessian_write);
  new GLUI_Column( panel00, false );
  spinner_evec_num = new GLUI_Spinner( panel00, "Mode", &ievec_num, CB_INST, control_cb );
  //spinner_evec_num->set_int_limits( 1, 20 );  //  spinner_evec_num->set_speed( 0.01 );
  //  new GLUI_Column( panel00, false);
  spinner_evec_len = new GLUI_Spinner( panel00, "Length", &evec_len, CB_INST, control_cb );
  new GLUI_Spinner(panel00, "Eig val", &eigval, CB_INST, control_cb);
  inst_btn = new GLUI_Button( panel00, "Inst analysis", INST_ID, pointer_cb );

  GLUI_Panel *panel1 = new GLUI_Panel( glui, "", GLUI_PANEL_NONE );
  GLUI_Rotation *view_rot = new GLUI_Rotation( panel1, "Rotation", rotate);
  //  GLUI_Rotation *view_rot = glui->add_rotation( "Rotation", rotate);
  float array[16] = {1.0, 0.0, 0.0, 0.0, 
		     0.0, 0.0,-1.0, 0.0, 
		     0.0, 1.0, 0.0, 0.0, 
		     0.0, 0.0, 0.0, 1.0};
  view_rot->set_float_array_val(array);
  new GLUI_Column( panel1, false );
  GLUI_Translation *trans_xy = new GLUI_Translation( panel1, "Objects XY", GLUI_TRANSLATION_XY, obj_pos );
  trans_xy->set_speed(0.005);
  new GLUI_Column( panel1, false );
  GLUI_Translation *trans_z =  new GLUI_Translation( panel1, "Objects Z", GLUI_TRANSLATION_Z, &scl );
  trans_z->set_speed(0.005);
  new GLUI_Column( panel1, false );
  GLUI_Button *capture_cutton = new GLUI_Button(panel1, "Capture", CAPTURE_ID, pointer_cb);
  GLUI_EditText *capture_counter = new GLUI_EditText(panel1, "Cap-file#", &capture_count);
  new GLUI_Button(panel1, "Reset view", CB_RESET_VIEW, control_cb);

  GLUI_Rollout *panel11 = new GLUI_Rollout(glui, "Status", true);
  //GLUI_Panel *panel110 = new GLUI_Panel(panel11, "", GLUI_PANEL_NONE);
  //panel110->set_alignment(GLUI_ALIGN_RIGHT);
  //GLUI_EditText *counter_edittext = new GLUI_EditText(panel110, "Step", &istep);
  //counter_edittext->disable();
  GLUI_Panel *panel111 = new GLUI_Panel(panel11, "", GLUI_PANEL_NONE);
  panel111->set_alignment(GLUI_ALIGN_RIGHT);
  GLUI_EditText *step_monitor = new GLUI_EditText(panel111, "Step", &istep);
  step_monitor->set_alignment(GLUI_ALIGN_RIGHT);
  GLUI_EditText *epot_monitor = new GLUI_EditText(panel111, "E_p(eV/atm)", &epotatom);
  epot_monitor->set_alignment(GLUI_ALIGN_RIGHT);
  epot_monitor->set_w(150);
  //GLUI_EditText *cnt_text1 = new GLUI_EditText(panel111, "CNTload(nN)", &cnt_pressure_ftot);
  //cnt_text1->set_alignment(GLUI_ALIGN_RIGHT);
  //cnt_text1->set_w(160);
  GLUI_EditText *temp_monitor = new GLUI_EditText(panel111, "Temp(K)", &tempc);
  temp_monitor->set_alignment(GLUI_ALIGN_RIGHT);
  GLUI_EditText *fmax_monitor = new GLUI_EditText(panel111, "Fmax(ev/A)", &f_max);
  fmax_monitor->set_alignment(GLUI_ALIGN_RIGHT);
  fmax_monitor->set_w(160);
  new GLUI_Column( panel11, false );
  GLUI_Panel *panel112 = new GLUI_Panel(panel11, "", GLUI_PANEL_NONE);
  panel112->set_alignment(GLUI_ALIGN_RIGHT);
  status_lx = new GLUI_EditText(panel112, "Cell (A)", &cellx,0,pointer_cb); 
  status_ly = new GLUI_EditText(panel112, "        ", &celly,0,pointer_cb);
  status_lz = new GLUI_EditText(panel112, "        ", &cellz,0,pointer_cb);
  status_dt = new GLUI_EditText(panel112, "dt (fs)", &dtm,0,pointer_cb);
  

  GLUI_Rollout *panel2 = new GLUI_Rollout( glui, "File control and config creation", false );
  GLUI_Panel *panel21 = new GLUI_Panel(panel2, "", GLUI_PANEL_NONE);
  open_file_btn = new GLUI_Button( panel21, "Read Config", OPEN_FILE_ID, pointer_cb );
  new GLUI_Column( panel21, false );
  writeconfig_btn = new GLUI_Button( panel21, "Write Config", WRITECONFIG_ID, pointer_cb );
  new GLUI_Column( panel21, false );
  createconfig_btn = new GLUI_Button( panel21, "Create Config", CREATECONFIG_ID, pointer_cb );

  GLUI_Panel *panel3 = new GLUI_Panel( glui, "", GLUI_PANEL_NONE );
  reset_btn = new GLUI_Button( panel3, "Reset", RESET_ID, pointer_cb );
  new GLUI_Column( panel3, false );
  new GLUI_Button( panel3, "Quit", 0, (GLUI_Update_CB)exit );
  glui->set_main_gfx_window( main_window );

  //GLUI_Master.set_glutIdleFunc( NULL );
  GLUI_Master.set_glutIdleFunc( myGlutIdle );

  //  writedata_initialize();

  // Loop

  glutMainLoop();



  //  fouttraj.close();
  //  fout.close();

  return EXIT_SUCCESS;
}

