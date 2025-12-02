# Choose fortran compiler.
# Normal compile
CC = g++ -O3 -ffast-math
#CC = g++ 
# Cuda compile
#CC = nvcc -arch=sm_13


MAIN = main.o
OBJS_MD = md.o writedata.o readconfig.o bookkeep.o vscale.o readqcelement.o \
	instability.o writeconfig.o stretch.o instability_QC.o e_force_geam.o \
	atom_number.o readsetdat.o writeqcelement.o createconfig.o \
	e_force_brenner.o cnt_wall.o myglutdisplay.o bond_set.o loading.o \
        e_force_tersoff.o out_of_cell.o velocity_random.o \
	set_atom_color.o relax.o capture.o e_force_mishin.o e_force_dipole.o \
	e_force_adp.o stress.o strain.o util.o change_atom_color.o \
	solve_eom.o qc_functions.o array_set.o multiply_cell.o forceconst.o \
	instability_atomcell.o instability_dsy.o instability_atom.o \
	instability_atom_noremove.o phonon_calc.o divide_cell.o neb_calc.o \
	remove_atom.o add_atom.o measure.o move_atom.o airebo.o recipe.o
OBJS_MD_NOCUDA = e_force_morse.o potential.o
#OBJS_MD_NOCUDA = e_force_morse_float.o potential.o
#OBJS_MD_CUDA = e_force_morse_cu.o potential_cu.o
PROGRAM = mdspass2

## Comment out when using g95.
#CFLAGS = -Vaxlib

#For lapack and blas
ifeq ($(PNG),no)
XXX = -DNOPNG
CFLAGS = -llapack  -lblas -lGL -lGLU -lglut -lglui -lm -L/usr/X11R6/lib -I/usr/X11R6/include -lX11 -lXmu
else
XXX =
CFLAGS = -llapack  -lblas -lGL -lGLU -lglut -lglui -lm -L/usr/X11R6/lib -I/usr/X11R6/include -lX11 -lXmu -lpng
endif
# For intel FC 10 + mkl10 (not tested)
#CFLAGS = -lmkl_lapack -lmkl_ia32 -lmkl_core -lguide -lpthread -lmkl_blacs -lGL -lGLU -lglut -lglui -lm -L/usr/X11R6/lib -I/usr/X11R6/include -lX11 -lXmu
#For intel FC 11 + mkl
#CFLAGS = -lmkl_lapack -lmkl_intel -lmkl_core -lmkl_sequential -lguide -lpthread -lmkl_blacs -lGL -lGLU -lglut -lglui -lm -L/usr/X11R6/lib -I/usr/X11R6/include -lX11 -lXmu

.SUFFIXES: .cpp .o .cu

.cpp.o:
	$(CC) $(XXX) -c $< -o $@
.cu.o:
	$(CC) -c $< -o $@

# Default
#default: cuda
mdspass2 : $(MAIN) $(OBJS_MD) $(OBJS_MD_NOCUDA)
	@echo "Linking $(PROGRAM) ..."
	$(CC) $(XXX) $(MAIN) $(OBJS_MD) $(OBJS_MD_NOCUDA) -o $(PROGRAM) $(CFLAGS)
cuda : $(MAIN) $(OBJS_MD) $(OBJS_MD_CUDA)
	@echo "Linking $(PROGRAM) ..."
	$(CC) $(XXX) $(MAIN) $(OBJS_MD) $(OBJS_MD_CUDA) -o $(PROGRAM) $(CFLAGS)

$(MAIN): myheader.h myclass.h
$(OBJS_MD): myheader.h myclass.h
$(OBJS_MD_NOCUDA): myheader.h myclass.h
$(OBJS_MD_CUDA): myheader.h myclass.h
#.o: myheader.h myclass.cpp
#.o: myheader.h

clean :; 
	rm -rf $(MAIN) $(OBJS_MD) $(OBJS_MD_CUDA) $(OBJS_MD_NOCUDA)
	rm -rf $(PROGRAM)
	rm -rf *~
