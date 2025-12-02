FC = ifort -static
OBJS = main.o bookkeep.o bookkeep_rp.o initconfsc.o \
       scalevelo.o enerkin.o initial_velocity.o \
       inverse.o calcvolume.o mkq.o mkenk.o \
       mkhmat.o mkhvmat.o out.o initmat.o \
       morseparam.o \
       search_neighbor.o kinetic.o var_initialize.o \
       forcecheck.o stresscheck.o writedata.o velocity_verlet.o gear_pc.o \
       gear_pc_nvt.o gear_pc_nph.o gear_pc_npt.o gear_pc_nptpr.o \
       calcstresstensor.o recipro.o outproduct.o gear_pc_nphpr.o strain.o \
       calcelastic.o cellrlx.o velocity_scaling.o atom_damping.o cellconstraint.o \
       emt_St_JCP1990.o emt_HM_JP1989.o forcecheck2.o potential.o calc_replica.o \
       spline.o eam_Mishin_PRB65_224114.o eam_Mishin.o emt_JSN_SS1996.o sw.o tersoff.o tersoff_terpotn.o \
       readcont.o readcondition.o readconfig.o atom_weight.o eam_FS.o tersoff_multi.o \
       instability.o eispc.o atom_number.o voronoi.o \
       totaldrift.o set_atom_species.o eam_geam_AM.o eam_geam_PRB.o tersoff_nomat.o \
       tersoff_multi_nomat.o
PROGRAM = ymd

CFLAGS = -Vaxlib

.f.o:
	$(FC) -c $< -o $@
.F.o:
	$(FC) -FR -c $< -o $@

# Default
$(PROGRAM) : $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(FC) $(CFLAGS) $(OBJS) -o $(PROGRAM)
configmake : configmake.o
	@echo "Linking $(PROGRAM) ..."
	$(FC) $(CFLAGS) configmake.o -o configmake

clean :; 
	rm -rf $(OBJS)
	rm -rf $(PROGRAM)
	rm -rf configmake configmake.o
	rm -rf *.mod
	rm -rf *~
