# Makefile: smmp 

#.SILENT:
%_p.o : %_p.f
	$(MPIF90) $(F_FLAGS) $<

.SUFFIXES:	.o .f
.f.o:
	$(F90) $(F_FLAGS) $<

.SUFFIXES:	.o .f90
.f90.o:
	$(F90) $(F_FLAGS) $<

# ==================================== Variables for compiling and linking
# L_FLAGS=-O0 -g
# Linker flags
# Default flags, e.g., for gfortran, g77
L_FLAGS=-O2 -g
# Flags for Intel's ifort
# L_FLAGS=-O3 -axW -g
# FLAGS for Cray
# L_FLAGS=-fastsse -Mipa
# Flags for xlf
# L_FLAGS=-O3 -qhot -q64 -qipa -qextname=flush
# L_FLAGS=-fPIC -O3 -funroll-loops -mmmx -msse2 -msse -march=i686 -malign-double -fomit-frame-pointer
# Compiler flags
# Default flags, e.g., for gfortran, g77
F_FLAGS=-c -O2 -g
# Flags for Intel's ifort
# F_FLAGS=-c -O3 -axW -g
# Flags for debugging build
# F_FLAGS=-c -O0 -g
# FLAGS for Cray
# F_FLAGS=-c -fastsse -Mipa
# Flags for xlf
# F_FLAGS=-c -O3 -qhot -q64 -qipa -qextname=flush
# F_FLAGS=-c -fPIC -O3 -funroll-loops -mmmx -msse2 -msse -march=i686 -malign-double -fomit-frame-pointer
# Select your compiler
F90=gfortran
#F90=ifort
#F90=pgf90
#F90=xlf
# F90=nagf90
MPIF90 = mpif90
# ______________________________________________ Name of program
BASEOBJ = ../redseq.o  ../bldmol.o ../getmol.o ../redvar.o ../setvar.o \
../difang.o ../setmvs.o ../mklist.o ../redstr.o  ../dihedr.o ../enyflx.o \
../addend.o ../opeflx.o ../opeshe.o \
../contacts.o ../hbond.o ../helix.o ../outvar.o ../outpdb.o  \
../pdbread.o  ../rmsdfun.o ../enyreg.o ../opereg.o ../cnteny.o \
../init_energy.o ../init_molecule.o ../gradient.o ../energy.o \
../regul.o ../nursvr.o ../twister.o ../eninteract.o ../eyabgn.o ../enylun.o\
../bgs.o ../minim.o ../minqsn.o ../mincjg.o ../opesol.o ../esolan.o
ENOBJ = ../enyshe.o ../enysol.o
PENOBJ = ../enyshe_p.o ../enysol_p.o
MCOBJ = ../canon.o ../metropolis.o ../rgyr.o

all: minimization annealing multicanonical regularization parallel_tempering_s # rootmeansquaredev

minimization:	$(BASEOBJ) $(ENOBJ) minimization.o 
	$(F90) -o minimization $(L_FLAGS) $(BASEOBJ) $(ENOBJ) minimization.o

annealing:	$(BASEOBJ) $(MCOBJS) $(ENOBJ) annealing.o
	$(F90) -o annealing $(L_FLAGS) $(BASEOBJ) $(MCOBJ) $(ENOBJ) ../anneal.o ../zimmer.o annealing.o

multicanonical:	$(BASEOBJ) $(MCOBJS) $(ENOBJ) multicanonical.o ../mulcan_par_mod.o 
	$(F90) -o multicanonical $(L_FLAGS) $(BASEOBJ) $(MCOBJ) $(ENOBJ) multicanonical.o ../mulcan_par_mod.o

regularization:	$(BASEOBJ) $(ENOBJ) regularization.o 
	$(F90) -o regularization $(L_FLAGS) $(BASEOBJ) $(ENOBJ) regularization.o

# rootmeansquaredev:	$(BASEOBJ) $(ENOBJ) rootmeansquaredev.o
#	$(F90) -o rootmeansquaredev $(L_FLAGS) $(BASEOBJ) $(ENOBJ) rootmeansquaredev.o

parallel_tempering: $(BASEOBJ) $(PENOBJ) $(MCOBJS) ../utilities.o partem_p.o parallel_tempering_p.o
	$(MPIF90) -o parallel_tempering $(L_FLAGS) $(BASEOBJ) $(PENOBJ) $(MCOBJ) ../utilities.o partem_p.o parallel_tempering_p.o

parallel_tempering_s: $(BASEOBJ) $(ENOBJ) $(MCOBJ) ../partem_s.o ../utilities.o parallel_tempering_s.o
	$(F90) -o parallel_tempering_s $(L_FLAGS) $(BASEOBJ) $(ENOBJ) $(MCOBJ) ../utilities.o ../partem_s.o parallel_tempering_s.o

.PHONY:	clean, restore
clean:
	rm -f  minimization annealing multicanonical regularization rootmeansquaredev parallel_tempering_s parallel_tempering core* *.o *~
# end
