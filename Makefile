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
# L_FLAGS=-g
# Linker flags
# Default flags, e.g., for gfortran, g77
#L_FLAGS=-c -O2

# Flags for Intel's ifort
# L_FLAGS=-O3 -axW -g

# FLAGS for Cray
# L_FLAGS=-fastsse -Mipa

# Flags for xlf
# L_FLAGS=-O3 -qhot -q64 -qipa -qextname=flush

# Optimized flag for gfortran on an i686 platform
# L_FLAGS=-fPIC -O3 -funroll-loops -mmmx -msse2 -msse -march=i686 -malign-double -fomit-frame-pointer

# Compiler flags
# Default flags, e.g., for gfortran, g77
#F_FLAGS=-c -O2 -g
F_FLAGS=-c -static -fopenmp
# Flags for Intel's ifort
# F_FLAGS=-c -O3 -axW -g
# Flags for debugging build
# F_FLAGS=-c -O0 -g
# FLAGS for Cray
# F_FLAGS=-c -fastsse -Mipa
# Flags for xlf
#F_FLAGS=-c -O3 -qhot -q64 -qipa -qextname=flush

# Optimized flag for gfortran on an i686 platform
# F_FLAGS=-c -fPIC -O3 -funroll-loops -mmmx -msse2 -msse -march=i686 -malign-double -fomit-frame-pointer

# Select your compiler
F90=gfortran
#F90=ifort
#F90=pgf90
#F90=xlf
MPIF90=mpif90


# ______________________________________________ Name of program
PROG=smmp
# ______________________________________________ Objects 
OBJ=redseq.o  bldmol.o getmol.o redvar.o setvar.o \
difang.o setmvs.o mklist.o redstr.o  dihedr.o enyflx.o addend.o opeflx.o opeshe.o minim.o minqsn.o\
contacts.o hbond.o helix.o sa.o metropolis.o rgyr.o zimmer.o \
canon.o mulcan_par_mod.o  outvar.o outpdb.o  \
pdbread.o  rmsdfun.o enyreg.o opereg.o mincjg.o cnteny.o \
init_energy.o init_molecule.o gradient.o energy.o \
regul.o nursvr.o twister.o eninteract.o bgs.o eyabgn.o enylun.o utilities.o\
partem_s.o esolan.o opesol.o
#_________Serial implementation of ECEPP/3 and solvent
SOBJ=enyshe.o enysol.o
#_________Parallel implementation of ECEPP/3 and solvent
POBJ=enyshe_p.o enysol_p.o partem_p.o

FILES=metropolis.f energy.f enyflx.f enyreg.f enyshe.f enysol.f redseq.f bldmol.f \
getmol.f redvar.f setvar.f difang.f setmvs.f mklist.f redstr.f dihedr.f addend.f \
opeflx.f opeshe.f minim.f minqsn.f contacts.f hbond.f helix.f sa.f90  rgyr.f \
zimmer.f canon.f outvar.f outpdb.f \
pdbread.f  rmsdfun.f opereg.f mincjg.f cnteny.f  init_energy.f init_molecule.f \
gradient.f regul.f nursvr.f twister.f eninteract.f bgs.f eyabgn.f enylun.f \
mulcan_par_mod.f90 esolan.f opesol.f
# partem_s.f


pyInterface=init_molecule init_energy mulcan_par mulcan_sim \
init_lund sgrnd energy helix zimmer rgyr rmsdfun \
sa canon  outpdb  minim  regul contacts interhbond hbond
# metropolis
# ============================================== Linking

$(PROG):	$(OBJ) $(SOBJ) main.o
	$(F90) -o $(PROG) $(L_FLAGS) main.o $(OBJ) $(SOBJ)

# Build parallel version of SMMP.
parallel: $(OBJ) $(POBJ) main_p.o 
	$(MPIF90) -o $(PROG) $(L_FLAGS) main_p.o $(OBJ) $(POBJ)

# Cross compile for BlueGene/P 
bgl: BGL_L_FLAGS = -L$(BGLSYS)/lib -lmpich.rts -lfmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts -qextname=flush
bgl: F_FLAGS = -c -O5 -qhot -g -I$(BGLSYS)/include -L$(BGLSYS)/lib -qarch=440 -qtune=440 -qextname=flush
bgl: BGL_F90 = blrts_xlf
bgl: F90 = $(BGL_F90)
bgl: MPIF90 = $(BGL_F90)
bgl: CC = blrts_xlc
bgl: LIBSF_MPI=-lmpich.rts -lfmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts

bgl: $(OBJ) $(POBJ) main_bgl_p.o
	$(BGL_F90) -o $(PROG).rts $(BGL_L_FLAGS) main_bgl_p.o $(OBJ) $(POBJ) $(LIBSF_MPI)

# Build the python bindings
pybind: $(FILES)
	#./rmexclpoint.py $(FILES)
	f2py -c -m smmp smmp.pyf $(FILES) --f90exec=$(F90) --f77exec=$(F90)
	#./restoreexclpoint.py $(FILES)

newpybind: $(FILES)
	./rmexclpoint.py $(FILES)
	mv smmp.pyf smmp.pyf.bak
	f2py -h smmp.pyf -m smmp only: ${pyInterface} : $(FILES)
	./restoreexclpoint.py $(FILES)

#doc:
#	$(MAKE) -w -C doc/
#	doxygen pySMMP.doxygen
	
examples: $(PROG)
	$(MAKE) -w -C EXAMPLES/
	
# ______________________________________________ Dependancies
main.o: multicanonical.mod

multicanonical.mod: mulcan_par_mod.f90
	$(F90) $(F_FLAGS) mulcan_par_mod.f90

# enyshe_p.o:
# 	$(MPIF90) $(F_FLAGS) enyshe_p.f
# enysol_p.o:
# 	$(MPIF90) $(F_FLAGS) enysol_p.f
# partem_p.o:
# 	$(MPIF90) $(F_FLAGS) partem_p.f
# main_p.o: 
# 	$(MPIF90) $(F_FLAGS) main_p.f

main.o redseq.o eyring.o bldmol.o getmol.o redvar.o setvar.o \
difang.o setmvs.o mklist.o redstr.o  dihedr.o enyflx.o\
enyshe.o addend.o opeflx.o opeshe.o minim.o minqsn.o enysol.o esolan.o\
contacts.o hbond.o helix.o sa.o metropolis.o rgyr.o zimmer.o\
pdbvars.o rmsdfun.o enyreg.o opereg.o cnteny.o opesol.o\
init_energy.o init_molecule.o gradient.o energy.o nursvr.o\
regul.o eninteract.o bgs.o eyabgn.o enylun.o\
canon.o mulcan_par.o mulcan_sim.o outvar.o outpdb.o partem_s.o : INCL.H 

pdbread.o  enyreg.o opereg.o init_molecule.o regul.o: INCP.H

.PHONY:	clean, restore, doc
clean:
	$(MAKE) -w -C EXAMPLES/ clean
	rm -f  smmp core* *.o *~ *.mod
restore:
	./restoreexclpoint.py $(FILES)

# end
