#      Makefile for moon shower program casc   (R.Engel)
#
#
# options:
#
# LINUX (Intel PC)
# OPT = -c -C  -Wunused -Wuninitialized -malign-double -mpentium -O
# OPT = -c -C  -Wunused -malign-double -m -g
# OPT = -c -C -Wunused -Wuninitialized -malign-double -fno-automatic \
#       -finit-local-zeros -mpentium -O
#OPT = -c -Wunused -Wuninitialized -fno-automatic -finit-local-zeros \
#-march=pentium -O -g -pg -fbounds-check -ffortran-bounds-check 

OPT = -c -Wunused -Wuninitialized -fno-automatic -O -g -pg -fbounds-check  

DEB = -pg
#OPT = -c -Wunused -Wuninitialized -fno-automatic -finit-local-zeros \
#-march=pentium -O -g  -fbounds-check -ffortran-bounds-check 
#DEB = 
#
# Libraries:
#

#LIB =  -L/home/cern/2005/lib -lkernlib_noshift -lmathlib -lpacklib_noshift
LIB =  
#
# fortran compiler, linker/loader:
#
FC = gfortran
LO = gfortran
#

.f.o:
	$(FC) $(OPT) $<

all: ZHS_time_freq_Fresnel_beam

clean:
	rm -f *.o *.chk *.prj ZHS_time_freq_Fresnel_beam core


ZHS_time_freq_Fresnel_beam:  ZHS_time_freq_Fresnel_beam.o ZHS_xsections_t.o greisen.o s_rndm.o trapfpe.o \
		tracklengths.o AMJ_medium_init_mod.o gauss.o dgausspkg.o CWJ_thinned_simult_hybrid.o divdif.o
	@echo " linking programs..."
	$(FC) -o ZHS_time_freq_Fresnel_beam ZHS_time_freq_Fresnel_beam.o ZHS_xsections_t.o greisen.o \
                 tracklengths.o s_rndm.o trapfpe.o AMJ_medium_init_mod.o gauss.o \
		 dgausspkg.o CWJ_thinned_simult_hybrid.o divdif.o $(DEB) $(LIB)
	@echo " run program with 'ZHS_time_freq_Fresnel_beam'"
	
ZHS_time_freq_Fresnel_beam.o: ZHS_time_freq_Fresnel_beam.f 
	$(FC) $(OPT) ZHS_time_freq_Fresnel_beam.f 

CWJ_thinned_simult_hybrid.o: CWJ_thinned_simult_hybrid.f
	$(FC) $(OPT) CWJ_thinned_simult_hybrid.f
	
ZHS_xsections_t.o: ZHS_xsections_t.f
	$(FC) $(OPT) ZHS_xsections_t.f

greisen.o: greisen.f
	$(FC) $(OPT) greisen.f

tracklengths.o: tracklengths.f
	$(FC) $(OPT) tracklengths.f

s_rndm.o: s_rndm.f
	$(FC) $(OPT) s_rndm.f

trapfpe.o: trapfpe.c
	gcc -c trapfpe.c 

AMJ_medium_init_mod.o: AMJ_medium_init_mod.f
	$(FC) $(OPT) AMJ_medium_init_mod.f

gauss.o: gauss.f
	$(FC) $(OPT) gauss.f
	
dgausspkg.o: dgausspkg.f
	$(FC) $(OPT) dgausspkg.f

divdif.o: divdif.f
	$(FC) $(OPT) divdif.f

