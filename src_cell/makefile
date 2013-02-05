###################################################
# MAKEFILE FOR MDFF ON LISA AMSTERDAM JUNE 2011
################################################### 

# The compiler
FC = mpif90 
# installation directory
INSTALL_DIR=../bin

# =================
#      on lisa 
# =================
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -fpp -check bounds -warn all
#FCFLAGS = -O3 -pg -warn all
#FCFLAGS = -O5 -ip -xW -warn all
#FCFLAGS = -fpp -O3 -ip -xW 

# =================
#   on superbock
# =================
#FCFLAGS = -cpp -O3 -ffree-line-length-0 -I/usr/include 

# =================
#   on mario 
# =================
FCFLAGS = -fpp -O3 -ip -xSSE4.2
#FCFLAGS = -fpp -O3 -pg
#FCFLAGS = -fpp -O3 -ip -xW
#FCFLAGS = -fpp -O0 -check all -warn all 

EXE= mdff_cell.x 

# libraries needed 

# on mario


# =================
#      on lisa 
# =================
#LDFLAGS = -lm -lfftw3 -lmkl_lapack -lmkl -lguide -lpthread 

# =================
#   on superbock
# =================
#LDFLAGS= -lm -llapack -lfftw3

# =================
#   on mario 
# =================
# mkl path
#MKL_PATH=/opt/intel/Compiler/11.1/064/mkl/lib/em64t
#LDFLAGS = $(MKL_PATH)/libmkl_intel_lp64.a -Wl,--start-group $(MKL_PATH)/libmkl_intel_lp64.a $(MKL_PATH)/libmkl_intel_thread.a $(MKL_PATH)/libmkl_core.a -Wl,--end-group -lguide
#INC=-I/opt/intel/Compiler/11.1/064/mkl/include/fftw
#LDFLAGS=-Wl,-rpath,/home/shared_apps/intel/composerxe-2011.0.084/mkl/lib/intel64/ -Wl,--start-group /home/shared_apps/intel/composerxe-2011.0.084/mkl/lib/intel64/libmkl_intel_lp64.a /home/shared_apps/intel/composerxe-2011.0.084/mkl/lib/intel64/libmkl_sequential.a /home/shared_apps/intel/composerxe-2011.0.084/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread  -limf -lm

# =================
#   on T5600
# =================

MKL_PATH=/opt/intel/composerxe/mkl/lib/intel64/
LDFLAGS =$(MKL_PATH)/libmkl_intel_lp64.a -Wl,--start-group $(MKL_PATH)/libmkl_intel_lp64.a $(MKL_PATH)/libmkl_intel_thread.a $(MKL_PATH)/libmkl_core.a -Wl,--end-group -lguide
INC=-I/opt/intel/composerxe/mkl/include/fftw/
LDFLAGS=-Wl,-rpath,/opt/intel/composerxe/mkl/lib/intel64/ -Wl,--start-group /opt/intel/composerxe/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_sequential.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread  -limf -lm
 


# List of executables to be built within the package
OBJ= constants.o io_file.o control.o prop.o time.o functions.o \
     fft_t.o rand.o lattice.o config.o md.o lbfgs.o m1qn3.o \
     block_average.o thermo.o kspace.o rspace.o kinetic.o field.o \
     tools.o read_pos.o integration.o opt.o vib.o msd.o multi_tau.o \
     radial_distrib.o efg.o NL_field.o vacf.o runmd.o 

#OBJ_= constants.o io_file.o control.o prop.o time.o functions.o \
     fft_t.o rand.o lattice.o config.o md.o lbfgs.o m1qn3.o \
     block_average.o thermo.o kspace.o rspace.o kinetic.o field.o \
     tools.o integration.o opt.o vib.o msd.o multi_tau.o \
     radial_distrib.o efg_new.o vacf.o runmd.o

#test: $(EXE_)
#	cp $(EXE_) $(INSTALL_DIR)

all: $(EXE)

mdff_cell.x: $(OBJ) mdff.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

#mdff_.x: $(OBJ_) mdff.o
#	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

fft_t.o:
	$(FC) $(FCFLAGS) $(INC) -c fft_t.f90

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< 

.PHONY: clean veryclean


clean:
	rm -f *.o *.mod *.MOD *__genmod.f90 *__genmod.mod

veryclean: clean
	rm -f $(EXE)
	rm -f *~ $(INSTALL_DIR)/$(EXE)

install: all
	cp $(EXE) $(INSTALL_DIR)

tar :
	@if test -f mdff.tar.gz ; then /bin/rm mdff.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v -e -e'/\.' -e'\.o$$' \
	     -e'\.mod$$'  -e'\.x$$' \
             -e'~$$'  | xargs tar rvf mdff.tar
	gzip mdff.tar

depend:
	@echo 'Checking dependencies...'


