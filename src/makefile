.SUFFIXES: .h .F90 .f90

###################################################
# MAKEFILE FOR MDFF ON LISA AMSTERDAM JUNE 2011
# module load fortran/intel
# module load openmpi/intel
# module load mkl
# module load lapack/intel
#
# #-lmkl_lapack -lmkl -lsvml -lguide -lpthread
# #-lscalapack -lfblacs -lcblacs -lblacs -lmkl_lapack -lmkl -lguide -lpthread
#
################################################### 
$(info ==============  MDFF compiling info ===============)
UNAME := $(shell uname -a)
HOST  := $(shell hostname -s)
#SVNREV := $(shell svn info -rHEAD | grep vision | head -n 1 | cut -d' ' -f2 )
$(info compiling on : $(HOST) )
#$(info svn revision : $(SVNREV) )
# suffix of preprocessed source file
SUFFIX=.F90

# The compiler
FC = mpif90 
#FC = gfortran 
$(info FC : $(shell which $(FC) ) )
FCL=$(FC)
# the pre-compiler
CP = cpp
CPPFLAGS=-P -C -CC -traditional

# installation directory
INSTALL_DIR=../bin
$(info instalation dir: : $(shell readlink -f $(INSTALL_DIR)) )
# List of executable to be built within the package
EXE= mdff.x 

# =====================================================
#                FORTRAN COMPILER FLAGS
# =====================================================
#
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
#   on mario / babar / charrue
# =================
FCFLAGS = -O3 -xSSE4.2 -align all 
#FCFLAGS = -O3 -xHost -align all 
#FCFLAGS = -O3 -ip -static -align all 
#FCFLAGS = -O3 -static -x  -align all
#FCFLAGS = -O3 -ip -xW
#FCFLAGS = -O0 -check all -warn all 
ifeq ($(HOST),debian-charrue)
FCFLAGS = -O3 -ip -align all
endif
$(info FCFLAGS :  $(FCFLAGS))
# =======================================
#               LIBRARIES 
# =======================================

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
#   on fer 
# =================
ifeq ($(HOST),iram-ca-004157)
	LDFLAGS= -Wl,-rpath,/opt/intel/composerxe/mkl/lib/intel64/ -Wl,--start-group /opt/intel/composerxe/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -limf -lm
	INC=-I/opt/intel/composerxe/mkl/include/fftw/
endif
# =================
#   on T5600
# =================
ifeq ($(HOST),IRAM-CE-004776)
	LDFLAGS= -Wl,-rpath,/opt/intel/composerxe/mkl/lib/intel64/ -Wl,--start-group /opt/intel/composerxe/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_sequential.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -limf -lm
	INC=-I/opt/intel/composerxe/mkl/include/fftw/
endif

# =================
#   on charrue
# =================
ifeq ($(HOST),debian-charrue)
	LDFLAGS=-Wl,-rpath,/opt/intel/composerxe/mkl/lib/ia32/ -Wl,--start-group /opt/intel/composerxe/mkl/lib/ia32/libmkl_intel.a /opt/intel/composerxe/mkl/lib/ia32/libmkl_sequential.a /opt/intel/composerxe/mkl/lib/ia32/libmkl_core.a -Wl,--end-group -limf -lm
	INC=-I/opt/intel/composerxe/mkl/include/fftw/
endif

# =================
#   on curie
# =================
ifeq ($(HOST),curie71)
	LDFLAGS=${MKL_LIBS}
	INC=-I/usr/local/fftw3-3.3.1/include/
endif

# =================
#   on drama 
# =================
ifeq ($(HOST),iram-ce-004819)
	LDFLAGS=-lmkl_intel_lp64  -lmkl_sequential -lmkl_core -lfftw3
	INC=-I/usr/local/fftw3-3.3.1/include/ -I$(MPI_INCLUDE)
endif


$(info LDFLAGS:  $(LDFLAGS))
$(info INC:  $(INC))
$(info ===================================================)
# ==================
#  target files
# ==================
OBJ= constants.o dumb.o io.o tt_damp.o tensors.o oxyde.o mpi_mdff.o control.o stochio.o time.o functions.o \
     fft_t.o rand.o lattice.o config.o md.o lbfgs.o m1qn3.o m1cg1.o \
     block_average.o kspace.o rspace.o thermo.o kinetic.o field.o harm.o \
     tools.o read_pos.o integration.o opt.o vib.o msd.o multi_tau.o \
     radial_distrib.o polarizability.o efg.o rmc.o vacf.o voisin.o runmd.o

# ================================================================
#   prepocessor :
#   -DHOST     : hostname
#   -Ddebug    : general debug flag !! too much outputs
#   -Dtiming   : timing info 
#   -DGFORTRAN : gfortran compiler 
# ================================================================
CPPDEFS = -DHOST=\"$(HOST)\"   \
          -Dtiming             \
          -DSVNREV=\"$(SVNREV)\" \
#          -Ddebug=debug

CPP =  $(CP) $(CPPFLAGS) $(CPPDEFS)  $*.f90 > $*$(SUFFIX)

all: $(OBJ) $(EXE)

$(EXE): $(OBJ) mdff.o
	$(FCL) -o $@ $^ $(LDFLAGS)

fft_t.o: fft_t.f90
	$(CPP)
	$(FC) $(FCFLAGS) $(INC) -c $*$(SUFFIX) 

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.F90: %.f90
	$(CPP)

%.o: %.f90
	$(CPP)
	$(FC) $(FCFLAGS) -c $*$(SUFFIX) 

.PHONY: clean veryclean doc

doc:
	doxygen Doxyfile


clean:
	rm -f *.F90 *.o *.mod *.MOD *__genmod.f90 *__genmod.mod 

veryclean: clean
	rm -f $(EXE)
	rm -f *~ $(INSTALL_DIR)/$(EXE)

install: all
	cp $(EXE) $(INSTALL_DIR)

tar :
	@if test -f mdff.tar.gz ; then /bin/rm mdff.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v "\.svn" | grep -v -e -e'/\./' -e'\.o$$' \
	     -e'\.mod$$'  -e'\.x$$' -e'\.F90$$' \
             -e'~$$'  | xargs tar rvf mdff.tar
	gzip mdff.tar



