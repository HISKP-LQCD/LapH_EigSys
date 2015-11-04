### PETSc and SLEPc variables
#PETSC_DIR  = /hiskp2/helmes/programm_libs/petsc-3.6.2
#PETSC_ARCH = arch-linux-complex-scalar-noMPI 
#SLEPC_DIR  = /hiskp2/helmes/programm_libs/slepc-3.6.1

include /hiskp2/helmes/programm_libs/slepc-3.6.1/lib/slepc/conf/slepc_common
### paths to different included libraries
MDIR        = ../modules
VPATH       = .:$(MDIR)
LIMEPATH    = $(WORK)/helmes/programm_libs/lime-1.3.2
EIGENPATH   = $(WORK)/helmes/programm_libs/eigen-3.2.0

INCPATH     = ../include $(EIGENPATH) $(LIMEPATH)/include /usr/local/include\
							$(PETSC_DIR)/include $(PETSC_DIR)/$(PETSC_ARCH)/include
							#$(PETSC_DIR)/include/mpiuni 
							#$(SLEPC_DIR)/include $(SLEPC_DIR)/arch-linux-complex-scalar-noMPI/include\
							#$(SLEPC_DIR)/$(PETSC_ARCH)/include
							#$(PETSC_DIR)/include 
               
LIBPATH     = $(LIMEPATH)/lib 
LIBS        = m lime 

### C and C++ options
CPPFLAGS    = -std=c++0x -O2 -Wall -fopenmp

### 
MAIN    = ev_ts
GENERAL = config_utils eigensystem io io_utils navigation par_io read_write recover_spec shell_matop timeslice
MODULES = $(GENERAL)

### programs to use
CC	   = mpicc.openmpi 
CPP	   = g++
CLINKER    = $(CC)


### begin of rules
all: rmexeq mkdep mkexeq

.PHONY: all



PGMS = $(MAIN) $(MODULES)

-include $(addsuffix .d,$(PGMS))

# rule to make dependencies
$(addsuffix .d, $(PGMS)): %.d: %.cpp Makefile
	@ $(CC) -ansi $< -MM $(CPPFLAGS) $(addprefix -I, $(INCPATH)) -o $@

# rule to compile source programs
$(addsuffix .o, $(PGMS)): %.o: %.cpp Makefile
	$(CC) $< -c $(CPPFLAGS) $(addprefix -I, $(INCPATH)) ${SLEPC_INCLUDE}
 
# rule to link main programs
$(MAIN): %: %.o $(addsuffix .o, $(MODULES)) Makefile
	$(CLINKER) $(CPPFLAGS) $< $(addsuffix .o, $(MODULES)) \
	$(addprefix -L, $(LIBPATH)) $(addprefix -l,$(LIBS)) $(addprefix -L, ${SLEPC_LIB_DIR}) ${SLEPC_EPS_LIB} -o $@

# produce executables
mkexeq: $(MAIN)

# remove executables
rmexeq:
	@ echo "delete old executables"
	@ -rm -f $(MAIN)

# make dependencies
mkdep: $(addsuffix .d,$(PGMS))
	@ echo "generate table of dependencies"

# clean directory
# clean is already predefined in the includes above

