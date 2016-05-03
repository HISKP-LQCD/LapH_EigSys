### PETSc and SLEPc variables
#PETSC_DIR  = $(HOME)/code/petsc
#PETSC_ARCH = arch-linux2-cxx-opt
#SLEPC_DIR  = $(HOME)/code/slepc

PETSC_INC = $(PETSC_DIR)/$(PETSC_ARCH)/include $(PETSC_DIR)/include $(PETSC_DIR)/include/mpiuni
SLEPC_INC = $(SLEPC_DIR)/$(PETSC_ARCH)/include $(SLEPC_DIR)/include

PETSC_LIB = $(PETSC_DIR)/$(PETSC_ARCH)/lib
SLEPC_LIB = $(SLEPC_DIR)/$(PETSC_ARCH)/lib

### paths to different included libraries
MDIR        = ../modules
VPATH       = .:$(MDIR)
LIMEPATH    = $(HOME)/code/lime-1.3.2
EIGENPATH   = $(HOME)/code/eigen

INCPATH     = ../include $(EIGENPATH) $(LIMEPATH)/include $(PETSC_INC) $(SLEPC_INC)
LIBPATH     = $(LIMEPATH)/lib $(PETSC_LIB) $(SLEPC_LIB)
LIBS        = m lime petsc slepc mpicxx

### C and C++ options
CPPFLAGS    = -std=c++0x -Ofast -Wall -pedantic -openmp

### 
MAIN    = test_eigensystems
GENERAL = timeslice config_utils eigensystem io io_utils navigation par_io read_write recover_spec shell_matop
MODULES = $(GENERAL)

### programs to use
CC	   = g++
CPP	   = gcc
#CC	   = mpic++
#CPP	   = mpic++
CLINKER    = $(CC)


### begin of rules
all: rmexeq mkdep mkexeq

.PHONY: all

#include $(PETSC_DIR)/conf/petscvariables
#include $(SLEPC_DIR)/conf/slepc_common

PGMS = $(MAIN) $(MODULES)

-include $(addsuffix .d,$(PGMS))

# rule to make dependencies
$(addsuffix .d, $(PGMS)): %.d: %.cpp Makefile
	@ $(CC) -ansi $< -MM $(CPPFLAGS) $(addprefix -I, $(INCPATH)) $(CC_INCLUDES) -o $@

# rule to compile source programs
$(addsuffix .o, $(PGMS)): %.o: %.cpp Makefile
	$(CC) $< -c $(CPPFLAGS) $(addprefix -I, $(INCPATH)) $(CC_INCLUDES)

# rule to link main programs
$(MAIN): %: %.o $(addsuffix .o, $(MODULES)) Makefile
	$(CLINKER) $(CPPFLAGS) $< $(addsuffix .o, $(MODULES)) \
	$(addprefix -L, $(LIBPATH)) $(addprefix -l,$(LIBS)) -o $@

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

myclean:
	@ -rm -rf *.d *.o *.alog *.clog *.slog $(MAIN)
.PHONY: myclean

