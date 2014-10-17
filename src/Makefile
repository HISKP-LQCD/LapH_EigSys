### PETSc and SLEPc variables
PETSC_DIR  = $(HOME)/jost/petsc
#PETSC_ARCH = arch-linux2-cxx-opt
SLEPC_DIR  = $(HOME)/jost/slepc-3.4.4

### paths to different included libraries
MDIR        = ../modules
VPATH       = .:$(MDIR)
LIMEPATH    = $(HOME)/jost/lime-1.3.2
EIGENPATH   = $(HOME)/jost/eigen-eigen-6b38706d90a9

INCPATH     = ../include $(EIGENPATH) $(LIMEPATH)/include
LIBPATH     = $(LIMEPATH)/lib 
LIBS        = m lime

### C and C++ options
CPPFLAGS    = -std=c++0x -O2 -Wall

### 
MAIN    = ev_ts
GENERAL = config_utils eigensystem io io_utils navigation par_io read_write recover_spec shell_matop timeslice
MODULES = $(GENERAL)

### programs to use
CC	   = gcc
CPP	   = g++
CLINKER    = $(CC)


### begin of rules
all: rmexeq mkdep mkexeq

.PHONY: all

include $(PETSC_DIR)/conf/petscvariables
include	$(SLEPC_DIR)/conf/slepc_common

PGMS = $(MAIN) $(MODULES)

-include $(addsuffix .d,$(PGMS))

# rule to make dependencies
$(addsuffix .d, $(PGMS)): %.d: %.cpp Makefile
	@ $(CC) -ansi $< -MM $(CPPFLAGS) $(addprefix -I, $(INCPATH)) $(CC_INCLUDES) -o $@

# rule to compile source programs
$(addsuffix .o, $(PGMS)): %.o: %.cpp Makefile
	$(CC) $< -c $(CPPFLAGS) $(addprefix -I, $(INCPATH)) $(CC_INCLUDES)

# rule to link main programs
$(MAIN): %: %.o $(addsuffix .o, $(MODULES)) Makefile chkopts
	$(CLINKER) $(CPPFLAGS) $< $(addsuffix .o, $(MODULES)) \
	$(addprefix -L, $(LIBPATH)) $(addprefix -l,$(LIBS)) $(SLEPC_LIB) -o $@

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

