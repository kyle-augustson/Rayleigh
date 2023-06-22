# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name
NAME := ULTRA

# Configuration settings
MKLROOT = /opt/intel/mkl
INCL = -I$(MKLROOT)/include
OPTS := -O0 -g
FC := gfortran $(OPTS) $(INCL)
AR := ar rcs
RM := rm -f

LIBS = -L$(MKLROOT)/lib/intel64/
LD = $(FC) $(LIBS) -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# List of all source files
SRCS := Ultraspherical_Polynomials_NoSparse.F90
TEST_SRCS := Test_Ultraspherical_NoSparse.F90

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
TEST_OBJS := $(addsuffix .o, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
BIN := $(patsubst %.F90, %.bin, $(TEST_SRCS))

# Declare all public targets
.PHONY: all clean
all: $(LIB) $(BIN)

# Create the static library from the object files
$(LIB): $(OBJS)
	$(AR) $@ $^

# Link the test executables
$(BIN): %.bin: %.F90.o $(LIB)
	$(LD) -o $@ $^

# Create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: %
	$(FC) -c -o $@ $<

# Define all module interdependencies
Ultraspherical_Polynomials_NoSparse.mod := Ultraspherical_Polynomials_NoSparse.F90.o
Test_Ultraspherical_NoSparse.F90.o: $(Ultraspherical_Polynomials_NoSparse.mod)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(filter %.bin, $(BIN)) $(filter %.a, $(LIB)) $(wildcard *.mod)
