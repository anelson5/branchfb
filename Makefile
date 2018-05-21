OS := $(shell uname -s)
FCSUN = f95
FCLIN = ifort
#FCLIN = f95
#FFLAGSSUN = -xtypemap=real:64
FFLAGSSUN = -xtarget=ultra -xtypemap=real:64 -O4
#FFLAGSSUN = -fast -xtarget=ultra3 -m64 -xarch=sparcvis
FFLAGSLIN = -r8  -O3 -p
#LFLAGSSUN = -03
#LFLAGSLIN = -i_dynamic
LDFLAGS = -Wl,-rpath,/usr/local/lib64 -L/usr/local/lib64
LIBS    = -llapack-$(FCLIN) -lblas-$(FCLIN) -llinpack
OBJS  = grid_mod.o bobconst_mod.o param_mod.o gc_mod.o mg_mod.o \
        diffusion_mod.o fluidsolve_mod.o  \
	gel_parm.o driver.o

.SUFFIXES: .f90
.SILENT:

#
# use different compilers depending on the machine
#
ifeq ($(OS),SunOS)
  FC     = $(FCSUN)
  FFLAGS = $(FFLAGSSUN)
  LFLAGS = $(LFLAGSSUN)
  XTARG := -$(shell fpversion -foption)
endif

ifeq ($(OS),Linux) 
  FC     = $(FCLIN)
  FFLAGS = $(FFLAGSLIN)
  LFLAGS = $(LFLAGSLIN)
endif

%.mod : %.o
	@if [! -f $@ ]; then \
	  rm $< \
	  $(MAKE) $< \
	fi

%.o : %.f90
	echo "Compiling $<"
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS) $(LIBS)

%.o : %.f
	echo "Compiling $<"
	$(FC) $(FFLAGS) -c -o $@ $< 

clot: $(OBJS) 
	echo "Creating Program $@"
	$(FC) $(FFLAGS) $(XTARG) $(LFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)


clean:
	rm -f *.o *.mod *~ 


