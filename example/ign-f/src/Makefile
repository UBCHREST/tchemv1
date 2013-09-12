include ../../../MakeIncl.Vars

TCHEM_DIR = ../../..
TCHEM_LIB = $(TCHEM_DIR)/lib

EXEC1 = ign
EXEC2 = ign-mod
SRCS  = ign.f rhsjac.f

DFLAGS = -DNO_OUTPUT
MAKECLEAN = clean

LIBS = -L$(DVODE_INSTALL) -ldvode -L$(TCHEM_LIB) -ltchem

ECHO = echo
LINKFFLS = ln -fs ign.f ign.F

ifeq "$(output)" "yes"
DFLAGS = 
endif

ifeq "$(OS_NAME)" "Darwin"
LINKFFLS = 
endif

default : execgri3 exec1step

execgri3 : $(MAKECLEAN) 
	$(LINKFFLS)
	$(F77) -c $(FFLAGS) $(DFLAGS) ign.F
	$(F77) -o $(EXEC1) ign.o $(LIBS)

exec1step : $(MAKECLEAN) 
	$(LINKFFLS)
	$(F77) -c $(FFLAGS) $(DFLAGS) -DMODKINMECH ign.F
	$(F77) -o $(EXEC2) ign.o $(LIBS)

clean:
	/bin/rm -rf *.o *~ ../include/*~

distclean: clean
	/bin/rm -rf $(EXEC1) $(EXEC2) 
