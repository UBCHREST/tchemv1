include ../../../MakeIncl.Vars

TCHEM_DIR = ../../../
TCHEM_LIB = $(TCHEM_DIR)/lib

EXEC = ign
SRCS = ign.cpp StiffInteg.cpp

OBJS = $(SRCS:.cpp=.o)
INCS = StiffInteg.h

CVODE_INC = $(CVODE_INSTALL)/include
CVODE_LIB = -L$(CVODE_INSTALL)/lib  -lsundials_cvode -lsundials_nvecserial 

TCHEM_INC = $(TCHEM_DIR)/include
INCLUDES = -I$(CVODE_INC) -I$(TCHEM_INC)

DFLAGS = -DUSEJAC -DALLSPEC

LIBS = $(CVODE_LIB) -L$(TCHEM_LIB) -ltchem

ifeq "$(output)" "no"
DFLAGS += -DNO_OUTPUT
endif

ECHO = echo

default : exec

exec : $(OBJS) $(INCS)
	$(CPP) -o $(EXEC) $(OBJS) $(LIBS)

classes : $(OBJS)

clean:
	/bin/rm -rf *.o *~ ../include/*~

distclean: clean
	/bin/rm -rf $(ARCH) $(EXEC)
	cd ../lib; make clean; cd ../src

.cpp.o : $(SRCS) $(INCS)
	@$(ECHO) "----------------------------------------------------------"
	@$(ECHO) "Compiling CPP Source File --" $<
	@$(ECHO) "----------------------------------------------------------"
	@echo "$(CPP) -c $(INCLUDES) $(CPPFLAGS) $(DFLAGS) -o $@ $<" ; \
	$(CPP) -c $(INCLUDES) $(CPPFLAGS) $(DFLAGS) -o $@ $<
