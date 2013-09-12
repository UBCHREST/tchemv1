include ../../../MakeIncl.Vars

TCHEM_DIR = ../../..
TCHEM_LIB = $(TCHEM_DIR)/lib

EXEC = ign
SRCS = ign.c setup.c rhsjac.c doign.c doignReinit.c doignDtSeq.c output.c initcvode.c

OBJS = $(SRCS:.c=.o)
INCS = ign.h TC_interface.h

CVODE_INC = $(CVODE_INSTALL)/include
CVODE_LIB = -L$(CVODE_INSTALL)/lib  -lsundials_cvode -lsundials_nvecserial 

TCHEM_INC = $(TCHEM_DIR)/include
INCLUDES = -I. -I$(CVODE_INC) -I$(TCHEM_INC)

DFLAGS = -DUSEJAC -DALLSPEC
#DFLAGS = -DALLSPEC

LIBS = $(CVODE_LIB) -L$(TCHEM_LIB) -ltchem -lm

ifeq "$(output)" "no"
DFLAGS += -DNO_OUTPUT
endif
ifeq "$(chemmod)" "yes"
DFLAGS += -DCHEMMODRUN
endif

ifeq "$(gettig)" "yes"
DFLAGS += -DGETTIG
endif

ECHO = echo

default : exec

exec : $(OBJS)
	$(CC) -o $(EXEC) $(OBJS) $(LIBS)

classes : $(OBJS)

clean:
	/bin/rm -rf *.o *~ ../include/*~

distclean: clean
	/bin/rm -rf $(EXEC)

.c.o :  $(SRCS) $(INCS)
	@$(ECHO) "----------------------------------------------------------"
	@$(ECHO) "Compiling C Source File --" $<
	@$(ECHO) "----------------------------------------------------------"
	@echo "$(CC) -c $(INCLUDES) $(CFLAGS) $(DFLAGS) -o $@ $<" ; \
	$(CC) -c $(INCLUDES) $(CFLAGS) $(DFLAGS) -o $@ $<
