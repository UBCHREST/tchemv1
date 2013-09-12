SHELL = /bin/sh

include ../MakeIncl.Vars

TOPDIR   = ..
INCLUDES = $(TOPDIR)/include
SRCS     = $(TOPDIR)/src
LIBDIR   = $(TOPDIR)/lib

TCHEMLIB = $(LIBDIR)/libtchem.a
TCHEMOBJ = $(SRCS)/TC_main.o $(SRCS)/TC_kmodint.o
TCHEMSUB = $(SRCS)/TC_chg.c    \
           $(SRCS)/TC_edit.c   \
           $(SRCS)/TC_mlms.c   \
           $(SRCS)/TC_rr.c     \
           $(SRCS)/TC_spec.c   \
           $(SRCS)/TC_src.c    \
           $(SRCS)/TC_thermo.c \
           $(SRCS)/TC_init.c   \
           $(SRCS)/TC_utils.c  \
           $(SRCS)/TC_for.c 
TCHEMINC = $(INCLUDES)/TC_defs.h $(INCLUDES)/TC_interface.h \
           $(INCLUDES)/TC_params.h $(INCLUDES)/TC_kmodint.h

all: $(TCHEMLIB) 
lib: $(TCHEMLIB)
obj: $(TCHEMOBJ)

$(TCHEMLIB): $(TCHEMOBJ)
	if [ ! -d "$(LIBDIR)" ]; then\
	  mkdir $(LIBDIR);\
	fi
	ar cr $(TCHEMLIB) $(TCHEMOBJ) 

clean:
	/bin/rm -f $(TCHEMLIB) $(TCHEMOBJ) *~

distclean: clean

$(SRCS)/TC_main.o: $(SRCS)/TC_main.c $(TCHEMSUB) $(TCHEMINC)
	$(CC) -I$(INCLUDES) $(CFLAGS) $(DFLAGS) -o $*.o -c $(SRCS)/TC_main.c

$(SRCS)/TC_kmodint.o: $(SRCS)/TC_kmodint.c $(TCHEMINC)
	$(CC) -I$(INCLUDES) $(CFLAGS) $(DFLAGS) -o $*.o -c $(SRCS)/TC_kmodint.c
