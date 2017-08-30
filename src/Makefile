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

TCHEMU_SRC = $(SRCS)/TC_reduce.cpp
TCHEMU_OBJ = $(SRCS)/TC_reduce.o
TCHEMU_LIB = $(LIBDIR)/libtchemutil.a

all: $(TCHEMLIB) $(TCHEMU_LIB)
lib: $(TCHEMLIB) $(TCHEMU_LIB)
obj: $(TCHEMOBJ) $(TCHEMU_OBJ)
tchemlib: $(TCHEMLIB)

$(TCHEMLIB): $(TCHEMOBJ)
	if [ ! -d "$(LIBDIR)" ]; then\
	  mkdir $(LIBDIR);\
	fi
	ar cr $(TCHEMLIB) $(TCHEMOBJ) 

$(TCHEMU_LIB): $(TCHEMU_OBJ)
	if [ ! -d "$(LIBDIR)" ]; then\
	  mkdir $(LIBDIR);\
	fi
	ar cr $(TCHEMU_LIB) $(TCHEMU_OBJ) 

clean:
	/bin/rm -f $(TCHEMLIB) $(TCHEMOBJ) $(TCHEMU_LIB) $(TCHEMU_OBJ) *~

distclean: clean

$(SRCS)/TC_main.o: $(SRCS)/TC_main.c $(TCHEMSUB) $(TCHEMINC)
	$(CC) -I$(INCLUDES) $(CFLAGS) $(DFLAGS) -o $*.o -c $(SRCS)/TC_main.c

$(SRCS)/TC_kmodint.o: $(SRCS)/TC_kmodint.c $(TCHEMINC)
	$(CC) -I$(INCLUDES) $(CFLAGS) $(DFLAGS) -o $*.o -c $(SRCS)/TC_kmodint.c

$(SRCS)/TC_reduce.o: $(SRCS)/TC_reduce.cpp $(TCHEMINC)
	$(CXX) -I$(INCLUDES) $(CXXFLAGS) $(DFLAGS) -o $*.o -c $(SRCS)/TC_reduce.cpp
