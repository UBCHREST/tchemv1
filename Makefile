SHELL = /bin/sh

include MakeIncl.Vars

all: library $(EXAMPLEDIR) $(DOCDIR)

library:
	cd src; make

tchemlib:
	cd src; make tchemlib

ex: 
	cd example; make; cd ..

doxy: Doxyfile
	$(DOXY)

doc:
	cd doc; make; cd ..

clean:
	cd src;     make clean; cd ..
	cd doc;     make clean; cd ..
	cd example; make clean; cd ..

distclean: clean
	/bin/rm -rf lib/libtchem.a
	cd example; make distclean; cd .. 


