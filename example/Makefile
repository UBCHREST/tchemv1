SHELL = /bin/sh

all: subdirs 
default: subdirs

subdirs: 
	cd ign-c;   make; cd ..
	cd ign-cpp; make; cd ..
	cd ign-f;   make; cd ..

clean:
	cd ign-c;   make clean; cd ..
	cd ign-cpp; make clean; cd ..
	cd ign-f;   make clean; cd ..

distclean: clean
	cd ign-c;   make distclean; cd ..
	cd ign-cpp; make distclean; cd ..
	cd ign-f;   make distclean; cd ..
