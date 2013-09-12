SHELL = /bin/sh

all: subdirs 
default: subdirs

subdirs: 
	cd src; make; cd ..

clean:
	cd src; make clean; cd ..

distclean: clean
	cd src; make distclean; cd ..
