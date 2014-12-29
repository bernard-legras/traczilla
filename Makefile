#$Id: 
#
#=======================================================================
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TRACZILLA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#=====|==1=========2=========3=========4=========5=========6=========7==
#
TARGET=X64-par-ng
HOME=/home/legras

SHELL = /bin/bash
FC       = pgf95
CC	= pgcc

ifeq ($(TARGET),AMD64-88)
  MAIN = TRACZILLA-D
  CFLAGS = -msse -m64 -DAMD64
  FFLAGS  = -fast -Mvect=sse -i8 -r8 -tp k8-64 -Minfo -Mextend -byteswapio
  LDFLAGS1  = -L/usr/local64/lib -lgrib
endif
ifeq ($(TARGET),AMD64-84)
  MAIN = TRACZILLA-T
  CC = pgcc
  CFLAGS = -fastsse -tp k8-64 -DAMD64
  FFLAGS  = -fastsse -i8 -r4 -tp k8-64 -Minfo -Mextend -byteswapio
  LDFLAGS1  = -lgrib84
endif
ifeq ($(TARGET),AMD64-44)
  MAIN = TRACZILLA-SS
  CC = pgcc
  #CC = pgcc
  CFLAGS = -fastsse -std=c89 -tp amd64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Bstatic 
  FFLAGS = -Mextend -byteswapio -tp amd64 -fastsse -O2 -Minfo -Minline -Bstatic -Ktrap=inv,divz,ovf
  LDFLAGS1 = -L$(HOME)/lib.amd64 -lgribex -lcom
endif
ifeq ($(TARGET),X64-par)
  MAIN = TRACZILLA-TT-par
  CC = pgcc
  CFLAGS = -mp -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf  -Mcache_align -Msmartalloc
  #FFLAGS = -I$(HOME)/local/include -gopt -mp -Mbounds -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  # include needed for netcdf.mod
  FFLAGS = -I$(HOME)/local/include -Mpreprocess -DPAR_RUN -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  #FFLAGS = -I$(HOME)/local/include -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Mbounds
  LDFLAGS1 = -rpath $(HOME)/local/lib.x64 -L$(HOME)/local/lib.x64 -lgribex -lcom -lnetcdf -lnetcdff
endif
ifeq ($(TARGET),X64-par-ng)
  MAIN = TRACZILLA-TT-par-ng
  CC = pgcc
  CFLAGS = -mp -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf  -Mcache_align -Msmartalloc
  # include needed for netcdf.mod
  FFLAGS = -I$(HOME)/local-ng/include -Mpreprocess -DPAR_RUN -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  #FFLAGS = -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Mbounds 
  LDFLAGS1 = -rpath $(HOME)/local-ng/lib -L$(HOME)/local-ng/lib -lgribex -lcom -lnetcdf -lnetcdff
endif
ifeq ($(TARGET),X64-ng)
  MAIN = TRACZILLA-TT-ng
  CC = pgcc
  CFLAGS = -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf  -Mcache_align -Msmartalloc
  FFLAGS = -I$(HOME)/local/include -Mpreprocess -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  LDFLAGS1 = -rpath $(HOME)/local-ng/lib.x64 -L$(HOME)/local-ng/lib.x64 -lgribex -lcom -lnetcdff -lnetcdf -L${HOME}/local-ng/grib_api/lib -lgrib_api_f90 -lgrib_api 
endif
ifeq ($(TARGET),win)
  MAIN = TRACZILLA-win
  CC = gcc
  FC=gfortran
  FFLAGS = -O2 -fopenmp -pedantic -std=f2003 -Wall -fconvert=swap -fall-intrinsics -fmax-errors=0
  CFLAGS = -O2 
  LDFLAGS1 =
endif

.SUFFIXES :	.f90

.f90.mod:
	        $(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

#
OBJECTS0 = commons.mod date.mod coord.mod thermo.mod io.mod lyapunov.mod\
isentrop_h.mod isentrop_m.mod ecmwf_diab.mod readinterp.mod\
sphereharmspe.mod mass_iso.mod ecmwf_inct.mod merra.mod combin.mod \
demar.mod advect.mod

OBJECTS1 =  TRACZILLA.o

OBJECTS2 = misc.o uvip3p.o random.o 

$(MAIN): $(OBJECTS0) $(OBJECTS1) $(OBJECTS2)
	$(FC) *.o -o $(MAIN) $(FFLAGS) $(LDFLAGS1)
$(OBJECTS1): $(OBJECTS0)

.PHONY : clean
clean:
	\rm *.o *.mod  
	ln -s ~/local/grib_api/include/grib_api.mod
#=====|==1=========2=========3=========4=========5=========6=========7==
#
#$Log: 
