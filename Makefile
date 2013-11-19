#$Id: 
#
#=======================================================================
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TRACZILLA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#=====|==1=========2=========3=========4=========5=========6=========7==
#
TARGET=X64
#TARGET=P7
HOME=/home/legras
#TARGET=XP

SHELL = /bin/bash
FC       = pgf95
CC	= pgcc

ifeq ($(TARGET),XP)
  MAIN = TRACZILLA-XP
  FFLAGS  = -tp athlonxp -fast -Mvect=sse -Minfo -Mextend -byteswapio -Bstatic
  CFLAGS  = -tp athlonxp -fast -Mvect=sse
#  LDFLAGS1 = -L/usr/local32/lib -lgrib
  LDFLAGS1 = -L$(HOME)/lib.i586 -lgribex-px -lcom
endif
ifeq ($(TARGET),AMD64-88)
  MAIN = TRACZILLA-D
  CFLAGS = -msse -m64 -DAMD64
  FFLAGS  = -fast -Mvect=sse -i8 -r8 -tp k8-64 -Minfo -Mextend -byteswapio
#  CFLAGS = -g -m64 -DAMD64
#  FFLAGS  = -g -i8 -r8 -tp k8-64 -Minfo -Mextend -byteswapio
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
  #CFLAGS = -fastsse -tp amd64 -Minfo -Mbounds -Ktrap=inv,divz,ovf -Bstatic 
  #CFLAGS = -fastsse -tp k8-64 
  #FFLAGS  = -fastsse -tp k8-64 -Minfo -Mextend -byteswapio
  FFLAGS = -Mextend -byteswapio -tp amd64 -fastsse -O2 -Minfo -Minline -Bstatic -Ktrap=inv,divz,ovf
  #FFLAGS = -Mextend -byteswapio -tp amd64 -Minfo -Bstatic -Mbounds -Ktrap=inv,divz,ovf
  #FFLAGS = -Mextend -byteswapio -tp amd64 -fastsse -O2 -Minfo -Minline  -Bstatic
  #FFLAGS = -Mextend -byteswapio -tp amd64 -Mbounds -g -Minfo -Bstatic
  #LDFLAGS1  = -L/net/gluck/home/priv/legras/install64/src/NgribexLinux -lgrib44
  LDFLAGS1 = -L$(HOME)/lib.amd64 -lgribex -lcom
endif
ifeq ($(TARGET),X64-par)
  MAIN = TRACZILLA-TT-par-test
  CC = pgcc
  #CFLAGS = -mp -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf  -Mcache_align -Msmartalloc
  CFLAGS = -mp -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf  -Mcache_align -Msmartalloc
  #FFLAGS = -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  #FFLAGS = -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,ovf -Mcache_align -Msmartalloc
  #FFLAGS = -I$(HOME)/local/include -gopt -mp -Mbounds -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  FFLAGS = -I$(HOME)/local/include -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  #FFLAGS = -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Mbounds 
  LDFLAGS1 = -rpath $(HOME)/local/lib.x64 -L$(HOME)/local/lib.x64 -lgribex -lcom -lnetcdf -lnetcdff
endif
ifeq ($(TARGET),X64)
  MAIN = TRACZILLA-TT
  CC = pgcc
  CFLAGS = -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf  -Mcache_align -Msmartalloc
  FFLAGS = -I$(HOME)/local/include -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  #FFLAGS = -I$(HOME)/local/include -Mextend -byteswapio -tp x64 -g -Ktrap=inv,divz,ovf
  #CFLAGS = -fastsse -tp x64 -O2 -Minfo -Mbounds  
  #FFLAGS = -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Mbounds -I$(HOME)/local/include
  #FFLAGS = -Mextend -byteswapio -tp x64 -Minfo -Mbounds -Bstatic -Ktrap=inv,divz,ovf 
  #FFLAGS = -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Mcache_align -Msmartalloc
  LDFLAGS1 = -rpath $(HOME)/local/lib.x64 -L$(HOME)/local/lib.x64 -lgribex -lcom -lnetcdff -lnetcdf -L${HOME}/local/grib_api/lib -lgrib_api_f90 -lgrib_api 
endif
ifeq ($(TARGET),win)
  MAIN = TRACZILLA-win
  CC = gcc
  FC=gfortran
  FFLAGS = -O2 -fopenmp -pedantic -std=f2003 -Wall -fconvert=swap -fall-intrinsics -fmax-errors=0
  CFLAGS = -O2 
  LDFLAGS1 =
endif
ifeq ($(TARGET),IFC)
  FC = ifc
  MAIN = TRACZILLA
  FFLAGS  = -132 -Vaxlib
  LDFLAGS1  = -L/net/parker/do2/privat/install-src1/gribexLinux-ifc -lgrib
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
