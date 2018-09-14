#$Id:
#
#=======================================================================
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TRACZILLA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#=====|==1=========2=========3=========4=========5=========6=========7==
#
TARGET=X64-par-ng
#TARGET=win
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
  LDFLAGS1 = -rpath $(HOME)/local/lib.x64 -L$(HOME)/local/lib.x64 -lcom -lnetcdf -lnetcdff
endif
ifeq ($(TARGET),X64-par-ng)
  MAIN = TRACZILLA-TT-par-ng
  CC = pgcc
  CFLAGS = -mp -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Msmartalloc
  # include needed for netcdf.mod
  FFLAGS = -I$(HOME)/local-ng/include -Mpreprocess -DPAR_RUN -mp -Munixlogical -Mextend -byteswapio -tp x64 -O2 -fastsse -Minfo -Mcache_align -Msmartalloc -Ktrap=inv,divz,ovf
  #FFLAGS = -I$(HOME)/local-ng/include -Mbounds -Mpreprocess -DPAR_RUN -mp -Mextend -byteswapio -tp x64 -O2 -Minfo -Ktrap=inv,divz,ovf
  #FFLAGS = -mp -Mextend -byteswapio -tp x64 -fastsse -O2 -Minfo -Mbounds 
  #LDFLAGS = -rpath $(HOME)/local-ng/lib.x64 -L$(HOME)/local-ng/lib.x64 -lgribex -lcom -lnetcdff -lnetcdf -L${HOME}/local-ng/grib_api/lib -lgrib_api_f90 -lgrib_api 
  LDFLAGS1 = -rpath $(HOME)/local-ng/lib.x64 -L$(HOME)/local-ng/lib.x64 -lcom -lnetcdff -lnetcdf -leccodes_f90 -leccodes 
endif
ifeq ($(TARGET),X64-ng)
  MAIN = TRACZILLA-TT-ng
  CC = pgcc
  CFLAGS = -fastsse -tp x64 -O2 -Minfo -Minline -Ktrap=inv,divz,ovf  -Mcache_align -Msmartalloc
  FFLAGS = -I$(HOME)/local/include -Mpreprocess -Mextend -Munixlogical -byteswapio -tp x64 -fastsse -O2 -Minfo -Minline -Ktrap=inv,divz,ovf -Mcache_align -Msmartalloc
  #FFLAGS = -I$(HOME)/local/include -g -Ktrap=inv,divz,ovf -Mbounds -Mpreprocess -Mextend -byteswapio -tp x64 -O2 -Minfo
  #LDFLAGS1 = -rpath $(HOME)/local-ng/lib.x64 -L$(HOME)/local-ng/lib.x64 -lgribex -lcom -lnetcdff -lnetcdf -L${HOME}/local-ng/grib_api/lib -lgrib_api_f90 -lgrib_api 
  LDFLAGS1 = -rpath $(HOME)/local-ng/lib.x64 -L$(HOME)/local-ng/lib.x64 -lgribex -lcom -lnetcdff -lnetcdf -leccodes_f90 -leccodes 
endif
ifeq ($(TARGET),win)
  MAIN = TRACZILLA-win
  CC = pgcc
  FC=gfortran
  FFLAGS = -O2 -fopenmp -pedantic -std=f2003 -cpp -Wall -fconvert=swap -fall-intrinsics -fmax-errors=20 -I /opt/netcdf42/gfortran/include -I /usr/lib64/gfortran/modules
  CFLAGS = -O2 
  LDFLAGS1 =
endif

.SUFFIXES :	.f90

.f90.mod:
	        $(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

#
OBJECTS0 = commons.mod date.mod coord.mod thermo.mod randme.mod io.mod interpol.mod\
isentrop_h.mod isentrop_m.mod ecmwf_diab.mod readinterpN.mod\
sphereharmspe.mod mass_iso.mod ecmwf_inct.mod merra.mod jra55.mod combin.mod \
demar.mod advect.mod

OBJECTS1 =  TRACZILLA.o

OBJECTS2 = misc.o 

$(MAIN): $(OBJECTS0) $(OBJECTS1) $(OBJECTS2)
	$(FC) *.o -o $(MAIN) $(FFLAGS) $(LDFLAGS1)
$(OBJECTS1): $(OBJECTS0) $(OBJECTS2)

.PHONY : clean
clean:
	\rm *.o 
	\rm *.mod  
	ln -s $(HOME)/local/include/grib_api.mod
	ln -s $(HOME)/local/include/netcdf.mod
lnmod:
	ln -s $(HOME)/local/include/grib_api.mod
	ln -s $(HOME)/local/include/netcdf.mod
#=====|==1=========2=========3=========4=========5=========6=========7==
#
#$Log: 
