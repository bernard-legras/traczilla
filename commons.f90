
!**********************************************************************
! Copyright 1996, 1997, 2001, 2002, 2006, 2007, 2012, 2013           *
! Andreas Stohl, Bernard Legras                                       *
!                                                                     *
! This file is part of TRACZILLA which is derived from FLEXPART V6    *
!                                                                     *
! TRACZILLA is free software: you can redistribute it and/or modify   *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! TRACZILLA is distributed in the hope that it will be useful,        *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with TRACZILLA.  If not, see <http://www.gnu.org/licenses/>.  *
!**********************************************************************

!=======================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TRACZILLA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7==
!
module commons

!******************************************************************************
!   Include file for calculation of particle trajectories (Program FLEXPART)   *
!        This file contains the parameter statements used in FLEXPART          *
!                                                                              *
!        Author: A. Stohl                                                      *
!                                                                              *
!        1997                                                                  *
!                                                                              *
!        Stripped down: 5/03/2013 by B. Legras
!                                                                              *
!*******************************************************************************

      implicit none
 
      integer, parameter:: dp=kind(0.0),sp=kind(0.0)
      
!********************************************************
! Number of threads used in parallel sections of the code
!********************************************************

      integer :: num_threads

!***********************************************************
! Number of directories/files used for FLEXPART input/output
!***********************************************************

      integer, parameter :: numpath=4

! numpath       L         Number of different pathnames for input/output files


!*****************************
! Physical and other constants
!*****************************

      real(dp), parameter :: pi=3.14159265_dp,r_earth=6.371e6_dp,r_air=287.05_dp,ga=9.81_dp
      real(dp), parameter :: cpa=1004.6_dp,kappa=0.2857_dp

! pi            L          number "pi"
! r_earth       L          radius of earth [m]
! r_air         L          individual gas constant for dry air [J/kg/K]
! ga            L          gravity acceleration of earth [m/s**2]
! cpa           L          specific heat for dry air
! kappa         L          exponent of formula for potential temperature


!********************
! Some time constants
!********************

      integer, parameter :: idiffnorm=10800,idiffmax=2*idiffnorm,minstep=1

! idiffnorm [s]           normal time interval between two wind fields
! idiffmax [s]            maximum time interval between two wind fields
! minstep [s]             minimum time step to be used within FLEXPART


!*****************************************************************
! Parameters for polar stereographic projection close to the poles
!*****************************************************************

      real, parameter :: switchnorth=75.,switchsouth=-75.

! switchnorth    use polar stereographic grid north of switchnorth
! switchsouth    use polar stereographic grid south of switchsouth


!*********************************************
! Maximum dimensions of the input mother grids
!*********************************************

!     parameter(nxmax=361,nymax=181,nuvzmax=51,nwzmax=51) ! 1 degree 51 levels
      integer, parameter :: nxmax=361,nymax=181,nuvzmax=92,nwzmax=92 ! 1 degree 61 levels
!     parameter(nxmax=721,nymax=181,nuvzmax=61,nwzmax=61) ! 1/2 degre 61 levels


! nxmax,nymax             maximum dimension of wind fields in x and y
!                         direction, respectively
! nuvzmax,nwzmax          maximum dimension of (u,v) and (w) wind fields in z
!                         direction (for fields on eta levels)
! nzmax                   maximum dimension of wind fields in z direction
!                         for the transformed Cartesian coordinates


!*********************************
! Parameters for GRIB file decoding
!*********************************

      integer, parameter :: jpack=4*nxmax*nymax,jpunp=4*jpack

! jpack,jpunp  L           maximum dimensions needed for GRIB file decoding


!**************************************************************
! Maximum number of particles, species, wind fields and similar
!**************************************************************

!     integer, parameter :: maxpart=60000000,maxpoint=50,maxspec=1
      integer, parameter :: maxpart=600000,maxpoint=50,maxspec=1
      integer, parameter :: maxwf=150000,maxtable=1000,numclass=9,ni=11

! maxpart                 Maximum number of particles
! maxpoint                Maximum number of release locations
! maxspec                 Maximum number of chemical species per release
! maxwf                   maximum number of wind fields to be used for simulation
! maxtable                Maximum number of chemical species that can be
!                         tabulated for FLEXPART
! numclass                Number of landuse classes available to FLEXPART
! ni                      Number of diameter classes of particles


!*********************************
! Dimension of random number field
!*********************************

      integer, parameter :: maxrand=2000000

! maxrand                 number of random numbers used


!************************************
! Unit numbers for input/output files
!************************************

      integer, parameter :: unitpath=1,unitcommand=1,unitageclasses=1,unitgrid=1
      integer, parameter :: unitmozaic=1
      integer, parameter :: unitavailab=1,unitreleases=1
      integer, parameter :: unitpartout=95,unitpartout2=96,unitpartout3=97,unitpartout4=98
      integer, parameter :: saveunit=101,saveunit_tmp=102
      integer, parameter :: unitpartin=93,unitflux=94,unitflux2=92
      integer, parameter :: unitvert=1,unitoro=1,unitpoin=1,unitreceptor=1
      integer, parameter :: unitoutgrid=97,unitoutgridppt=99,unitoutinfo=1
      integer, parameter :: unitspecies=1,unitoutrecept=91,unitoutreceptppt=92
      integer, parameter :: unitER2_O3=81,unitER2_MG=82,unitER2_MM=83,unitER2_FP=84
      integer, parameter :: unitlsm=1,unitsurfdata=1,unitland=1,unitwesely=1

!**********************************************
! constant pressure to define vertical velocity
!**********************************************

      real(dp), parameter :: p0=100000._dp

!*******************************************************************************
!        Include file for particle diffusion model FLEXPART                    *
!        This file contains a global common block used by FLEXPART             *
!                                                                              *
!        Author: A. Stohl                                                      *
!                                                                              *
!        June 1996                                                             *
!                                                                              *
!        Last update: 9 August 2000  
!       
!        Changes : B. Legras, Apr 2002 
!                  fields on z grid cancelled    
!                  fields on nested grid cancelled                             *
!                                                                              *
!*******************************************************************************

!****************************************************************
! Variables defining where FLEXPART input/output files are stored
!****************************************************************

      character(len=128):: path(numpath)
      integer len_path(numpath)

! path                    path names needed for trajectory model
! len_path                length of path names needed for trajectory model


!********************************************************
! Variables defining the general model run specifications
!********************************************************

      integer ibdate,ibtime,iedate,ietime,hrstart,itime0
      double precision bdate,edate
      logical command_old,releases_old,restart,perpetual,shuffling
      character(len=4):: vert_interpol


! ibdate                  beginning date (YYYYMMDD)
! ibtime                  beginning time (HHMISS)
! iedate                  ending date (YYYYMMDD)
! ietime                  ending time (HHMISS)
! bdate                   beginning date of simulation (julian date)
! edate                   ending date of simulation (julian date)
! vert_interpol           type of vertical interpolation
! restart                 restart run if .true.
! hrstart                 restart hour for restart run
! perpetual               perpetual run for age of air calculations
! itime0                  starting time, non zero for restart run
! shuffling               controls whether the parcels are treated in order or shuffled


      integer ldirect,ideltas

! ldirect                 1 for forward, -1 for backward simulation
! ideltas                 length of trajectory loop from beginning to
!                         ending date (s)

      integer loutstep,loutaver,loutsample,method,lsynctime,loffset,loffset2
      logical savfull

! loffset [s]             offset time for the first output (useful if it is
!                         not a multiple of loutstep)
! loffset2 [s]		  second offset of the first output if not in sequence
!                         with the followers
! loutstep [s]            gridded concentration output every loutstep seconds
! loutaver [s]            concentration output is an average over [s] seconds
! loutsample [s]          sampling interval of gridded concentration output
! lsynctime [s]           synchronisation time of all particles
! method                  indicator which dispersion method is to be used
! savfull		  makes backups of the trajectories in full resolution

      real ctl,fine
      integer ifine,iout,ipout,ipin,iflux
      logical turbswitch

! ctl      factor, by which time step must be smaller than Lagrangian time scale
! ifine    reduction factor for time step used for vertical wind
!          Langevin equation for the vertical wind component
! iflux    flux calculation options: 1 calculation of fluxes, 2 no fluxes
! iout     output options: 1 conc. output (ng/m3), 2 mixing ratio (pptv), 3 both
! ipout    particle dump options: 0 no, 1 every output interval, 2 only at end
! ipin     read in particle positions from dumped file from a previous run
! fine     float(ifine)
! turbswitch              determines how the Markov chain is formulated

      integer mintime,itsplit

! mintime                 minimum time step to be used by FLEXPART
! itsplit                 time constant for splitting particles

      integer lsubgrid,lconvection,lagespectra

! lsubgrid     1 if subgrid topography parameterization switched on, 2 if not
! lconvection  1 if convection parameterization switched on, 0 if not
! lagespectra  1 if age spectra calculation switched on, 2 if not


!*********************************************************************
! Variables defining the release locations, released species and their
! properties, etc.
!*********************************************************************

      character(len=45):: compoint(maxpoint)
      integer numpoint,ireleasestart(maxpoint),ireleaseend(maxpoint)
      integer ireleaseinterval(maxpoint)
      real xpoint1(maxpoint),ypoint1(maxpoint)
      real xpoint2(maxpoint),ypoint2(maxpoint)
      real zpoint1(maxpoint),zpoint2(maxpoint)
      real mesh_size_lat, mesh_size_long
      logical switch_diff_off
      integer delay_switch_diff_off
      logical TTLactiv, AccurateTemp, AGEFactiv, TTLFILLactiv
      logical delayed_initialization, press2theta, theta2press
      logical uniform_spread
      real, allocatable :: lat_list(:)
      logical make_curtain, make_layer, make_uni3D, uniform_mesh
      integer curtain_type
      real lev1,lev2,inclev
      real uppertheta
      real pcut,thetacut, thetalowcut
      
      integer npart(maxpoint)
      character(len=12):: release_plan

      


! compoint                comment, also "name" of each starting point
! numpoint                actual number of trajectory starting/ending points
! ireleasestart,ireleaseend [s] starting and ending time of each release
! ireleaseinterval        interval between successive release [s]
! release_plan	L	  type of release method
! xmass                   total mass emitted
! xpoint1,ypoint1         lower left coordinates of release area
! xpoint2,ypoint2         upper right coordinates of release area
! zpoint1,zpoint2         min./max. z-coordinates of release points
! npart                   number of particles per release point
! nspec                   number of different species allowed for one release
! species                 name of species
! link                    index referring each species to the species table


!******************************************************************************
! Variables associated with the ECMWF meteorological input data ("wind fields")
!******************************************************************************

      integer numbwf,wftime(maxwf),wfldat(maxwf),wfltim(maxwf),lwindinterv
      character(len=16):: wfname(maxwf),wfspec(maxwf)

! lwindinterv [s]         Interval between wind fields currently in memory
! numbwf                  actual number of wind fields
! wftime(maxwf) [s]       times relative to beginning time of wind fields
! wfname(maxwf)           file names of wind fields
! wfspec(maxwf)           specifications of wind field file, e.g. if on hard 
!                         disc or on tape

      integer memtime(2),memind(2)

! memtime [s]             validation times of wind fields in memory
! memind                  pointer to wind field, in order to avoid shuffling
!                         of wind fields


!***********************************************************
! Variable associated with the usage of z velocities 
! or isentropic motion (see also headings of isentrop.f90)
!***********************************************************

     logical z_motion

! z_motion                use z=ln(p0/p) as vertical coordinate
!                         and interpolate wind to z level
	
!***********************************************************
! Variable associated with the usage of diab velocities 
! (see also headings of ecmwf_diab, ecmwf_inct, merra
!***********************************************************
		
logical clear_sky, cloud_sky
logical mass_correction, mean_diab_output

! clear_sky		  Use of clear sky tendency (LW+SW)
! cloud_sky               Use of cloud sky radiative tendency (LW+SW) no latent
!                         heat
! mass_correction         apply a correction to heating on isobaric surfaces
!                         such that mean mass flux is zero across this surface
! mean_diab_output	  output of the mean mass flux before correction	

integer numbwf_diab,wftime_diab(maxwf),lwindinterv_diab
character(len=11):: wfname_diab(maxwf),wfspec_diab(maxwf)

! lwindinterv_diab [s]    Interval between wind fields currently in memory
! numbwf_diab             actual number of wind fields
! wftime_diab(maxwf) [s]  times relative to beginning time of wind fields
! wfname_diab(maxwf)      file names of wind fields
! wfspec_diab(maxwf)      specifications of wind field file, e.g. if on hard 
!                         disc or on tape

integer memtime_diab(2),memind_diab(2)

! memtime_diab [s]        validation times of wind fields in memory
! memind_diab             pointer to wind field, in order to avoid shuffling
!                         of wind fields

real, allocatable :: w_diab(:,:,:,:)

! w_diab [K/s]            vertical velocity as tendency in potential temperature
 
	
character(len=128):: path_diab(2)
integer len_diab(2)

! path_diab               paths for diabatic directory and associated AVAILABLE file
! len_diab                lengths of the two previous strings
		
!****************************************************************************
! Variables defining actual size and geographical location of the wind fields
!****************************************************************************
      integer nx,ny,nxfield,nuvz,nwz,nz,nmixz,nlev_ec
      real dx,dy,xlon0,ylat0,dxconst,dyconst,zmax

! nx,ny,nz  L             actual dimensions of wind fields in x,y and z
!                         direction, respectively
! nuvz,nwz  L             vertical dimension of original ECMWF data
! nxfield                 same as nx for limited area fields,
!                         but for global fields nx=nxfield+1
! nmixz                   number of levels up to maximum PBL height (3500 m)

! nuvz is used for u,v components
! nwz is used for w components (staggered grid)
! nz is used for the levels in transformed coordinates (terrain-following Cartesian
! coordinates)

! nlev_ec  L              number of levels ECMWF model
! dx       L              grid distance in x direction
! dy       L              grid distance in y direction
! dxconst,dyconst         auxiliary variables for utransform,vtransform
! xlon0    L              geographical longitude and
! ylat0    L              geographical latitude
!                         of lower left grid point of nested wind fields



!*************************************************
! Variables used for vertical model discretization
!*************************************************

      real akm(nwzmax),bkm(nwzmax)
      real akz(nuvzmax),bkz(nuvzmax)

! akm,bkm: coeffizients which regulate vertical discretization of ecmwf model
!          (at the border of model layers)
! akz,bkz: model discretization coeffizients at the centre of the layers


!*****************************************
! Variables associated with vertical noise
!*****************************************

! diffus: m^2/s (diftype=1) or K^2/s (diftype=2)
! diftype: diffusion type. 1: in z; 2: in pot. temp.
! nsample

       real diffus
       logical verdiff,hordiff
       integer diftype
       integer nsample(maxpoint),numlevel(maxpoint)
       
!*****************************************
! Variables associated with MOZAIC and ER2  flights
!*****************************************

       integer start_index, end_index, n_sample, n_loc
       character(len=48)::  MOZAIC_filename, ER2_day
       character(len=128):: MOZAIC_dir, ER2_dir
       character(len=12):: index_type, campaign
       logical instant_release,interp_release
       real start_ER2time, end_ER2time

!*********************************
! Variables associated with CLAUS 
!*********************************
   integer,parameter      ::  unitCLAUS = 85                
   character(len=100)     ::  Claus_dir
   logical                ::  diabatic_Claus, CLAUSactiv
   integer                ::  TB_max
   integer                ::  latmin_Claus
   integer                ::  latmax_Claus
   


! Fixed fields, unchangeable with time
!*************************************

      real, allocatable :: oro(:,:), lsm(:,:), excessoro(:,:)

! oro [m]              orography of the ECMWF model
! excessoro            excess orography mother domain
! lsm                  land sea mask of the ECMWF model


! 3d fields
!**********

      real, allocatable :: uupol(:,:,:,:),vvpol(:,:,:,:)
      real, allocatable :: uuh(:,:,:,:),vvh(:,:,:,:),wwh(:,:,:,:)
      real, allocatable :: tth(:,:,:,:),qvh(:,:,:,:)
 
! uuh,vvh,wwh [m/2] L   wind components in x,y and z direction
! uupol,vvpol [m/s] L   hor. wind components in polar stereographic projection
! rho [kg/m3]          air density
! drhodz [kg/m2]       vertical air density gradient
! tth [K]           L   temperature on origibal eta levels
! qvh               L   water vapor mixing ratio on eta levels

! 2d fields
!**********

      real, allocatable:: ps(:,:,:,:),u10(:,:,:,:),v10(:,:,:,:),tt2(:,:,:,:)

! ps             L     surface pressure
! sd                   snow depth
! msl                  mean sea level pressure
! tcc                  total cloud cover
! u10                  10 meter u
! v10                  10 meter v
! tt2                  2 meter temperature
! td2                  2 meter dew point


!******************************************************
! Variables defining the polar stereographic projection
!******************************************************

      logical xglobal,sglobal,nglobal
      real switchnorthg,switchsouthg

!     xglobal      L       T for global fields, F for limited area fields
!     sglobal      L       T if domain extends towards south pole
!     nglobal      L       T if domain extends towards north pole
!     switchnorthg,switchsouthg L  same as parameters switchnorth,
!                                  switchsouth, but in grid units

      real southpolemap(9),northpolemap(9)

!     southpolemap,northpolemap L  define stereographic projections
!                         at the two poles


!***************************************
! Variables characterizing each particle
!***************************************

      integer :: numpart,itra1(maxpart)
!      integer nclass(maxpart)
      integer :: itramem(maxpart)
      real(dp) ::  xtra1(maxpart),ytra1(maxpart),ztra1(maxpart)
      integer, allocatable :: itra0(:)
      real(dp), allocatable :: ttra1(:), qtra1(:)

! numpart                 actual number of particles in memory
! itra1 (maxpart) [s]     temporal positions of the particles
! npoint/(maxpart)         indicates the release point of each particle
! nclass (maxpart)        one of nclassunc classes to which the particle is attributed
! itramem (maxpart) [s]   memorized release times of the particles
! itrasplit (maxpart) [s] next time when particle is to be split into two
! idt(maxpart) [s]        time step to be used for next integration
! xtra1,ytra1,ztra1       spatial positions of the particles
! xmass1 [kg]             particle masses

!**************************************
! Correcting vertical wind
!**************************************

      logical correct_vertwind
 

!********************
! Random number field
!********************

      real rannumb(maxrand)

! rannumb                 field of normally distributed random numbers


!*******************
! Debug
!******************

     logical debug_out


end module commons

!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log: 
