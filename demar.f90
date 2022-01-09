! *********************************************************************
! Copyright  1996,2002,2004, 2007, 2012, 2013, 2014, 2016, 2020             *
! Andreas Stohl, Bernard Legras, Ann'Sophie Tissier                   *
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

!$Id:
!
!###############################################################################
!----------------------------------- DEMAR -------------------------------------
!###############################################################################

module demar
use commons
use date
use isentrop_m
use io
use interpol
implicit none

! variables local to this module
! external source of initial positions
logical :: external_pos0
integer :: year_b, year_e
character(len=256) :: external_directory
character(len=16) :: external_type

private check_launch_time, check_numpart
private check_input_date_2, check_sample
private check_input_date_1, check_output_options, check_synchro, check_consist

contains

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ READCOMMANDB @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

 subroutine readcommandB(error)

!***************************************************************************
!    This routine reads the user specifications for the current model run
!
!    Replaces readcommand.f by using namelist, providing more flexibility
!    for extensions
!
!    Author: B. Legras
!    1 June 2002
!
!***************************************************************************
!    Variables (not in namelists)
!
!  bdate                beginning date as Julian date 
!  edate                ending date as Julian date 
!  hhh                  hour
!  mi                   minute
!  ss                   second
!  error                .true., if error ocurred in subprogram, else .false.
!  ideltas [s]          modelling period
!  method               method used to compute the particle pseudovelocities  
 
!  unitcommand          unit connected to file COMMAND 
!***************************************************************************

 use isentrop_h
 use ecmwf_diab
 use ecmwf_inct
 use mass_iso
 use merra
 use jra55
 use era5
 use readinterpN, only : ecmwf_data 
 
 character(len=72):: line
 logical error
 real epsil
 
 namelist /COMMAND/     &
   ldirect,             & ! direction of integration
   ibdate, ibtime,      & ! beginnning date and time (YYYYMMDD, HHMISS)
   iedate, ietime,      & ! ending date and time (YYYYMMDD, HHMISS)
   loutstep,            & ! time interval of position output [s]
   loutprint,           & ! time interval of printed diagnostics
   loutsample,          & ! average is computed from samples taken every 
                          ! loutsample [s]
   itsplit,             & ! time constant for particle splitting [s]
   loutsav ,            & ! time interval of saving full fields
   lsynctime,           & ! synchronisation time interval for all particles[s]
   loffset,             & ! offset time for the first output [s]
                          ! (allows to start at any time)
   loffset2,            & ! second offset of the first output when out of
                          ! sequence with the following ones [s]
   ctl,                 & ! factor by which time step must be smaller than *
                          ! Lagrangian time scale
   ifine,               & ! reduction factor for vertical wind time step
   iout,                & ! 1 for conc.output, 2 for mixing ratio output, 3 both 
   ipout,               & ! 0 no part. dump, 1 every output time, 3 only at end
   lsubgrid,            & ! switch to turn on/off subgrid topography param.
   lagespectra,         & ! switch to turn on (1)/off (2) calc. of age spectra
   ipin,                & ! 1 continue simul. with dumped particle data, 2 no 
   iflux,               & ! switch to turn on (1)/off (2) flux calculations
   lconvection,         & ! switch to turn on (1)/off (0)convection param.
   diffus,              & ! vertical diffusivity  [unit depending on diftype]
   diftype,             & ! type of vertical diffusivity, 1:z, 2:theta
   verdiff,             & ! vertical diffusion
   hordiff,             & ! horizontal diffusion
   release_plan,        & ! type of release method, 'o3sonde', 'MOZAIC'
   vert_interpol,       & ! type of vertical interpolation in p, 'log', 'lin'
   restart,hrstart,     & ! restart run, hour labeling the restart file
   correct_vertwind,    & ! correcting vertical wind
   delay_switch_diff_off,& ! diffusion is switched off after this delay (in s)
   data_source,         & ! data source (for compatibility version, set after reading)
   ecmwf_data ,         & ! use ECMWF data
   diabatic_w,          & ! diabatic velocities are used
   ecmwf_diabatic_w,    & ! diabatic velocities from ecmwf archive
   clear_sky,           & ! clear sky version of diabatic velocities
   cloud_sky,           & ! cloud sky version of diabatic velocities (no latent heat)
   ecmwf_inct_w,        & ! temperature increment from ecmwf
   z_motion,            & ! motion is in z coordinate
   isentropic_motion,   & ! isentropic motion
   theta_bounds,        & ! theta excursion is bounded
   upper_theta_level,   & ! upper level on which theta is calculated  
   lower_theta_level,   & ! lower level on which theta is calculated
   debug_out,           & ! test prints to check diabatic versus mass velocities
   perpetual,           & ! perpetual run
   iso_mass,            & ! data on isentropic levels, mass equilibrated
   mass_diabat,         & ! use heating rates on isentropic levels
   mass_isent,          & ! isentropic motion from data on isentropic levels
   mass_correction,     & ! heating rate correction to preserve mass conservation in the stratosphere
   mean_diab_output,    & ! Output of the mass averaged heating rate
   merra_data,          & ! merra winds
   merra_diab,          & ! merra heating rates
   jra55_data,          & ! jra55 winds
   jra55_diab,          & ! jra55 heating rates
   era5_data,           & ! era5 winds
   era5_diab,           & ! era5 heating rates
   ecmwf_netcdf,        & ! ECMWF data are in netcdf4 format
   hour_accu              ! number of hours in the ECMWF accumulations
   
!Open the command file and read user options
!-------------------------------------------

 error=.false.
 print *,trim(path(1))//'COMMAND'
 open(unitcommand,file=trim(path(1))//'COMMAND',status='old',err=999)

! Check the format of the COMMAND file and call old readcommand
! if new heading with COMMAND in the first line is not detected
!--------------------------------------------------------------

 epsil=1.e-3_dp  ! must be consistent with the value in demar
 read (unitcommand,'(a)') line
 if (index(line,'COMMAND') == 0 ) then
   close(unitcommand)
   command_old=.true.
!   call readcommand(error)
   print *,' OLD COMMAND FORMAT NOW OBSOLETE'
   print *,' USE NEW FORMAT'
   error = .true.
   return
 endif

! Default values for some parameters
 ldirect = -1
 loutstep = 43200
 loutaver = 3600
 loutsample = 1800
 itsplit = BIG_INT
 loutsav = BIG_INT
 lsynctime = 900
 loutprint = 96*lsynctime
 loffset = 0
 loffset2 = 0
 ctl = 5.
 ifine = 4
 iout = 2
 ipout = 1
 lsubgrid = 2
 lagespectra = 2
 ipin = 2
 iflux = 2
 lconvection = 0
 diftype = 1
 diffus = 0.
 verdiff=.false.
 hordiff=.false.
 release_plan='o3sonde'
 vert_interpol='log'
 restart=.false.
 hrstart=0
 correct_vertwind = .true.
 delay_switch_diff_off=0
 switch_diff_off=.false.
 diabatic_w=.false.
 ecmwf_diabatic_w=.false.
 clear_sky=.false.
 cloud_sky=.false.
 ecmwf_inct_w=.false.
 z_motion=.true.
 isentropic_motion=.false.
 theta_bounds=.true.
 lower_theta_level=2
 upper_theta_level=1000
 debug_out=.false.
 perpetual=.false.
 iso_mass=.false.
 mass_diabat=.false.
 mass_isent=.false.
 mass_correction=.false.
 mean_diab_output=.false.
 merra_data=.false.
 merra_diab=.false.
 jra55_data=.false.
 jra55_diab=.false.
 era5_data=.false.
 era5_diab=.false.
 ecmwf_data=.true.
 data_source='EI'
 ecmwf_netcdf=.false.
 hour_accu = 3
  
 read(unitcommand,NML=COMMAND)
 write(6,COMMAND)
 
 ifine=max(ifine,1)

! Determine how Markov chain is formulated (for w or for w/sigw)
!---------------------------------------------------------------
 if (ctl.ge.0.1) then
   turbswitch=.true.
 else
   turbswitch=.false.
   ifine=1
 endif
 fine=1./float(ifine)
 ctl=1./ctl
 
! Fix data source
!---------------- 
 if(ecmwf_data) then
   data_source='EI'
 else if (jra55_data) then
   data_source='JRA55'
 else if (merra_data) then
   data_source='MERRA'
 else if (era5_data) then
   data_source='ERA5'
 endif
 
 select case (release_plan)
   case('StratoCim')
     parcel_dyn_alloc=.false.
   case default
     parcel_dyn_alloc=.true.
 end select
      

! Check input dates
!------------------
 error = (check_input_date_1() == 1).or.error
      
! Determine kind of dispersion method
!------------------------------------
 if (ctl > 0.) then    
   method=1
   mintime=minstep
 else
   method=0
   mintime=lsynctime
 endif

! Checks 
!    output options
!    whether synchronisation interval is sufficiently short
!    consist of gridded concentration output
!----------------------------------------------------------
 error = (check_output_options() == 1).or.error
 error = (check_synchro() == 1).or.error
 error = (check_consist() == 1).or.error

! Compute modeling time in seconds and beginning date in Julian date
!-------------------------------------------------------------------
! Notice that this estimate is wrong when 29 Febs are removed from 
! the calendar in perpetuals or CO2 runs. It follows that the number 
! of steps is slightly in excess. Hence the CO2 run encounter the 1st 
! Jan 1979 boundary and terinate in error. Pertual are just overdoing 
! by a few days.
! In both direct and indirect mode, iedate is posterior to ibdate
! In direct mode, edate (iedate) is posterior to bdate (ibdate)
! In indirect mode, edate (ibdate) is anterior to bdate (iedate)
! In direct and indirect mode, time is counted respectiveley forward 
! and backward from 0 at bdate
! ideltas is positive for a direct run and negative for an indirect run
 if (ldirect.eq.1) then
   bdate=juldate(ibdate,ibtime)
   edate=juldate(iedate,ietime)
   ideltas=nint((edate-bdate)*86400._dbl)
 else if (ldirect.eq.-1) then
   loutaver=-loutaver
   loutstep=-loutstep
   loutsample=-loutsample
   lsynctime=-lsynctime
   bdate=juldate(iedate,ietime)
   edate=juldate(ibdate,ibtime)
   ideltas=nint((edate-bdate)*86400._dbl)
 else
     write(*,*) ' #### DIRECTION" MUST BE EITHER -1 OR 1  #### '
     error=.true.
     return
 endif
 
 print *,'readcommand > direction ',ldirect
 print *,'readcommand > ',ibdate, ibtime
 print *,'readcommand > ',iedate, ietime
 print *,'readcommand > diffusion ',diftype,diffus
 print *,'readcommand > ',release_plan
 print *,'readcommand > ',vert_interpol
 print *,'readcommand > restart ',restart, hrstart
 print *,'readcommand > z_motion ',z_motion
 print *,'readcommand > isentropic_motion ',isentropic_motion
 print *,'readcommand > diabatic_w ',diabatic_w
 if (merra_data) print *,'readcommand > merra_data n _diab ',merra_data,merra_diab
 if (jra55_data) print *,'readcommand > jra55_data n _diab ',jra55_data,jra55_diab 
 if (era5_data) print *,'readcommand > era5_data n _diab ',era5_data,era5_diab 
 print *,'readcommand > ecmwf_data _diabatic_w ',ecmwf_data,ecmwf_diabatic_w 
 print *,'readcommand > correct_vertwind ',correct_vertwind
 print *,'readcommand > clear_sky ',clear_sky
 print *,'readcommand > cloud_sky ',cloud_sky
 print *,'readcommand > ecmwf_inct_w ',ecmwf_inct_w
 print *,'readcommand > theta_bounds ',theta_bounds,lower_theta_level,&
                                                   upper_theta_level
 print *,'readcommand > debug_out ',debug_out
 print *,'readcommand > iso_mass mass_diabat ',iso_mass,mass_diabat
 print *,'readcommand > mass_correction mean_diab_output ',mass_correction,mean_diab_output
 
! Conversion of format HHHMISS to seconds
!----------------------------------------
! hhh=ietime/10000
! mi=(ietime-10000*hhh)/100
! ss=ietime-10000*hhh-100*mi

! reset loffset2 if restart run
! so that one does not need to rest it manually in the auto script 
 if(restart) loffset2=0

 return    

 999  write(*,*) ' #### FLEXPART MODEL ERROR! FILE "COMMAND"    #### ' 
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      write(*,*) ' #### xxx/flexpart/options                    #### '
      print *,path(1)(1:len_path(1))//'COMMAND'
      error=.true.

 return
      
 end subroutine readcommandB
      
      function check_input_date_1 ()
      integer :: check_input_date_1
      check_input_date_1=0
      if (iedate.lt.ibdate) then
        write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING DATE    #### ' 
        write(*,*) ' #### IS LARGER THAN ENDING DATE. CHANGE      #### '
        write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
        write(*,*) ' #### "COMMAND".                              #### '
        check_input_date_1=1
      else if (iedate.eq.ibdate) then
        if (ietime.lt.ibtime) then
          write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING TIME    #### ' 
          write(*,*) ' #### IS LARGER THAN ENDING TIME. CHANGE      #### '
          write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
          write(*,*) ' #### "COMMAND".                              #### '
          check_input_date_1=1
        endif
      endif
      end function check_input_date_1

      function check_output_options()
      integer :: check_output_options
      check_output_options = 0
      if ((iout.ne.1).and.(iout.ne.2).and.(iout.ne.3)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! NEITHER CONCENTRAT#### ' 
        write(*,*) ' #### ION NOR MIXING RATIO OUTPUT SELECTED.   #### '
        check_output_options = 1
      endif
      if ((ipout.ne.0).and.(ipout.ne.1).and.(ipout.ne.2)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### ' 
        write(*,*) ' #### IPOUT MUST BE 1, 2 OR 3!                #### '
        check_output_options = 1
      endif
      if(lsubgrid.ne.1) then
        write(*,*) '             ----------------               '
        write(*,*) ' INFORMATION: SUBGRIDSCALE TERRAIN EFFECT IS'
        write(*,*) ' NOT PARAMETERISED DURING THIS SIMULATION...'
        write(*,*) '             ----------------               '
      endif
      end function check_output_options
      
      function check_synchro()
      integer :: check_synchro
      check_synchro = 0      
      if (lsynctime.gt.(idiffnorm/2)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SYNCHRONISATION   #### ' 
        write(*,*) ' #### TIME IS TOO LONG. MAKE IT SHORTER.      #### '
        check_synchro = 1
      endif
      end function check_synchro
  
      function check_consist()
      integer :: check_consist
      check_consist = 0
      if (loutaver.eq.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### ' 
        write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
        write(*,*) ' #### ZERO.                                   #### '
        write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
        check_consist = 1
      endif
      if (loutaver.gt.loutstep) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### ' 
        write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
        write(*,*) ' #### GREATER THAN INTERVAL OF OUTPUT.        #### '
        write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
        check_consist = 1
      endif
      if (loutsample.gt.loutaver) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### ' 
        write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
        write(*,*) ' #### GREATER THAN TIME AVERAGE OF OUTPUT.    #### '
        write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
        check_consist = 1
      endif
      if (mod(loutaver,lsynctime).ne.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### ' 
        write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
        write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
        check_consist = 1
      endif
      if (mod(loutstep,lsynctime).ne.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### ' 
        write(*,*) ' #### CONCENTRATION FIELDS MUST BE A MULTIPLE #### '
        write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
        check_consist = 1
      endif
      if (mod(loutsample,lsynctime).ne.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### ' 
        write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
        write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
        check_consist = 1
      endif
      if (itsplit.lt.loutaver) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SPLITTING TIME FOR#### ' 
        write(*,*) ' #### PARTICLES IS TOO SHORT. PLEASE INCREASE #### '
        write(*,*) ' #### SPLITTING TIME CONSTANT.                #### '
        check_consist = 1
      endif
      if ((lconvection.ne.0).and.(lconvection.ne.1)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! LCONVECTION MUST  #### ' 
        write(*,*) ' #### BE SET TO EITHER 1 OR 0 IN FILE COMMAND.#### ' 
        check_consist = 1
      endif
      end function check_consist

 !===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ READRELEASESB2 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

 subroutine readreleasesB2(error)

!*******************************************************************************
! 
!     This routine reads the release point specifications for the current 
!     model run. Several release points can be used at the same time.
!     meant to replace readreleasesB, and to have original readreleases in place 
!
!     Author: B. Legras
!     
!     1 June 2002
!     updated: 2004, 2005, 2007, 2010, 2012, 2013 : B. Legras, A.S. Tissier
! 
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! decay               decay constant of species                                *
! dquer [um]          mean particle diameters                                  *
! dsigma              e.g. dsigma=10 or dsigma=0.1 means that 68% of the mass  *
! npart               number of particles to be released                       *
! density [kg/m3]     density of the particles                                 *
! rm [s/m]            Mesophyll resistance                                     *
! species             name of species                                          *
! xmass               total mass of each species                               *
! xpoint1,ypoint1     geograf. coordinates of lower left corner of release area*
! xpoint2,ypoint2     geograf. coordinates of upper right corner of release are*
! weta, wetb          parameters to determine the wet scavenging coefficient   *
! zpoint1,zpoint2     height range, over which release takes place             *
!                                                                              *
!*******************************************************************************

 !use lyapunov, only: activ_Lyapunov, height_scale, sheet_slope, delta_ver, &
 !   delta_hor, tau_orthog
 use isentrop_h
 use mass_iso
 implicit none

 integer i,k,indz
 integer launch_date,launch_time,nb_level,launch_int,launch_date_end,launch_time_end
 real lat_station,long_station,p_min,p_max
 character(len=16) :: name_station
 real long_min,long_max,lat_min,lat_max,press_l,theta_l,mesh_size
 real list_l(maxpoint)
 real (dbl) :: jul1,jul2
 logical, intent(out):: error
 real (dp) :: epsil

! All the namelists which can used for specific case

 namelist /CO2/        &
    launch_date,       & ! launch date (yyyymmdd)
    launch_time,       & ! launch time (hhmmss)
    p_min, p_max,      & ! lower and upper p level (Pa above ground)
    lev1, lev2, inclev,& ! lev1, lev2, inclev 
    nb_level, n_sample,& ! number of levels to define, sammple size per level
    long_min,          & ! minimum longitude 
    long_max,          & ! maximum longitude 
    lat_min ,          & ! minimum latitude 
    lat_max ,          & ! maximum latitude 
    press_l,           & ! launch pressure
    theta_l,           & ! launch theta 
    mesh_size,         & ! size of the mesh in degree
    n_sample,          & ! sample size per point 
    uniform_spread,    & ! 
    make_curtain,      &
    make_layer,        &
    curtain_type,      &  
    make_uni3D,        &
    uniform_mesh
    
 namelist /O3SONDE/    &
    name_station,      & ! name of launching station
    long_station,      & ! longitude of launching station (deg)
    lat_station,       & ! latitude of launching station  (deg)
    p_min, p_max,      & ! lower and upper p level (Pa above ground)
    nb_level, n_sample,& ! number of levels to define, sammple size per level
    launch_date,       & ! launch date (yyyymmdd)
    launch_time          ! launch time (hhmmss)
 namelist /TTL/         & 
    long_min,           & ! minimum longitude
    long_max,           & ! maximum longitude
    lat_min ,           & ! minimum latitude
    lat_max ,           & ! maximum latitude
    press_l,            & ! launch pressure
    theta_l,            & ! launch theta 
    mesh_size,          & ! size of the mesh in degree
    launch_date,        & ! launch date (yyyymmdd)
    launch_time,        & ! launch time (hhmmss)
    n_sample,           & ! sample size per point
    tropodir,           & ! directory where to find the surf380K file
    TTLactiv              ! activation of humidity saturation calculations
    
 namelist /AGEF/        &
    long_min,           & ! minimum longitude
    long_max,           & ! maximum longitude
    lat_min ,           & ! minimum latitude
    lat_max ,           & ! maximum latitude
    nb_level ,          & ! number of levels  
    list_l,             & ! list of levels
    mesh_size_lat,      & ! lat size of the mesh in degree
    mesh_size_long,     & ! long size of the mesh in degree
    launch_date,        & ! start launch date (yyyymmdd)
    launch_time,        & ! start launch time (hhmmss)
    launch_date_end,    & ! end launch date (yyyymmdd)
    launch_time_end,    & ! end launch time (hhmmss)
    launch_int,         & ! launch interval (s)
    pcut,               & ! max pressure threshold to stop trajectories (Pa)
    plowcut,            & ! min pressure threshold to stop trajectories (Pa)
    thetacut,           & ! theta top boundary to stop trajs (K)   
    TTLFILLactiv,       & ! activate low theta cut and deactivate pcut
    thetalowcut,        & ! low theta cut
    n_sample,           & ! sample size per point
    uniform_spread,     & ! activate uniform spread
    uniform_mesh,       & ! activate the uniform mesh in longitude
    external_pos0,      & ! activate reading of external files for initial position
    external_type,      & ! type of the external initialization
    external_directory, & ! directory in which external files must be found
    shuffling             ! activate shuffling of parcels in timemanager

 namelist /AGEB/        &
    long_min,           & ! minimum longitude
    long_max,           & ! maximum longitude
    lat_min ,           & ! minimum latitude
    lat_max ,           & ! maximum latitude
    max_life_time,      & ! maximum life time of parcels (in y)
    nb_level ,          & ! number of levels  
    list_l,             & ! list of levels
    mesh_size_lat,      & ! lat size of the mesh in degree
    mesh_size_long,     & ! long size of the mesh in degree
    launch_date,        & ! start launch date (yyyymmdd)
    launch_date_end,    & ! end launch date (yyyymmdd)
    thetacut,           & ! theta top boundary to stop trajs (K)   
    n_sample,           & ! sample size per point
    uniform_spread,     & ! activate uniform spread
    uniform_mesh,       & ! activate the uniform mesh in longitude
    tropodir,           & ! directory where to find the tropopause files
    track_kill,         & ! track the killing of parcels
    track_cross,        & ! track the crossing of 1800K and 2300K surfaces
    shuffling             ! activate shuffling of parcels in timemanager
    
 namelist /AGE/         & 
    long_min,           & ! minimum longitude
    long_max,           & ! maximum longitude
    lat_min ,           & ! minimum latitude
    lat_max ,           & ! maximum latitude
    press_l,            & ! launch pressure
    theta_l,            & ! launch theta 
    mesh_size,          & ! size of the mesh in degree
    launch_date,        & ! launch date (yyyymmdd)
    launch_time,        & ! launch time (hhmmss)
    n_sample,           & ! sample size per point
    uniform_spread,     & !
    make_curtain,       &
    make_layer,         &
    curtain_type,       &
    uniform_mesh
    
 namelist /MOZAIC/      &
    MOZAIC_dir,         & ! MOZAIC directory
    MOZAIC_filename,    & ! name of the file containing flight data
    start_index,        & ! beginning of the analysed segment
    end_index,          & ! end of the analysed segment
    index_type,         & ! index type, 'rank' or 'time'
    n_sample,           & ! number of particles in each sample
    instant_release,    & ! if .true., release at observed time with partial step
    interp_release,     & ! if .true., interpolated release
    n_loc                 ! number of points (if interp_release=.true.)
    
 namelist /ER2/         &
    ER2_dir,            & ! ER2 directory
    ER2_day,            & ! flight day	
    start_ER2time,      & ! start time of the analysed segment
    end_ER2time,        & ! end time of the analysed segment
    n_sample,           & ! number of particles in each sample
    instant_release,    & ! if .true., release at observed time with partial step
    interp_release,     & ! if .true., interpolated release
    n_loc,              & ! number of points (if interp_release=.true.)
    campaign,           & ! name of the campaign
    instrument,         & ! instrument
    pcut,               & ! max pressure threshold to stop trajectories (Pa)
    plowcut,            & ! min pressure threshold to stop trajectories (Pa)
    thetacut,           & ! potential temperature killing boundary
    setxylim,           & ! determine whether the box boundaries differ from default
    correct_vertwind,   & ! should not be here as it is already in the command file, maintained for compatibility
    xlmin, xlmax, ylmin, ylmax, & ! box boundaries (in degree)
    startfrom0,         & ! define that starting is made from part_000 file in 107 format
    flight                ! flight number to remove ambiguity when several flights in a day          
    
 namelist /CLAUS/       &
    Claus_dir,          & ! CLAUS directory
    iedate_Claus, ietime_Claus, & ! ending date and time (YYYYMMDD, HHMISS) for the released parcels
    diabatic_Claus,     & ! if .true., theta levels; if .false. z levels
    TB_max,             & ! maximum brightness temperature (Kelvin) of the released parcel (top cloud)
    latmin_Claus,       & ! minimum latitude where it can exist a released parcel at the initial time
    latmax_Claus,       & ! maximum latitude where it can exist a released parcel at the initial time
    thetacut              ! potential temperature killing boundary 
  
 namelist /StratoClim/  &
    TBStratoClim_dir,   & ! StratoClim directory
    diabatic_StratoClim,& ! if .true., theta levels; if .false. z levels
    TB_max,             & ! maximum brightness temperature (Kelvin) of the released parcels
    pcut,               & ! max pressure threshold to stop trajectories (Pa)
    plowcut,            & ! min pressure threshold to stop trajectories (Pa)
    thetacut,           & ! potential temperature killing boundary
    fulldate_cloudtop,  & ! date of the first launch (integer yyyymmddhh format)
    NearRealTime,       & ! mode NearRealTime that shifts time0 by -6h or -18h
    NearRealTimeShiftHours, & ! number of hours of the shift in the near realtime case
    nb_halfdays_in_batch, &  ! number of half days to be processed in the batch sequence
    setxylim,           & ! determine whether the box boundaries deffer from default
    correct_vertwind,   & ! should not be here as it is already in the command file, maintained for compatibility
    xlmin, xlmax, ylmin, ylmax, &  ! box boundaries (in degree)
    startfrom0,         & ! define that starting is made from part_000 file in 107 format
    sigma_cut,          & ! defines the sigma value to cut the trajectory at the bottom (default 1)
    shuffling,          & ! activate shuffling of parcels in timemanager
    max_life_time         ! maximum life time in days

! Initialize logical variables
!-----------------------------
 delayed_initialization=.false.
 press2theta=.false.
 TTLactiv=.false.
 AGEFactiv=.false.
 shuffling=.false.
 external_pos0=.false.
 CLAUSactiv=.false.
 StratoClimactiv =.false.
 track_kill=.false.
 track_cross=.false.
 startfrom0 = .false.           ! only used in StratoClim plan
 sigma_cut = 1.
 
! Initialize max_life_time to a large value but not too large to avoid 
! overflow when converted to seconds
max_life_time=BIG_INT

! Initialize box boundaries
!--------------------------
 xlmin = 0._dp
 xlmax = float(nx-1)
 ylmin = 0._dp
 ylmax = float(ny-1)
 setxylim = .false.

! Initialization default index origin
!------------------------------------
 idx_orgn = 1
   
! Open the releases file and read user options
!---------------------------------------------

 error=.false.
 open(unitreleases,file=path(1)(1:len_path(1))//'RELEASES',status='old',err=999)

! Depreciated: old RELEASE format

! Check the format of the RELEASES file and call old readcommand
! if new heading with RELEASES in the first line is not detected
!--------------------------------------------------------------

! read (unitreleases,'(a)') line
! if (index(line,'RELEASES') == 0 ) then
!   close(unitreleases)
!   releases_old=.true.
!   call readreleasesB(error)
!   return
! endif

 epsil=1.e-3_dp

 select case (release_plan)

 ! Launch particles at a given time to match an                                                                                                              
!  CO2 simulation                                                                                                                                         
  case ('CO2')
    print *,'readreleases> CO2 MOTION IN AIR'
    numpoint=1
    theta_l=0.
    make_curtain=.false.
    make_layer=.false.
    make_uni3D=.true.
    uniform_spread=.true.
    uniform_mesh=.true.
    AccurateTemp=.false.
    savfull=.false.
    lev1 = 10._dp
    lev2 = 50._dp
    inclev = 1._dp
    read(unitreleases,CO2)
    write(6,CO2)
    indz=floor(lev1)
    p_min=akz(indz) + bkz(indz)*p0 &
         + (lev1-indz)*(akz(indz+1)-akz(indz)+p0*(bkz(indz+1)-bkz(indz)))
    indz=floor(lev2)    
    p_max=akz(indz) + bkz(indz)*p0 &
         + (lev2-indz)*(akz(indz+1)-akz(indz)+p0*(bkz(indz+1)-bkz(indz)))
    zpoint1(1) = p_min ; zpoint2(1) = p_max
    jul1=juldate(launch_date,launch_time)
    numlevel(1) = nb_level
    nsample(1)  = n_sample
    xpoint1(1)=long_min; xpoint2(1)=long_max
    ypoint1(1)=lat_min;  ypoint2(1)=lat_max
    mesh_size_lat=mesh_size
    mesh_size_long=mesh_size
    uppertheta=2485.30_dp
    if(iso_mass .or. diabatic_w .or. isentropic_motion) then
        delayed_initialization=.true.
        press2theta=.true.
    endif
    if(make_uni3D) delayed_initialization=.true.

 case ('AGEF')
    print *,'readreleases> AGEF'
    AGEFactiv=.true.
    TTLFILLactiv=.false.
    numpoint=1
    theta_l=0.
    pcut=25000._dp
    plowcut=0._dp
    thetacut=1800._dp
    thetalowcut=100._dp
    make_uni3D=.false.
    uniform_spread=.false.
    uniform_mesh=.false.
    AccurateTemp=.true.
    savfull=.true.
    nb_level=1
    long_min = xlon0
    long_max = xlon0 + (nx-1)*dx ! global assumed
    lat_min = -60._dp
    lat_max =  60._dp
    mesh_size_lat=2 
    mesh_size_long=2
    n_sample=1
    read(unitreleases,AGEF)
    ! if external_pos0 is true, the release parameters
    ! are not used since setpos0fromext is called
    ! instead of fixmultilayer
    write(6,AGEF)
    numpoint=nb_level
    jul1=juldate(launch_date,launch_time)
    jul2=juldate(launch_date_end,launch_time_end)
    if (nb_level > maxpoint) then
      go to 997
    else
      do k=1,nb_level
        zpoint1(k)=list_l(k)
        if(z_motion) zpoint1(k)=log(p0/zpoint1(k))
        nsample(k)  = n_sample
        xpoint1(k)=long_min; xpoint2(k)=long_max
        ypoint1(k)=lat_min;  ypoint2(k)=lat_max
        ireleasestart(k)=nint((jul1-bdate)*86400._dbl)
        ireleaseend(k)  =nint((jul2-bdate)*86400._dbl)
        ireleaseinterval(k)=launch_int
        nsample(k)  = n_sample
      enddo
    endif
    uppertheta=2485.30_dp
    
 case ('AGEB')
    print *,'readreleases> AGEB'
    AGEBactiv=.true.
    track_kill=.true.
    track_cross=.true.
    numpoint=1
    theta_l=0._dp
    thetacut=10000._dp
    max_life_time=10
    make_uni3D=.false.
    uniform_spread=.true.
    uniform_mesh=.true.
    AccurateTemp=.true.
    savfull=.false.
    nb_level=1
    long_min = xlon0
    long_max = xlon0 + (nx-1)*dx ! global assumed
    lat_min = -90._dp
    lat_max =  90._dp
    mesh_size_lat=2._dp 
    mesh_size_long=3._dp
    n_sample=1
    max_life_time = 1000
    read(unitreleases,AGEB)
    write(6,AGEB)
    ! convert max_life_time from year to s
    max_life_time=floor(365.25_dp*86400*max_life_time)
    numpoint=nb_level
    jul1=juldate(launch_date,120000)
    jul2=juldate(launch_date_end,120000)
    if (nb_level > maxpoint) then
      go to 997
    else
      do k=1,nb_level
        zpoint1(k)=list_l(k)
        !if(z_motion) zpoint1(k)=log(p0/zpoint1(k))
        !nsample(k)  = n_sample
        xpoint1(k)=long_min; xpoint2(k)=long_max
        ypoint1(k)=lat_min;  ypoint2(k)=lat_max
      enddo
    endif
    uppertheta=2485.30_dp 
    year_b=floor(launch_date/10000._dp)
    year_e=floor(launch_date_end/10000._dp)

 case ('o3sonde')
    numpoint=1
    name_station='noname'
    read(unitreleases,O3SONDE)
    xpoint1(1) = long_station ; xpoint2(1)=xpoint1(1)
    ypoint1(1) = lat_station  ; ypoint2(1)=ypoint1(1)
    zpoint1(1) = p_min ; zpoint2(1) = p_max
    jul1=juldate(launch_date,launch_time)
    ireleasestart(1) = nint((jul1-bdate)*86400._dbl)
    ireleaseend (1) = ireleasestart(1)
    numlevel(1) = nb_level
    nsample(1)  = n_sample
    npart(1) = nb_level * n_sample
    if(npart(1) > maxpart) go to 996
    error=(check_input_date_2(jul1)==1).or.error
    print *,'readreleases> O3SONDE'
    print *,'readreleases> ',name_station
    print *,'readreleases> ',launch_date,launch_time
    print *,'readreleases> ',xpoint1(1),ypoint1(1)
    print *,'readreleases> ',zpoint1(1),zpoint2(1)
    print *,'readreleases> ',numlevel(1),nsample(1)

 ! Releases particles in the TTL
 ! Type a: particles are released within a lat-long box
 ! at a fixed pressure or theta level 
    case('TTL')
    print *,'readreleases> TTL'
    numpoint=1
    n_sample=1
    long_min = xlon0
    long_max = xlon0 + (nx-2)*dx ! global assumed
    press_l=0._dp
    theta_l=0._dp
    n_loc=0
    TTLactiv=.true.
    uniform_spread=.false.
    AccurateTemp=.true.
    savfull=.true.
    read(unitreleases,TTL)
    write(6,TTL)
    xpoint1(1)=long_min; 
    xpoint2(1)=long_max
    ypoint1(1)=lat_min;  
    ypoint2(1)=lat_max
    delayed_initialization=.true.
    press2theta=.false.
    theta2press=.false.
    if(press_l > 0.) then      
      if(diabatic_w .or. isentropic_motion) then
        press2theta=.true.
      endif
      zpoint1(1)=-log(press_l/p0)
    else if (theta_l > 0.) then
      if(z_motion) then
        theta2press=.true.
      endif
      zpoint1(1)=theta_l
    else  
      go to 998
    endif
    nsample(1)=n_sample
   ! calculation of n_loc delayed to fixlayerpart
    jul1=juldate(launch_date,launch_time)
    ireleasestart(1) = nint((jul1-bdate)*86400._dbl)
    ireleaseend (1) = ireleasestart(1)
    error=(check_input_date_2(jul1)==1).or.error
    mesh_size_lat=mesh_size
    mesh_size_long=mesh_size
    print *,'readreleases> min long lat (deg) ',xpoint1(1),ypoint1(1)
    print *,'readreleases> max long lat (deg) ',xpoint2(1),ypoint2(1)

! Age of air calculation (akin to TTL)
  case('AGE')
    print *,'readreleases> AGE OF AIR'
    numpoint=1
    n_sample=1
    long_min = xlon0
    long_max = xlon0 + (nx-1)*dx ! global assumed
    lat_min = -90._dp
    lat_max =  90._dp
    press_l=0._dp
    theta_l=0._dp
    n_loc=0
    AccurateTemp=.true.
    uniform_spread=.true.
    savfull=.false.
    make_curtain=.false.
    curtain_type=1
    uniform_mesh=.false.
    make_layer=.true.
    read(unitreleases,AGE)
    write(6,AGE)
    xpoint1(1)=long_min; xpoint2(1)=long_max
    ypoint1(1)=lat_min;  ypoint2(1)=lat_max
    nsample(1)=n_sample
    if(press_l > 0.) then
      if(iso_mass .or. diabatic_w .or. isentropic_motion) then
        delayed_initialization=.true.
        press2theta=.true.
      endif
      zpoint1(1)=-log(press_l/p0)
    else if (theta_l > 0.) then
      delayed_initialization=.true.
      if(z_motion) then   
        theta2press=.true.
      endif
      zpoint1(1)=theta_l
    else
      go to 998
    endif
    nsample(1)=n_sample
    ! calculation of n_loc delayed to fixlayerpart
    jul1=juldate(launch_date,launch_time)
    ireleasestart(1) = nint((jul1-bdate)*86400._dbl)
    ireleaseend (1) = ireleasestart(1)
    error=(check_input_date_2(jul1)==1).or.error
    mesh_size_lat=mesh_size
    mesh_size_long=mesh_size
    if (make_curtain) then
      select case (curtain_type)
      case (1)  
        allocate(lat_list(15))
        lat_list = (/ -70._dp,-60._dp,-50._dp,-40._dp,-30._dp,-20._dp, &
         -10._dp,0._dp,10._dp,20._dp,30._dp,40._dp,50._dp,60._dp,70._dp/)
        lev1 = 39._dp
        lev2 = 57._dp
        inclev= 1._dp
        uppertheta=2485.30_dp
      case (2)
        print *,'readreleasesB2> WARNING: specific option will not work for all cases'
        make_layer=.false.
        allocate(lat_list(31))
        do i=1,31
          lat_list(i)=2*i-32._dp
        enddo
        lev1 = 31._dp
        lev2 = 40._dp
        inclev = 0.5_dp
        uppertheta=2485.30_dp
        numpoint = 8
        xpoint1(1:8)=long_min
        xpoint2(1:8)=long_max
        nsample(1:8)=n_sample
        zpoint1(2:8)=zpoint1(1) ! just to avoid WARNING in fixlayerpart
        launch_time=120000
        ireleasestart(1) = nint((juldate(20001231,launch_time)-bdate)*86400._dbl)
        ireleasestart(2) = nint((juldate(20001221,launch_time)-bdate)*86400._dbl)
        ireleasestart(3) = nint((juldate(20001211,launch_time)-bdate)*86400._dbl)
        ireleasestart(4) = nint((juldate(20001201,launch_time)-bdate)*86400._dbl)
        ireleasestart(5) = nint((juldate(20000731,launch_time)-bdate)*86400._dbl)
        ireleasestart(6) = nint((juldate(20000721,launch_time)-bdate)*86400._dbl)
        ireleasestart(7) = nint((juldate(20000711,launch_time)-bdate)*86400._dbl)
        ireleasestart(8) = nint((juldate(20000701,launch_time)-bdate)*86400._dbl)
      end select        
    endif 
    if(make_curtain) print *,'lev1 lev2 inclev ',lev1,lev2,inclev
    if(make_curtain) print *,'lat_list',lat_list(1),lat_list(size(lat_list)), &
             lat_list(2)-lat_list(1)
    print *,'readreleases> min long lat (deg) ',xpoint1(1),ypoint1(1)
    print *,'readreleases> max long lat (deg) ',xpoint2(1),ypoint2(1)

 ! Releases particles along a MOZAIC flight
 case('MOZAIC')
    MOZAIC_dir='/net/parker/d502/legras/MOZAIC/'
    index_type='time'
    n_sample=1
    instant_release=.false.
    interp_release =.false.
    numpart=0
    n_loc=0
    read(unitreleases,MOZAIC)
    nsample(1) = n_sample
    if(numpart > 0) n_loc=numpart ! to keep valid old RELEASE files
    numlevel(1) = n_loc           ! to perform correct output of n_loc
    print *,'readreleases> MOZAIC'
    print *,'readreleases> ',MOZAIC_filename
    print *,'readreleases> ',index_type
    print *,'readreleases> ',start_index,end_index
    print *,'readreleases> n_sample ',n_sample
    print *,'readreleases> n_loc    ',n_loc
    print *,'readreleases> instant_release ',instant_release
    print *,'readreleases> interp_release  ',interp_release
 ! Releases particles along an ER2 flight assuming
 ! standardized NASA-Ames format of meteorological files
 case('ER2')
    campaign='SOLVE'
    instrument='UCSE'
    flight = ''
    ER2_dir='/user/legras/h4/NASA/SOLVE/er2'
    n_sample=1
    instant_release=.false.
    interp_release =.true.
    AccurateTemp=.true.
    numpart=0
    n_loc=0
    savfull=.true.
    thetacut = 800._dp
    pcut=100000._dp
    plowcut=0._dp
    setxylim = .false.
    read(unitreleases,ER2)
    nsample(1) = n_sample
    if(numpart > 0) n_loc=numpart ! to keep valid old RELEASE files valid
    numlevel(1) = n_loc           ! to perform correct output of n_loc
    print *,'readreleases> ER2'
    print *,'readreleases> ',campaign
    print *,'readreleases> ',instrument, flight
    print *,'readreleases> ',ER2_day
    print *,'readreleases> ',start_ER2time,end_ER2time
    print *,'readreleases> n_sample ',n_sample
    print *,'readreleases> n_loc    ',n_loc
    print *,'readreleases> instant_release ',instant_release
    print *,'readreleases> interp_release  ',interp_release
    if (campaign == 'STC') then
      ! dynamical allocation in fixM55 
      parcel_dyn_alloc = .true.
      ! output with 107 format
      stratoClimActiv = .true.
    endif  
! case('ER2lyapou')
!    activ_Lyapunov=.true.
!    height_scale = 5000.
!    sheet_slope = 250.
!    delta_hor = 2500.
!    tau_orthog = -43200
!    campaign='SOLVE'
!    ER2_dir='/user/legras/h4/NASA/SOLVE/er2'
!    n_sample=1
!    instant_release=.false.
!    interp_release =.true.
!    numpart=0
!    n_loc=0
!    read(unitreleases,ER2)
!    call check_sample()
!    nsample(1) = n_sample
!    numlevel(1) = n_loc           ! to perform correct output of n_loc
!    delta_ver = delta_hor / sheet_slope
!    print *,'readreleases> ER2lyapou'
!    print *,'readreleases> ',ER2_day
!    print *,'readreleases> ',start_ER2time,end_ER2time
!    print *,'readreleases> n_sample ',n_sample
!    print *,'readreleases> n_loc    ',n_loc
!    print *,'readreleases> instant_release ',instant_release
!    print *,'readreleases> interp_release  ',interp_release 
! 
 ! Releases particles at the top of clouds in the tropics
 case('CLAUS')
    CLAUS_dir='/data/atissier/flexpart_in/CLAUS_TBmax/'
    iedate_Claus = 20050707
    ietime_Claus = 235959
    diabatic_Claus=.true.
    TB_max=230
    latmin_Claus=-20
    latmax_Claus=40
    thetacut=380._dp
    delayed_initialization=.true.
    AccurateTemp=.true.

    read(unitreleases,CLAUS)
    CLAUSactiv = .true.
    write(6,CLAUS)
    print *,'readreleases> CLAUS'
    print *,'readreleases> iedate_Claus ', iedate_Claus
    print *,'readreleases> ietime_Claus ', ietime_Claus
    print *,'readreleases>',Claus_dir
    print *,'readreleases> diabatic_Claus ',diabatic_Claus
    print *,'readreleases> TB_max ',TB_max
    print *,'readreleases> latmin_Claus ',latmin_Claus
    print *,'readreleases> latmax_Claus ',latmax_Claus
    print *,'readreleases> thetacut ',thetacut
 
 ! Releases parcels at the top of clouds from meteosat/himawari
 case('StratoClim')
    TBStratoClim_dir = "/data/atissier/flexpart_in/StratoClim+1km/"
    diabatic_StratoClim = .not. z_motion
    TB_max = 230
    fulldate_cloudtop=2016032900
    thetacut = 380._dp
    pcut=100000._dp
    plowcut=0._dp
    NearRealTime=.false.
    NearRealTimeShiftHours = -18
    nb_halfdays_in_batch = 0 
    M10 = .false.
    !debug_out=.true.
    
    delayed_initialization = .false.
    AccurateTemp = .true.
    ! correct_vertwind is redefined here while it belongs to the 
    ! parameters of the command file. This is not good practice. 
    ! We keep it as it is used in a number of runs with StratoClim
    ! release but this should be corrected
    correct_vertwind=.false.
    
    max_life_time = 1000
    
    read(unitreleases, StratoClim)
    !try to guess this is an M10 run
    !allows to set M10 as input when nb_halfdays_in_batch > 0 
    ! M10 should not be set on a restart run as it is not needed and 
    ! this generates an error in the io code which is looking for idxcross
    if(.not.NearRealTime .and. nb_halfdays_in_batch==0 .and. ldirect==1 .and. .not.restart) M10=.true. 

    if(.not.z_motion .and. (sigma_cut<1.)) then
      print *,'INCONSISTENT CHOICE OF SIGMA_CUT WITH DIABATIC MOTION'
      stop
    endif

    if (M10) then
      print *,'BEWARE THAT M10 RUNS CANNNOT BE RESTARTED '
    endif
    
    StratoClimactiv = .true.
    write(6,StratoClim)
    print *,'readreleases> StratoClim'
    print *,'readreleases> TBStratoClim_dir ', TBStratoClim_dir
    print *,'readreleases> fulldate_cloudtop ',fulldate_cloudtop 
    print *,'readreleases> diabatic_StratoClim ', diabatic_StratoClim
    print *,'readreleases> TB_max ',TB_max
    print *,'readreleases> thetacut ',thetacut
    print *,'readreleases> sigma_cut ',sigma_cut
    print *,'readreleases> pcut plowcut ',pcut,plowcut
    print *,'readreleases> NearRealTime ',NearRealTime
    print *,'readreleases> NearRealTimeShiftHours ',NearRealTimeShiftHours
    print *,'readreleases> M10 ',M10
    print *,'readreleases> nb_halfdays_in_batch ',nb_halfdays_in_batch
    print *,'readreleases> correct_vertwind ',correct_vertwind
    
    
    ! convert max_life_time from day to s
    max_life_time=floor(86400._dp*max_life_time)
    
 case default
   error=.true.
   write(*,*)  
   write(*,*) '#####################################################'
   write(*,*) '   FLEXPART MODEL SUBROUTINE READRELEASES2B: '
   write(*,*)
   write(*,*) 'FATAL ERROR - RELEASES PLAN NON AVAILABLE'
   write(*,*) '#####################################################'
   return
 end  select
! Common processing
 if (setxylim) then
     xlmin = (xlmin-xlon0)/dx
     xlmax = (xlmax-xlon0)/dx
     ylmin = (ylmin-ylat0)/dy
     ylmax = (ylmax-ylat0)/dy
 endif

 close(unitreleases) 
 
 call check_sav

 return

999   error=.true.
      write(*,*)  
      write(*,*) '#####################################################'
      write(*,*) '   FLEXPART MODEL SUBROUTINE READRELEASES: '
      write(*,*)
      write(*,*) 'FATAL ERROR - FILE CONTAINING PARTICLE RELEASE POINTS'
      write(*,*) 'POINTS IS NOT AVAILABLE OR YOU ARE NOT'
      write(*,*) 'PERMITTED FOR ANY ACCESS'
      write(*,*) '#####################################################'
      return

996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" EXCEEDS THE MAXIMUM      ####'
      write(*,*) '#### ALLOWED NUMBER. REDUCE EITHER NUMBER OF     ####'
      write(*,*) '#### PARTICLES PER RELEASE POINT OR REDUCE NUMBER####'
      write(*,*) '#### OF RELEASE POINTS.                          ####'
      write(*,*) '#####################################################'
      return

998   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - MUST PROVIDE A VALUE FOR A PRESSURE ####'
      write(*,*) '#### OR THETA LEVEL ON WHICH TRACER ARE          ####'
      write(*,*) '#### INITIALIZED                                 ####'
      write(*,*) '#####################################################'
      return

997   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOO MUCH LEVELS REQUIRED            ####'
      write(*,*) '#### REDUCE OR CHANGE MAXPOINT                   ####'
      write(*,*) '#####################################################'
      return
end subroutine readreleasesB2

      function check_input_date_2 (jul1)
      integer :: check_input_date_2
      real (dbl), intent(in) :: jul1
      check_input_date_2 = 0
      if (ldirect.eq.1) then
        if ((jul1.lt.bdate).or.(jul1.gt.edate)) then
          write(*,*) 'TRACZILLA MODEL ERROR'
          write(*,*) 'Release before simulation begins or '
          write(*,*) 'after simulation stops.'
          write(*,*) 'Make files COMMAND and RELEASES consistent.'
          check_input_date_2=1
          return
        endif        
      else if (ldirect.eq.-1) then
        if ((jul1.lt.edate).or.(jul1.gt.bdate)) then
          write(*,*) 'TRACZILLA MODEL ERROR'
          write(*,*) 'Release after simulation begins or '
          write(*,*) 'before simulation stops.'
          write(*,*) 'Make files COMMAND and RELEASES consistent.'
          check_input_date_2=1 
          return
        endif
      endif
      return
      end function check_input_date_2
      
      subroutine check_sample ()
      if (n_sample > 1) then
        write(*,*) 'TRACZILLA WARNING'
        write(*,*) 'n_sample forced to 1'
        n_sample = 1
      endif
      return
      end subroutine check_sample
      
      subroutine check_sav
      if(savfull .and. loutsav<BIG_INT) then
        write(*,*) 'TRACZILLA WARNING'
        write(*,*) 'savfull and loutsav both operate'
        write(*,*)'savfull forced to .false.'
        savfull=.false.
      endif
      return
      end subroutine check_sav

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXPARTICLESB @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
	
      subroutine fixparticlesB(error)
      
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles. It also attributes the mass to individual particles.          *
!                                                                              *
!     Author: A. Stohl                                                         
!                                                                              *
!     18 May 1996                                                              *
!
!     Changes: B. Legras, Apr. 2002
!              version for distribution along a profile
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]          simulation period                                       *
!                                                                              *
!*******************************************************************************

      real thelp,xaux,yhelp,zhelp
      integer idummy,i,j1,j2,ih
      logical error
      data idummy/-7/

      error=.false.
      
      numpart=0
      do i=1,numpoint
        thelp=float(ireleaseend(i)-ireleasestart(i))/float(npart(i))
        xaux=xpoint2(i)-xpoint1(i)
        yhelp=ypoint2(i)-ypoint1(i)
        zhelp=log(zpoint2(i)/zpoint1(i))/(numlevel(i)-1)
        do j1=1,numlevel(i)
          do j2=1,nsample(i)
            numpart=numpart+1
            if (numpart.gt.maxpart) goto 996
!fossil            nclass(numpart)=min(int(ran1(idummy)*float(nclassunc))+1, &
!fossil            nclassunc)
            xtra1(numpart)=xpoint1(i)
            ytra1(numpart)=ypoint1(i)
            ztra1(numpart)=-(log(zpoint1(i)/p0)+j1*zhelp)
!           print *,numpart,xtra1(numpart),ytra1(numpart),ztra1(numpart)
            itra1(numpart)=ireleasestart(i)+int((0.5_dp+float(j1-1))*thelp)
            ih=mod(itra1(numpart),abs(lsynctime)) ! only intervals of lsynctime
            itra1(numpart)=itra1(numpart)-ih    ! seconds are allowed
            !itramem(numpart)=itra1(numpart)
!           itrasplit(numpart)=itra1(numpart)+ldirect*itsplit
          enddo
        enddo
      enddo 
      print *,'fixparticlesB > completed, numpart ',numpart
      return

996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXPARTICLES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
      write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
      write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
      write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBE OF        ####'
      write(*,*) '#### RELEASE POINTS.                             ####'
      write(*,*) '#####################################################'
      return

      end subroutine fixparticlesB
      
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIX_SIMPLE_LAYER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
	
      subroutine fix_simple_layer(error)
      
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles within several pressure or potential temperature layers.
!     A simple code extracted from fixlayerpart without the complicated
!     and confusing options of that code 
!                                                                              *
!     Author: B. Legras      
!                                                                              *
!     07/07/2015               
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      real xc,yc, zc, epsil, epsil2, ylat, pint, lev 
      real corrected_mesh_size
      integer idummy,i,j,it,numparti,count_lon,jl,jy,indz,ix,jyps
      integer nb_points_latcircle
      real ddx, ddy, rddx, rddy, p1, p2, p3, p4, psurf
      logical error
      data idummy/-7/

      error = .false.
      epsil =  1.e-4 ; epsil2=1.e-7
      
      numpart=0 
      print *,'numpoint ',numpoint
      open(unitpartout,file=trim(path(2))//'release_table')
      print *,'path',path(2)
      layerloop: do i=1,numpoint
        n_loc = 0
        numparti=numpart+1
        xpoint1(i) = (xpoint1(i)-xlon0)/dx
        xpoint2(i) = (xpoint2(i)-xlon0)/dx
        ypoint1(i) = (ypoint1(i)-ylat0)/dy
        ypoint2(i) = (ypoint2(i)-ylat0)/dy
        mesh_size_lat  = mesh_size_lat/dy
        mesh_size_long = mesh_size_long/dx
        print *,'fix_simple_layer > sizes  ',mesh_size_long,mesh_size_lat
        it = ireleasestart(i) - mod(ireleasestart(i),abs(lsynctime)) 
        print *,'fix_simple_layer > xpoint ',xpoint1(i),xpoint2(i)
        print *,'fix_simple_layer > ypoint ',ypoint1(i),ypoint2(i)
        print *,'fix_simple_layer > dx, dy  ',dx,dy
        if ((xpoint1(i) < 0.).or.(xpoint1(i) > nx+epsil2)) go to 999
        if ((xpoint2(i) < 0.).or.(xpoint2(i) > nx+epsil2)) go to 999
        if ((ypoint1(i) < 0.).or.(ypoint1(i) > ny+epsil2)) go to 999
        if ((ypoint2(i) < 0.).or.(ypoint2(i) > ny+epsil2)) go to 999  

! set horizontal locations
! if horizontal spread, number of points along a contour depends on cos of lat
! otherwise it is constant over each latitude circle and latitude is the fast
! variable in the array (unlike usual conventions)
        if (uniform_spread) then     
            yc = ypoint1(i)
            do while (yc <= ypoint2(i)+epsil)
               ylat = (ylat0 + yc*dy)*pi/180.
               xc = xpoint1(i)
               count_lon=0
               if (uniform_mesh) then
                 nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
                 if(nb_points_latcircle>0) then
                   corrected_mesh_size = float(nx)/nb_points_latcircle
                 else
!                  set a large value to get only 1 pt on circle               
                   corrected_mesh_size = 2*float(nx) 
                 endif
               else
                 corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
               endif
               do while (xc <= xpoint2(i)+epsil)
                 n_loc=n_loc+1
                 count_lon=count_lon+1  
                 do j=1,nsample(i)
                   numpart = numpart+1 
                   if(numpart > maxpart) go to 996
                   xtra1(numpart) = xc
                   ytra1(numpart) = yc
                   ztra1(numpart) = zpoint1(i)
                 enddo
                 xc = xc + corrected_mesh_size
               enddo
               write(unitpartout,'(4G15.6)') ylat0+yc*dx,xlon0, &
                  corrected_mesh_size*dx,count_lon
               yc = yc+mesh_size_lat
            enddo
            print *,'fix_simple_layer > numpart after uniform spread ',numpart
        else  
          xc = xpoint1(i)
          do while (xc <= xpoint2(i)+epsil)
             yc = ypoint1(i)
             do while (yc <= ypoint2(i)+epsil)
                n_loc = n_loc+1
                do j=1,nsample(i)
                   numpart = numpart+1 
                   if(numpart > maxpart) go to 996
                   xtra1(numpart) = xc
                   ytra1(numpart) = yc
                   ztra1(numpart) = zpoint1(i)
                enddo
                yc = yc+mesh_size_lat
             enddo
             xc = xc+mesh_size_long
             !----------- temporary test
             !if((xc.lt.0).or.(xc.gt.nx-1)) then
             !   print *,'fixlayerpart > ALARM!!: xc ',xc
             !endif
             !-----------        
          enddo
          print *,'fix_simple_layer > numpart for reg grid ',numpart
        endif
! set vertical positions and initialise time
        do j=numparti,numpart          
          itra1(j) = it
          !itramem(j) = it
        enddo
      enddo layerloop
      close(unitpartout)
      print *,'fix_simple_layer > completed, numpart ',numpart
      return

996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### TRACZILLA MODEL SUBROUTINE FIXPARTICLES:    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
      write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
      write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
      write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBE OF        ####'
      write(*,*) '#### RELEASE POINTS.                             ####'
      write(*,*) '#####################################################'
      return
      
     
999   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '###  ERROR IN THE XPOINT YPOINT INITIALIZATION    ###'
      write(*,*) '#####################################################'
      return            
  
      end subroutine fix_simple_layer
    
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SET_TEMP_FOR_THETA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

subroutine set_temp_for_theta(error)
!*******************************************************************************
!
!    This routine defines the temperature after fix_simple_layer and can be used to
!    define the pressure in a z_motion run when the initialization is required
!    in theta.  
!    Should be used in delayed initialization mode when the temperature of the 
!    first field has been read.
!    An extension to be made can apply this code in AGEF runs when the temperature
!    needs to be initialized.  
!
!    Author: B. Legras
!    
!    07/07/2015
!
!*******************************************************************************
!
!    Variables:
!
!*******************************************************************************   
  
  use isentrop_h
  use readinterpN, only : ecmwf_data
  use jra55, only : jra55_data
  use merra, only : merra_data,merra_diab 
  use era5, only : era5_data
  
  logical, intent(out):: error
  real(dp), allocatable :: logzp(:),logpzp(:),mtab(:)
  real(dp), allocatable :: log_press_prof(:), log_theta_prof(:), temp_prof(:)
  real(dp):: psaver 
  integer :: errm
  integer :: i,ix,jy,nbtheta
  real(dp) :: p1,p2,p3,p4,ddx,ddy,rddx,rddy
  
  call alloc_isentrop_perm
  allocate (log_theta_prof(lower_theta_level:upper_theta_level))
  allocate (log_press_prof(lower_theta_level:upper_theta_level))
  allocate (temp_prof(lower_theta_level:upper_theta_level))
  allocate (mtab(lower_theta_level:upper_theta_level))
  nbtheta=upper_theta_level-lower_theta_level+1
  allocate (logzp(1),logpzp(1))
  error=.false.
  
  if (.not.allocated(ttra1)) then
     allocate(ttra1(maxpart))
     print *,'allocated ttra1 in set_temp_for_theta '
  endif  
  
! Loop on the particles
! This can be made more selective for enlarged applications (e.g. AGEB)
! where this routine is called several times on a selected range of parcels
! This code can be slow as the procedure is repeated for each parcel
! one per one without taking vertical distribution into account when
! there is more than one layer 
  dopart: do i=1,numpart
     if(xglobal) ix=modulo(floor(xtra1(i)),nx-1)       
      jy=max(min(floor(ytra1(i)),ny-1),-1)
      ddx=modulo(xtra1(i)-float(ix),1.)
      ! accounts for the two polar regions in JRA55
      ! the case should not occur for EI since 0 <= yt <= ny-1
      ! check change for MERRA
      if (ytra1(i)<0) then
         ddy=ytra1(i)/(90._dp/dy-ny/2)+1._dp
      else if (ytra1(i)>ny-1) then
         ddy=(ytra1(i)+1._dp-ny)/(90._dp/dy-ny/2)
      else
         ddy=ytra1(i)-float(jy)
      endif
      rddx=1.-ddx      ;  rddy=1.-ddy
      p1=rddx*rddy     ;  p2=ddx*rddy
      p3=rddx*ddy      ;  p4=ddx*ddy
      ! Calculate theta column on adjacent points
      ! Works for ERA5, JRA55 and ERA-Int
      ! For MERRA, replace by calc_col_theta_merra (CHECK)
      if (era5_data.or.jra55_data.or.ecmwf_data) then      
         call calc_col_theta(ix,jy,1)
         call calc_col_theta(ix+1,jy,1)
         call calc_col_theta(ix,jy+1,1)
         call calc_col_theta(ix+1,jy+1,1)
      else if (merra_diab) then
         print *,'ERROR set_temp_for_theta'
         error=.true.
         return  
      endif   
      ! Calculate pressure column 
      ! Works for ERA5, JRA55 and ERA-Int
      ! For MERRA, just use the fixed pressure levels
      ! psaver in Pascal
      psaver=ps(ix,jy,1,1)*(1-ddx)*(1-ddy) + ps(ix+1,jy,1,1)*(1-ddx)*ddy &
            +ps(ix,jy+1,1,1)*ddx*(1-ddy) + ps(ix+1,jy+1,1,1)*ddx*ddy
      log_theta_prof=log(theta_g(:,ix,jy,1)*(1-ddx)*(1-ddy) &
                        +theta_g(:,ix+1,jy,1)*(1-ddx)*ddy &
                        +theta_g(:,ix,jy+1,1)*ddx*(1-ddy) &
                        +theta_g(:,ix+1,jy+1,1)*ddx*ddy)
      log_press_prof=log(akz(lower_theta_level:upper_theta_level) &
                        +bkz(lower_theta_level:upper_theta_level)*psaver)
      ! Calculate temperature column
      temp_prof=tth(ix,jy,lower_theta_level:upper_theta_level,1)*(1-ddx)*(1-ddy) &
               +tth(ix+1,jy,lower_theta_level:upper_theta_level,1)*ddx*(1-ddy)   &
               +tth(ix,jy+1,lower_theta_level:upper_theta_level,1)*(1-ddx)*ddy   &
               +tth(ix+1,jy+1,lower_theta_level:upper_theta_level,1)*ddx*ddy
      ! Perform the interpolation
      ! Calculate slopes
      call slopes(log_theta_prof,log_press_prof,mtab,nbtheta)
      logzp(1)=log(ztra1(i))
      ! Find pressure for the required theta level
      call meval_spe(logzp(1:1),logpzp(1:1),& 
        log_theta_prof,log_press_prof,mtab,nbtheta,1,0.,errm)
      if (errm>0) then
         print *,'error in getting pressure from theta ',errm
         error=.true.
         return
      endif
      ! Interpolate temperature to the pressure level 
      ! (inverting pressure to get increasing values)      
      call uvip3p(3,nbtheta,-log_press_prof,temp_prof,1,-logpzp(1:1),ttra1(i:i))
      ! Set the z coordinate if pressure levels are used
      if (z_motion.and.theta2press) then
         ztra1(i)= log(p0)-logpzp(1)
      endif              
  enddo dopart

end subroutine set_temp_for_theta      

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SET_TEMP_FOR_PRESS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

subroutine set_temp_for_press(error)
!*******************************************************************************
!
!    This routine defines the temperature and can be used to
!    define the pressure in a z_motion run when the initialization is required
!    in theta.  
!    Should be used in delayed initialization mode when the temperature of the 
!    first field has been read.
!    An extension to be made can apply this code in AGEF runs when the temperature
!    needs to be initialized.  
!
!    Author: B. Legras
!    
!    07/07/2015
!
!*******************************************************************************
!
!    Variables:
!
!*******************************************************************************   
  
  use isentrop_h
  use readinterpN, only : ecmwf_data
  use jra55, only : jra55_data
  use merra, only : merra_data 
  use era5, only:era5_data
  
  logical, intent(out) :: error
  real(dp), allocatable :: logzp(:)
  real(dp), allocatable :: log_z_prof(:), temp_prof(:)
  real(dp):: psaver, errm
  integer :: i,ix,jy,nblevel
  real(dp) :: p1,p2,p3,p4,ddx,ddy,rddx,rddy
    
  allocate (log_z_prof(lower_theta_level:upper_theta_level))
  ! notice the non-orthodox usage of lower_theta_level and upper_theta_level
  ! used as generic lower_level and upper_level
  nblevel=upper_theta_level-lower_theta_level+1
  allocate (logzp(1))
  error=.false.
  
  if (.not.allocated(ttra1)) then
     allocate(ttra1(maxpart))
     print *,'allocated ttra1 in set_temp_for_press'
  endif
  
! Loop on the particles
! This can be made more selective for enlarged applications (e.g. AGEB)
! where this routine is called several times on a selected range of parcels
! This code can be slow as the procedure is repeated for each parcel
! one per one without taking vertical distribution into account when
! there is more than one layer 
  dopart: do i=1,numpart
     if(xglobal) ix=modulo(floor(xtra1(i)),nx-1)       
      jy=max(min(floor(ytra1(i)),ny-1),-1)
      ddx=modulo(xtra1(i)-float(ix),1.)
      ! accounts for the two polar regions in JRA55
      ! the case should not occur for EI since 0 <= yt <= ny-1
      ! check change for MERRA
      if (ytra1(i)<0) then
         ddy=ytra1(i)/(90._dp/dy-ny/2)+1._dp
      else if (ytra1(i)>ny-1) then
         ddy=(ytra1(i)+1._dp-ny)/(90._dp/dy-ny/2)
      else
         ddy=ytra1(i)-float(jy)
      endif
      rddx=1.-ddx      ;  rddy=1.-ddy
      p1=rddx*rddy     ;  p2=ddx*rddy
      p3=rddx*ddy      ;  p4=ddx*ddy
      
      ! Calculate pressure column 
      ! Works for ERA5, JRA55 and ERA-Int
      ! For MERRA, just use the fixed pressure levels
      ! psaver in Pascal
      psaver=ps(ix,jy,1,1)*(1-ddx)*(1-ddy) + ps(ix+1,jy,1,1)*(1-ddx)*ddy &
            +ps(ix,jy+1,1,1)*ddx*(1-ddy) + ps(ix+1,jy+1,1,1)*ddx*ddy
      log_z_prof=log(p0)-log(akz(lower_theta_level:upper_theta_level) &
                        +bkz(lower_theta_level:upper_theta_level)*psaver)
      ! Calculate temperature column
      temp_prof=tth(ix,jy,lower_theta_level:upper_theta_level,1)*(1-ddx)*(1-ddy) &
               +tth(ix+1,jy,lower_theta_level:upper_theta_level,1)*ddx*(1-ddy)   &
               +tth(ix,jy+1,lower_theta_level:upper_theta_level,1)*(1-ddx)*ddy   &
               +tth(ix+1,jy+1,lower_theta_level:upper_theta_level,1)*ddx*ddy
      ! Perform the interpolation
      ! Calculate slopes
      logzp(1)=ztra1(i)
      ! Interpolate temperature to the pressure level 
      ! (inverting pressure to get increasing values)      
      call uvip3p(3,nblevel,log_z_prof,temp_prof,1,logzp(1:1),ttra1(i:i))
      ! Set the z coordinate if pressure levels are used
      if (diabatic_w.and.press2theta) then
         ztra1(i)= ttra1(i)*exp(kappa*ztra1(i))
      endif              
  enddo dopart

end subroutine set_temp_for_press 

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXLAYERPART @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
	
      subroutine fixlayerpart(error)
      
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles within a pressure or potential temperature layer.              *
!                                                                              *
!     Author: B. Legras      
!                                                                              *
!     01 February 2004                                                         *
!
!     mod: 17/02/04: add tolerance to boundaries
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      real xc,yc, zc, epsil, epsil2, ylat, pint, lev 
      real corrected_mesh_size
      integer idummy,i,j,it,numparti,count_lon,jl,jy,indz,ix,jyps
      integer nb_points_latcircle
      real ddx, ddy, rddx, rddy, p1, p2, p3, p4, psurf
      logical error
      data idummy/-7/

      error = .false.
      epsil =  1.e-4 ; epsil2=1.e-7
      
      numpart=0 
      print *,'numpoint ',numpoint
      if(make_layer) open(unitpartout,file=trim(path(2))//'release_table')
      if(make_curtain) open(unitpartout2,file=trim(path(2))//'release_table_2')
      if(make_uni3D) open(unitpartout3,file=trim(path(2))//'release_table_3')          !!
      if(make_uni3D) open(unitpartout4,file=trim(path(2))//'release_table_4') 
      print *,'path',path(2)
      do i=1,numpoint
        n_loc = 0
        numparti=numpart+1
        xpoint1(i) = (xpoint1(i)-xlon0)/dx
        xpoint2(i) = (xpoint2(i)-xlon0)/dx
        ypoint1(i) = (ypoint1(i)-ylat0)/dy
        ypoint2(i) = (ypoint2(i)-ylat0)/dy
        mesh_size_lat  = mesh_size_lat/dy
        mesh_size_long = mesh_size_long/dx
        print *,'fixlayerpart > sizes  ',mesh_size_long,mesh_size_lat
        it = ireleasestart(i) - mod(ireleasestart(i),abs(lsynctime)) 
        print *,'fixlayerpart > xpoint ',xpoint1(i),xpoint2(i)
        print *,'fixlayerpart > ypoint ',ypoint1(i),ypoint2(i)
        print *,'fixlayerpart > dx, dy  ',dx,dy
        if ((xpoint1(i) < 0.).or.(xpoint1(i) > nx+epsil2)) go to 999
        if ((xpoint2(i) < 0.).or.(xpoint2(i) > nx+epsil2)) go to 999
        if ((ypoint1(i) < 0.).or.(ypoint1(i) > ny+epsil2)) go to 999
        if ((ypoint2(i) < 0.).or.(ypoint2(i) > ny+epsil2)) go to 999  

! set horizontal locations
! if horizontal spread, number of points along a contour depends on cos of lat
! otherwise it is constant over each latitude circle and latitude is the fast
! variable in the array (unlike usual conventions)
        if (uniform_spread) then
! The following option works for pressure, not for diabatic runs 
! TODO: reorganize into a more general 3D setting           
          if(make_uni3D) then
            if (diabatic_w) go to 997
            yc = ypoint1(i)
            do while (yc <= ypoint2(i)+epsil)
               jy=floor(yc)
               jyps=jy
!              Avoid calling psurf beyond the pole
               if(jy==ny-1)jyps=jy-1
               ylat = (ylat0 + yc*dy)*pi/180._dp
!              print *, 'fixlayerpartM > ylat ',ylat
               xc = xpoint1(i)
               count_lon=0
               if (uniform_mesh) then
                 nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
                 if(nb_points_latcircle>0) then
                   corrected_mesh_size = float(nx)/nb_points_latcircle
                 else 
    !              set a large value to get only 1 pt on circle             	 
                   corrected_mesh_size = 2*float(nx) 
                 endif
               else
                 corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
               endif
               do while (xc < xpoint2(i)+epsil)
                 n_loc=n_loc+1
                 ix=floor(xc)
                 count_lon=count_lon+1
!                Avoid to call non existent ps(.,jy+1,.,.) when at the NP
                 ddx=modulo(xc-float(ix),1._dp) ;  ddy=yc-float(jyps)
                 rddx=1.-ddx      ;  rddy=1.-ddy
                 p1=rddx*rddy     ;  p2=ddx*rddy
                 p3=rddx*ddy      ;  p4=ddx*ddy
                 psurf=p1*ps(ix,jyps,1,1)  +p2*ps(ix+1,jyps,1,1)+ &
                       p3*ps(ix,jyps+1,1,1)+p4*ps(ix+1,jyps+1,1,1)
                 lev=lev1
                 do while (lev<lev2+inclev)
                    indz=floor(lev)
                    pint=akz(indz) + bkz(indz)*psurf &
                       + (lev-indz)*(akz(indz+1)-akz(indz)+ &
                          psurf*(bkz(indz+1)-bkz(indz)))
                    zc=log(p0/pint)

                    do j=1,nsample(i)
                      numpart = numpart+1
                      if(numpart > maxpart) go to 996
                      xtra1(numpart) = xc
                      ytra1(numpart) = yc
                      ztra1(numpart) = zc
                    enddo
                 !write for keeping zc and the other parameters
                    write(unitpartout4,'(8G15.6)') ylat0+yc*dx,xlon0, &
                       corrected_mesh_size*dx,count_lon,zc,lev1,lev2,inclev
                    lev=lev+inclev
                 enddo
                 xc = xc + corrected_mesh_size
               enddo
               write(unitpartout3,'(8G15.6)') ylat0+yc*dx,xlon0, &
                  corrected_mesh_size*dx,count_lon,zc,lev1,lev2,inclev
               yc = yc+mesh_size_lat
            enddo
            print *,'fixlayerpart > numpart after make_uni3D ',numpart
          endif

          if(make_layer) then
            yc = ypoint1(i)
            do while (yc <= ypoint2(i)+epsil)
               ylat = (ylat0 + yc*dy)*pi/180.
               xc = xpoint1(i)
               count_lon=0
               if (uniform_mesh) then
                 nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
                 if(nb_points_latcircle>0) then
                   corrected_mesh_size = float(nx)/nb_points_latcircle
                 else
!                  set a large value to get only 1 pt on circle               
                   corrected_mesh_size = 2*float(nx) 
                 endif
               else
                 corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
               endif
               do while (xc <= xpoint2(i)+epsil)
                 n_loc=n_loc+1
                 count_lon=count_lon+1  
                 do j=1,nsample(i)
                   numpart = numpart+1 
                   if(numpart > maxpart) go to 996
                   xtra1(numpart) = xc
                   ytra1(numpart) = yc
                   ztra1(numpart) = zpoint1(i)
                 enddo
                 xc = xc + corrected_mesh_size
               enddo
               write(unitpartout,'(4G15.6)') ylat0+yc*dx,xlon0, &
                  corrected_mesh_size*dx,count_lon
               yc = yc+mesh_size_lat
            enddo
            print *,'fixlayerpart > numpart after make_layer ',numpart
          endif
          
!         Adding vertical curtains of parcel at selected latitudes 
!         (only available under uniform spread

! This mode also restrained to initialize pressure values      
          if (make_curtain) then
             do jl=1,size(lat_list)
               ylat = lat_list(jl)*pi/180._dp  
               yc = (lat_list(jl)-ylat0)/dy
               jy = floor(yc)
               if (uniform_mesh) then
                 nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
                 if(nb_points_latcircle>0) then
                   corrected_mesh_size = float(nx)/nb_points_latcircle
                 else
!                  set a large value to get only 1 pt on circle             	 
                   corrected_mesh_size = 2*float(nx) 
                 endif
               else
                 corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
               endif
               lev=lev1
               do while (lev<lev2+inclev)
                 indz=floor(lev)
                 pint=akz(indz) + bkz(indz)*p0 &
                     + (lev-indz)*(akz(indz+1)-akz(indz)+p0*(bkz(indz+1)-bkz(indz)))
                 zc=log(p0/pint)
                 xc=xpoint1(i)
                 count_lon=0
                 do while (xc <= xpoint2(i)+epsil)
                   count_lon=count_lon+1
                   do j=1,nsample(i)
                     numpart = numpart+1
                     if(numpart > maxpart) go to 996
                     xtra1(numpart) = xc
                     ytra1(numpart) = yc
                     ztra1(numpart) = zc
                   enddo
                   xc = xc + corrected_mesh_size
                 enddo
                 write(unitpartout2,'(5G15.6)') ylat0+yc*dx,xlon0, &
                    corrected_mesh_size*dx,count_lon,zc
                 lev=lev+inclev
               enddo
             enddo
             print *,'fixlayerpart > numpart after make_curtain ',numpart
          endif 

        else  
          xc = xpoint1(i)
          do while (xc <= xpoint2(i)+epsil)
             yc = ypoint1(i)
             do while (yc <= ypoint2(i)+epsil)
                n_loc = n_loc+1
                do j=1,nsample(i)
                   numpart = numpart+1 
                   if(numpart > maxpart) go to 996
                   xtra1(numpart) = xc
                   ytra1(numpart) = yc
                   ztra1(numpart) = zpoint1(i)
                enddo
                yc = yc+mesh_size_lat
             enddo
             xc = xc+mesh_size_long
             !----------- temporary test
             if((xc.lt.0).or.(xc.gt.nx-1)) then
                print *,'fixlayerpart > ALARM!!: xc ',xc
             endif
             !-----------        
          enddo
          print *,'fixlayerpart > numpart for reg grid ',numpart
        endif
! set vertical positions and initialise time
        do j=numparti,numpart          
          itra1(j) = it
          !itramem(j) = it
        enddo
      enddo
      if(make_layer) close(unitpartout)
      if(make_curtain) close(unitpartout2)
      if(make_uni3D) close(unitpartout3)
      print *,'fixlayerpart > completed, numpart ',numpart
      return

996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### TRACZILLA MODEL SUBROUTINE FIXPARTICLES:    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
      write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
      write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
      write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBE OF        ####'
      write(*,*) '#### RELEASE POINTS.                             ####'
      write(*,*) '#####################################################'
      return
      
997   error=.true.
      write(*,*) '### DO NOT USE MAKE_UNI3D IN FIXLAYERPART        ####'
      write(*,*) '### FOR DIABATIC RUNS                            ####'
      return      

999   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '###  ERROR IN THE XPOINT YPOINT INITIALIZATION    ###'
      write(*,*) '#####################################################'
      return
         
      end subroutine fixlayerpart
      
     
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXPRESSLAYERPART @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

 subroutine fixpresslayerpart(error)
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles within a pressure layer for the case of potential T used       *
!     as vertical coordinate
!                                                                              *
!     Author: B. Legras      
!                                                                              *
!     from fixlayerpart                                                        *
!     04/03/06
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************
 
      use mass_iso

      real xc,yc, epsil, epsil2, ylat
      real theta, thetamax,thetamin,thetainf,thetasup
      real psurf, corrected_mesh_size
      real pint,coef_interp,lev
      real ddx,ddy,rddx,rddy,p1,p2,p3,p4
      integer idummy,i,j,it,numparti,indz,ix,jy,count_lon,jl,jyps
      integer nb_points_latcircle
      logical error
      data idummy/-7/

      error = .false.
      epsil =  1.e-4 ; epsil2 = 1.e-7
      thetamin=10000._dp ; thetamax=0._dp

      if(iso_mass) then
        print *,'fixpresslayerpart does not allow to set pressure'
        print *,'in a iso mass run' 
        stop
      endif

      numpart=0 
      print *,'numpoint ',numpoint
      if(make_layer) open(unitpartout,file=trim(path(2))//'release_table')
      if(make_curtain) open(unitpartout2,file=trim(path(2))//'release_table_2')
      if(make_uni3D) open(unitpartout3,file=trim(path(2))//'release_table_3')
      if(make_uni3D) open(unitpartout4,file=trim(path(2))//'release_table_4')
      do i=1,numpoint
        n_loc = 0
        numparti=numpart+1
        xpoint1(i) = (xpoint1(i)-xlon0)/dx
        xpoint2(i) = (xpoint2(i)-xlon0)/dx
        ypoint1(i) = (ypoint1(i)-ylat0)/dy
        ypoint2(i) = (ypoint2(i)-ylat0)/dy
        mesh_size_lat  = mesh_size_lat/dy
        mesh_size_long = mesh_size_long/dx
        print *,'fixpresslayerpart > sizes  ',mesh_size_long,mesh_size_lat
        it = ireleasestart(i) - mod(ireleasestart(i),abs(lsynctime))
        print *,'fixpresslayerpart > xpoint ',xpoint1(i),xpoint2(i)
        print *,'fixpresslayerpart > ypoint ',ypoint1(i),ypoint2(i)
        print *,'fixpresslayerpart > dx, dy  ',dx,dy
        if ((xpoint1(i) < 0.).or.(xpoint1(i) > nx+epsil2)) go to 999
        if ((xpoint2(i) < 0.).or.(xpoint2(i) > nx+epsil2)) go to 999
        if ((ypoint1(i) < 0.).or.(ypoint1(i) > ny+epsil2)) go to 999
        if ((ypoint2(i) < 0.).or.(ypoint2(i) > ny+epsil2)) go to 999  

! set horizontal locations
! if horizontal spread, number of points along a contour depends on cos of lat
! otherwise it is constant over each latitude circle and latitude is the fast
! variable in the array (unlike usual conventions)
        if (uniform_spread) then
          if(make_uni3D) then
            yc = ypoint1(i)
            do while (yc <= ypoint2(i)+epsil)
               ylat = (ylat0 + yc*dy)*pi/180._dp
               jy = floor(yc)
!              avoid calling psurf beyond North pole
               jyps=jy
               if(jy==ny-1) jyps=jy-1
!              print *, 'fixpresslayerpart > ylat ',ylat
               xc = xpoint1(i)
               count_lon=0
               if (uniform_mesh) then
                 nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
                 if(nb_points_latcircle>0) then
                   corrected_mesh_size = float(nx)/nb_points_latcircle
                 else 
                   corrected_mesh_size = 2*float(nx) ! set a large value to get only 1 pt on circle
                 endif
               else
                 corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
               endif
               do while (xc <= xpoint2(i)+epsil)
                 n_loc=n_loc+1
                 count_lon=count_lon+1
                 ix = floor(xc)
                 ddx=modulo(xc-float(ix),1.) ;  ddy=yc-float(jyps)
                 rddx=1.-ddx      ;  rddy=1.-ddy
                 p1=rddx*rddy     ;  p2=ddx*rddy
                 p3=rddx*ddy      ;  p4=ddx*ddy
                 psurf=p1*ps(ix,jyps,1,1)+p2*ps(ix+1,jyps,1,1)+p3*ps(ix,jyps+1,1,1)+p4*ps(ix+1,jyps+1,1,1)
                 lev=lev1
                 do while (lev<lev2+inclev)
                    indz=floor(lev)
                    pint=akz(indz) + bkz(indz)*psurf&
                       + (lev-indz)*(akz(indz+1)-akz(indz)+psurf*(bkz(indz+1)-bkz(indz)))
                    thetasup=tth(ix,jy,indz+1,1)*(p0/(akz(indz+1) + bkz(indz+1)*psurf))**kappa
                    thetainf=tth(ix,jy,indz  ,1)*(p0/(akz(indz  ) + bkz(indz  )*psurf))**kappa
                    theta=thetainf+(lev-indz)*(thetasup-thetainf)
                    theta=min(theta,uppertheta) 
                    thetamin=min(theta,thetamin)
                    thetamax=max(theta,thetamax)
                    do j=1,nsample(i)
                      numpart = numpart+1
                      if(numpart > maxpart) go to 996
                       xtra1(numpart) = xc
                       ytra1(numpart) = yc
                       ztra1(numpart) = theta
                    enddo
                 ! write for keeping pint,and the others parameters
                 write(unitpartout4,'(9G15.6)') ylat0+yc*dx,xlon0,corrected_mesh_size*dx,count_lon,theta,lev1,lev2,inclev,pint
                 lev=lev+inclev
                 enddo
                 xc = xc + corrected_mesh_size
               enddo
               write(unitpartout3,'(8G15.6)') ylat0+yc*dx,xlon0,corrected_mesh_size*dx,count_lon,theta,lev1,lev2,inclev
               yc = yc+mesh_size_lat
            enddo
            print *,'fixpresslayerpart > numpart after make_uni3D ',numpart
          endif

          if(make_layer) then
!       find the ecmwf levels bracketting the chosen pressure level
!       we assume that initialization is done at an altitude for which
!       levels are pure pressure levels
!       If this is not the case, calculate indz inside the loop on xc, yc
!       with pint calculated
            pint=p0*exp(-zpoint1(i))
            indz=locuv(nuvz_b,nuvz,p0)
            coef_interp=(pint- akz(indz+1) - bkz(indz+1)*p0) &
                   /(akz(indz)-akz(indz+1)+p0*(bkz(indz)-bkz(indz+1)))
            print *,'fixpresslayerpart > coef_interp, ',coef_interp
            yc = ypoint1(i)
            do while (yc <= ypoint2(i)+epsil)
               ylat = (ylat0 + yc*dx)*pi/180._dp
               xc = xpoint1(i)
               count_lon=0
               if (uniform_mesh) then
                 nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
                 if(nb_points_latcircle>0) then
                   corrected_mesh_size = float(nx)/nb_points_latcircle
                 else 
                   corrected_mesh_size = 2*float(nx) ! set a large value to get only 1 pt on circle
                 endif
               else
                 corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
               endif
               do while (xc <= xpoint2(i)+epsil)
                 n_loc=n_loc+1
                 count_lon=count_lon+1
                 ix=floor(xc); jy=floor(yc)
                 thetasup=tth(ix,jy,indz+1,1)*(p0/(akz(indz+1) + bkz(indz+1)*p0))**kappa
                 thetainf=tth(ix,jy,indz  ,1)*(p0/(akz(indz  ) + bkz(indz  )*p0))**kappa
                 theta=thetasup+(thetainf-thetasup)*coef_interp
                 thetamax=max(theta,thetamax) ; thetamin=min(theta,thetamin)      
                 do j=1,nsample(i)
                   numpart = numpart+1 
                   if(numpart > maxpart) go to 996
                   xtra1(numpart) = xc
                   ytra1(numpart) = yc
                   ztra1(numpart) = theta
                 enddo
                 xc = xc + corrected_mesh_size
               enddo
               write(unitpartout,'(4G15.6)') ylat0+yc*dx,xlon0,corrected_mesh_size*dx,count_lon
               yc = yc+mesh_size_lat
            enddo
            print *,'fixpresslayerpart > numpart after make_layer ',numpart
          endif

!         Adding vertical curtains of parcel at selected latitudes (only available under uniform spread)
          if (make_curtain) then
             do jl=1,size(lat_list)
               ylat = lat_list(jl)*pi/180.   
               yc = (lat_list(jl)-ylat0)/dy
               jy = floor(yc)
               if (uniform_mesh) then
                 nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
                 if(nb_points_latcircle>0) then
                   corrected_mesh_size = float(nx)/nb_points_latcircle
                 else 
                   corrected_mesh_size = 2*float(nx) ! set a large value to get only 1 pt on circle
                 endif
               else
                 corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
               endif
               lev=lev1
               do while (lev<lev2+inclev)
                 indz=floor(lev)
                 pint=akz(indz) + bkz(indz)*p0 &
                     + (lev-indz)*(akz(indz+1)-akz(indz)+p0*(bkz(indz+1)-bkz(indz)))
                 xc=xpoint1(i)
                 count_lon=0
                 do while (xc <= xpoint2(i)+epsil)
                   count_lon=count_lon+1
                   ix = floor(xc)
                   thetasup=tth(ix,jy,indz+1,1)*(p0/(akz(indz+1) + bkz(indz+1)*p0))**kappa
                   thetainf=tth(ix,jy,indz  ,1)*(p0/(akz(indz  ) + bkz(indz  )*p0))**kappa
                   theta=thetainf+(lev-indz)*(thetasup-thetainf)
                   theta=min(theta,uppertheta) ! bounded by the max value of theta in iso_mas
                   do j=1,nsample(i)
                     numpart = numpart+1
                     if(numpart > maxpart) go to 996
                     xtra1(numpart) = xc
                     ytra1(numpart) = yc
                     ztra1(numpart) = theta
                   enddo
                   xc = xc + corrected_mesh_size
                 enddo
                 write(unitpartout2,'(5G15.6)') ylat0+yc*dx,xlon0,corrected_mesh_size*dx,count_lon,pint
                 lev=lev+inclev
               enddo
             enddo
          endif 

        else  
          pint=p0*exp(-zpoint1(i))
          indz=locuv(nuvz_b,nuvz,p0)
          coef_interp=(pint- akz(indz+1) - bkz(indz+1)*p0) &
                   /(akz(indz)-akz(indz+1)+p0*(bkz(indz)-bkz(indz+1)))
          xc = xpoint1(i) 
          print *,'fixpresslayerpart > xpoint ',xpoint1(i),xpoint2(i)
          do while (xc <= xpoint2(i)+epsil)
             yc = ypoint1(i) 
             do while (yc <= ypoint2(i)+epsil)
                n_loc = n_loc+1
                ix=int(xc); jy=int(yc)
                thetasup=tth(ix,jy,indz+1,1)*(p0/(akz(indz+1) + bkz(indz+1)*p0))**kappa
                thetainf=tth(ix,jy,indz  ,1)*(p0/(akz(indz  ) + bkz(indz  )*p0))**kappa
                theta=thetasup+(thetainf-thetasup)*coef_interp
                thetamax=max(theta,thetamax) ; thetamin=min(theta,thetamin)
                do j=1,nsample(i)
                   numpart = numpart+1 
                   if(numpart > maxpart) go to 996
                   xtra1(numpart) = xc
                   ytra1(numpart) = yc
                   ztra1(numpart) = theta
                enddo
                yc = yc+mesh_size_lat
             enddo
             !----------- temporary test
             if((xc.lt.0).or.(xc.gt.nx-1)) then
                print *,'fixpresslayerpart > ALARM!!: xc ',xc
             endif
             !-----------
             xc = xc+mesh_size_long
          enddo
        endif      
        do j=numparti,numpart
          itra1(j) = it
          !itramem(j) = it
        enddo
      enddo
      if(make_layer) close(unitpartout)
      if(make_curtain) close(unitpartout2)
      if(make_uni3D) close(unitpartout3)
      if(make_uni3D) close(unitpartout4)
      print *,'fixpresslayerpart > completed, numpart ',numpart
      print *,'fixpresslayerpart > thetamin, thetamax ',thetamin,thetamax
      return

996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXPARTICLES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
      write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
      write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
      write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBE OF        ####'
      write(*,*) '#### RELEASE POINTS.                             ####'
      write(*,*) '#####################################################'
      return

999   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '###   ERROR IN THE XPOINT YPOINT INITIALIZATION   ###'
      write(*,*) '#####################################################'
      return

!************************************
      contains
!************************************     

      function locuv(ib1,iu1,psu)
      integer, intent(in) :: ib1,iu1
      real, intent(in) :: psu
      integer :: locuv,ib,iu,im
      real :: pm
      ib=ib1 ; iu=iu1
      do while(iu-ib>1)
        im = (iu+ib)/2
        pm = akz(im) + bkz(im)*psu
        if( pm >= pint ) then
          ib=im
        else
          iu=im
        endif
      enddo
      if (bkz(ib) > 0.) then
         print *,'WARNING: initialization within sigma range'
         print *,'inaccurate conversion of theta to p'
      endif
      locuv = ib
      end function locuv

      end subroutine fixpresslayerpart
      
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXMULTILAYER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

      subroutine fixmultilayer(error)
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles a fixed number of potential temperature or pressure layers
!     and over a time interval                                                 *
!                                                                              *
!     Author: B. Legras
!     2012                                                                     *
!    
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      real xc,yc,zc, epsil,epsil2, ylat
      real corrected_mesh_size
      integer i,ix,jy,it,it_start,it_end,count_lon,count_lat,count_layer
      integer nb_points_latcircle
      logical error

      error = .false.
      epsil =  1.e-4 ; epsil2=1.e-7

      allocate (itra0(maxpart))
      if (.not.allocated(xtra1)) then
         allocate(xtra1(maxpart),ytra1(maxpart),ztra1(maxpart))
         allocate(itra1(maxpart))
      endif

      numpart=0 
      print *,'numpoint ',numpoint
      print *,'path',path(2)
      mesh_size_lat  = mesh_size_lat/dy
      mesh_size_long = mesh_size_long/dx
      print *,'fixmultilayer > meshs ',mesh_size_long,mesh_size_lat
      print *,'fixmultilayer > xlon0,ylat0,dx,dy ',xlon0,ylat0,dx,dy
!     loop on the number of layers
      dolayer: do i=1,numpoint
        print *,'fixmultilayer  layer ',i
        zc=zpoint1(i)
        xpoint1(i) = (xpoint1(i)-xlon0)/dx
        xpoint2(i) = (xpoint2(i)-xlon0)/dx
        ypoint1(i) = (ypoint1(i)-ylat0)/dy
        ypoint2(i) = (ypoint2(i)-ylat0)/dy
!       define time to be on a time-step boundary       
        it = ireleasestart(i) - mod(ireleasestart(i),abs(lsynctime))
        print *,'fixmultilayer > zc     ',zc
        print *,'fixmultilayer > xpoint ',xpoint1(i),xpoint2(i)
        print *,'fixmultilayer > ypoint ',ypoint1(i),ypoint2(i)
        if ((xpoint1(i) < 0.).or.(xpoint1(i) > nx-1+epsil2)) go to 999
        if ((xpoint2(i) < 0.).or.(xpoint2(i) > nx-1+epsil2)) go to 999
        if ((ypoint1(i) < 0.).or.(ypoint1(i) > ny-1+epsil2)) go to 999
        if ((ypoint2(i) < 0.).or.(ypoint2(i) > ny-1+epsil2)) go to 999
!       define start and end times to be on a time-step boundary       
        it_start = ireleasestart(i) - mod(ireleasestart(i),abs(lsynctime))
        it_end   = ireleaseend(i)   - mod(ireleaseend(i),  abs(lsynctime))
!       check that the release interval is a divider of end - start
        if (mod(it_start-it_end,ireleaseinterval(i)).ne.0) go to 995
!       loop on time
        dotime: do it=it_start,it_end,ireleaseinterval(i)
!         loop on latitude
          yc = ypoint1(i)
          dolat: do while (yc <= ypoint2(i)+epsil)
            jy=floor(yc)
            ylat = (ylat0 + yc*dy)*pi/180._dp
!           defined corrected mesh if needed           
            if (uniform_mesh.and.uniform_spread) then
              nb_points_latcircle = int((nx-1)*abs(cos(ylat))/mesh_size_long)
              if(nb_points_latcircle>0) then
                corrected_mesh_size = float(nx-1)/nb_points_latcircle
              else 
                corrected_mesh_size = 2*float(nx-1) ! set a large value to get only 1 pt on circle
              endif
            else if (uniform_spread) then
              corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
            else
              corrected_mesh_size = mesh_size_long
            endif
!           loop in longitude
            xc = xpoint1(i)
            dolong: do while (xc <= xpoint2(i)+epsil)
              ix=floor(xc)
              numpart=numpart+1
              if(numpart > maxpart) go to 996
              xtra1(numpart) = xc
              ytra1(numpart) = yc
              ztra1(numpart) = zc
              itra1(numpart) = it
              itra0(numpart) = it
              xc = xc + corrected_mesh_size
            enddo dolong
            yc = yc + mesh_size_lat
          enddo dolat
        enddo dotime
      enddo dolayer
      open(unitpartout4,file=trim(path(2))//'release_table_5')
!     loop to make a release table (the same for each time)
      write(unitpartout4,*)'plan: AGEF'
      write(unitpartout4,'(A,I3)')'nb_layer: ',numpoint
      dolayer2: do i=1,numpoint
        count_layer=0
        count_lat=0
        yc = ypoint1(i)
        do while (yc <= ypoint2(i)+epsil)
          count_lat=count_lat+1
          yc=yc+mesh_size_lat
        enddo
        write(unitpartout4,'(A,I3,I4)')'layer,nb_lat: ',i,count_lat
        yc = ypoint1(i)
        dolat2: do while (yc <= ypoint2(i)+epsil)
          jy=floor(yc)
          ylat = (ylat0 + yc*dy)*pi/180._dp
          if (uniform_mesh.and.uniform_spread) then
             nb_points_latcircle = int((nx-1)*abs(cos(ylat))/mesh_size_long)
             if(nb_points_latcircle>0) then
               corrected_mesh_size = float(nx-1)/nb_points_latcircle
             else 
               corrected_mesh_size = 2*float(nx-1) ! set a large value to get only 1 pt on circle
             endif
          else if (uniform_spread) then
              corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
          else
              corrected_mesh_size = mesh_size_long
          endif
          xc = xpoint1(i)
          count_lon=0
          dolong2: do while (xc <= xpoint2(i)+epsil)
            ix=floor(xc)
            count_lon=count_lon+1           
            xc = xc + corrected_mesh_size
            count_layer=count_layer+1
          enddo dolong2
          write(unitpartout4,'(3G15.6,I6)') ylat0+yc*dx,xlon0,corrected_mesh_size*dx,count_lon
          yc = yc + mesh_size_lat
        enddo dolat2
        write(unitpartout4,'(A,I3,I6)')'layer,count_layer: ',i,count_layer
      enddo dolayer2
      close(unitpartout4)
!     Save the initial positions
      open(unitpartout4,file=trim(path(2))//'sav_init',form='unformatted')
      write(unitpartout4) xtra1(1:numpart)
      write(unitpartout4) ytra1(1:numpart)
      write(unitpartout4) ztra1(1:numpart)
      write(unitpartout4) itra0(1:numpart)
      close(unitpartout4)
      print *,'fixmultilayer > completed, count_layer, numpart ',count_layer, numpart
      return

996   error=.true.
      print *,'numpart ',numpart
      print *,'xc yc ',xc,yc
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMULTILAYER:    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
      write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
      write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
      write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBE OF        ####'
      write(*,*) '#### RELEASE POINTS.                             ####'
      write(*,*) '#####################################################'
      return

999   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '###  ERROR IN THE XPOINT YPOINT INITIALIZATION    ###'
      write(*,*) '#####################################################'
      return

995   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '###  ERROR THE LAUNCH INTERVAL IS NOT A DIVIDER   ###'
      write(*,*) '###  OF THE START - END PERIOD                    ###'
      write(*,*) '#####################################################'
      return

      end subroutine fixmultilayer
      
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXMULTILAYER_TROPOMASK @@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

 subroutine fixmultilayer_tropomask(error)
!*******************************************************************************b
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles a fixed number of potential temperature layers 
!     and over a interval of time.
!     In addition this version masks particles located under the tropopause    *
!     Used in AGEB 
!                                                                              *
!     Author: B. Legras
!     25 October 2014
!     from the previous version of fixmultilayer                               *
!    
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

 real (dp) :: xc,yc,zc, epsil,epsil2, ylat, xlm
 real (dp) :: corrected_mesh_size, ddx, ddy
 real (dbl) :: datef
 real (dp), allocatable :: tropoz(:,:,:),buff(:)
 integer :: k,ix,jy,it,tt,yy,mm,dd,ll,dyr
 integer :: year_length,ldayscum(12)
 character (len=4) :: yyyy
 character (len=256) :: tropofile
 integer nb_points_latcircle
 logical error
 
 ldayscum=(/0,31,59,90,120,151,181,212,243,273,304,334/)

 error = .false.
 epsil =  1.e-3 ; epsil2=1.e-7

 allocate (itra0(maxpart))
      
 numpart=0 
 print *,'numpoint ',numpoint
 print *,'path',path(2)
 mesh_size_lat  = mesh_size_lat/dy
 mesh_size_long = mesh_size_long/dx
 print *,'fixmultilayer_tm > meshs ',mesh_size_long,mesh_size_lat
 print *,'fixmultilayer_tm > xlon0,ylat0,dx,dy ',xlon0,ylat0,dx,dy
 print *,'fixmultilayer_tm > year_b year_e ',year_b,year_e
! loops on the years, months and days
! note that time note used
 tt=0
 allocate(buff(ny*nx))
 ! calculate the size of packet_size
 ll=0
 do yy=year_b,year_e
   if(mod(yy,4)==0) then
     ll=366+ll
   else
     ll=365+ll
   endif
 enddo
 allocate (packet_len(ll))
 
! Shift of xpoint and ypoint 
 do k=1,numpoint
   xpoint1(k) = (xpoint1(k)-xlon0)/dx
   xpoint2(k) = (xpoint2(k)-xlon0)/dx
   ypoint1(k) = (ypoint1(k)-ylat0)/dy
   ypoint2(k) = (ypoint2(k)-ylat0)/dy
   if ((xpoint1(k) < 0.).or.(xpoint1(k) > nx-1+epsil2)) go to 999
   if ((xpoint2(k) < 0.).or.(xpoint2(k) > nx-1+epsil2)) go to 999
!  Deactivated because it fails for MERRA due to the poles not belonging
!  to the list of latitudes
 !  if ((ypoint1(k) < 0.).or.(ypoint1(k) > ny+epsil2)) go to 999
   if ((ypoint2(k) < 0.).or.(ypoint2(k) > ny+epsil2)) go to 999
 enddo
 print *,'ypoint ',ypoint1(1),ypoint2(1)
 print *,'ypoint ',ypoint1(numpoint),ypoint2(numpoint)
 
 doyy: do yy=year_b,year_e  
!  load the tropopause height in theta
!  should have been generated with same grid, that is with an exyra column
   write(yyyy,'(I4)')yy
   select case (data_source)
     case ('EI')
       tropofile=trim(tropodir)//'/tropo-theta-EI-'//yyyy//'.bin'
     case ('MERRA')
       tropofile=trim(tropodir)//'/tropo-theta-MERRA-'//yyyy//'.bin'
     case ('JRA55')
       tropofile=trim(tropodir)//'/tropo-theta-JRA-'//yyyy//'.bin'
     case ('ERA5')
       tropofile=trim(tropodir)//'/tropo-theta-ERA5-'//yyyy//'.bin'
   end select  
   open(tropunit,file=tropofile,access='direct',status='old',recl=4*ny*nx)
   if(mod(yy,4)==0) then
     year_length=366
   else
     year_length=365
   endif
   allocate (tropoz(0:nx-1,0:ny-1,year_length))
   do dyr=1,year_length
     read(tropunit,rec=dyr) buff
     tropoz(:,:,dyr)=reshape(buff,(/nx,ny/))
     !if((yy==year_b).and.(dyr==180)) then
     !  print *,'tropo meridional section'
     !  print *,sum(tropoz(:,:,dyr),1)/nx
     !endif
   enddo
   print *,'tropopause loaded ',yy
   domm: do mm=1,12
!    days 10 and 20 of each month are used
     dodd: do dd=10,20,10
       tt=tt+1
       datef=juldate(10000*yy+100*mm+dd,120000)
       it=nint((datef-bdate)*86400._dp)
!      define time to be on a time-step boundary
       it=it-mod(it,abs(lsynctime)) 
       dyr=dd+ldayscum(mm)
       if (mod(yy,4)==0 .and. mm>2) dyr=dyr+1
!      loop on the number of layers
       packet_len(tt)=0
       
       dolayer: do k=1,numpoint
         zc=zpoint1(k)
         if (xglobal) then
           xlm=xpoint2(k)
         else
           xlm=xpoint2(2)+epsil
         endif
         !print *,xpoint1(k),xpoint2(k)
         !print *,ypoint1(k),ypoint2(k)

         yc = ypoint1(k)
         dolat: do while (yc < ypoint2(k)+epsil)
           ! minmax for the case of the pole off the grid to avoid getting -1
           jy=min(ny-2,max(0,floor(yc)))
           ddy=yc-jy
           ylat = (ylat0 + yc*dy)*pi/180._dp
!          defined corrected longitude mesh if needed           
           if (uniform_mesh.and.uniform_spread) then
             nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
             if(nb_points_latcircle>0) then
               corrected_mesh_size = float(nx-1)/nb_points_latcircle
             else 
               corrected_mesh_size = 2*float(nx-1) ! set a large value to get only 1 pt on circle
             endif
           else if (uniform_spread) then
             corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
           else
             corrected_mesh_size = mesh_size_long
           endif
!          loop in longitude
           xc = xpoint1(k)
           dolong: do while (xc < xlm)
             ix=floor(xc)
             ddx=xc-ix
             ! test the tropopause
             ! it is assumed that the tropopause has same grid as the wind
             ! Only true for the ERA-Interim at the moment
             if (zc < 380) then
               if (zc < tropoz(ix,jy,dyr)*(1-ddx)*(1-ddy) &
                       +tropoz(ix+1,jy,dyr)*(1-ddx)*ddy &
                       +tropoz(ix,jy+1,dyr)*ddx*(1-ddy) &
                       +tropoz(ix+1,jy+1,dyr)*ddx*ddy) then
                 xc = xc + corrected_mesh_size
                 cycle dolong
               endif     
             endif
             numpart=numpart+1
             packet_len(tt)=packet_len(tt)+1
             if(numpart > maxpart) go to 996
             xtra1(numpart) = xc
             ytra1(numpart) = yc
             ztra1(numpart) = zc
             itra1(numpart) = it
             itra0(numpart) = it
             xc = xc + corrected_mesh_size
           enddo dolong
           yc = yc + mesh_size_lat
         enddo dolat
       enddo dolayer
       write(*,'("tt it yyyymmdd packet_len",I5,I12,3X,I4,2I2.2, I10)') &
             tt,it,yy,mm,dd,packet_len(tt)
     enddo dodd
   enddo domm
   deallocate(tropoz)
   close(tropunit)
 enddo doyy
 flush 6
 deallocate(buff)
 open(unitpartout4,file=trim(path(2))//'sav_init',form='unformatted')
 write(unitpartout4) xtra1(1:numpart)
 write(unitpartout4) ytra1(1:numpart)
 write(unitpartout4) ztra1(1:numpart)
 write(unitpartout4) itra0(1:numpart)
 close(unitpartout4)


 print *,'fixmultilayer_tm > completed, layers, numpart ',numpoint, numpart
 return

996   error=.true.
      print *,'numpart ',numpart
      print *,'xc yc ',xc,yc
      write(*,*) '#####################################################'
      write(*,*) '#### SUBROUTINE FIXMULTILAYER_TROPOMASK          ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
      write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
      write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
      write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBE OF        ####'
      write(*,*) '#### RELEASE POINTS.                             ####'
      write(*,*) '#####################################################'
      return

999   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '###  ERROR IN THE XPOINT YPOINT INITIALIZATION    ###'
      write(*,*) '#####################################################'
      print *,k,xpoint1(k),xpoint2(k),ypoint1(k),ypoint2(k)
      return

      end subroutine fixmultilayer_tropomask
      
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@ FIXMULTILAYER_TROPOMASK_TEMP @@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

 subroutine fixmultilayer_tropomask_temp(error)
!*******************************************************************************b
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles a fixed number of potential temperature layers.
!     In addition calculate the temperature at each launch location            *
!                                                                              *
!     Author: B. Legras
!     23 March 2015
!     from the previous version of fixmultilayer_tropomask
!     designed for jra-55, cannot be used for other reanalysis 
!     at the moment 
!     Could be easily adapted for ERA5                                                                      *
!    
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************
 use jra55
 real (dp) :: xc,yc,zc, epsil,epsil2, ylat, xlm
 real (dp) :: corrected_mesh_size, ddx, ddy
 real (dbl) :: datef
 real (dp), allocatable :: tropoz(:,:,:),buff(:)
 integer :: k,ix,jy,it,tt,yy,mm,dd,ll,dyr
 integer :: year_length,ldayscum(12)
 character (len=4) :: year
 character (len=2) :: month,day
 character (len=256) :: tropofile
 integer nb_points_latcircle
 logical error
 integer :: idxt,idxps,lat,nbtheta,iret,igrib,l, errm
 real(dp), allocatable :: values(:),val_gauss(:,:),logzp(:),logpzp(:),mtab(:)
 real(dp), allocatable :: log_press_prof(:), log_theta_prof(:), temp_prof(:)
 real(dp):: psaver,tropoxy
 
 ldayscum=(/0,31,59,90,120,151,181,212,243,273,304,334/)

 error = .false.
 epsil =  1.e-3 ; epsil2=1.e-7
 
 if (.not.allocated(itra0)) allocate (itra0(maxpart))
 if (.not.allocated(ttra1)) allocate (ttra1(maxpart))
 allocate (values((nx-1)*ny),val_gauss(nx-1,ny))
 call alloc_isentrop_perm
 allocate (log_theta_prof(lower_theta_level:upper_theta_level))
 allocate (log_press_prof(lower_theta_level:upper_theta_level))
 allocate (temp_prof(lower_theta_level:upper_theta_level))
 allocate (mtab(lower_theta_level:upper_theta_level))
 nbtheta=upper_theta_level-lower_theta_level+1
 allocate (logzp(numpoint),logpzp(numpoint))
 logzp=log(zpoint1(1:numpoint))
      
 numpart=0 
 print *,'numpoint ',numpoint
 print *,'path',path(2)
 mesh_size_lat  = mesh_size_lat/dy
 mesh_size_long = mesh_size_long/dx
 print *,'fixmultilayer_temp > meshs ',mesh_size_long,mesh_size_lat
 print *,'fixmultilayer_temp > xlon0,ylat0,dx,dy ',xlon0,ylat0,dx,dy
 print *,'fixmultilayer_temp > year_b year_e ',year_b,year_e
! loops on the years, months and days
! note that time note used
 tt=0
 allocate(buff(ny*nx))
 ! calculate the size of packet_size
 ll=0
 do yy=year_b,year_e
   if(mod(yy,4)==0) then
     ll=366+ll
   else
     ll=365+ll
   endif
 enddo
 allocate (packet_len(ll))
 
! Shift of xpoint and ypoint 
 do k=1,numpoint
   xpoint1(k) = (xpoint1(k)-xlon0)/dx
   xpoint2(k) = (xpoint2(k)-xlon0)/dx
   ypoint1(k) = (ypoint1(k)-ylat0)/dy
   ypoint2(k) = (ypoint2(k)-ylat0)/dy
   if ((xpoint1(k) < 0.).or.(xpoint1(k) > nx-1+epsil2)) go to 999
   if ((xpoint2(k) < 0.).or.(xpoint2(k) > nx-1+epsil2)) go to 999
!  Deactivated because it fails for MERRA due to the poles not belonging
!  to the list of latitudes
 !  if ((ypoint1(k) < 0.).or.(ypoint1(k) > ny+epsil2)) go to 999
   if ((ypoint2(k) < 0.).or.(ypoint2(k) > ny+epsil2)) go to 999
 enddo
 print *,'ypoint ',ypoint1(1),ypoint2(1)
 print *,'ypoint ',ypoint1(numpoint),ypoint2(numpoint)
 
 
 !#####################################
 ! TEMPORAY TEST START
! numpart=26214494
! if(numpart>0) then
!   open(unitpartout4,file=trim(path(2))//'sav_init',form='unformatted')
!   read(unitpartout4) xtra1(1:numpart)
!   read(unitpartout4) ytra1(1:numpart)
!   read(unitpartout4) ztra1(1:numpart)
!   read(unitpartout4) ttra1(1:numpart)
!   read(unitpartout4) itra0(1:numpart)
!   close(unitpartout4)
!   itra1(1:numpart)=itra0(1:numpart)
!   print *,'special restart from sav_init'
!   return
! endif
 !##################################### 
 
 doyy: do yy=year_b,year_e  
!  load the tropopause height in theta
!  should have been generated with same grid, that is with an exyra column
   write(year,'(I4)')yy
   select case (data_source)
     case ('EI')
       tropofile=trim(tropodir)//'/tropo-theta-EI-'//year//'.bin'
     case ('MERRA')
       tropofile=trim(tropodir)//'/tropo-theta-MERRA-'//year//'.bin'
     case ('JRA55')
       tropofile=trim(tropodir)//'/tropo-theta-JRA-'//year//'.bin'
     case ('ERA5')
       tropofile=trim(tropodir)//'/tropo-theta-ERA5-'//year//'.bin'
   end select
   open(tropunit,file=tropofile,access='direct',status='old',recl=4*ny*nx)
   if(mod(yy,4)==0) then
     year_length=366
   else
     year_length=365
   endif
   allocate (tropoz(0:nx-1,0:ny-1,year_length))
   do dyr=1,year_length
     read(tropunit,rec=dyr) buff
     tropoz(:,:,dyr)=reshape(buff,(/nx,ny/))
     !if((yy==year_b).and.(dyr==180)) then
     !  print *,'tropo meridional section'
     !  print *,sum(tropoz(:,:,dyr),1)/nx
     !endif
   enddo
   print *,'tropopause loaded ',yy
   domm: do mm=1,12
!  Reading section depending on JRA-55
!    read indexes here
     write(month,'(I2.2)')mm
     call grib_index_read(idxt, &
     path(3)(1:len_path(3))//year//'/'//'uvwt-'//year//'-'//month//'.gribidx');
     call grib_index_read(idxps, &
     path(3)(1:len_path(3))//year//'/'//'ps-'//year//'-'//month//'.gribidx');
!    days 10 and 20 of each month are used
     dodd: do dd=10,20,10
       write(day,'(I2.2)')dd
       tt=tt+1
       datef=juldate(10000*yy+100*mm+dd,120000)
       it=nint((datef-bdate)*86400._dp)
!      define time to be on a time-step boundary
       it=it-mod(it,abs(lsynctime)) 
       dyr=dd+ldayscum(mm)
       if (mod(yy,4)==0 .and. mm>2) dyr=dyr+1
!      loop on the number of layers
       packet_len(tt)=0
       
!      here read the temperature from the required time
       call grib_index_select(idxt,'mars.date',year//month//day)
       call grib_index_select(idxt,'mars.time','1200')
       call grib_index_select(idxt,'mars.param','11.200')
       do l=1,nuvz
         call grib_index_select(idxt,'mars.levelist',l)
         call grib_new_from_index(idxt,igrib, iret)
         call grib_get_real4_array(igrib,'values',values,iret)      
         val_gauss=reshape(values,(/nx-1,ny/))
!        Interpolation to regular grid
!        Non polar latitudes only in JRA55 (that is 0:ny-1 in lat index)
!        Polar ring of no usage here
         do lat=1,ny
           tth(0:nx-2,ny-lat,l,1) = w1(lat)*val_gauss(1:nx-1,g1(lat)) &
                                   +w2(lat)*val_gauss(1:nx-1,g2(lat))
         enddo
         if(xglobal) tth(nx-1,0:ny-1,l,1) = tth(0,0:ny-1,l,1)
       enddo
!      now read surface pressure and interpolate 
       call grib_index_select(idxps,'mars.date',year//month//day)
       call grib_index_select(idxps,'mars.time','1200')
       call grib_new_from_index(idxps,igrib, iret)
       call grib_get_real4_array(igrib,'values',values)
       val_gauss=reshape(values,(/nx-1,ny/))      
       do lat=1,ny
         ps(0:nx-2,ny-lat,1,1) = w1(lat)*val_gauss(1:nx-1,g1(lat)) &
                                +w2(lat)*val_gauss(1:nx-1,g2(lat))
       enddo
       if(xglobal) ps(nx-1,0:ny-1,1,1) = ps(0,0:ny-1,1,1)
!      Assume that xpoint and ypoint do not depend on the level
!      End of JRA-55 reading section
 
       if (xglobal) then
         xlm=xpoint2(1)
       else
         xlm=xpoint2(1)+epsil
       endif

       yc = ypoint1(1)
       dolat: do while (yc < ypoint2(1)+epsil)
           ! minmax for the case of the pole off the grid to avoid getting -1
           jy=min(ny-2,max(0,floor(yc)))
           ddy=yc-jy
           ylat = (ylat0 + yc*dy)*pi/180._dp
!          defined corrected longitude mesh if needed           
           if (uniform_mesh.and.uniform_spread) then
             nb_points_latcircle = int(nx*abs(cos(ylat))/mesh_size_long)
             if(nb_points_latcircle>0) then
               corrected_mesh_size = float(nx-1)/nb_points_latcircle
             else 
               corrected_mesh_size = 2*float(nx-1) ! set a large value to get only 1 pt on circle
             endif
           else if (uniform_spread) then
             corrected_mesh_size = mesh_size_long/(abs(cos(ylat))+epsil2)
           else
             corrected_mesh_size = mesh_size_long
           endif
!          loop in longitude
           xc = xpoint1(1)
           dolong: do while (xc < xlm)
             ix=floor(xc)
             ddx=xc-ix
             ! test the tropopause
             ! it is assumed that the tropopause has same grid as the wind
             ! make sure to read the file that corresponds to the reanalysis in use
             tropoxy  = tropoz(ix,jy,dyr)*(1-ddx)*(1-ddy) &
                       +tropoz(ix+1,jy,dyr)*(1-ddx)*ddy &
                       +tropoz(ix,jy+1,dyr)*ddx*(1-ddy) &
                       +tropoz(ix+1,jy+1,dyr)*ddx*ddy
             ! calculate theta column on adjacent points
             call calc_col_theta(ix,jy,1)
             call calc_col_theta(ix+1,jy,1)
             call calc_col_theta(ix,jy+1,1)
             call calc_col_theta(ix+1,jy+1,1)
             ! psaver in Pascal
             psaver=ps(ix,jy,1,1)*(1-ddx)*(1-ddy) + ps(ix+1,jy,1,1)*(1-ddx)*ddy &
                   +ps(ix,jy+1,1,1)*ddx*(1-ddy) + ps(ix+1,jy+1,1,1)*ddx*ddy
             log_theta_prof=log(theta_g(:,ix,jy,1)*(1-ddx)*(1-ddy) &
                               +theta_g(:,ix+1,jy,1)*(1-ddx)*ddy &
                               +theta_g(:,ix,jy+1,1)*ddx*(1-ddy) &
                               +theta_g(:,ix+1,jy+1,1)*ddx*ddy)
             log_press_prof=log(akz(lower_theta_level:upper_theta_level) &
                           +bkz(lower_theta_level:upper_theta_level)*psaver)
             temp_prof=tth(ix,jy,lower_theta_level:upper_theta_level,1)*(1-ddx)*(1-ddy) &
                      +tth(ix+1,jy,lower_theta_level:upper_theta_level,1)*ddx*(1-ddy)   &
                      +tth(ix,jy+1,lower_theta_level:upper_theta_level,1)*(1-ddx)*ddy   &
                      +tth(ix+1,jy+1,lower_theta_level:upper_theta_level,1)*ddx*ddy
             k=1
             do while (zpoint1(k)<tropoxy)
                k=k+1
             enddo
             packet_len(tt)=packet_len(tt)+numpoint-k+1
             xtra1(numpart+1:numpart+numpoint-k+1)=xc
             ytra1(numpart+1:numpart+numpoint-k+1)=yc
             if (diabatic_w) then
                ztra1(numpart+1:numpart+numpoint-k+1)=zpoint1(k:numpoint)
             endif
             itra1(numpart+1:numpart+numpoint-k+1)=it
             itra0(numpart+1:numpart+numpoint-k+1)=it
             ! check that the interpolation range is Ok
             if (zpoint1(k) < exp(log_theta_prof(lower_theta_level))) then
                print *,'lower theta exceeded'
                print *,k,zpoint1(k),exp(log_theta_prof(lower_theta_level)),psaver
                print *,ix,jy,xlon0+ix*dx,ylat0+jy*dy
             endif
             if (zpoint1(numpoint) > exp(log_theta_prof(upper_theta_level))) then
                print *,'upper theta exceeded'
                print *,k,zpoint1(numpoint),exp(log_theta_prof(upper_theta_level)),psaver
                print *,ix,jy,xlon0+ix*dx,ylat0+jy*dy
             endif
             call slopes(log_theta_prof,log_press_prof,mtab,nbtheta)
             call meval_spe(logzp(k:numpoint),logpzp(k:numpoint),& 
                 log_theta_prof,log_press_prof,mtab,nbtheta,numpoint-k+1,0.,errm)
             if (errm>0) then
               print *,'error in getting pressure from theta ',errm
               error=.true.
             endif
             ! failing uvip3p due to large variations in the gradient of theta
             ! replaced by slopes/meval, see interpol in sandbox
             !call uvip3p(1,nbtheta,log_theta_prof,log_press_prof, &
             !           numpoint-k+1,logzp(k:numpoint),logpzp(k:numpoint))
             call uvip3p(3,nbtheta,-log_press_prof,temp_prof, & 
                        numpoint-k+1,-logpzp(k:numpoint),&
                        ttra1(numpart+1:numpart+numpoint-k+1))
             ! check interpolation errors generating negative temperatures
             if (z_motion) then
                ztra1(numpart+1:numpart+numpoint-k+1)= log(p0)-logpzp(k:numpoint)
             endif
             if(minval(ttra1(numpart+1:numpart+numpoint-k+1))<=0) then
               print *,'ttra1 ',ttra1(numpart+1:numpart+numpoint-k+1)
               print *,'press_prof ',exp(log_press_prof)
               print *,'temp_prof ',temp_prof
               print *,'pzp ',exp(logpzp(k:numpoint))
               print *,'theta_prof ',exp(log_theta_prof)
               print *,'zp ',exp(logzp(k:numpoint))
             endif
             numpart=numpart+numpoint-k+1
             if(numpart>maxpart) go to 996      
             xc = xc + corrected_mesh_size
           enddo dolong
           yc = yc + mesh_size_lat
         enddo dolat
       write(*,'("tt it yyyymmdd packet_len",I5,I12,3X,I4,2I2.2, I10)') &
             tt,it,yy,mm,dd,packet_len(tt)
     enddo dodd
     call grib_index_release(idxt)
     call grib_index_release(idxps)
   enddo domm
   deallocate(tropoz)
   close(tropunit)
   flush 6
 enddo doyy
 deallocate(buff,values,val_gauss)
 open(unitpartout4,file=trim(path(2))//'sav_init',form='unformatted')
 write(unitpartout4) xtra1(1:numpart)
 write(unitpartout4) ytra1(1:numpart)
 write(unitpartout4) ztra1(1:numpart)
 write(unitpartout4) ttra1(1:numpart)
 write(unitpartout4) itra0(1:numpart)
 close(unitpartout4)
 call savsav(0)
 
 print *,'fixmultilayer_temp > completed, layers, numpart ',numpoint, numpart
 return

996   error=.true.
      print *,'numpart ',numpart
      print *,'xc yc ',xc,yc
      write(*,*) '#####################################################'
      write(*,*) '#### SUBROUTINE FIXMULTILAYER_TROPOMASK_TEMP     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
      write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
      write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
      write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
      write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBER OF       ####'
      write(*,*) '#### RELEASE POINTS.                             ####'
      write(*,*) '#####################################################'
      return
      
999   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '###  ERROR IN THE XPOINT YPOINT INITIALIZATION    ###'
      write(*,*) '#####################################################'
      print *,k,xpoint1(k),xpoint2(k),ypoint1(k),ypoint2(k)
      return

      end subroutine fixmultilayer_tropomask_temp
      
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SETPOS0FROMEXT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine setpos0fromext
!*******************************************************************************
!
!    This routine defines the release times and locations from external files
!    as a replacement of fixmultilayer
!
!    First version: B. Legras, Aug. 2014
!                   Quick-'n-dirty version to be parameterized in a clean way
!                   in future applications.
!                   Turn this mess into a clean setup reading a netcdf file.
!                   Presently read the files generated by the matlab script
!                   external_pos0_generator.m in AGEF-annex using tropopause 
!                   mat files.
!
!*******************************************************************************
  use ecmwf_diab
 
! Files
  character(len=256) :: sizes_file,itra_file,ztra_file,xytra_file
! Experiment type
  character(len=8) :: exp_type
! Local sizes
  integer :: packet_size, count_t
! Local itra, xtra, ytra, ztra
  integer(kind=4), allocatable, dimension(:) :: itra_l
  real(sp), allocatable, dimension(:) :: xtra_l, ytra_l
  real(sp), allocatable, dimension(:,:) :: ztra_l
! do index
  integer :: i
  

! Determination of exp type from COMMAND parameters
  if (diabatic_w) then
    exp_type='DIA'
  else if (z_motion) then
    exp_type='Z'
  else
    exp_type=''
  endif

! process differently according the the type of initialization
  select case (external_type)
  
! First case: initialisation of the tropause on a fixed rectangular grid
  case ('TROPOPAUSE')
  
!   File names
    sizes_file=trim(external_directory)//'extern_pos0_'//trim(exp_type)//'_sizes.dat'
    itra_file =trim(external_directory)//'extern_pos0_'//trim(exp_type)//'_itra.dat'
    ztra_file =trim(external_directory)//'extern_pos0_'//trim(exp_type)//'_ztra.dat'
    xytra_file=trim(external_directory)//'extern_pos0_'//trim(exp_type)//'_xytra.dat'
  
!   Read what is needed from the size file
    open(701,FILE=sizes_file,STATUS='OLD',ACCESS='DIRECT',RECL=4)
    read(701,REC=1) packet_size
    read(701,REC=2) count_t
    close(701)
 
!   Check the dimension against maxpart
    if (packet_size*count_t > maxpart) then
       print *,'[ACHTUNG ALARM] setpos0from ext'
       print *,'packet_size,count_t,maxpart ',packet_size,count_t,maxpart
       stop 750
    endif
  
!   Read the itra file
    allocate (itra_l(count_t))
    open(702,FILE=itra_file,STATUS='OLD',ACCESS='DIRECT',RECL=4*count_t)
    read(702,REC=1) itra_l
    close(702)
  
!   Read the xytra file
    allocate (xtra_l(packet_size),ytra_l(packet_size))
    open(703,FILE=xytra_file,STATUS='OLD',ACCESS='DIRECT',RECL=4*packet_size)
    read(703,REC=1) xtra_l
    read(703,REC=2) ytra_l
    close(703)
  
!   Read the ztra file
    allocate (ztra_l(packet_size,count_t))
    open(704,FILE=ztra_file,STATUS='OLD',ACCESS='DIRECT',RECL=4*packet_size)
    do i=1,count_t
       read(704,REC=i) ztra_l(:,i)
    enddo
    close(704)
  
!   Print a few lines
    print *,'setpos0fromext > packet_size, count_t ',packet_size,count_t
    print *,'setpos0fromext > ',maxpart,packet_size*count_t
    print *,'setpos0fromext > min, max x ',minval(xtra_l),maxval(xtra_l)
    print *,'setpos0fromext > min, max y ',minval(ytra_l),maxval(ytra_l)
    print *,'setpos0fromext > min, max z ',minval(ztra_l),maxval(ztra_l)
  
!   Proceed with the initialization

!   Allocate itra0
    allocate(itra0(maxpart))
  
!   Copy files
    numpart=0
    do i=1,count_t
      xtra1(numpart+1:numpart+packet_size)=xtra_l
      ytra1(numpart+1:numpart+packet_size)=ytra_l
      ztra1(numpart+1:numpart+packet_size)=ztra_l(:,i)
      itra0(numpart+1:numpart+packet_size)=itra_l(i)
      itra1(numpart+1:numpart+packet_size)=itra_l(i)
      numpart=numpart+packet_size 
    enddo
    print *,'setpos0fromext > numpart ',numpart
    flush 6
    deallocate (itra_l,xtra_l,ytra_l,ztra_l)
    
  case default
    print *,'This case of external input is not valid ',external_type
    stop
  
  end select
  
! Diagnostic output
! release_table_5 not done
  open(unitpartout4,file=trim(path(2))//'sav_init',form='unformatted')
  write(unitpartout4) xtra1(1:numpart)
  write(unitpartout4) ytra1(1:numpart)
  write(unitpartout4) ztra1(1:numpart)
  write(unitpartout4) itra0(1:numpart)
  close(unitpartout4)
  
  print *,'setpos0fromext completed'

  
  end subroutine setpos0fromext

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXER2INTERP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine fixER2interp(error)
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles.           						       *
!     Interpolate selected part of the flight path to n_loc points regularly
!     distributed in time. Put n_sample particles at each point.
!
!     Changes: B. Legras, Nov. 2002
!              version for distribution along an aircraft transect 
!              with temporal interpolation
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]          simulation period                                       *
!                                                                              *
!*******************************************************************************

! use lyapunov, only: numpart_np, activ_Lyapunov
       
 real, allocatable, dimension(:) :: timeO3, UTGPS, GMTMM,  &
      off_jul,interp_time, &
      interp_xtra1, interp_ytra1, interp_ztra1, &
      F_z, F_lat, F_long
 
 integer, allocatable, dimension(:) :: O3, &
      DaltGPS, DlatGPS, DLongGPS, QGPS, TedrGPS, ReynGPS, &
      PStaMM, TstaMM, ThtaMM, UMM, VMM, WMM
 
 character(len=96):: line
 logical, intent(out):: error
 logical do_rel,stop_loop
 integer lhead, FDate, YYYY, MM, DD
 integer i,j,j1,nobs,nO3obs,nGPSobs,nMMobs
 real diff_time 
 
 error=.false.
 open(unitER2_O3,file=trim(ER2_dir)//'O3'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 open(unitER2_MG,file=trim(ER2_dir)//'MG'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 open(unitER2_MM,file=trim(ER2_dir)//'MM'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 read(unitER2_O3,*) lhead
 do i=2,lhead
   read(unitER2_O3,*)
 enddo   
 read(unitER2_MG,*) lhead
 do i=2,lhead
   read(unitER2_MG,*)
 enddo
 read(unitER2_MM,*) lhead
 do i=2,6
   read(unitER2_MM,*)
 enddo
 read(unitER2_MM,*) YYYY, MM, DD
 do i=8,lhead
   read(unitER2_MM,*)
 enddo 
 FDate=10000*YYYY + 100*MM + DD
 diff_time = end_ER2time - start_ER2time
 nobs=int(1.05_dp*diff_time)
 allocate (timeO3(nobs),O3(nobs),UTGPS(nobs),GMTMM(nobs))
 allocate (DaltGPS(nobs),DlatGPS(nobs),DlongGPS(nobs))
 allocate (QGPS(nobs),TedrGPS(nobs),ReynGPS(nobs))
 allocate (PStaMM(nobs),TstaMM(nobs),ThtaMM(nobs),UMM(nobs))
 allocate (VMM(nobs),WMM(nobs))

 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_o3,*,err=997,end=120) timeO3(j), O3(j)
   if(.not.do_rel) then
     if((timeO3(j)>(start_ER2time-0.9999_dp)) &
     .and.(timeO3(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((timeO3(j)<(end_ER2time+0.9999_dp)) &
     .and.(timeO3(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nO3obs=j
 print *,'fixER2interp > O3 obs read ',nO3obs
 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_MG,*,err=997,end=120) UTGPS(j), &
     DaltGPS(j), DlatGPS(j), DLongGPS(j), QGPS(j), &
     TedrGPS(j), ReynGPS(j)
   if(.not.do_rel) then
     if((UTGPS(j)>(start_ER2time-0.9999_dp)) &
     .and.(UTGPS(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((UTGPS(j)<(end_ER2time+0.9999_dp)) &
     .and.(UTGPS(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nGPSobs=j
 print *,'fixER2interp > GPS obs read ',nGPSobs
 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_MM,*,err=997,end=120) GMTMM(j), &
     PStaMM(j), TstaMM(j), ThtaMM(j), UMM(j), VMM(j), WMM(j)
   if(.not.do_rel) then
     if((GMTMM(j)>(start_ER2time-0.9999_dp)) &
     .and.(GMTMM(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((GMTMM(j)<(end_ER2time+0.9999_dp)) &
     .and.(GMTMM(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nMMobs=j
 print *,'fixER2interp > Meteo obs read ',nMMobs

 allocate (F_z(nMMobs),F_lat(nGPSobs),F_long(nGPSobs))  
 do j=1,nMMobs
!   julian_date(j)=juldate(F_date(j),F_time(j))
   F_z(j) = -log(10.*real(PStaMM(j))/p0) ! p0 in Pascal
 enddo
 do j=1,nGPSobs
   F_lat(j)  = real(DlatGPS(j)) *1.e-5
   F_long(j) = real(DlongGPS(j))*1.e-5
 enddo
 
! if( activ_Lyapunov ) then
!   numpart_np = n_loc
!   numpart = numpart_np * 7
! else
   numpart = n_loc * n_sample
! endif
 error = (check_numpart() == 1).or.error
 
!  Interpolation
 print *,'fixER2interp > interpolating and filling positions'
 allocate(interp_time(n_loc))
 allocate(interp_xtra1(n_loc),interp_ytra1(n_loc),interp_ztra1(n_loc))
 do i=1,n_loc
   interp_time(i) = Start_ER2time+(i-1)*(End_ER2time-Start_ER2time)/(n_loc-1)
 enddo
! print *,julian_date(1),julian_date(2),julian_date(3)
! print *,off_jul(1),off_jul(nobs)
! write(*,'(a,2g12.5)') 'fixER2interp> interp_time ',interp_time(1),interp_time(n_loc)
 call uvip3p(3,nGPSobs,UTGPS,F_long,n_loc,interp_time,interp_xtra1)
 call uvip3p(3,nGPSobs,UTGPS,F_lat ,n_loc,interp_time,interp_ytra1)
 call uvip3p(3,nMMobs, GMTMM,F_z   ,n_loc,interp_time,interp_ztra1)
 j1=1   
 do i=1,n_loc
   itra1(j1) = nint((juldate(FDate,000000)-bdate)*86400._dbl + interp_time(i))
   error = (check_launch_time(i,j1) == 1).or.error
   do j=j1,j1+n_sample-1
     itra1(j) = nint((juldate(FDate,000000)-bdate)*86400._dbl + interp_time(i))
     xtra1(j) = (interp_xtra1(i) - xlon0)/dx
     ytra1(j) = (interp_ytra1(i) - ylat0)/dy
     ztra1(j) = interp_ztra1(i)
     !itramem(j) = itra1(j)
   enddo
   j1=j1+n_sample
 enddo
 
 print *,'fixER2interp> release end'
 print *,'fixER2interp> ',numpart
 write(*,'(a,2g12.5)') ' fixER2interp> itra1 ',itra1(1),itra1(numpart)
 write(*,'(a,2g12.5)') ' fixER2interp> xtra1 ',xtra1(1),xtra1(numpart)
 write(*,'(a,2g12.5)') ' fixER2interp> ytra1 ',ytra1(1),ytra1(numpart)  
 write(*,'(a,2g12.5)') ' fixER2interp> ztra1 ',ztra1(1),ztra1(numpart) 
 
 close(unitER2_O3)
 close(unitER2_MG)
 close(unitER2_MM)
 deallocate(timeO3, UTGPS, GMTMM,interp_time,F_z, F_lat, F_long, &
   O3, DaltGPS, DlatGPS, DLongGPS, QGPS, TedrGPS, ReynGPS, &
   PStaMM, TstaMM, ThtaMM, UMM, VMM, WMM, &
   interp_ztra1,interp_xtra1,interp_ytra1)
 return
                   

120   print *, 'fixER2interp> unexpected EOF'
      error=.true.
      print *, 'special diag ',j,GMTMM(j-1),GMTMM(j-2)
      return

995   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2INTERP :    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - NOBS TOO SMALL                      ####'
      write(*,*) '#####################################################'      
      return
996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2INTERP :    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
      write(*,*) '#####################################################'      
      return
997   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2INTERP:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR WHILE READING THE FILE                ####'
      write(*,*) '#####################################################'      
      return     

      return    
      
  end subroutine fixER2interp
  
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXER2SOLVEINTERPO @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine fixER2SOLVEinterp(error)
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles.           						       *
!     Interpolate selected part of the flight path to n_loc points regularly
!     distributed in time. Put n_sample particles at each point.
!
!     Changes: B. Legras, Nov. 2002
!              version for distribution along an aircraft transect 
!              with temporal interpolation
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]          simulation period                                       *
!                                                                              *
!*******************************************************************************

! use lyapunov, only: numpart_np, activ_Lyapunov
 use isentrop_m, only: isentropic_motion
 use ecmwf_diab
       
 real, allocatable, dimension(:) :: timeO3, UTGPS, GMTMM,  &
      interp_time, &
      interp_xtra1, interp_ytra1, interp_ztra1, &
      F_z, F_lat, F_long
 
 integer, allocatable, dimension(:) :: O3, &
      DaltGPS, DlatGPS, DLongGPS, QGPS, TedrGPS, ReynGPS, &
      PStaMM, TstaMM, ThtaMM, UMM, VMM, WMM
 
 logical, intent(out):: error
 logical do_rel,stop_loop
 integer lhead, FDate, YYYY, MM, DD
 integer i,j,j1,nobs,nO3obs,nGPSobs,nMMobs
 real diff_time 
 
 error=.false.
 print * ,trim(ER2_dir)//'O3'//trim(ER2_day)//'.ER2'
 open(unitER2_O3,file=trim(ER2_dir)//'O3'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 print * ,trim(ER2_dir)//'MG'//trim(ER2_day)//'.ER2'
 open(unitER2_MG,file=trim(ER2_dir)//'MG'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 print * ,trim(ER2_dir)//'MM'//trim(ER2_day)//'.ER2'
 open(unitER2_MM,file=trim(ER2_dir)//'MM'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 read(unitER2_O3,*) lhead
 do i=2,lhead
   read(unitER2_O3,*)
 enddo   
 read(unitER2_MG,*) lhead
 do i=2,lhead
   read(unitER2_MG,*)
 enddo
 read(unitER2_MM,*) lhead
 do i=2,6
   read(unitER2_MM,*)
 enddo
 read(unitER2_MM,*) YYYY, MM, DD
 do i=8,lhead
   read(unitER2_MM,*)
 enddo 
 FDate=10000*YYYY + 100*MM + DD
 diff_time = end_ER2time - start_ER2time
 nobs=int(1.05_dp*diff_time) ! rough estimate for less than 1Hz signal
 allocate (timeO3(nobs),O3(nobs),UTGPS(nobs),GMTMM(nobs))
 allocate (DaltGPS(nobs),DlatGPS(nobs),DlongGPS(nobs))
 allocate (QGPS(nobs),TedrGPS(nobs),ReynGPS(nobs))
 allocate (PStaMM(nobs),TstaMM(nobs),ThtaMM(nobs),UMM(nobs))
 allocate (VMM(nobs),WMM(nobs))

 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_o3,*,err=997,end=120) timeO3(j), O3(j)
   if(.not.do_rel) then
     if((timeO3(j)>(start_ER2time-0.9999_dp)) &
     .and.(timeO3(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((timeO3(j)<(end_ER2time+0.9999_dp)) &
     .and.(timeO3(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nO3obs=j
 print *,'fixER2interp > O3 obs read ',nO3obs
 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_MG,*,err=997,end=120) UTGPS(j), &
     DaltGPS(j), DlatGPS(j), DLongGPS(j), QGPS(j), &
     TedrGPS(j), ReynGPS(j)
   if(.not.do_rel) then
     if((UTGPS(j)>(start_ER2time-0.9999_dp)) &
     .and.(UTGPS(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((UTGPS(j)<(end_ER2time+0.9999_dp)) &
     .and.(UTGPS(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nGPSobs=j
 print *,'fixER2interp > GPS obs read ',nGPSobs
 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_MM,*,err=997,end=120) GMTMM(j), &
     PStaMM(j), TstaMM(j), ThtaMM(j), UMM(j), VMM(j), WMM(j)
   if(.not.do_rel) then
     if((GMTMM(j)>(start_ER2time-0.9999_dp)) &
     .and.(GMTMM(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((GMTMM(j)<(end_ER2time+0.9999)) &
     .and.(GMTMM(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nMMobs=j
 print *,'fixER2interp > Meteo obs read ',nMMobs

 allocate (F_z(nMMobs),F_lat(nGPSobs),F_long(nGPSobs))
 if(z_motion) & 
   F_z(1:nMMobs) = -log(10.*real(PStaMM(1:nMMobs))/p0) ! p0 in Pascal
 if(diabatic_w.or.isentropic_motion) &
   F_z(1:nMMobs) = ThtaMM(1:nMMobs) * 0.1_dp
 F_lat(1:nGPSobs)  = real(DlatGPS(1:nGPSobs)) *1.e-5
 F_long(1:nGPSobs) = real(DlongGPS(1:nGPSobs)) *1.e-5
 
! if( activ_Lyapunov ) then
!   numpart_np = n_loc
!   numpart = numpart_np * 7
! else
   numpart = n_loc * n_sample
! endif
 error = (check_numpart() == 1).or.error
 
!  Interpolation
 print *,'fixER2interp > interpolating and filling positions'
 allocate(interp_time(n_loc))
 allocate(interp_xtra1(n_loc),interp_ytra1(n_loc),interp_ztra1(n_loc))
 do i=1,n_loc
   interp_time(i) = Start_ER2time+(i-1)*(End_ER2time-Start_ER2time)/(n_loc-1)
 enddo
! print *,julian_date(1),julian_date(2),julian_date(3)
! print *,off_jul(1),off_jul(nobs)
! write(*,'(a,2g12.5)') 'fixER2interp> interp_time ',interp_time(1),interp_time(n_loc)
 call uvip3p(3,nGPSobs,UTGPS,F_long,n_loc,interp_time,interp_xtra1)
 call uvip3p(3,nGPSobs,UTGPS,F_lat ,n_loc,interp_time,interp_ytra1)
 call uvip3p(3,nMMobs, GMTMM,F_z   ,n_loc,interp_time,interp_ztra1)
 j1=1   
 do i=1,n_loc
   itra1(j1) = nint((juldate(FDate,000000)-bdate)*86400._dbl + interp_time(i))
   error = (check_launch_time(i,j1) == 1).or.error
   do j=j1,j1+n_sample-1
     itra1(j) = nint((juldate(FDate,000000)-bdate)*86400._dbl + interp_time(i))
     xtra1(j) = (interp_xtra1(i) - xlon0)/dx
     ytra1(j) = (interp_ytra1(i) - ylat0)/dy
     ztra1(j) = interp_ztra1(i)
     !itramem(j) = itra1(j)
   enddo
   j1=j1+n_sample
 enddo
 
 print *,'fixER2interp> release end'
 print *,'fixER2interp> ',numpart
 write(*,'(a,2g12.5)') ' fixER2interp> itra1 ',itra1(1),itra1(numpart)
 write(*,'(a,2g12.5)') ' fixER2interp> xtra1 ',xtra1(1),xtra1(numpart)
 write(*,'(a,2g12.5)') ' fixER2interp> ytra1 ',ytra1(1),ytra1(numpart)  
 write(*,'(a,2g12.5)') ' fixER2interp> ztra1 ',ztra1(1),ztra1(numpart) 
 
 close(unitER2_O3)
 close(unitER2_MG)
 close(unitER2_MM)
 deallocate(timeO3, UTGPS, GMTMM,interp_time,F_z, F_lat, F_long, &
   O3, DaltGPS, DlatGPS, DLongGPS, QGPS, TedrGPS, ReynGPS, &
   PStaMM, TstaMM, ThtaMM, UMM, VMM, WMM, &
   interp_ztra1,interp_xtra1,interp_ytra1)
 return
                   

120   print *, 'fixER2interp> unexpected EOF'
      error=.true.
      print *, 'special diag ',j,GMTMM(j-1),GMTMM(j-2)
      return

995   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2SOLVEINTERP ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - NOBS TOO SMALL                      ####'
      write(*,*) '#####################################################'      
      return
996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2SOLVEINTERP ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
      write(*,*) '#####################################################'      
      return
997   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2SOLVEINTERP ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR WHILE READING THE FILE                ####'
      write(*,*) '#####################################################'      
      return     

      return    

 end subroutine fixER2SOLVEinterp
 
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXER2SPADEINTERP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine fixER2SPADEinterp(error)
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles.           						       *
!     Interpolate selected part of the flight path to n_loc points regularly
!     distributed in time. Put n_sample particles at each point.
!
!     Changes: B. Legras, Nov. 2002
!              version for distribution along an aircraft transect 
!              with temporal interpolation
!              B. Legras, Jan 2005
!              Version specialized for Spade campaign
!
!     In order to design files for new campaign, pay attention to 
!     following points:
!     - how to get pressure, lat and long
!     - time intervals of the data
!     - format of files with same name may change from campaign to campaign
!     - how to find date
!     - check scale factors (not read but hard-coded here, not good)
!
!     to do: a fortran routine reading the NASA Ames files, at least 
!            with format 1020
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]          simulation period                                       *
!                                                                              *
!*******************************************************************************

! use lyapunov, only: numpart_np, activ_Lyapunov 
 
 real, allocatable, dimension(:) :: UTFP, GMTMM,  &
      interp_time, &
      interp_xtra1, interp_ytra1, interp_ztra1, &
      F_z, F_lat, F_long
 
 integer, allocatable, dimension(:) ::  &
      PaltFP, DlatGPS, DLongGPS, TASFP,  &
      PStaMM, TstaMM, ThtaMM, UMM, VMM, WMM
 
 logical, intent(out):: error
 logical do_rel,stop_loop
 integer lhead, FDate, YYYY, MM, DD
 integer i,j,j1,nobs,nGPSobs,nMMobs
 real diff_time
 
 error=.false.
 print * ,trim(ER2_dir)//'MM'//trim(ER2_day)
 open(unitER2_MM,file=trim(ER2_dir)//'MM'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 print * ,trim(ER2_dir)//'FP'//trim(ER2_day)
 open(unitER2_FP,file=trim(ER2_dir)//'FP'//trim(ER2_day)//'.ER2', &
      status='OLD',err=996,FORM='FORMATTED')
 read(unitER2_MM,*) lhead
 do i=2,lhead
   read(unitER2_MM,*)
 enddo
 read(unitER2_FP,*) lhead
 do i=2,6
   read(unitER2_FP,*)
 enddo
 read(unitER2_FP,*) YYYY, MM, DD
 do i=8,lhead
   read(unitER2_FP,*)
 enddo
 FDate=10000*YYYY + 100*MM + DD
 diff_time = end_ER2time - start_ER2time
 nobs=int(1.05_dp*diff_time)
 allocate (UTFP(nobs),GMTMM(nobs))
 allocate (PaltFP(nobs),DlatGPS(nobs),DlongGPS(nobs),TASFP(nobs))
 allocate (PStaMM(nobs),TstaMM(nobs),ThtaMM(nobs),UMM(nobs))
 allocate (VMM(nobs),WMM(nobs))

 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_FP,*,err=997,end=120) UTFP(j), &
     PaltFP(j), DlatGPS(j), DLongGPS(j), TASFP(j)
   if(.not.do_rel) then
     if((UTFP(j)>(start_ER2time-4.9999_dp)) &
     .and.(UTFP(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((UTFP(j)<(end_ER2time+4.9999_dp)) &
     .and.(UTFP(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nGPSobs=j
 print *,'fixER2SPADEinterp > GPS obs read ',nGPSobs
 j=1 ; do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitER2_MM,*,err=997,end=120) GMTMM(j), &
     PStaMM(j), TstaMM(j), ThtaMM(j), UMM(j), VMM(j), WMM(j)
   if(.not.do_rel) then
     if((GMTMM(j)>(start_ER2time-0.9999_dp)) &
     .and.(GMTMM(j).le.start_ER2time)) then
        do_rel=.true. ; j=j+1
     else
        cycle
     endif
   else   
     if((GMTMM(j)<(end_ER2time+0.9999_dp)) &
     .and.(GMTMM(j).ge.end_ER2time)) then
        stop_loop=.true.
     else
       j=j+1 ; if (j.gt.nobs) go to 995
     endif
   endif
 enddo
 nMMobs=j
 print *,'fixER2SPADEinterp > Meteo obs read ',nMMobs

 allocate (F_z(nMMobs),F_lat(nGPSobs),F_long(nGPSobs))  
 do j=1,nMMobs
!   julian_date(j)=juldate(F_date(j),F_time(j))
   F_z(j) = -log(10.*real(PStaMM(j))/p0) ! p0 in Pascal
 enddo
 do j=1,nGPSobs
   F_lat(j)  = real(DlatGPS(j)) *1.e-3
   F_long(j) = real(DlongGPS(j))*1.e-3
 enddo
 

   numpart = n_loc * n_sample
! endif
 error = (check_numpart() == 1).or.error
 
!  Interpolation
 print *,'fixER2SPADEinterp > interpolating and filling positions'
 allocate(interp_time(n_loc))
 allocate(interp_xtra1(n_loc),interp_ytra1(n_loc),interp_ztra1(n_loc))
 do i=1,n_loc
   interp_time(i) = Start_ER2time+(i-1)*(End_ER2time-Start_ER2time)/(n_loc-1)
 enddo
! print *,julian_date(1),julian_date(2),julian_date(3)
! print *,off_jul(1),off_jul(nobs)
! write(*,'(a,2g12.5)') 'fixER2interp> interp_time ',interp_time(1),interp_time(n_loc)
 call uvip3p(3,nGPSobs,UTFP,F_long,n_loc,interp_time,interp_xtra1)
 call uvip3p(3,nGPSobs,UTFP,F_lat ,n_loc,interp_time,interp_ytra1)
 call uvip3p(3,nMMobs, GMTMM,F_z   ,n_loc,interp_time,interp_ztra1)
 j1=1   
 do i=1,n_loc
   itra1(j1) = nint((juldate(FDate,000000)-bdate)*86400._dbl + interp_time(i))
   error = (check_launch_time(i,j1) == 1).or.error
   do j=j1,j1+n_sample-1
     itra1(j) = itra1(j1)
     xtra1(j) = (interp_xtra1(i) - xlon0)/dx
     ytra1(j) = (interp_ytra1(i) - ylat0)/dy
     ztra1(j) = interp_ztra1(i)
     !itramem(j) = itra1(j1)
   enddo
   j1=j1+n_sample
 enddo
 
 print *,'fixER2SPADEinterp> release end'
 print *,'fixER2SPADEinterp> ',numpart
 
 close(unitER2_FP)
 close(unitER2_MM)
 deallocate(UTFP, GMTMM,interp_time,F_z, F_lat, F_long, &
   PaltFP, DlatGPS, DLongGPS, TASFP, &
   PStaMM, TstaMM, ThtaMM, UMM, VMM, WMM, &
   interp_ztra1,interp_xtra1,interp_ytra1)
 return
                   

120   print *, 'fixER2SPADEinterp> unexpected EOF'
      error=.true.
      print *, 'special diag ',j,GMTMM(j-1),GMTMM(j-2)
      return

995   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2INTERP :    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - NOBS TOO SMALL                      ####'
      write(*,*) '#####################################################'      
      return
996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2INTERP :    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
      write(*,*) '#####################################################'      
      return
997   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXER2INTERP:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR WHILE READING THE FILE                ####'
      write(*,*) '#####################################################'      
      return     

      return    
      

  end subroutine fixER2SPADEinterp
  
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXM55 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine fixM55(error)
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles for a M55 flight.           
!     Interpolate selected part of the flight path to n_loc points regularly
!     distributed in time. Put n_sample particles at each point.
!     Possibility to turn off the interpolation and using all recorded points
!
!     January 2017: Bernard Legras
!
!     In order to design files for new campaign, pay attention to 
!     following points:
!     - how to get pressure, lat and long
!     - time intervals of the data
!     - format of files with same name may change from campaign to campaign
!     - how to find date
!     - check scale factors (not read but hard-coded here, not good)
!
!     to do: a fortran routine reading the NASA Ames files, at least 
!            with format 1020
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]          simulation period                                       *
!                                                                              *
!*******************************************************************************
 
 real, allocatable, dimension(:) :: UTC,  &
      interp_time, &
      interp_xtra1, interp_ytra1, interp_ztra1, interp_ttra1, &
      F_Pstat, F_Tout, F_Alt, F_Lat, F_Lon, F_Halt, F_z
 
 logical, intent(out):: error
 integer lhead, FDate, YYYY, MM, DD
 integer i,j,j1,nobs_prior,nGPSobs,maxj, iday
 real diff_time, time
 real, allocatable :: buffer(:)
 
 integer (kind=4) :: outfmt,lheadidx,mode,itime
 
 print *,'fixM55 start'
 error=.false.
 select case (instrument)
   case ('UCSE') 
     print * ,trim(ER2_dir)//trim(ER2_day)//'_1_ucse.nas'
     open(unitER2_MM,file=trim(ER2_dir)//trim(ER2_day)//'_1_ucse.nas', &
         status='OLD',err=997,FORM='FORMATTED')
     read(unitER2_MM,*) lhead
     ! make buffer the size of the line to be read
     allocate (buffer(17))
     do i=2,6
       read(unitER2_MM,*)
     enddo
     read(unitER2_MM,*) YYYY, MM, DD
     do i=8,lhead
       read(unitER2_MM,*)
     enddo
     FDate=10000*YYYY + 100*MM + DD
     print *,'fixM55> ',FDate,juldate(FDate,000000),bdate
     diff_time = end_ER2time - start_ER2time
     ! allocate for 1Hz observations with a 5% margin  
     nobs_prior=int(1.05_dp*diff_time)
     allocate (UTC(nobs_prior),F_Pstat(nobs_prior),F_Tout(nobs_prior))
     allocate (F_Lat(nobs_prior),F_Lon(nobs_prior))
     j=1
     maxj=100000 ! to avoid an infinite loop
     do while (j<maxj)
       read(unitER2_MM,*,err=997,end=120) time, buffer
       if(time .le. start_ER2time-0.9999_dp) then
         cycle
       else if (time .ge. end_ER2time+0.9999_dp) then
         exit
       else
         UTC(j)   = time
         F_Tout(j)  = buffer(1)+273.15
         F_Pstat(j) = buffer(2) ! in hPa
         !F_Halt(j)  = buffer(4) 
         F_Lat(j)   = buffer(13)
         F_Lon(j)   = buffer(14)
         !F_Alt(j)   = buffer(17)    
         j=j+1
       endif
     enddo
 120 continue         
     
   case ('GLORIA')
     print * ,trim(ER2_dir)//trim(ER2_day)//'_geolocations.dat-T'
     open(unitER2_MM,file=trim(ER2_dir)//trim(ER2_day)//'_geolocations.dat-T', &
         status='OLD',err=997,FORM='FORMATTED')
     ! make buffer the size of the line to be read
     allocate (buffer(5))
     read(unitER2_MM,*) nobs_prior
     read(unitER2_MM,*) YYYY, MM, DD
     do i=1,2
       read(unitER2_MM,*)
     enddo
     FDate=10000*YYYY + 100*MM + DD
     print *,'fixM55> ',FDate,juldate(FDate,000000),bdate
     allocate (UTC(nobs_prior),F_Pstat(nobs_prior),F_Tout(nobs_prior))
     allocate (F_Alt(nobs_prior),F_Lat(nobs_prior),F_Lon(nobs_prior))
     j=1
     maxj=nobs_prior ! to avoid an infinite loop
     do while (j<maxj+1)
       read(unitER2_MM,*,err=997,end=130) iday,itime, buffer
       time = 3600*int(itime/10000)+60*modulo(int(itime/100),100)+modulo(itime,100)
       if(buffer(4) > 450) then
         cycle
       else
         UTC(j)   = time
         F_Pstat(j) = buffer(4)  ! in hPa
         F_Lat(j)   = buffer(2) 
         F_Lon(j)   = buffer(3)
         F_Tout(j)  = buffer(5)
         j=j+1 
       endif 
     enddo
 130 continue    

case ('BAL')
     print * ,trim(ER2_dir)//trim(ER2_day)//'_1_balloon'//trim(flight)//'.nas'
     open(unitER2_MM,file=trim(ER2_dir)//trim(ER2_day)//'_1_balloon'//trim(flight)//'.nas', &
         status='OLD',err=997,FORM='FORMATTED')
     read(unitER2_MM,*) lhead
     ! make buffer the size of the line to be read
     allocate (buffer(17))
     do i=2,6
       read(unitER2_MM,*)
     enddo
     read(unitER2_MM,*) YYYY, MM, DD
     do i=8,lhead
       read(unitER2_MM,*)
     enddo
     FDate=10000*YYYY + 100*MM + DD
     print *,'fixM55> ',FDate,juldate(FDate,000000),bdate
     diff_time = end_ER2time - start_ER2time
     ! allocate for 1Hz observations with a 5% margin  
     nobs_prior=int(1.05_dp*diff_time)
     allocate (UTC(nobs_prior),F_Pstat(nobs_prior),F_Tout(nobs_prior))
     allocate (F_Lat(nobs_prior),F_Lon(nobs_prior))
     j=1
     maxj=100000 ! to avoid an infinite loop
     do while (j<maxj)
       read(unitER2_MM,*,err=997,end=140) time, buffer
       if(time .le. start_ER2time-0.9999_dp) then
         cycle
       else if (time .ge. end_ER2time+0.9999_dp) then
         exit
       else
         UTC(j)   = time
         ! eliminate points above the pressure low cut
         if(buffer(2) <= 0.01*plowcut) cycle
         ! replace missing PTU data
         if (buffer(1)>990) then
           F_Tout(j) = F_tout(j-1)
           F_Pstat(j) = F_Pstat(j-1)           
         else
           F_Tout(j)  = buffer(1)+273.15 ! conversion from Celsius
           F_Pstat(j) = buffer(2) ! in hPa
         endif
         ! replace missing GPS data
         if ((buffer(13)>990).or.(buffer(14)>990)) then
           F_lat(j) = F_lat(j-1)
           F_lon(j) = F_lon(j-1)
         else
           F_Lat(j)   = buffer(13)
           F_Lon(j)   = buffer(14)
         endif
         j=j+1
       endif
     enddo
 140 continue

 end select
 
 close(unitER2_MM)
 nGPSobs=j-1
 print *,'fixM55> obs read ',nGPSobs
 ! print a couple of lines of dignostics
 print *,'min max F_lat ',minval(F_lat(1:nGPSobs)),maxval(F_lat(1:nGPSobs))
 print *,'min max F_lon ',minval(F_lon(1:nGPSobs)),maxval(F_lon(1:nGPSobs))
 print *,'min max F_Pstat ',minval(F_Pstat(1:nGPSobs)),maxval(F_Pstat(1:nGPSobs))
 
 allocate (F_z(nGPSobs))
 if (z_motion) then
   F_z(1:nGPSobs) = -log((100.*F_PStat(1:nGPSobs))/p0) ! p0 in Pascal
 else
   ! This will not work with GLORIA !!!
   F_z(1:nGPSobs) = F_Tout(1:nGPSobs) * (p0/(100.*F_Pstat(1:nGPSobs)))**kappa 
 endif

 if (interp_release) then  
   !  Interpolation
   print *,'fixM55> interpolating and filling positions'
   allocate(interp_time(n_loc))
   allocate(interp_xtra1(n_loc),interp_ytra1(n_loc))
   allocate(interp_ttra1(n_loc),interp_ztra1(n_loc))
   do i=1,n_loc
     interp_time(i) = Start_ER2time+(i-1)*(End_ER2time-Start_ER2time)/(n_loc-1)
   enddo
   call uvip3p(3,nGPSobs,UTC,F_Lon,n_loc,interp_time,interp_xtra1)
   call uvip3p(3,nGPSobs,UTC,F_Lat,n_loc,interp_time,interp_ytra1)
   call uvip3p(3,nGPSobs,UTC,F_z,n_loc,interp_time,interp_ztra1)
   call uvip3p(3,nGPSobs,UTC,F_Tout,n_loc,interp_time,interp_ttra1)
 else
   n_loc =  nGPSobs
   print *,'fixM55> filling positions n_loc ',n_loc 
   allocate(interp_time(n_loc))
   allocate(interp_xtra1(n_loc),interp_ytra1(n_loc))
   allocate(interp_ttra1(n_loc),interp_ztra1(n_loc))
   interp_time  = UTC
   interp_xtra1 = F_Lon
   interp_ytra1 = F_lat
   interp_ztra1 = F_z
   interp_ttra1 = F_Tout
 endif
 
 deallocate(UTC,F_lon,F_lat,F_z,F_Tout)
 
 numpart = n_loc*n_sample
 print *,'fixM55> numpart ',numpart
 !error = (check_numpart() == 1).or.error
  
 if(.not.allocated(itra1)) allocate(itra1(numpart))
 if(.not.allocated(xtra1)) allocate(xtra1(numpart))
 if(.not.allocated(ytra1)) allocate(ytra1(numpart))
 if(.not.allocated(ztra1)) allocate(ztra1(numpart))
 if(.not.allocated(ttra1)) allocate(ttra1(numpart))
 if(.not.allocated(itra0)) allocate(itra0(numpart))
 if(.not.allocated(flagTB)) allocate(flagTB(numpart)) 
 j1=1   
 do i=1,n_loc
   itra1(j1) = nint((juldate(FDate,000000)-bdate)*86400._dbl + interp_time(i))
   error = (check_launch_time(i,j1) == 1).or.error
   do j=j1,j1+n_sample-1
     itra1(j) = itra1(j1)
     xtra1(j) = (interp_xtra1(i) - xlon0)/dx
     ytra1(j) = (interp_ytra1(i) - ylat0)/dy
     ztra1(j) = interp_ztra1(i)
     ttra1(j) = interp_ttra1(i)
   enddo
   j1=j1+n_sample
 enddo
 itra0(:) = itra1(:)
 if (z_motion) then
   !flagTB(:) = Z'3E'
   flagTB(:) = 62
 else
   !flagTB(:) = Z'7E'
   flagTB(:) = 126
 endif 
 print *,minval(itra1),maxval(itra1) 
 print *,minval(xtra1),maxval(xtra1)
 print *,minval(ztra1),maxval(ztra1)
 print *,minval(ttra1),maxval(ttra1) 
 
 print *,'fixM55> release end'
 
 deallocate(interp_time,interp_ztra1,interp_xtra1,interp_ytra1,interp_ttra1)
   
! Output section to generate an index file that can be used for the    
! analysis

 print *,'fixM55> store index' 
 outfmt=107       ! Output format
 lheadidx = 3        ! Header length (# of lines)
 mode = 1
 open(unitpartout,file=trim(path(2))//'index_old',  &
          form='unformatted')
 write(unitpartout) lheadidx,outfmt,mode  
 itime=0
 write(unitpartout) stamp_date, itime, int(lsynctime,4) 
 write(unitpartout) int(numpart,4), int(numpart,4), int(idx_orgn,4), &
                    int(0,4), int(0,4), int(0,4)   
 write(unitpartout) int(flagTB(:),4)
 write(unitpartout) int(itra0(:),4)
 write(unitpartout) real(xlon0+xtra1(:)*dx,4)
 write(unitpartout) real(ylat0+ytra1(:)*dy,4)
 if (z_motion) then
   write(unitpartout) real(p0*exp(-ztra1(:)),4)
 else
   write(unitpartout) real(p0*(ttra1(:)/ztra1(:))**(1/kappa),4)
 endif
 write(unitpartout) real(ttra1(:),4)   
 write(unitpartout) (int(i,4),i=idx_orgn,idx_orgn+numpart-1)
 close(unitpartout)
   
 return
                   

997 error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### TRACZILLA MODEL SUBROUTINE FIXM55      :    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
      write(*,*) '#####################################################'      
      return    

  end subroutine fixM55 

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXMOZAIC INTERP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine fixmozaicinterp(error)
!*******************************************************************************
!
!     This routine fixes the release times and release locations of all
!     particles.
!     Interpolate selected part of the flight path to n_loc points regularly
!     distributed in time. Put n_sample particles at each point.
!
!     Changes: B. Legras, Jun. 2002
!              version for distribution along an aircraft transect 
!              with temporal interpolation
!
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]          simulation period                                       *
!                                                                              *
!*******************************************************************************

 real, allocatable, dimension(:) :: F_lat,F_long,F_baroalt,F_radioalt, &
      F_press, A340_statictemp,A340_airspeed,F_groundspeed,F_zonalwind, &
      F_meridionalwind,F_ozone,F_statictemp,F_relhumid,F_RHvalid,F_RHaccur, &
      F_water,F_z,off_jul,interp_time, &
      interp_xtra1, interp_ytra1, interp_ztra1
 integer, allocatable, dimension (:) :: F_date, F_time 
 character(len=96):: line
 logical error,do_rel,stop_loop
 integer hh,mm,ss,hh_start,mm_start,ss_start,hh_end,mm_end,ss_end
 integer i,j,j1,rank_obs,nobs
 real diff_time
 real (dbl), allocatable :: julian_date(:)
 
 error=.false.
 open(unitmozaic,file=trim(MOZAIC_dir)//trim(MOZAIC_filename), &
      status='OLD',err=996,FORM='FORMATTED')
 
! read first 3 lines
 read(unitmozaic,*)
 read(unitmozaic,'(a)') line
 read(unitmozaic,*)
 print *,'fixmozaicinterp >',line
 
 select case (index_type)
 case ('rank')
   nobs = end_index - start_index +1
 case ('time')
   hh_start = start_index/10000
   mm_start = (start_index-10000*hh_start)/100
   ss_start = (start_index-10000*hh_start-100*mm_start)
   hh_end = end_index/10000
   mm_end = (end_index-10000*hh_end)/100
   ss_end = end_index-10000*hh_end-100*mm_end
   if(hh_end > hh_start) then 
      diff_time = (hh_end - hh_start)*3600
   else
      diff_time = (24+hh_end - hh_start)*3600
   endif
   diff_time=diff_time+(mm_end - mm_start)*60 + &
                        ss_end - ss_start
   nobs=int(diff_time/4 * 1.05_dp)
   allocate (F_lat(nobs),F_long(nobs),F_baroalt(nobs),F_radioalt(nobs))
   allocate (F_press(nobs),A340_statictemp(nobs),A340_airspeed(nobs))
   allocate (F_groundspeed(nobs),F_zonalwind(nobs),F_meridionalwind(nobs))
   allocate (F_ozone(nobs),F_statictemp(nobs),F_relhumid(nobs),F_RHvalid(nobs))
   allocate (F_RHaccur(nobs),F_water(nobs),F_date(nobs),F_time(nobs))
   allocate (F_z(nobs),julian_date(nobs),off_jul(nobs))
 end select
 
 rank_obs=0
 j=1
 do_rel=.false. ; stop_loop=.false.
 do while(.not. stop_loop)
   read(unitmozaic,*,err=997,end=120) &
      F_date(j),F_time(j),F_lat(j),F_long(j),F_baroalt(j),F_radioalt(j),&
      F_press(j),A340_statictemp(j),A340_airspeed(j),F_groundspeed(j),  &
      F_zonalwind(j),F_meridionalwind(j),F_ozone(j),F_statictemp(j),    &
      F_relhumid(j),F_RHvalid(j),F_RHaccur(j),F_water(j)
   rank_obs=rank_obs+1
   select case (index_type)
   case ('rank')
      if(rank_obs < start_index)cycle
      if(rank_obs == end_index) then
        stop_loop=.true. ; cycle
      endif
      do_rel=.true.
      j=j+1
   case ('time')
      hh = F_time(j)/10000
      mm = (F_time(j)-10000*hh)/100
      ss = F_time(j)-10000*hh-100*mm
      if(.not.do_rel) then
        if((hh==hh_start).and.(mm==mm_start).and.(ss==ss_start)) then
          do_rel=.true.
          print *,'fixmozaicinterp> release start'
          j=j+1
        else
          cycle
        endif
      else   
        if((hh==hh_end).and.(mm==mm_end).and.(ss==ss_end)) then
          stop_loop=.true.
        else
          j=j+1
          if (j.gt.nobs) go to 995
      endif
      endif
   end select
 enddo
 
 nobs=j
 print *,'fixmozaicinterp > nobs ',nobs
 do j=1,nobs
   julian_date(j)=juldate(F_date(j),F_time(j))
   F_z(j) = -log(F_press(j)/p0)
 enddo
 do j=1,nobs
   off_jul(j)=julian_date(j)-julian_date(1)
 enddo
 
 numpart = n_loc * n_sample
 error = (check_numpart() == 1).or.error
 
!  Interpolation
 print *,'fixmozaicinterp > interpolating and filling positions' 
 allocate(interp_time(n_loc))
 allocate(interp_xtra1(n_loc),interp_ytra1(n_loc),interp_ztra1(n_loc))
 do i=1,n_loc
   interp_time(i) = (i-1)*(julian_date(nobs)-julian_date(1))/(n_loc-1)
 enddo
 print *,julian_date(1),julian_date(2),julian_date(3)
 print *,off_jul(1),off_jul(nobs)
 print *,interp_time(1),interp_time(n_loc)
 call uvip3p(3,nobs,off_jul,F_long,n_loc,interp_time,interp_xtra1)
 call uvip3p(3,nobs,off_jul,F_lat ,n_loc,interp_time,interp_ytra1)
 call uvip3p(3,nobs,off_jul,F_z   ,n_loc,interp_time,interp_ztra1)
 j1=1
 do i=1,n_loc
   itra1(j1) = nint((julian_date(1)-bdate & 
            +(i-1)*(julian_date(nobs)-julian_date(1))/(n_loc-1))*86400)
   error = (check_launch_time(i,j1) == 1).or.error
   do j=j1,j1+n_sample-1
     itra1(j) = nint((julian_date(1)-bdate & 
            +(i-1)*(julian_date(nobs)-julian_date(1))/(n_loc-1))*86400)
     xtra1(j) = (interp_xtra1(i) - xlon0)/dx
     ytra1(j) = (interp_ytra1(i) - ylat0)/dy
     ztra1(j) = interp_ztra1(i)
     !itramem(j) = itra1(j)
   enddo
   j1=j1+n_sample
 enddo
 
 print *,'fixmozaicinterp> release end'
 print *,'fixmozaicinterp> ',numpart
 print *,'fixmozaicinterp> F_long(obs) ',F_long(1),F_long(nobs)
 print *,'fixmozaicinterp>             ',xtra1(1)*dx + xlon0,xtra1(numpart)*dx + xlon0
 print *,'fixmozaicinterp> F_lat(obs)  ',F_lat(1),F_lat(nobs)
 print *,'fixmozaicinterp>             ',ytra1(1)*dy + ylat0,ytra1(numpart)*dy + ylat0
 print *,'fixmozaicinterp> F_z  (obs)  ',F_z(1),F_z(nobs)
 print *,'fixmozaicinterp>             ',ztra1(1),ztra1(numpart)
 
 close(unitmozaic)
 deallocate(F_lat,F_long,F_baroalt,F_radioalt, &
      F_press, A340_statictemp,A340_airspeed,F_groundspeed,F_zonalwind, &
      F_meridionalwind,F_ozone,F_statictemp,F_relhumid,F_RHvalid,F_RHaccur, &
      F_water,F_time,F_date,F_z,interp_time,julian_date,off_jul, &
      interp_ztra1,interp_xtra1,interp_ytra1)
 return                  

120   print *, 'fixmozaicpart> unexpected EOF'
      error=.true.
      return

995   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMOZAICPART:    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - NOBS TOO SMALL                      ####'
      write(*,*) '#####################################################'      
      return
996   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMOZAICPART:    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
      write(*,*) '#####################################################'      
      return
997   error=.true.
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMOZAICPART:    ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR WHILE READING THE FILE                ####'
      write(*,*) '#####################################################'      
      return 
      
      end subroutine fixmozaicinterp
      
!===============================================================================
      
     function check_launch_time(i,j1)
     integer, intent(in):: i,j1
     integer :: check_launch_time
     
     check_launch_time=0
     if(itra1(j1)>0) then
       write(*,*) ' #### FLEXPART MODEL ERROR! LAUNCH TIME >0    #### ' 
       write(*,*) i,itra1(j1)
       check_launch_time=1
     endif
     end function check_launch_time
     
     function check_numpart()
     integer :: check_numpart
     check_numpart=0
     if(numpart > maxpart) then
       write(*,*) '#####################################################'
       write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMOZAICINTERP : ####'
       write(*,*) '####                                             ####'
       write(*,*) '#### NB OF REQUIRED PARTICLES EXCEEDS MAXPART    ####'
       write(*,*) '#####################################################'             
       check_numpart=1
     endif
     end function check_numpart

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXMOZAIC PART @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
  
 subroutine fixmozaicpart(error)
!*******************************************************************************
!                                                                              *
!     This routine fixes the release times and release locations of all        *
!     particles.           						       *
!     Put n_sample particles at each measurement point of the selected     
!     part of the flight
!
!     Changes: B. Legras, Apr. 2002
!              version for distribution along a profile
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]          simulation period                                       *
!                                                                              *
!*******************************************************************************

 real F_lat,F_long,F_baroalt,F_radioalt,F_press, &
      A340_statictemp,A340_airspeed,F_groundspeed,F_zonalwind, &
      F_meridionalwind,F_ozone,F_statictemp,F_relhumid,F_RHvalid,F_RHaccur, &
      F_water
 integer F_date,F_time 
 character(len=96):: line
 logical error,do_rel,stop_loop
 integer hh,mm,ss,hh_start,mm_start,ss_start,hh_end,mm_end,ss_end
 integer i,rank_obs,itra_part
 real x_part,y_part
 real (dbl) :: julian_date
 
 error=.false.
 open(unitmozaic,file=trim(MOZAIC_dir)//trim(MOZAIC_filename), &
      status='OLD',err=996,FORM='FORMATTED')
 
! read first 3 lines
 read(unitmozaic,*)
 read(unitmozaic,'(a)') line
 read(unitmozaic,*)
 print *,'fixmozaicpart >',line
 
 if(index_type=='time') then
   hh_start = start_index/10000
   mm_start = (start_index-10000*hh_start)/100
   ss_start = (start_index-10000*hh_start-100*mm_start)
   hh_end = end_index/10000
   mm_end = (end_index-10000*hh_end)/100
   ss_end = end_index-10000*hh_end-100*mm_end
 endif  
 
 rank_obs=0
 do_rel=.false. ; stop_loop=.false.
 numpart=0
 do while(.not. stop_loop)
   read(unitmozaic,*,err=997,end=120) &
      F_date,F_time,F_lat,F_long,F_baroalt,F_radioalt,F_press,&
      A340_statictemp,A340_airspeed,F_groundspeed,F_zonalwind,F_meridionalwind,&
      F_ozone,F_statictemp,F_relhumid,F_RHvalid,F_RHaccur,F_water
   rank_obs=rank_obs+1
   select case (index_type)
   case ('rank')
      if(rank_obs < start_index)cycle
      if(rank_obs > end_index) then
        stop_loop=.true. ; cycle
      endif
      do_rel=.true.
   case ('time')
      hh = F_time/10000
      mm = (F_time-10000*hh)/100
      ss = F_time-10000*hh-100*mm
      if(.not.do_rel) then
        if((hh==hh_start).and.(mm==mm_start).and.(ss==ss_start)) then
          do_rel=.true.
          print *,'fixmozaicpart> release start'
        else
          cycle
        endif
      else if((hh==hh_end).and.(mm==mm_end).and.(ss==ss_end)) then
        stop_loop=.true.
      endif
   end select
   
   if(.not.do_rel) cycle
   julian_date=juldate(F_date,F_time)
   if(instant_release) then
      itra_part=nint((julian_date-bdate)*86400)
   else
      itra_part=sign(nint(abs(julian_date-bdate)*86400/abs(lsynctime)+0.5) &
         *lsynctime,nint((julian_date-bdate)*86400))
   endif
   x_part=(F_long-xlon0)/dx
   y_part=(F_lat-ylat0)/dy
!  print *,'fixmozaic> ',F_date,F_time,itra_part
   if(numpart+n_sample > maxpart) go to 998
   do i=numpart+1,numpart+n_sample
      itra1(i) = itra_part
      xtra1(i) = x_part
      ytra1(i) = y_part
      ztra1(i) = -log(F_press/p0)
      !itramem(i) = itra_part                       
   enddo
   numpart=numpart+n_sample
 enddo
 print *,'fixmozaicpart> release end'
 print *,'fixmozaicpart> ',numpart
 close(unitmozaic)
 return                  

120 print *, 'fixmozaicpart> unexpected EOF'
    error=.true.
    return

996 error=.true.
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMOZAICPART:    ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
    write(*,*) '#####################################################'      
    return
997 error=.true.
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMOZAICPART:    ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### ERROR WHILE READING THE FILE                ####'
    write(*,*) '#####################################################'      
    return     
998 error=.true.
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXMOZAICPART:    ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### NUMBER OF REQUIRED POINTS EXCEEDS MAXPART   ####'
    write(*,*) '#####################################################'      
    return          

  end subroutine fixmozaicpart



  
  !===============================================================================
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIXPARTICLESCLAUS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !=====|==1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine fixparticlesClaus(error)
    !*******************************************************************************
    !                                                                              *
    !     This routine fixes the release times and release locations of all        *
    !     particles.           						           *
    !                                                                              *
    !     Released particles : at the top of clouds in the tropics                 *
    !                                                                              *
    !     Changes: A.-S. Tissier, may 2013                                         *
    !                                                                              *
    !*******************************************************************************
    
    use commons
    use mass_iso

    implicit none

    logical                                   :: error
    character(3)                              :: Tbmax_string
    character(2)                              :: latmin_string_short, latmax_string_short
    character(3)                              :: latmin_string, latmax_string
    character(2)                              :: yearstring_short
    character(10)                             :: fulldate_Claus_string
    character(26)                             :: filename

    integer                                   :: ok_Claus
    
    integer                                   :: ibhourtrac, ibmintrac, ibsectrac
    integer                                   :: ibyeartrac, ibmonthtrac, ibdaytrac
    integer                                   :: hourtrac, mintrac, sectrac
    integer                                   :: yeartrac, monthtrac, daytrac
    integer                                   :: iehourtrac
    integer                                   :: ieyeartrac, iemonthtrac, iedaytrac
    integer                                   :: fulldate_Claus
    integer                                   :: fulldate_Claus_end
    integer,dimension(12)                     :: tab_day_month

    integer                                   :: ios
    integer                                   :: numpart_date
    logical                                   :: test_exist
    integer                                   :: numpart_prev

    real, dimension(:), allocatable           :: presspart, temppart
    real, dimension(:), allocatable           :: xpart, ypart, zpart, thetapart
    real, dimension(:), allocatable           :: ir_start

    integer                                   :: start_date, start_time
    integer                                   :: ipart

    error=.false.

    ok_Claus = 1

    ! Check input date :
    call check_input_date_Claus(error)
    if(error) return

    ! Determination de la date du premier lancer :
    call date_first_release()

    ! Recherche de la date de fin : fulldate_Claus_end
    ! Determination de la date du dernier lancer :
    call date_last_release()


    ! Conversion Tb, latitude, longitude, etc. en chaine de caracteres
    ! pour les titres de fichier
    call conversion_string()

    ! Lecture et stockage des donnees :
    ok_claus = 1
    numpart_prev = 0

    ! Premier temps :
    yeartrac = ibyeartrac
    monthtrac = ibmonthtrac
    daytrac = ibdaytrac
    hourtrac = ibhourtrac

    fulldate_Claus = yeartrac*10**(6) + monthtrac*10**(4) + daytrac*10**(2) + hourtrac
    write(fulldate_Claus_string,'(i10)') fulldate_Claus
    filename = fulldate_Claus_string//'_Tb'//Tbmax_string//'_lat'//latmin_string//latmax_string

    do while(ok_claus.EQ.1)
       ! Test d'existence du fichier a lire :
       ! on utilise une variable logique
       inquire( file=trim(Claus_dir)//yearstring_short//'/'//filename, exist=test_exist)
       if( .NOT. test_exist) go to 998

       ! Lecture du fichier : juste la premiere variable pour connaitre le nombre de particules
       open(UNIT=unitCLAUS, FILE=trim(Claus_dir)//yearstring_short//'/'//filename, &
            FORM="unformatted", ACCESS="sequential", &
            ACTION="read", POSITION="rewind", &
            IOSTAT=ios )
       write(*,*) trim(Claus_dir)//yearstring_short//'/'//filename
       if(ios /=0) go to 996
       read(UNIT=unitCLAUS, IOSTAT=ios) numpart_date
       if(ios /=0) go to 997

       ! Allocation des tableaux intermediaires :
       allocate(xpart(numpart_date), ypart(numpart_date), presspart(numpart_date), temppart(numpart_date))
       allocate(thetapart(numpart_date), zpart(numpart_date))
       allocate(ir_start(numpart_date))

       read(UNIT=unitCLAUS, IOSTAT=ios) xpart
       if(ios /=0) go to 997
       read(UNIT=unitCLAUS, IOSTAT=ios) ypart
       if(ios /=0) go to 997
       read(UNIT=unitCLAUS, IOSTAT=ios) presspart
       if(ios /=0) go to 997
       read(UNIT=unitCLAUS, IOSTAT=ios) temppart
       if(ios /=0) go to 997

       close(UNIT=unitCLAUS)

       ! Computing the time to release the particle
       do ipart = 1,numpart_date ! in forward calculation 
          start_date = yeartrac*10**(4) + monthtrac*10**(2)+daytrac
          start_time = hourtrac*10**(4)
          ir_start(ipart) = int((juldate(start_date,start_time)-bdate)*86400.)
       enddo

       ! Check initialisation :
       do ipart = 1,numpart_date
         if(temppart(ipart)>TB_max) go to 1000
       enddo

       ! Calcul de la coordonnee verticale :
       if(diabatic_Claus) then
          ! coordonnees theta
          thetapart = temppart*(p0/presspart)**kappa
       else
          ! coordonnees z
          zpart = log(p0/presspart)
       endif

       ! Copie dans les tableaux finaux :
       do ipart = 1,numpart_date
          numpart_prev = numpart_prev +1
          xtra1(numpart_prev) = modulo((xpart(ipart)-xlon0)/dx,real(nx-1))
          ytra1(numpart_prev) = (ypart(ipart)-ylat0)/real(dy)
          ttra1(numpart_prev) = temppart(ipart)
          if(diabatic_Claus) then
             ztra1(numpart_prev) = thetapart(ipart)
          else
             ztra1(numpart_prev) = zpart(ipart)
          endif
          itra1(numpart_prev) = ir_start(ipart)
          !if(numpart_prev==1148648) then
          !  write(*,*) 'fixparticlesCLAUS'
          !  write(*,*) 'numpart_prev = 1148648'
          !  write(*,*) 'ttra1=',thetapart(ipart)*(p0/presspart(ipart))**(-kappa)
          !  write(*,*) 'xtra1=',xtra1(numpart_prev)
          !  write(*,*) 'ytra1=',ytra1(numpart_prev)
          !  write(*,*) 'ztra1=',thetapart(ipart)
          !  write(*,*) 'presspart=',presspart(ipart)
          !endif
       enddo
       ! Deallocation des tableaux intermediaires
       deallocate(xpart,ypart,presspart,temppart)
       deallocate(thetapart,zpart)
       deallocate(ir_start)


       if(fulldate_Claus .EQ. fulldate_Claus_end) then
          ! il n'y a plus de donnees Claus a lire
          ok_claus=0
       else
          ! Passage a l'heure suivante :
          hourtrac = hourtrac + 3 ! data every 3 hours for Claus

          call newdate()

          ! Prochain filename a lire :
          fulldate_Claus = yeartrac*10**(6) + monthtrac*10**(4) + daytrac*10**(2) + hourtrac
          write(fulldate_Claus_string,'(i10)') fulldate_Claus
          filename = fulldate_Claus_string//'_Tb'//Tbmax_string//'_lat'//latmin_string//latmax_string
       endif ! if(fulldate_Claus_string.EQ.fulldate_Claus_string_end)

    enddo

    numpart = numpart_prev
    write(*,*)
    write(*,*) 'number of particles numpart=',numpart

    if(numpart>maxpart) goto 994

    return   ! subroutine principale : fixparticlesClaus

    ! ---------------------------------
    ! ----- Affichage des erreurs -----
    ! ---------------------------------

994 error=.true.
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXPARTICLES:     ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES SPECIFIED ####'
    write(*,*) '#### IN FILE "RELEASES" OR CARRIED FORWARD FROM  ####'
    write(*,*) '#### PREVIOUS RUN EXCEEDS THE MAXIMUM ALLOWED    ####'
    write(*,*) '#### NUMBER. REDUCE EITHER NUMBER OF PARTICLES   ####'
    write(*,*) '#### PER RELEASE POINT OR REDUCE NUMBER OF       ####'
    write(*,*) '#### RELEASE POINTS.                             ####'
    write(*,*) '#####################################################'
    return

996 error=.true.
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXPARTICLES      ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
    write(*,*) '#### '//trim(Claus_dir)//yearstring_short//'/'//filename//' ####'
    write(*,*) '#####################################################'
    return

997 error=.true.
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXPARTICLES      ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### ERROR WHILE READING THE FILE                ####'
    write(*,*) '#### '//trim(Claus_dir)//yearstring_short//'/'//filename//' ####'
    write(*,*) '#####################################################'
    return

998 error=.true.
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE FIXPARTICLES      ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### ERROR THIS FILE DOESN''T EXIST :            ####'
    write(*,*) '#### ' //trim(Claus_dir)//yearstring_short//'/'//filename//' ####'
    write(*,*) '#####################################################'
    return

1000 error=.true.
     write(*,*) '#####################################################'
     write(*,*) '#### FLEXPART MODEL SURBROUTINE FIXPARTICLES     ####'
     write(*,*) '####                                             ####'
     write(*,*) '#### ERROR : temperature > TB_max                ####'
     write(*,*) '#####################################################'

    ! ----------------------------------
    ! ---------- Contains --------------
    !-----------------------------------

  contains

    subroutine check_input_date_Claus(error)
      logical, intent(inout) :: error

      if(iedate < iedate_Claus) go to 9998
      if(iedate.EQ.iedate_Claus) then
         if(ietime<=ietime_Claus) go to 9998
      endif

      if(ibdate > iedate_Claus) go to 9999
      if(ibdate.EQ.iedate_Claus)then
         if(ibtime>=ietime_Claus) go to 9999
      endif

      return ! subroutine check_input_date_Claus

9998  error=.true.
      write(*,*) ' ################################################# '
      write(*,*) ' #### FLEXPART MODEL ERROR! CLAUS ENDING DATE #### '
      write(*,*) ' #### IS LARGER THAN ENDING DATE. CHANGE      #### '
      write(*,*) ' #### CLAUS ENDING DATE IN FILE "RELEASES"    #### '
      write(*,*) ' ################################################# '
      return

9999  error=.true.
      write(*,*) ' ################################################# '
      write(*,*) ' #### FLEXPART MODEL ERROR! CLAUS ENDING DATE #### '
      write(*,*) ' #### IS SMALLER THAN BEGINNIG DATE. CHANGE   #### '
      write(*,*) ' #### CLAUS ENDING DATE IN FILE "RELEASES"    #### '
      write(*,*) ' ################################################# '
      return

    end subroutine check_input_date_Claus

    ! -------------------------
    ! -------------------------

    subroutine date_first_release()

      ibyeartrac = int(aint(real(ibdate)*10.0**(-4)))
      ibmonthtrac = int(aint(real(ibdate-ibyeartrac*10**(4))*10.0**(-2)))
      ibdaytrac = int(aint(real(ibdate-ibyeartrac*10**(4)-ibmonthtrac*10**(2))))

      ibhourtrac = int(aint(real(ibtime)*10.0**(-4)))
      ibmintrac = int(aint(real(ibtime-ibhourtrac*10**(4))*10.0**(-2)))
      ibsectrac = int(aint(real(ibtime-ibhourtrac*10**(4)-ibmintrac*10**(2))))

      if(ibsectrac>0) then
         ibsectrac=0
         ibmintrac = ibmintrac+1
      endif
      if(ibmintrac>0) then
         ibhourtrac = ibhourtrac + 1
         ibmintrac = 0
      endif
      do while(mod(ibhourtrac,3).NE. 0)  ! data every 3 hours for Claus
         ibhourtrac = ibhourtrac + 1
      enddo
      if(ibhourtrac ==24) then
         ibhourtrac =0
         ibdaytrac = ibdaytrac + 1
      endif

      if((ibyeartrac<1997).OR.(ibyeartrac>2009)) then
         stop 'check ibdate to use Claus data'
      endif

      tab_day_month = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      if((ibyeartrac.EQ.2000) .OR. (ibyeartrac.EQ.2004) .OR. &
         (ibyeartrac.EQ.2008) .or. (ibyeartrac.eq.2012)) then
         tab_day_month(2) = 29
      endif
      if(ibdaytrac>tab_day_month(ibmonthtrac)) then
         ibdaytrac = 1
         ibmonthtrac = ibmonthtrac + 1
      endif
      if(ibmonthtrac==13) then
         ibmonthtrac=1
         ibyeartrac = ibyeartrac+1
      endif

      if(ibyeartrac>=2000) then
         write(yearstring_short,'(a1,i1)') '0',ibyeartrac-2000
      else
         write(yearstring_short,'(i2)') ibyeartrac-1900
      endif

    end subroutine date_first_release

    ! -------------------------
    ! -------------------------

    subroutine date_last_release()

      ieyeartrac = int(aint(real(iedate_Claus)*10.0**(-4)))
      iemonthtrac = int(aint(real(iedate_Claus-ieyeartrac*10**(4))*10.0**(-2)))
      iedaytrac = int(aint(real(iedate_Claus-ieyeartrac*10**(4)-iemonthtrac*10**(2))))
      ! Determination de l'heure du dernier lancer :
      iehourtrac = int(aint(real(ietime_Claus)*10.0**(-4)))

      do while(mod(iehourtrac,3).NE.0) ! data every 3 hours for Claus
         iehourtrac = iehourtrac-1
      enddo

      fulldate_Claus_end = ieyeartrac*10**(6) + iemonthtrac*10**(4) + iedaytrac*10**(2) + iehourtrac

    end subroutine date_last_release

    ! -------------------------
    ! -------------------------

    subroutine conversion_string()
      write(Tbmax_string,'(i3)') TB_max

      write(latmin_string_short,'(i2)') abs(latmin_Claus)
      write(latmax_string_short,'(i2)') abs(latmax_Claus)

      if(latmin_Claus<0) then
         if(latmin_Claus>-10) then
            write(latmin_string,'(a2,i1)') 'm0',abs(latmin_Claus)
         else
            write(latmin_string,'(a1,i2)') 'm',abs(latmin_Claus)
         endif
      else
         if(latmin_Claus<10) then
            write(latmin_string,'(a2,i1)') 'p0',abs(latmin_Claus)
         else
            write(latmin_string,'(a1,i2)') 'p',abs(latmin_Claus)
         endif
      endif

      if(latmax_Claus<0) then
         if(latmax_Claus>-10) then
            write(latmax_string,'(a2,i1)') 'm0',abs(latmax_Claus)
         else
            write(latmax_string,'(a1,i2)') 'm',abs(latmax_Claus)
         endif
      else
         if(latmax_Claus<10) then
            write(latmax_string,'(a2,i1)') 'p0',abs(latmax_Claus)
         else
            write(latmax_string,'(a1,i2)') 'p',abs(latmax_Claus)
         endif
      endif

    end subroutine conversion_string


    ! -------------------------
    ! -------------------------

    subroutine newdate()

      integer,dimension(12)                     :: tab_day_month

      ! Possible passage au jour suivant :
      if(hourtrac.EQ.24) then
         hourtrac =0
         daytrac = daytrac + 1

         ! Possible passage au mois suivant :
         tab_day_month = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         if((yeartrac.EQ.2000) .OR. (yeartrac.EQ.2004) .or. (yeartrac==2008)) then
            tab_day_month(2) = 29
         endif
         if(daytrac>tab_day_month(monthtrac)) then
            daytrac = 1
            monthtrac = monthtrac + 1

            ! Possible passage a l'annee suivante :
            if(monthtrac==13) then
               monthtrac=1
               yeartrac = yeartrac+1
               if(yeartrac<2000) then
                  write(yearstring_short,'(i2)') yeartrac-2000
               else
                  write(yearstring_short,'(i2)') yeartrac-1900
               endif
            endif
         endif
      endif ! if(hourtrac.EQ.24)

    end subroutine newdate

  end subroutine fixparticlesClaus

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@ fixparticlesTBStratoClim @@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine fixparticlesTBStratoClim(error)
  
! This routine releases the parcel defined by the Stratoclim processing 
! of cloud tops from geostationry satellites written in 107 format
! This is the single file version which is well adapted to the real time
! mode
! May 2016: B. Legras
      
  use commons
  use mass_iso

  implicit none

  !----- Declarations ---

  ! Input :
  !integer :: ibyeartrac, ibmonthtrac, ibdaytrac
  !integer :: ibhourtrac, ibmintrac, ibsectrac
  !integer :: Tbmax   ! brightness temperature threshold

  ! Output :
  logical :: error

  ! Read :
  !integer :: numpart

  ! Only here :
  integer (kind=4), allocatable :: ir_start(:), flagTB0(:)
  real (sp), allocatable :: xpart(:), ypart(:), presspart(:), temppart(:)
  real, allocatable :: thetapart(:)  
  character(26) :: filename
  character(3)  :: Tbmax_string
  integer :: ios, ipart
  logical :: test_exist 
  integer :: nss(4),outofrange, i, j, jc, irmin, irmax, numpart0 
  integer (kind=4) :: lhead, infmt, mode, naktiv, idx_orgn_s, &
                     numpart_s
  character(10) :: fulldate_cloudtop_string 
  real (kind=8) :: real_date 
  
  ! Switch to the multi case if required
  if (nb_halfdays_in_batch>0) then
     call fixparticlesMultiTBStratoClim(error)
     return
  endif                 
  
  ! --- Init ---
  error =.false.
  nss(1:4)=0
  outofrange=-2**30+1

  !------ Define ---
  write(Tbmax_string,'(i3)') TB_max

  write(fulldate_cloudtop_string, '(i10)') fulldate_cloudtop
  filename = trim(fulldate_cloudtop_string)//'_TB'//Tbmax_string

  ! --------------------------------------------
  ! -----       Reading the input file       ---
  ! --------------------------------------------

  write(*,*) trim(TBStratoClim_dir)//filename

  ! Check : file existence
  inquire(file=trim(TBStratoClim_dir)//filename, exist=test_exist)
  if(.NOT. test_exist) go to 9982

  ! Read the file
  open(UNIT=unitStratoClim, FILE=trim(TBStratoClim_dir)//filename, & 
       FORM='unformatted', ACCESS="sequential", &
       ACTION="read", POSITION="rewind", &
       IOSTAT=ios)
  if(ios /=0) go to 9962

  ! Read the number of cloud top :
  ! TODO ACCOUNT THAT THIS CAN BE COMBINED WITH PREVIOUS LIST OF PARCELS
  ! start to use first version index files
  read(unit=unitStratoClim)lhead
  if(lhead>3) then
    mode=1000
    print *,'temporary start'
    numpart=lhead
    idx_orgn=1
    stamp_date=fulldate_cloudtop
    stamp_date=10000*stamp_date  ! Achtung: cannot be done in one line
    ! test
    !write(*,'(i10,i16)')fulldate_cloudtop,stamp_date
  else
  rewind(unit=unitStratoClim)      
  
  read(UNIT=unitStratoClim, IOSTAT=ios)lhead,infmt,mode
  if ((infmt /=107).or.(mode/=0)) then
    print *,'fixparticlesTBStratoClim ERROR'
    print *,'wrong input file infmt mode ',infmt,mode
    stop
  endif
  ! ignores the time step
  read(UNIT=unitStratoClim, IOSTAT=ios) stamp_date 
  if(ios /=0) go to 9972
  ! check that the date is the same as fulldate_cloudtop
  real_date=dble(stamp_date)/10000._dbl+0.1_dbl
  if (floor(real_date)/=fulldate_cloudtop) then
    print *,'fixparticlesTBStratoClim ERROR'
    print *,'stamp_date does not match file name'
    print *,stamp_date,fulldate_cloudtop
    stop
  endif   
  read(UNIT=unitStratoClim, IOSTAT=ios) numpart_s,naktiv,idx_orgn_s
  if(ios /=0) go to 9972
  numpart=numpart_s
  ! overwrite the value read (it is 0 instead of 1)
  idx_orgn=1
  ! The three last words of the record are skipped 
  ! (nact_last0, nact_lastNM, nact_lastNH)    
  
  endif ! stop temporary start here
  
  print *, 'stamp_date numpart ',stamp_date,numpart
  ! do not pay attention to naktiv

  ! Allocation for the other variables
  ! Achtung: flagTB is global and permanent
  allocate(flagTB0(numpart), ir_start(numpart))
  allocate(xpart(numpart), ypart(numpart))
  allocate(presspart(numpart), temppart(numpart))
  
  ! Read the flag and the launching time
  read(UNIT=unitStratoClim, IOSTAT=ios) flagTB0  
  if(ios /=0) go to 9972
  read(UNIT=unitStratoClim, IOSTAT=ios) ir_start
  if(ios /=0) go to 9972
  ! Read the other variables :
  ! longitude, latitude, pressure and temperature
  ! of the cloud tops :
  read(UNIT=unitStratoClim, IOSTAT=ios) xpart
  if(ios /=0) go to 9972
  read(UNIT=unitStratoClim, IOSTAT=ios) ypart
  if(ios /=0) go to 9972
  read(UNIT=unitStratoClim, IOSTAT=ios) presspart
  if(ios /=0) go to 9972
  read(UNIT=unitStratoClim, IOSTAT=ios) temppart
  if(ios /=0) go to 9972
  ! Don't read the index which contains no information
  ! in a new parcels file  
  
  if(mode==1000) flagTB(:)=0
  
  close(UNIT=unitStratoClim)
  
  ! Computing the time to release the particle•
  !do ipart = 1,numpart  ! in forward calculation
  !        start_date = ibyeartrac*10**(4) + ibmonthtrac*10**(2)  ibdaytrac
  !        start_time = ibhourtrac*10**(4)
  !        ir_start(ipart) = int((juldate(start_date, start_time)- bdate) * 86400.)
  !enddo
  
  print *,'fixparticlesTBStrat:oClim stat lect'
  print *,xlon0, ylat0, xlon0+(nx-1)*dx, ylat0+(ny-1)*dy
  print *,minval(xpart),minval(ypart),maxval(xpart),maxval(ypart)
  print *,maxval(presspart),minval(presspart)
  allocate (thetapart(numpart))
  thetapart=temppart*(p0/presspart)**kappa
  print *,minval(thetapart),maxval(thetapart),sum(thetapart)/numpart

  ! Longitude and discretization : (180 to 179 -> 0 to 360)
  xpart = (xpart-xlon0)/dx

  ! Discretization in latitude :
  ypart = (ypart-ylat0)/dy

  if (M10) then
    ! Quick 'n dirty but easy version
    ! Determines the number of parcels
    irmin=minval(ir_start)
    irmax=maxval(ir_start)
    numpart0=numpart
    numpart=0
    do i=1,numpart0
      if(ir_start(i)==irmin) then
        numpart=numpart+4
      else if (ir_start(i)==irmax) then
        numpart=numpart+3
      else
        numpart=numpart+6
      endif
    enddo
    ! Allocate parcel fields
    allocate (xtra1(numpart),ytra1(numpart))
    allocate (itra0(numpart),itra1(numpart))
    allocate (ztra1(numpart),ttra1(numpart))
    allocate (flagTB(numpart),idxcross(numpart))
    ! Copy in the final arrays with duplication over 3h
    ! every 30' (from 1h before to 1h30 after) but for 
    ! the first and the last launch time of the sequence
    ! (this is intended to apply in the forecast mode) 
    ! idxcross contains the index in the origin file
    ! as each point is duplicated in time
    j=1 
    do i=1,numpart0
      if(ir_start(i)==irmin) then
        jc=3
        itra1(j)=ir_start(i)
        itra1(j+1)=ir_start(i)+1800
        itra1(j+2)=ir_start(i)+3600
        itra1(j+3)=ir_start(i)+5400
      else if (ir_start(i)==irmax) then
        jc=2
        itra1(j)=ir_start(i)-3600
        itra1(j+1)=ir_start(i)-1800
        itra1(j+2)=ir_start(i)
      else
        jc=5
        itra1(j)=ir_start(i)-3600
        itra1(j+1)=ir_start(i)-1800
        itra1(j+2)=ir_start(i)
        itra1(j+3)=ir_start(i)+1800
        itra1(j+4)=ir_start(i)+3600
        itra1(j+5)=ir_start(i)+5400     
      endif
      xtra1(j:j+jc)=xpart(i)
      ytra1(j:j+jc)=ypart(i)
      itra0(j:j+jc)=itra1(j:j+jc)
      ttra1(j:j+jc)=temppart(i)
      idxcross(j:j+jc)=i
      ztra1(j:j+jc)=log(p0/presspart(i)) ! no choice offered
      flagTB(j:j+jc)=ibclr(flagTB0(i),6)
      j=j+jc+1
    enddo
    write(*,'(A,3i10)')' M10 fixparticlesTBStratoClim numpart j numpart0',numpart,j,numpart0  
  else 
    ! Allocate parcel fields
    allocate (xtra1(numpart),ytra1(numpart))
    allocate (itra0(numpart),itra1(numpart))
    allocate (ztra1(numpart),ttra1(numpart))
    allocate (flagTB(numpart)) 
    xtra1(:)=xpart(:)
    ytra1(:)=ypart(:)
    itra0(:)=ir_start(:)
    ttra1(:)=temppart(:)
    itra1(:)=ir_start(:)
    do ipart = 1,numpart
      if (xtra1(ipart).lt.nearest(0.,-1.)) then 
        nss(1)=nss(1)+1
        itra1(ipart)=outofrange
      else if (xtra1(ipart).ge.nearest(float(nx-1),1.)) then
        nss(3)=nss(3)+1
        itra1(ipart)=outofrange
      else if (ytra1(ipart).lt.0.) then
        nss(4)=nss(4)+1
        itra1(ipart)=outofrange
      else if (ytra1(ipart).gt.float(ny-1)) then
        nss(2)=nss(2)+1
        itra1(ipart)=outofrange
      endif
    enddo
    write(*,'(A,i10,4i8)')' fixparticlesTBStratoClim numpart elim ',numpart,nss
    if(diabatic_StratoClim) then
      ztra1(:) = temppart(:)*(p0/presspart(:))**kappa
      flagTB(:)=ibset(flagTB0(:),6)
    else
      ztra1(:) = log(p0/presspart(:))
      flagTB(:)=ibclr(flagTB0(:),6)
    endif
  endif
  
  ! set new parcel and relative time flags (5 and 6)
  ! bits are numbered from zero, not one!!!
  ! This overwrite the value read  
  flagTB(:)=ibset(flagTB(:),4)
  flagTB(:)=ibset(flagTB(:),5)
   
  ! Deallocate :
  deallocate(ir_start)
  deallocate(xpart, ypart, presspart, temppart, flagTB0)
    
  return ! Main subroutine : fixparticlesTBStratoClim

  ! -------------------------------------------------
  ! --- Error display for fixparticlesTBStratoClim---
  ! -------------------------------------------------

9962 error=.true.
  write(*,*) '#####################################################'
  write(*,*) '#### TRACZILLA SUBROUTINE FIXPARTICLESTBSTRATOCLIM ##' 
  write(*,*) '####                                             ####'
  write(*,*) '#### ERROR - FILE CANNOT BE OPENED               ####'
  !write(*,*)
  !'####'//trim(TBClaus_dir)//yearstring_short//'/'//filename//'
  !####'
  write(*,*) '#####################################################'
  return

9972 error=.true.
  write(*,*) '#####################################################'
  write(*,*) '#### TRACZILLA SUBROUTINE FIXPARTICLESTBSTRATOCLIM ##' 
  write(*,*) '####                                             ####'
  write(*,*) '#### ERROR WHILE READING THE FILE                ####'
  ! write(*,*) '####
  ! '//trim(TBClaus_dir)//yearstring_short//'/'//filename//'
  ! ####'
  write(*,*) '#####################################################'
  return

9982 error=.true.
  write(*,*) '#####################################################'
  write(*,*) '#### TRACZILLA SUBROUTINE FIXPARTICLESTBSTRATOCLIM ##' 
  write(*,*) '####                                             ####'
  write(*,*) '#### ERROR THIS FILE DOESN''T EXIST :            ####'
  !write(*,*) '####
  !'//trim(TBClaus_dir)//yearstring_short//'/'//filename//'
  !####'
  write(*,*) '#####################################################'
  return
  end subroutine fixparticlesTBStratoClim
  
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@ fixparticlesMultiTBStratoClim @@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

  subroutine fixparticlesMultiTBStratoClim(error)
  
! This routine releases the parcel defined by the Stratoclim processing 
! of cloud tops from geostationry satellites written in 107 format
! This is the multi file version which is well adapted to the 
! reprocessing mode
! November 2016: B. Legras
      
  use commons
  use mass_iso
  use date

  implicit none

  !----- Declarations ---

  ! Input :
  !integer :: ibyeartrac, ibmonthtrac, ibdaytrac
  !integer :: ibhourtrac, ibmintrac, ibsectrac
  !integer :: Tbmax   ! brightness temperature threshold

  ! Output :
  logical :: error

  ! Read :
  !integer :: numpart

  ! Only here :
  integer (kind=4), allocatable :: ir_start(:), flagTB0(:)
  real (sp), allocatable :: xpart(:), ypart(:), presspart(:), temppart(:)
  real, allocatable :: thetapart(:)  
  character(26) :: filename
  character(3)  :: Tbmax_string
  integer :: ios, ipart
  logical :: test_exist 
  integer :: nss(4),outofrange, i, j, jc, irmin, irmax, numpart0 
  integer (kind=4) :: lhead, infmt, mode, naktiv, idx_orgn_s, &
                     numpart_s, numpart_s0
  character(10) :: fulldate_cloudtop_string 
  real (kind=8) :: real_date, &
                   julian_date, julian_date_beg, julian_date_end
  integer :: current_fulldate_cloudtop
  integer (kind=4) :: yyyymmdd, hh, hhmmss, offs, nbread, numx               
  
  ! --- Init ---
  error =.false.
  nss(1:4)=0
  outofrange=-2**30+1

  !------ Define ---
  write(Tbmax_string,'(i3)') TB_max
  
  ! --------------------------------------------
  !         First reading loop to check that
  !         all files are present and 
  !         determines the number of particles
  !         which are realeased
  ! --------------------------------------------
  
  ! fulldate is an integer input date in the format YYYYMMDDHH
  ! we need to extract YYYYMMDD and HH
  ! Calculate the Julian date of the first file and the Julian date
  ! of the last file 
  
  yyyymmdd = fulldate_cloudtop/100
  hh = fulldate_cloudtop - 100*yyyymmdd
  julian_date_beg = juldate(yyyymmdd, 10000*hh)
  julian_date = julian_date_beg
  julian_date_end = julian_date_beg + 0.5_dp*(nb_halfdays_in_batch-1)
  ! block counter
  nbread = 0
  ! parcel counter
  numpart = 0
  
  write(*,*)'fixparticlesmultiTBStratoClim begins'
  
  ! Generates a reference file to facilitate exploitation
  open(UNIT=unitpartout2, FILE=trim(path(2))//'multiparams', &
       STATUS='REPLACE',FORM='FORMATTED',ACTION='WRITE')
  write(unitpartout2,*) trim(TBStratoClim_dir)   
  
  ! First read to determine the size  
  do while (julian_date .le. julian_date_end)
     nbread = nbread+1
     call caldate(julian_date, yyyymmdd, hhmmss)
     write(fulldate_cloudtop_string, '(i10)') 100*yyyymmdd+hhmmss/10000
     filename = trim(fulldate_cloudtop_string)//'_TB'//Tbmax_string
     inquire(file=trim(TBStratoClim_dir)//filename, exist=test_exist)
     if(.not. test_exist) then
        write(*,*) 'fixparticlesmultiTBStratoClim: Missing file '
        write(*,*) trim(TBStratoClim_dir)//filename
        error = .true.
        return
     endif   
     open(UNIT=unitStratoClim, FILE=trim(TBStratoClim_dir)//filename, & 
       FORM='unformatted', ACCESS="sequential", &
       ACTION="read", POSITION="rewind", IOSTAT=ios)
     if(ios /=0) then
        write(*,*) 'fixparticlesmultiTBStratoClim: File cannot be opened'
        write(*,*) trim(TBStratoClim_dir)//filename
        error = .true.
        return
     endif  
     read(UNIT=unitStratoClim, IOSTAT=ios)lhead,infmt,mode
     if ((infmt /=107).or.(mode/=0)) then
        write(*,*) 'fixparticlesmultiTBStratoClim: wrong input file infmt mode'
        write(*,*) trim(TBStratoClim_dir)//filename
        write(*,*) 'fmt ',infmt,' mode ',mode
        error = .true.
        return
     endif
     if (nbread==1) then
        read(UNIT=unitStratoClim) stamp_date 
        ! check that the date is the same as fulldate_cloudtop
        real_date=dble(stamp_date)/10000._dbl+0.1_dbl
        if (floor(real_date)/=fulldate_cloudtop) then
           write(*,*) 'fixparticlesmultiTBStratoClim ERROR'
           write(*,*) 'stamp_date does not match name for first file'
           write(*,*) stamp_date,fulldate_cloudtop
           error=.true.
           return
        endif
        write(unitpartout2,*)'stamp date'
        write(unitpartout2,'(I8)')stamp_date
        write(unitpartout2,*)'nb halfdays in batch'
        write(unitpartout2,'(I3)')nb_halfdays_in_batch           
     else
        read(UNIT=unitStratoClim)
     endif  
     read(UNIT=unitStratoClim) numpart_s
     if (M10) numpart_s=6*numpart_s
     numpart = numpart+numpart_s
     write(*,*) 'filename numpart_s cumul ', trim(filename), numpart_s, numpart
     write(unitpartout2,*) '    filename'
     write(unitpartout2,'(4X,A)') trim(filename)
     write(unitpartout2,*) '    numpart_s cumul'
     write(unitpartout2,'(4X,I8,I11)') numpart_s, numpart 
     
     close(UNIT=unitStratoClim)
     julian_date = julian_date + 0.5 
  enddo ! end of the first loop on files
  write(*,*)'All files in batch checked'
  write(*,*)'stamp_date numpart', stamp_date, numpart
  write(unitpartout2,*) 'numpart'
  write(unitpartout2,'(I9)') numpart
  close(unitpartout2)
  ! No longer needed as all fields are dynamically allocated 
  !if (numpart > maxpart) then
  !    write(*,*) 'numpart > maxpart ',numpart, maxpart
  !    error = .true.
  !    return
  !endif
      
  allocate (xtra1(numpart),ytra1(numpart))
  allocate (itra0(numpart),itra1(numpart))
  allocate (ztra1(numpart),ttra1(numpart))
  allocate (flagTB(numpart))
  if (M10) allocate (idxcross(numpart)) 
  
  ! setting of the index origin
  idx_orgn=1
  
  ! --------------------------------------------
  ! -----       MAIN LOOP ON FILES           ---
  ! --------------------------------------------
  julian_date = julian_date_beg
  ! index first particles
  j=1
  ! offset first block (file)
  offs=0
  ! parcel counter accumulation
  numx=0
  do while (julian_date .le. julian_date_end)
     call caldate(julian_date, yyyymmdd, hhmmss)
     write(fulldate_cloudtop_string, '(i10)') 100*yyyymmdd+hhmmss/10000
     filename = trim(fulldate_cloudtop_string)//'_TB'//Tbmax_string
     open(UNIT=unitStratoClim, FILE=trim(TBStratoClim_dir)//filename, & 
       FORM='unformatted', ACCESS="sequential", &
       ACTION="read", POSITION="rewind", IOSTAT=ios)
     ! --------------------------------------------
     ! -----       Reading the input file       ---
     ! --------------------------------------------
     write(*,*) trim(TBStratoClim_dir)//filename
     ! skip headers
     read(UNIT=unitStratoClim) ! zap fmt+mode
     read(UNIT=unitStratoClim) ! zap stamp_date
     read(UNIT=unitStratoClim) numpart_s
     ! Allocation of temporary variables
     allocate(flagTB0(numpart_s), ir_start(numpart_s))
     allocate(xpart(numpart_s), ypart(numpart_s))
     allocate(presspart(numpart_s), temppart(numpart_s))
     read(UNIT=unitStratoClim, IOSTAT=ios) flagTB0
     read(UNIT=unitStratoClim, IOSTAT=ios) ir_start
     read(UNIT=unitStratoClim, IOSTAT=ios) xpart
     read(UNIT=unitStratoClim, IOSTAT=ios) ypart
     read(UNIT=unitStratoClim, IOSTAT=ios) presspart
     read(UNIT=unitStratoClim, IOSTAT=ios) temppart
     close(UNIT=unitStratoClim)
  
     print *,'fixparticlesMultiTBStratoClim stat lect'
     print *,xlon0, ylat0, xlon0+(nx-1)*dx, ylat0+(ny-1)*dy
     print *,'position',minval(xpart),minval(ypart),maxval(xpart),maxval(ypart)
     print *,'pressure',maxval(presspart),minval(presspart)
     allocate (thetapart(numpart_s))
     thetapart=temppart*(p0/presspart)**kappa
     print *,'theta---',minval(thetapart),maxval(thetapart),sum(thetapart)/numpart_s

     ! Longitude and discretization : (180 to 179 -> 0 to 360)
     xpart = (xpart-xlon0)/dx

     !  Discretization in latitude :
     ypart = (ypart-ylat0)/dy
     
     ! --------------------------------------------
     ! -----        Allocate the data           ---
     ! --------------------------------------------

     if (M10) then
        ! this mode never tested 
        irmin=minval(ir_start)
        irmax=maxval(ir_start)
        ! Copy in the final arrays with duplication over 3h
        ! every 30' (from 1h before to 1h30 after) including 
        ! the first and the last time of the sequence
        ! Beware: this differs from the single block case
        ! and implies to shift itime0 backward by 3h
        ! with respect to stamp_time (if M10 and 
        ! fixparticlesTBStratoClim(error)>0).
        do i=1,numpart_s
           itra1(j)=ir_start(i)-3600+offs
           itra1(j+1)=ir_start(i)-1800+offs
           itra1(j+2)=ir_start(i)+offs
           itra1(j+3)=ir_start(i)+1800+offs
           itra1(j+4)=ir_start(i)+3600+offs
           itra1(j+5)=ir_start(i)+5400+offs     
           xtra1(j:j+5)=xpart(i)
           ytra1(j:j+5)=ypart(i)
           itra0(j:j+5)=itra1(j:j+5)
           ttra1(j:j+5)=temppart(i)
           idxcross(j:j+5)=i+numx
           ! the coordinate is changed if required at the end of the routine 
           ztra1(j:j+5)=log(p0/presspart(i))
           flagTB(j:j+5)=flagTB0(i)
           j=j+6
        enddo
        write(*,'(A,3i10)')' M10 fixparticlesMultiTBStratoClim numpart j ',numpart_s,j-1  
     else 
        xtra1(j:j+numpart_s-1)=xpart(:)
        ytra1(j:j+numpart_s-1)=ypart(:)
        itra0(j:j+numpart_s-1)=ir_start(:)+offs
        ttra1(j:j+numpart_s-1)=temppart(:)
        itra1(j:j+numpart_s-1)=ir_start(:)+offs
        flagTB(j:j+numpart_s-1)=flagTB0(:)
        ! the coordinate is changed if required at the end of the routine 
        ztra1(j:j+numpart_s-1)=log(p0/presspart(:))
        ! Check for outlier (should not happen).
        ! This mechanism can be used as a quick way to input data prepaired for
        ! a bigger box but this is a waste of memory.
        do ipart = j,j+numpart_s-1 
           if (xtra1(ipart).lt.nearest(0.,-1.)) then 
              nss(1)=nss(1)+1
              itra1(ipart)=outofrange
           else if (xtra1(ipart).ge.nearest(float(nx-1),1.)) then
              nss(3)=nss(3)+1
              itra1(ipart)=outofrange
           else if (ytra1(ipart).lt.0.) then
              nss(4)=nss(4)+1
              itra1(ipart)=outofrange
           else if (ytra1(ipart).gt.float(ny-1)) then
              nss(2)=nss(2)+1
              itra1(ipart)=outofrange
           endif     
        enddo
        j=j+numpart_s
        write(*,'(A,2i11,4i4)')' fixparticlesMultiTBStratoClim numpart j elim ',numpart,j-1,nss              
     endif
     deallocate(flagTB0,xpart,ypart,ir_start,temppart,presspart)
     deallocate(thetapart)
     ! Update the time offset by adding 12 hours for next block 
     ! So we make the assumption here that date are available every 12h 
     offs = offs + 12*3600
     numx=numx+numpart_s
     julian_date = julian_date + 0.5
  enddo
  
  ! Switch the vertical coordinate to potential temperature if required
  if(diabatic_StratoClim) then
    ztra1(:) = ttra1(:)*exp(ztra1(:)*kappa)
    flagTB(:)=ibset(flagTB(:),6)
  else
    flagTB(:)=ibclr(flagTB(:),6)
  endif
      
  ! Set new parcel and relative time flags (5 and 6)
  ! bits are numbered from zero, not one!!!
  ! This overwrites the read value 
  flagTB(:)=ibset(flagTB(:),4)
  flagTB(:)=ibset(flagTB(:),5)
     
  return ! Main subroutine : fixparticlesTBStratoClim  
  
  end subroutine fixparticlesMultiTBStratoClim   

!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
  end module demar
