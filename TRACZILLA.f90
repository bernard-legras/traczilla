
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
!
!=======================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TRACZILLA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7==
!
      program traczilla
!*******************************************************************************
!                                                                              *
!     This is the Lagrangian Particle Dispersion Model TRACZILLA
!     derived from FLEXPART
!     The main program manages the reading of model run specifications, and 
!     initializations
!     All actual computing is done within subroutine timemanager.              *
!                                                                              *
!     Authors: A. Stohl, B. Legras                                                         *
!                                                                              *                                             *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! error                .true., if error ocurred in subprogram, else .false.    *
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!*******************************************************************************
       
      use commons
!      use lyapunov
      use isentrop_h
      use ecmwf_diab
      use mass_iso
      use ecmwf_inct
      use merra
      use jra55
      use thermo
      use advect
      use date
      use demar
      implicit none

      integer :: i
      logical :: error,oronew

      character(len=256):: pathfile
      
      integer :: OMP_GET_NUM_PROCS, OMP_GET_NUM_THREADS
      
      integer :: iargc
    
! Determines the number of threads for parallel run
!**************************************************
!$OMP PARALLEL
!$OMP MASTER
        num_threads=OMP_GET_NUM_THREADS()
        print *,'NUM_THREADS ',num_threads
!$OMP END MASTER
!$OMP END PARALLEL

        print *,'NUM_PROCS ',OMP_GET_NUM_PROCS()

! Read the pathnames where input/output files are stored
!*******************************************************

      select case (iargc())
      case (1)
        call getarg(1,pathfile)
      case (0)
        print *,'*** PATHNAME as argument is required ***'
        stop
      end select

      call readpaths(pathfile,error)
      if (error) goto 999

! Read the user specifications for the current model run
!*******************************************************
! based on backward compatible namelists unlike FLEXPART
! defines the release plan
      call readcommandB(error)
      
      if (error) goto 999
       
! Read the age classes to be used
!********************************
! deactivated feature
!      call readageclasses(error)
!      if (error) goto 999
       
! Read, which wind fields are available within the modelling period
!******************************************************************
! This is the main legacy from FLEXPART extended to include 
! other list descriptions of diabatic heating
! Notice that iso_mass runs do not use this procedure (not a good idea)

      call readavailable(error)
      if(error) goto 999
      if(diabatic_w .and. & 
        (ecmwf_diabatic_w.or.merra_diab.or.jra55_diab) ) then
        call readpaths_diab(pathfile,error)
        call readavailable_diab(error)
        if (error) goto 999
      endif
      if(diabatic_w .and. ecmwf_inct_w) then
        call readpaths_inct(pathfile,error)
        call readavailable_inct(error)
        if (error) goto 999
      endif
      
! Read the model grid specifications of the current model run
!************************************************************
! distinguish the type of wind format
      if(iso_mass) then
! binary input format made by genesis
        call gridcheck_iso(error)
      else if (merra_data) then
        call gridcheck_merra(error)
      else if (jra55_data) then
        call gridcheck_jra55(error)
      else
! standard FLEXPART wind grib files
        call gridcheck(oronew,error)
! ADD HERE A CALL TO gridcheck_merra
      endif
      if (error) goto 999

! Set the coefficients for mass correction whenever required
!***********************************************************
      if (ecmwf_diabatic_w.and.mass_correction) call diab_mass_init
      if (ecmwf_inct_w.and.mass_correction) call diab_mass_inct_init
      if (merra_diab.and.mass_correction) call diab_mass_merra_init
      if (jra55_diab.and.mass_correction) call diab_mass_jra55_init

! Read the parameters of the release plan defined in COMMAND
!***********************************************************
! Release plan is selected in readreleasesB2
! Initialize the random_seed if needed
      print *,'call readreleasesB2'
      call readreleasesB2(error)
      if (error) goto 999
      if(nsample(1) > 1) then
        print *,'randomly perturbed experiment'
        call random_seed()
      endif
      
! Allocate the 3D fields
!***********************
      call alloc_3D      
      

! Fix the release times and release positions of all particles
!*************************************************************
! delayed_initialization is required for fields which need 
! a preliminary acquisation of some fields (typically T)
! The delayed initialization is then done in timemanager
! TO DO: modify this to use the fields read by gridcheck
! and discard this option which complicates the beginning of the run
      if(restart) then
! call to restart routine, the branching between a restart from a part file 
! or from a full save is done within restartpart
        call restartpart(error)
        loffset=itime0+loutstep
!      Depreciated as no old initialization still in use
!      else if (releases_old) then
!        itime0=0
!        call coordtrafo(error) ! conversion of starting coordinates
!        if (error) goto 999
!        call fixparticlesB(error)
      else
        itime0=0
        select case (release_plan)
        case('CO2','TTL','AGE')
           if(delayed_initialization) then
             print *,'flexpart> delayed initialization'
           else
             call fixlayerpart(error)
           endif
        case('o3sonde')
!          conversion of initial coordinates
           call coordtrafo(error)
           if (error) goto 999
           call fixparticlesB(error)
        case('AGEF')
          if(external_pos0) then
            call setpos0fromext
          else
            call fixmultilayer(error)
          endif
        case('AGEB')
           if(jra55_data) then
             call fixmultilayer_tropomask_temp(error)
           else
             call fixmultilayer_tropomask(error)
           endif
           if (track_kill) then
              allocate(x_kill(numpart),y_kill(numpart),z_kill(numpart))
              allocate(it_kill(numpart),nstop_kill(numpart))
              it_kill(:)=ZERO_INIT
              nstop_kill(:)=ZERO_INIT
              x_kill(:)=MISSING
              y_kill(:)=MISSING
              z_kill(:)=MISSING
           endif
           if (track_cross) then
              allocate(it_1800(numpart),it_2300(numpart))
              it_1800(:)=ZERO_INIT
              it_2300(:)=ZERO_INIT
           endif
        case('MOZAIC')
           if(interp_release) then
              call fixmozaicinterp(error)
           else
              call fixmozaicpart(error)
           endif
        case('ER2')
           select case (campaign)
           case('SOLVE')
              call fixER2SOLVEinterp(error)
           case('SPADE')
              call fixER2SPADEinterp(error)
           case('AASEII')
              call fixER2SPADEinterp(error)
           end select
!        case('ER2lyapou')
!           print *,'FLEXPART> initialize Lyapunov'
!           call fixER2SOLVEinterp(error)
!           call fix_perturb
!           call alloc_Lyapunov
        case('CLAUS')
            if(delayed_initialization) then
              print *,'flexpart> delayed initialization'
           else
              call fixparticlesClaus(error)
           endif
        end select
      endif
      if (error) goto 999
      
!     tell the interpolator to calculate temperature with best accuracy
      AccurateTemp=.true.

! Put longitudes within the mother domain by shifting
! *****************************************************
! In principle, the initializations are already doing
! this step
      if (xglobal) then
        do i=1, numpart
          xtra1(i) = modulo(xtra1(i), nearest(float(nx-1),-1.))
        enddo
      endif

! Calculate particle trajectories
!********************************
! This does the real job of integrating the parcels in time      
      call timemanagerB()

      write(*,*) 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED' 
      write(*,*) 'A TRACZILLA MODEL RUN!'

      goto 1000

999   write(*,*) 'TRACZILLA MODEL ERROR: EXECUTION HAD TO BE TERMINATED'

1000  continue

      contains

      subroutine alloc_3D

      ! always
      if(merra_data) then
        call alloc_merra
      else if (jra55_data) then
        call alloc_jra55
      else
         print *,'alloc_3D uuh vvh tth uupol vvpol ps tt2 u10 v10'
         allocate (uupol(0:nx-1,0:ny-1,nuvz,2),vvpol(0:nx-1,0:ny-1,nuvz,2))
         allocate (uuh(0:nx-1,0:ny-1,nuvz,2),vvh(0:nx-1,0:ny-1,nuvz,2))
         allocate (tth(0:nx-1,0:ny-1,nuvz,2))
         allocate (ps(0:nx-1,0:ny-1,1,2),tt2(0:nx-1,0:ny-1,1,2))
         allocate (u10(0:nx-1,0:ny-1,1,2),v10(0:nx-1,0:ny-1,1,2))
         ! if needed according to the release plan
         if(TTLactiv .or. CLAUSactiv) then
            print *,'alloc_3D qvh'
            allocate (qvh(0:nx-1,0:ny-1,nuvz,2))
         endif
         if(diabatic_w .or. isentropic_motion) then
            call alloc_isentrop_perm
            if(diabatic_w .and. ecmwf_diabatic_w) &
               call alloc_ecmwf_diab
            if(diabatic_w .and. ecmwf_inct_w) &
               call alloc_ecmwf_inct
         elseif (z_motion) then
            print *,'alloc_3D wwh'
            allocate (wwh(0:nx-1,0:ny-1,nwz,2))
         elseif (mass_diabat) then
            call alloc_iso
         endif
      endif
      
      end subroutine alloc_3D

      end program traczilla
!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log: 
