!**********************************************************************
! Copyright 1996, i1997, 2001, 2002, 2006, 2007, 2012, 2013           *
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

!
!###############################################################################
!----------------------------------- ADVECT ------------------------------------
!###############################################################################

! This module of TRACZILLA contains the routines controlling time evolution

module advect
use commons
!use lyapunov, only: activ_Lyapunov, numpart_np, tau_orthog, &
!                          lambdaa, lambdab, &
!                          ortho, &
!                          mass_triad_b  !test----
use ecmwf_diab
use ecmwf_inct
use mass_iso
use merra
use jra55
use combin
use io
use demar
use thermo
use readinterp
implicit none
public :: timemanagerB
private :: advanceB

contains

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TIMEMANAGERB @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

      subroutine timemanagerB()                                  

!*******************************************************************************
!                                                                              *
! Handles the computation of trajectories, i.e. determines which               *
! trajectories have to be computed at what time.                               *
! Manages dry+wet deposition routines, radioactive decay and the computation   *
! of concentrations.                                                           *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     20 May 1996                                                              *
!                                                                              *
!*******************************************************************************
!  Changes
!           B. Legras, Apr. 2002:
!           simplified version for stratospheric studies
!           (...)
!           B. Legras, March 2006:
!           delayed initialization
!           B. Legras, September 2012
!           AGEF mode
!           // version
!           B. Legras 16/9/2012
!           stopped trajs removed from npproc
!           MERRA mode: B. Legras, March 2013
!           B. Legras, October 2014      
!           AGEB mode
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! ideltas [s]        modelling period                                          *
! itime [s]          actual temporal position of calculation                   *
! ldeltat [s]        time since computation of radioact. decay of depositions  *
! loutaver [s]       averaging period for concentration calculations           *
! loutend [s]        end of averaging for concentration calculations           *
! loutnext [s]       next time at which output fields shall be centered        *
! loutsample [s]     sampling interval for averaging of concentrations         *
! loutstart [s]      start of averaging for concentration calculations         *
! loutstep [s]       time interval for which concentrations shall be calculated*
! xtra1(maxpart), ytra1(maxpart), ztra1(maxpart) =                             *
!                    spatial positions of trajectories                         *
!                                                                              *
! Constants:                                                                   *
! maxpart            maximum number of trajectories                            *
!                                                                              *
!*******************************************************************************

      integer :: i,j,itime,nstop, ix, jy
      integer :: loutnext,loutnext0
      integer :: npproc,npstop(7),nb_bounce, nerr
      integer :: count_clock_1,count_clock_2,count_rate,count_max
      integer :: next_time_orthog, prev_time_orthog
      integer, allocatable :: shuffle(:)
      integer :: npproc_thread,nb_bounce_thread,npstop_thread(7),num_thread
      logical :: tropotest
      integer :: yyyymmdd,hhmmss,ytr,mtr,dtr,dyr,year_length
      integer :: ldayscum(12)
      character (4) :: yyyy
      character (len=256) :: tropofile
      real (dp), allocatable :: buff(:),tropozdd(:,:)
      real (dp) :: ddx,ddy
      real (dp) :: tint
      real (dp) :: dt
      real (dp) :: epsil=1.e-3
!     real lambdaa(3,maxlyap), lambdab(3,maxlyap)
!     real xm,xm1
      logical :: error,trace_time, tropopen
      logical, allocatable :: mask_part(:)
      integer :: n,nb_processors,OMP_GET_THREAD_NUM
!     external ctrper

      !test-----
      ! real cc
      !---------

      ldayscum=(/0,31,59,90,120,151,181,212,243,273,304,334/)

!------------------------------------
! Performs further initializations   
!------------------------------------
            
      if (shuffling) then
! Generate a random permutation among the numpart elements
! needs numpart to be defined (no delayed initialisation)
         print *,'shuffling over numpart ',numpart
         allocate(shuffle(numpart))
         call ranper(shuffle,.true.)
!        j=1
!        do n=1,num_threads
!           do i=n,numpart,num_threads
!              shuffle(j)=i
!              j=j+1
!           enddo
!        enddo       
!        print *,j-numpart-1,shuffle(1:5)
!        print *,minval(shuffle),maxval(shuffle)
      endif
      
      !print *,'jra logical variables ',jra55_data,jra55_diab
      
      if (TTLactiv.or.CLAUSactiv) then
! Allocate arrays for temperature and saturation mixing ratio
        print *,'timemanager> allocate ttra1 and qtra1'
        if (.not.allocated(ttra1)) then
          allocate(ttra1(maxpart))
          print *,'allocated ttra1'
        endif       
        if (.not.allocated(qtra1)) then 
          allocate (qtra1(maxpart))
          print *,'allocated qtra1'
        endif
! mixing ratio read if restart but here it is reinitalized
! mixing ratio initialized at the beginning, after each output or after restart
        qtra1=0.05_dp
      endif

      if (AGEBactiv.or.AGEFactiv.or.merra_diab.or.jra55_data) then
        if (.not.allocated(ttra1)) then 
          allocate(ttra1(maxpart))
          print *,'allocated ttra1'
        endif
      endif
      
      if(AGEBactiv) then
        allocate (buff(ny*nx),tropozdd(0:nx-1,0:ny-1))
      endif     

! Initialisation of the variables controling the output
      loutnext = loffset
!     loutnext0 allows first output out of sequence
!     and also to have output at time 0 if loffset2 = -loutnext
      loutnext0 = loutnext + loffset2
      
!! Lyapunov deactivated by default
!      if( activ_Lyapunov ) then
!        print *,'timemanager> activation of Lyapunov calculation'
!!       first orthogonalization performed at the same time as the first output
!        prev_time_orthog = 0
!        next_time_orthog = loutnext0
!        !test----
!        mass_triad_b(1:numpart_np) = exp(-ztra1(1:numpart_np))
!        !--------
!      endif
 
      flush 6

!***********************************************************************     
!     T H E   T I M E    L O O P
!
!**********************************************************************
! Loop over the whole modelling period in time steps of mintime seconds
!**********************************************************************

      npstop(:)=0
      npproc=numpart ! particle counter initialized to numpart (which may not
                     ! be defined at this stage)
      print *,'start time loop'  

      tropoyy=0
      trace_time=.false.
      dotime: do itime=itime0,ideltas,lsynctime
  
!       print *,'timemanager> itime, loutnext, npproc', itime, loutnext0, npproc 

! Get necessary wind fields if not available
!===========================================
        if(iso_mass) then
          call getfields_iso(itime,nerr)
        else if (merra_data) then     
          call getfields_merra(itime,nerr)
        else if (jra55_data) then
          call getfields_jra55(itime,nerr)
        else
          call getfields(itime,nerr)
        endif
        if (nerr > 1) stop 'NO METEO FIELDS AVAILABLE'
        if(diabatic_w .and. ecmwf_diabatic_w) then
          call getfields_diab(itime,nerr)
          if (nerr > 1) stop 'NO DIAB FIELDS AVAILABLE'
        endif
        if(diabatic_w .and. ecmwf_inct_w) then
          call getfields_inct(itime,nerr)
          if (nerr > 1) stop 'NO INCT FIELDS AVAILABLE'
        endif

! Switch off the diffusion if needed 
        if ((delay_switch_diff_off > 0).and. &
              (-itime > delay_switch_diff_off)) then
           switch_diff_off=.true.
           print *,'switch_diff_off changed to true'
        endif

! Delayed initialization if required
        if(itime==itime0.and.delayed_initialization.and.(.not.restart)) then
          if(press2theta) then
            call fixpresslayerpart(error)
            if(error) stop 'ERROR in fixpresslayerpart for' &
                           //' delayed initialization'
          else if(theta2press) then 
            stop 'THETA2PRESS: OPTION NOT IMPLEMENTED'
          else if(make_uni3D) then ! case of CO2 runs
            call fixlayerpart(error)
            if(error) stop 'ERROR in fixlayerpart for' &
                           //' delayed initialization'
          else if (CLAUSactiv) then
            call fixparticlesCLAUS(error)
            if(error) stop 'ERROR in fixparticlesCLAUS for' &
                           //' delayed initialization'
          else
            stop 'NO OPTION FOR DELAYED INITIALIZATION'
          endif
          if(TTLactiv.or.CLAUSactiv) then
            qtra1=0.05
          endif
        endif
        if(itime==itime0 .and. TTLactiv .and. .not.JRA55_data) then
           ! fake time step in order to calculate the temperature
           print *,'fake initial step to calculate T'
           do j=1,numpart
             call advanceB(j,itime,0.,nstop,xtra1(j),ytra1(j),ztra1(j),ttra1(j))
           enddo
        endif
        
! Calculate date and read a new tropopause file if needed
        if(AGEBactiv) then
          call caldate(bdate+itime/86400._dp,yyyymmdd,hhmmss)
          ytr=floor(yyyymmdd/10000._dp+epsil)
          if(.not.(ytr .eq. tropoyy)) then
            inquire(unit=tropunit,opened=tropopen)
            if (tropopen) close(tropunit)
            write(yyyy,'(I4)') ytr
            if (merra_data) then
               tropofile=trim(tropodir)//'tropo-theta-MERRA-'//yyyy//'.bin'
            else if (jra55_data) then
               tropofile=trim(tropodir)//'tropo-theta-JRA-'//yyyy//'.bin'
            else
               tropofile=trim(tropodir)//'tropo-theta-EI-'//yyyy//'.bin'          
            endif
            open(tropunit,file=tropofile,access='direct',status='old', &
              recl=4*ny*nx)
            print *,'opening tropopause file for ',yyyy
            tropoyy=ytr
            tropopen=.true.
          endif 
          mtr=floor(yyyymmdd/100._dp-100*ytr+epsil)
          dtr=yyyymmdd-10000*ytr-100*mtr
          dyr=dtr+ldayscum(mtr)
          if (mod(ytr,4)==0 .and. mtr>2) dyr=dyr+1
        endif       
 
! Set tropotest
        if(AGEBactiv) then
          tropotest=(mod(itime,loutprint)==0)
          if (tropotest) then
            read(tropunit,rec=dyr)buff
            tropozdd=reshape(buff,(/nx,ny/))
            write(*,'(" tropotest ",I4,2I2.2,I5)') ytr,mtr,dtr,dyr
          endif
        else
          tropotest=.false.
        endif

! Output of particle positions
        if((itime == loutnext0).and.(ipout == 1)) then
          if (TTLactiv.or.CLAUSactiv) then
            call partout_qv(itime)
       !     print *, &
       !       'output of particle positions, T and qv, format 103 ', &
       !       itime/86400,' days'
            qtra1= 0.05  ! reinitialize qtra1 after each output 
          else if (AGEFactiv) then
            call partout_agef(itime)
       !     print *,'output of particle positions and T, format 105 ', &
       !       itime/86400,' days'
          else if (AGEBactiv) then
            call partout_agef(itime)     
          else
            call partout_fast(itime)
       !     print *,'output of particle positions format 102 ', &
       !        itime/86400,' days' 
          endif  
          loutnext=loutnext+loutstep
          loutnext0 = loutnext
          if(savfull.or.(kind(1.0)==8)) then
            call savsav(itime)
            print *,'full backup of trajectories ',itime
          endif
          flush (6)
        endif
        
        if(mod(itime,loutsav)==0) then
          call savsav(itime)
          if (track_kill)  call savkill(itime)
          if (track_cross) call savcross(itime)
          flush (6)
        endif                           
 
! Calculate Lyapunov exponents: TO BE UPDATED WITH I. PISSO CODE
! deactivated by default
     
!        if( activ_Lyapunov .and. (itime == next_time_orthog) ) then
!           print *,'timemanager> Gram-Schmidt'
!           call ortho(itime)
!!	   if( prev_time_orthog == 0 ) then	  
!!	      do i = 1, numpart_np
!!	         lambda(:,i) = log / (itime-itramem(i))
!!             enddo              
!!	   else
!!	      do i = 1, numpart_np
!!	         lambda(:,i) = 0.5* (lambdaa(:,i)+lambdab(:,i)) / (itime - prev_time_orthog)
!!	      enddo
!!	   endif
!           !test----
!           do i=1,numpart_np
!             cc = 1./(itime - itramem(i))   
!             write(*,'(a,i8,i5,3G15.5)') ' Lyapunov> ',itime,i,  &
!                cc*lambdaa(1,i),cc*lambdaa(2,i),cc*lambdaa(3,i)
!             write(*,'(a,i8,i5,3G15.5)') ' Lyapunov> ',itime,i,  &
!                cc*lambdab(1,i),cc*lambdab(2,i),cc*lambdab(3,i)
!           enddo
!           mass_triad_b(1:numpart_np) = exp(-ztra1(1:numpart_np))
!           !--------
!           prev_time_orthog = itime
!           next_time_orthog = itime + tau_orthog
!           !call lambda_out(itime)
!        endif
   
        if (itime == ideltas) goto 99           ! do not advect beyond final time
        
 
!************************
! Loop over all particles
!************************
        npproc=0 ; nb_bounce=0
        call system_clock(count_clock_1,count_rate,count_max)
        if(trace_time) print *,'itime > ',itime

! Beginning of the parallel section
!$OMP PARALLEL DEFAULT(SHARED) SHARED(npproc,nb_bounce,npstop,shuffle,itime) &
!$OMP PRIVATE(nstop,npproc_thread,nb_bounce_thread,npstop_thread,num_thread)
        npproc_thread=0
        nb_bounce_thread=0
        npstop_thread(:)=0
        num_thread=OMP_GET_THREAD_NUM()
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,dt,tint,ix,jy,ddx,ddy)
        dopart: do i=1,numpart
          if (shuffling) then
             j=shuffle(i)
          else
             j=i
          endif

        ! nstop now controls if an individual trajectory needs to be terminated
        ! set to 0 in advance       

        ! TEST!
         !if(mod(j,833)==0 ) then
         !  debug_out=.true.
         !else 
         !  debug_out=.false.
         !endif

! If integration step is due, do it
!==================================
          ifdue: if ((abs(itime+lsynctime)> abs(itra1(j))) .and. &
              (abs(itra1(j))>= abs(itime))) then
            if(itra1(j) == itime) then
               dt = float(abs(lsynctime))
            else
               dt = float(abs(itime+lsynctime-itra1(j)))
            endif

! Integrate advection equation for lsynctime seconds
!===================================================
            if (merra_diab.or.jra55_data) tint=ttra1(j)
!           if(debug_out) &
!              print "(' timemanager> particle ',i6,2X,3G12.5)",&
!                 j,xtra1(j),ytra1(j),ztra1(j)
            call advanceB(j,itime,dt,nstop, &
              xtra1(j),ytra1(j),ztra1(j),tint)     
            npproc_thread=npproc_thread+1

! Test ------------------------------
! may generate large amount of output
!cc         if(nstop.gt.1) print *,'stop traj ',j,xtra1(j),ytra1(j),ztra1(j)
!------------------------------------

! Apply tropopause test when due to discard parcels having crossed
!=================================================================
! Apply simultaneously the duration test to discard parcels beyond maximum life
!============================================================================== 
! It is assumed that the tropopause data have same grid as the wind
! Management of the date is done locally without the weight of 
! the available/getfields mech
! This test has priority over those performed in advanceB
 
! This test is inexact at the pole but who cares?
           
            if (tropotest) then
              ix=floor(xtra1(j))
              jy=min(ny-2,max(0,floor(ytra1(j))))
              ddx= xtra1(j)-ix
              ddy= ytra1(j)-jy
              if(ztra1(j)<tropozdd(ix,jy)*(1-ddx)*(1-ddy) &
                         +tropozdd(ix+1,jy)*ddx*(1-ddy) &
                         +tropozdd(ix,jy+1)*(1-ddx)*ddy &
                         +tropozdd(ix+1,jy+1)*ddx*ddy) nstop=6
              if (itra0(j)-itime > max_life_time) nstop=7
              if (nstop >=6 .and. track_kill) then
                it_kill(j)=itime
                nstop_kill(j)=nstop
                x_kill(j)=xlon0+xtra1(j)*dx
                y_kill(j)=ylat0+ytra1(j)*dy
                z_kill(j)=ztra1(j)
              endif
              if (track_cross) then
                if(ztra1(j)>1800 .and. it_1800(j)==0) it_1800(j)=0
                if(ztra1(j)>2300 .and. it_2300(j)==0) it_2300(j)=0
              endif                       
            endif

! Calculate qv for final pressure and initial temperature
! Copy tint in ttra1
!========================================================
            if(TTLactiv) then
              ttra1(j)=tint
              qtra1(j)=min(qtra1(j),satratio(p0*exp(-ztra1(j)),tint))
            else if (AGEBactiv.or.AGEFactiv.or.merra_diab.or.jra55_data) then
              ttra1(j)=tint
            endif
           
! Determine, when next time step is due
! If trajectory is terminated, do not increment its time
!=======================================================

            if (nstop>1) then           
              npstop_thread(nstop)=npstop_thread(nstop)+1                                
!             TEST    
!              write(*,'(a,2i8,3f12.4)') ' timemanager > stopped parcel',&
!                 j,nstop,xtra1(j),ytra1(j),ztra1(j)
!             itra1(j)=-999999999
            elseif (nstop<0) then
              nb_bounce_thread=nb_bounce_thread+1
!             print *,'timemanager > bouncing parcel'
!             print *,j,nstop,xtra1(j),ytra1(j),ztra1(j)
              itra1(j)=itime+lsynctime
            else
              itra1(j)=itime+lsynctime
            endif

!          if(j>500000) print *,j
            
          endif ifdue
         
        enddo dopart                      ! end of loop over particles

!$OMP END DO

!$OMP CRITICAL
        npproc=npproc+npproc_thread
        nb_bounce=nb_bounce+nb_bounce_thread
        npstop(:)=npstop(:)+npstop_thread(:)
        npproc=npproc-sum(npstop_thread)
!$OMP END CRITICAL
!       print *,'THREAD ',num_thread,' completed'
!$OMP END PARALLEL

        if (mod(itime,loutprint)==0) then
          call system_clock(count_clock_2,count_rate,count_max)        
          write(*,'(a,i9,i6,f10.2,a)') ' timemanager> npproc*, bounce, cpu ',&
             npproc, nb_bounce, float(count_clock_2-count_clock_1)/       &
             float(count_rate), 's'
          write(*,'(a,7i8)') ' timemanager> npstop ',npstop
          allocate (mask_part(numpart))
          mask_part(:) = itra1(:)==itime+lsynctime
          print*,'ztra1 ', minval(ztra1(1:numpart),DIM=1,MASK=mask_part), &
                           maxval(ztra1(1:numpart),DIM=1,MASK=mask_part), &
                           minloc(ztra1(1:numpart),DIM=1,MASK=mask_part), &
                           maxloc(ztra1(1:numpart),DIM=1,MASK=mask_part)
          
          if(iso_mass) then                 
             mask_part(:) = mask_part(:) .and. (ztra1(:) .ge. ThetaLev(NbTheta))
             print*,'# hanged ',sum(transfer(mask_part,1,size=size(mask_part)))
          endif
          deallocate (mask_part)
          npstop(:)=0
        endif
!        if (nb_bounce >0) print *,'timemanager > nb bouncing ',nb_bounce
        
!        if(npstop.eq.1) stop
          
      enddo dotime                         ! end of loop over simulation time

      write(*,'(a,7i8)') 'timemanager > npstop ',npstop

99    continue

      if (ipout.eq.2) call partout_fast(itime)     ! dump particle positions
      if(shuffling) deallocate (shuffle)

      return
      end subroutine timemanagerB

!=======================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ADVANCEB @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7==

      subroutine advanceB(jp,itime,dt,nstop,x,y,z,tint)
            
!*******************************************************************************
!                                                                              *
!  Calculation of trajectories utilizing a zero-acceleration scheme.           *
!  The time step is determined by the Courant-Friedrichs-Lewy (CFL) criterion. *
!  This means that the time step must be so small that the displacement within *
!  this time step is smaller than 1 grid distance. Additionally, a temporal    *
!  CFL criterion is introduced: the time step must be smaller than the time    *
!  interval of the wind fields used for interpolation.                         *
!  For random walk simulations, these are the only time step criteria.         *
!  For the other options, the time step is also limited by the Lagrangian time *
!  scale.                                                                      *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     16 December 1997                                                         *
!     
!     Changes: B. Legras, April 2002
!            boundary layer calculation cancelled
!              simple advection with vertical noise
!              B. Legras, September 2007
!              option for horizontal diffusion
!              B. Legras, September 2012
!              tests for pcut and thetacut
!              MERRA mode, B. Legras, March 2013
!              test for thetalowcut
!              B. Legras, September 2013
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! h [m]              Mixing height                                             *
! ladvance [s]       Total integration time period                             *
! ngrid              index which grid is to be used                            *
! rannumb(maxrand)   normally distributed random variables                     *
! x,y,z              Next time step's spatial position of trajectory           *
!                                                                              *
!                                                                              *
!*******************************************************************************

      integer, parameter :: nss=50
      real, parameter :: pih=pi/180.,eps=1.e-4
      
      integer, intent(in):: jp,itime
      integer, intent(out):: nstop
      real(dp), intent(in) :: dt
      real(dp), intent(inout):: x,y,z
      real(dp), intent(inout):: tint

!      integer i,j,k,nrand,loop
!      integer memindnext,ixnext,jynext     
!      real rhog	
!      logical partial_step	
!      integer ix,jy,ixp,jyp

      integer :: ngrid, idummy
      real(dp) :: wdiff
      real(dp) :: xlon,ylat,xpol,ypol,gridsize
      real(dp) :: cosfact,psaver,theta_inf,theta_sup
      real(dp) :: u,v,w,dxsave,dysave,z_factor,rand_angle
      real(dp), dimension(nss) :: harvest
      save idummy

      data idummy/-7/

      dxsave=0.           ! reset position displacements
      dysave=0.           ! due to mean wind
      nstop=0

! Determine whether lat/long grid or polarstereographic projection
! is to be used
! Furthermore, determine which nesting level to be used
!=================================================================

      if (nglobal.and.(y.gt.switchnorthg)) then
        ngrid=-1
      else if (sglobal.and.(y.lt.switchsouthg)) then
        ngrid=-2
      else
        ngrid=0
      endif

! Interpolate the wind
!=====================

!      if (debug_out) print *,'advanceB call interpol'
      if(isentropic_motion) then
        call interpol_wind_theta (itime,x,y, z, & 
           u,v, ngrid, theta_inf, theta_sup, z_factor,tint)
!           if(debug_out) print *,'advanceB (u,v)> ',u,v             ! test
      elseif (diabatic_w .and. ecmwf_diabatic_w ) then
        if (ecmwf_inct_w) then
          call interpol_wind_theta_inct (itime,x,y, z, & 
            u,v,w, ngrid, theta_inf, theta_sup, z_factor,tint,nstop)
!            if(debug_out) print *,'advanceB (u,v,w)> ',u,v,w*86400  ! test
        else
          call interpol_wind_theta_diab (itime,x,y, z, & 
            u,v,w, ngrid, theta_inf, theta_sup, z_factor,tint,nstop)
!            if(debug_out) print *,'advanceB (u,v,w)> ',u,v,w*86400  ! test
        endif
      elseif (merra_data) then
        call interpol_wind_merra(itime,x,y,z,u,v,w,ngrid, &
           theta_inf,theta_sup,psaver,z_factor,tint,nstop)
        !if(debug_out) print *,'abvanceB (x,y,z)> ',x,y,z
        !if(debug_out) print *,'lat lon        > ',y*dy+ylat0,x*dx+xlon0
        !if(debug_out) print *,'(u,v,w,tint)    > ',u,v,w*86400,tint
      elseif (jra55_data) then
        call interpol_wind_jra55(itime,x,y,z,u,v,w,ngrid, &
           theta_inf,theta_sup,psaver,z_factor,tint,nstop)       
      elseif (z_motion) then
        call interpol_windB(itime,x,y,z,u,v,w,ngrid, &
           psaver,z_factor,tint)
        !if(debug_out) print *,'abvanceB (x,y,z)> ',x,y,z
        !if(debug_out) print *,'(u,v,w,tint)> ',u,v,w,tint        ! test
      elseif (iso_mass) then
        call interpol_wind_iso(itime, x, y, z, & 
           u, v, w, ngrid, z_factor, tint, nstop)
!           if(debug_out) print *,'advanceB (u,v,w)> ',u,v,w*86400  ! test
      endif

! Compute everything for above the PBL

! Calculate horizontal displacement (in m)
!=========================================
      dxsave=u*dt*float(ldirect)   
      dysave=v*dt*float(ldirect)
      
! Calculate vertical position at time step itime+lsynctime
!==========================================================
!     kinematic trajectories with w in Pa/s
      if(z_motion) z=z-w*dt*float(ldirect)/(p0*exp(-z))
      !if(z_motion.and.debug_out) then
      !  print*,p0*exp(-z),w*dt,z,psaver
      !endif
!     diabatic trajectories with w in K/s
      if(diabatic_w .and. (ecmwf_diabatic_w.or.merra_data.or.jra55_data)) &
         z=z+w*dt*float(ldirect)
      if(iso_mass .and. mass_diabat)  z=z+w*dt*float(ldirect)


! Add random vertical/horizontal component to the displacement
!=============================================================
       
      if((nsample(1) > 1).and.(.not.switch_diff_off)) then
        call random_number(harvest)
        wdiff=2*sum(harvest) - nss
        wdiff = wdiff*sqrt(6*diffus/(nss*dt))
        if(verdiff) z=z+wdiff*dt*float(ldirect)*z_factor
        if(hordiff) then
          call random_number(rand_angle)
          rand_angle=2*pi*rand_angle
          dxsave = dxsave + wdiff*dt*cos(rand_angle)
          dysave = dysave + wdiff*dt*sin(rand_angle)
        endif
      endif    
      
! Projects onto the frontier if going above or below
!====================================================
      if (z_motion) then 
         z=min(z,zmax)
         if (z.lt.log(p0/psaver)) z=2*log(p0/psaver)-z ! bouncing
      endif
      if (theta_bounds) then
         if (diabatic_w.and.(ecmwf_diabatic_w.or.merra_diab.or.jra55_diab)) &
            z=max(min(z,theta_sup),theta_inf)
         if (mass_diabat) &
            z=max(min(z,ThetaLev(NbTheta)),ThetaLev(1))
      endif
!     if (debug_out) print *,'advanceB z sup inf ',z,theta_sup,theta_inf
      
! Calculates new horizontal position in grid coordinates
!=======================================================
! Manage polar regions

      if (ngrid.ge.0) then
        cosfact=dxconst/cos((y*dy+ylat0)*pih)
        x=x+dxsave*cosfact
        y=y+dysave*dyconst
      else if (ngrid.eq.-1) then      ! around north pole
        xlon=xlon0+x*dx
        ylat=ylat0+y*dy
        call cll2xy(northpolemap,ylat,xlon,xpol,ypol)
        gridsize=1000.*cgszll(northpolemap,ylat,xlon)
        dxsave=dxsave/gridsize
        dysave=dysave/gridsize
        xpol=xpol+dxsave
        ypol=ypol+dysave
        call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)
        x=(xlon-xlon0)/dx
        y=(ylat-ylat0)/dy
      else if (ngrid.eq.-2) then    ! around south pole
        xlon=xlon0+x*dx
        ylat=ylat0+y*dy
        call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
        gridsize=1000.*cgszll(southpolemap,ylat,xlon)
        dxsave=dxsave/gridsize
        dysave=dysave/gridsize
        xpol=xpol+dxsave
        ypol=ypol+dysave
        call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
        x=(xlon-xlon0)/dx
        y=(ylat-ylat0)/dy
      endif


! If global data are available, use cyclic boundary condition
!------------------------------------------------------------
     
      if (xglobal) x = modulo(x,nearest(float(nx-1),-1.))
      
! Check position: If trajectory outside model domain, terminate it
! (should not happen)
!=================================================================

!     nstop must be less or equal to 7 or change bound in timemanager

      if(merra_data) then
          if ((x.lt.nearest(0.,-1.)).or.(x.ge.nearest(float(nx-1),1.)).or.(y.lt.nearest(-0.5,-1.)).or. &
	       (y.gt.nearest(float(ny)-0.5,1.)).or.(z_motion.and.(z>zmax))) then   !trajectory terminated
	        nstop=3
	        print *, 'stop parcel:trajectory terminated ',jp, ngrid
	        if (ngrid.ge.0) then
	          write(*,'(9G14.5)') x,y,z,dxsave,dysave,w,u,v,dt
	        else
	          write(*,'(6g14.5)') x,y,dxsave,dysave,xpol,ypol
	        endif
	      endif
	  else if (jra55_data)  then
	      ! (90-ylat0)/dy=0.7675033...     1-0.7675033...=0.23249669...
	      if ((x.lt.nearest(0.,-1.)).or.(x.ge.nearest(float(nx-1),1.)).or.(y.lt.-0.767504).or. &
	       (y.gt.float(ny)-0.232497).or.(z_motion.and.(z>zmax))) then   !trajectory terminated
	        nstop=3
	        print *, 'stop parcel:trajectory terminated ',jp, ngrid
	        if (ngrid.ge.0) then
	          write(*,'(9G14.5)') x,y,z,dxsave,dysave,w,u,v,dt
	        else
	          write(*,'(6g14.5)') x,y,dxsave,dysave,xpol,ypol
	        endif
	      endif
      else
          if ((x.lt.nearest(0.,-1.)).or.(x.ge.nearest(float(nx-1),1.)).or.(y.lt.0.).or. &
            (y.gt.float(ny-1)).or.(z_motion.and.(z>zmax))) then   !trajectory terminated
            nstop=3
            print *, 'stop parcel:trajectory terminated ',jp, ngrid
            if (ngrid.ge.0) then
              write(*,'(9G14.5)') x,y,z,dxsave,dysave,w,u,v,dt
!             print *,y-dysave*dyconst*float(ldirect),switchnorthg
            else
              write(*,'(6g14.5)') xlon,ylat,dxsave,dysave,xpol,ypol
            endif
!          write (*,'(a,3g12.3)')'stopped trajectory ',x,y,z
          endif
      endif

!     if (debug_out) print *,'abvanceB (x,y,z)> ',x,y,z

      if (AGEFactiv) then
        if (z_motion) then
          if (p0*exp(-z) > pcut) then
            nstop=4
          endif
        else
          if (TTLFILLactiv) then
            if (z<thetalowcut) then 
              nstop=4
            else if (z>thetacut) then
              nstop=5
            endif
          else 
            if (p0*(tint/z)**(1/kappa) > pcut) then
              nstop=4
            else if (z > thetacut) then
              nstop=5
            endif
          endif
        endif
      endif

      if (CLAUSactiv) then
        if (z_motion) then
           if (tint*exp(kappa*z) > thetacut) nstop=5
        else
           if (z > thetacut) nstop=5
        endif
      endif
 
      return
      end subroutine advanceB

end module advect
