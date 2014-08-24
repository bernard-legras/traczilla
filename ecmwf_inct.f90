
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

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TRACZILLA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

module ecmwf_inct

!*******************************************************************************
! This modules allow to use temperature tendencies archived at ECMWF as
! vertical velocities.
! Must be combined with isentropic interpolation
!*******************************************************************************

use commons
use isentrop_m
use ecmwf_diab
implicit none
save

integer idiffmax_inct,idiffnorm_inct
parameter(idiffmax_inct=43200,idiffnorm_inct=43200)

logical ecmwf_inct_w

! ecmwf_inct_w            Specify that temperature increments are added
!                         to ECMWF diabatic heating

integer numbwf_inct,wftime_inct(maxwf),lwindinterv_inct
character(len=32):: wfname_inct(maxwf),wfspec_inct(maxwf)

! lwindinterv_inct [s]    Interval between wind fields currently in memory
! numbwf_inct             actual number of wind fields
! wftime_inct(maxwf) [s]  times relative to beginning time of wind fields
! wfname_inct(maxwf)      file names of wind fields
! wfspec_inct(maxwf)      specifications of wind field file, e.g. if on hard 
!                         disc or on tape

      integer memtime_inct(2),memind_inct(2)

! memtime_inct [s]        validation times of wind fields in memory
! memind_inct             pointer to wind field, in order to avoid shuffling
!                         of wind fields

real, allocatable :: w_inct(:,:,:,:)

! w_inct [K/s]            vertical velocity as tendency in potential temperature
 
character(len=128):: path_inct(2)
integer len_inct(2)

! path_inct               paths for inct directory and associated AVAILABLE file
! len_inct                lengths of the two previous strings

contains

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine alloc_ecmwf_inct

allocate (w_inct(0:nx-1,0:ny-1,nuvz,2))
print *,'alloc_ecmwf_inct w_inct'

end subroutine alloc_ecmwf_inct


!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine getfields_inct(itime,nstop)

!******************************************************************************
!                                                                              *
!  This subroutine manages the 3 data fields to be kept in memory.             *
!  During the first time step of petterssen it has to be fulfilled that the    *
!  first data field must have |wftime|<itime, i.e. the absolute value of wftime*
!  must be smaller than the absolute value of the current time in [s].         *
!  The other 2 fields are the next in time after the first one.                *
!  Pointers (memind) are used, because otherwise one would have to resort the  *
!  wind fields, which costs some computing time. Here only the pointers are    *
!  resorted.                                                                   *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     29 April 1994                                                            *
!                                                                              *
!*******************************************************************************
!  Changes, Bernd C. Krueger, Feb. 2001:
!        Variables tth,qvh,tthn,qvhn (on eta coordinates) in common block.
!        Function of nstop extended.
!           B. Legras, April 2002:
!        Nested operations cancelled.
!        uuh, vvh, wwh on two times and in common
!        Message printed when reading a field
!           B. Legras, August 2005
!        diabatic version for ECMWF temperature tendencies
!           B. Legras, June 2009
!        version for temperature increment
!
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! lwindinterval [s]    time difference between the two wind fields read in     *
! indj                 indicates the number of the wind field to be read in    *
! indmin               remembers the number of wind fields already treated     *
! memind(2)            pointer, on which place the wind fields are stored      *
! memtime(2) [s]       times of the wind fields, which are kept in memory      *
! itime [s]            current time since start date of trajectory calculation *
! nstop                > 0, if trajectory has to be terminated                 *
! nx,ny,nuvz,nwz       field dimensions in x,y and z direction                 *
!                                                                              *
! Constants:                                                                   *
! idiffmax_inct        maximum allowable time difference between 2 wind fields *
!                                                                            *
!*******************************************************************************

      integer indj,indmin,itime,nstop,memaux
      save indmin
      data indmin/1/

! Check, if wind fields are available for the current time step
!**************************************************************

      nstop=0

      if ((ldirect*wftime_inct(1).gt.ldirect*itime).or.    &
        (ldirect*wftime_inct(numbwf_inct).lt.ldirect*itime)) then
        write(*,*) 'FLEXPART WARNING: NO TEMPERATURE INCREMENT ARE AVAILABLE.'
        write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
        nstop=4
        return
      endif

      if ((ldirect*memtime_inct(1).le.ldirect*itime).and.   &
        (ldirect*memtime_inct(2).gt.ldirect*itime)) then

! The right wind fields are already in memory -> don't do anything
!*****************************************************************

        continue

      else if ((ldirect*memtime_inct(2).le.ldirect*itime).and. &
        (memtime_inct(2).ne.999999999)) then
 
! Current time is after 2nd wind field
! -> Resort wind field pointers, so that current time is between 1st and 2nd
!***************************************************************************

        memaux=memind_inct(1)
        memind_inct(1)=memind_inct(2)
        memind_inct(2)=memaux
        memtime_inct(1)=memtime_inct(2)

! Read a new wind field and store it on place memind(2)
!******************************************************

        do indj=indmin,numbwf_inct-1
           if (ldirect*wftime_inct(indj+1).gt.ldirect*itime) then
              call read_inct(indj+1,memind_inct(2))
              memtime_inct(2)=wftime_inct(indj+1)
              call verttransform_inct(memind_inct(2),2)
!              print *, 'getfields_inct> reads new field at time ',memtime_inct(2)
              write(*,'(a,a,a,i11,a,i11)') &
                 ' getfields_inct> ',trim(wfname_inct(indj+1)), &
                 '  memtime ',memtime_inct(2),'  time ',itime
              nstop = 1
              goto 40
           endif
       enddo
 40    indmin=indj

      else

! No wind fields, which can be used, are currently in memory 
! -> read both wind fields
!***********************************************************

         do indj=indmin,numbwf_inct-1
            if ((ldirect*wftime_inct(indj).le.ldirect*itime).and.   &
                  (ldirect*wftime_inct(indj+1).gt.ldirect*itime)) then
               memind_inct(1)=1
               call read_inct(indj,memind_inct(1))
               memtime_inct(1)=wftime_inct(indj)
!               print*,'getfields_inct> reads first field:  time ',memtime_inct(1)
               write(*,'(a,a,a,i11,a,i11)') &
                 ' getfields_inct> ',trim(wfname_inct(indj)), &
                 '  memtime ',memtime_inct(1),'  time ',itime
               memind_inct(2)=2
               call read_inct(indj+1,memind_inct(2))
               memtime_inct(2)=wftime_inct(indj+1)
               call verttransform_inct(memind_inct(1),1)
               call verttransform_inct(memind_inct(2),2)
!               print*,'getfields_inct> reads second field: time ',memtime_inct(2)
               write(*,'(a,a,a,i11,a,i11)') &
                 ' getfields_inct> ',trim(wfname_inct(indj+1)), &
                 '  memtime ',memtime_inct(2),'  time ',itime
               nstop = 1
               goto 60
            endif
         enddo
 60      indmin=indj

      endif

      lwindinterv_inct=abs(memtime_inct(2)-memtime_inct(1))
 
      if (lwindinterv_inct.gt.idiffmax_inct) nstop=3

      return
      end subroutine getfields_inct

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine read_inct(indj,n)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READWIND                    *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      G. WOTAWA                                  *
!             DATE:        1997-08-05                                 *
!             LAST UPDATE: 2000-10-17, Andreas Stohl                  *
!                                                                     * 
!**********************************************************************
!  Changes, Bernd C. Krueger, Feb. 2001:
!        Variables tth and qvh (on eta coordinates) in common block
!           B. Legras, Apr. 2002:
!        Two times for uuh, vvh, wwh, in common
!        select case + loop inside
!           B. Legras, August 2005
!        reads temperature tendencies
!**********************************************************************
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL TEMPERATURE TENDENCIES. THE         *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!                                                                     *
! IMPORTANT VARIABLES FROM SHARED VARIABLES:                          *
!                                                                     *
! wfname_inct        File name of data to be read in                  *
! nx,ny,nuvz,nwz     expected field dimensions                        *
! nlev_ec            number of vertical levels ecmwf model            *
! w_inct             temperature tendency over 3h                     *
!                                                                     *
!**********************************************************************

      integer, intent(in):: indj,n
      integer i,j,k,levdiff2,ifield,iumax,iwmax,lunit

! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

! dimension of isec2 at least (22+n), where n is the number of parallels or
! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

! dimension of zsec2 at least (10+nn), where nn is the number of vertical
! coordinate parameters

      integer isec0(2),isec1(56),isec3(2)
      integer isec4(64),ilen,ierr,iword,count130
      integer, allocatable :: inbuff(:), isec2(:)
      real zsec3(2)
      real, allocatable :: zsec4(:), zsec2(:)
      real xaux,yaux,xaux0,yaux0

      character(len=1):: yoper
      logical hflswitch,strswitch

      data yoper/'D'/

      hflswitch=.false.
      strswitch=.false.
      levdiff2=nlev_ec-nwz+1
      iumax=0
      iwmax=0
      allocate (zsec4(jpunp),inbuff(jpack))
      allocate (isec2(22+nx+ny),zsec2(60+nuvz+nwz))

! OPENING OF DATA FILE (GRIB CODE)

      call pbopen(lunit,path_inct(1)(1:len_inct(1))//wfname_inct(indj),'r',ierr)
      if(ierr < 0) goto 999

      ifield=0
      count130=0

10    ifield=ifield+1
!
! GET NEXT FIELDS
!
      call pbgrib(lunit,inbuff,jpack,ilen,ierr)
      if(ierr.eq.-1) goto 50    ! EOF DETECTED
      if(ierr.lt.-1) goto 888   ! ERROR DETECTED

      ierr=1

      call gribex(isec0,isec1,isec2,zsec2,isec3,zsec3,isec4,  &
                  zsec4,jpunp,inbuff,jpack,iword,yoper,ierr)
      if (ierr.ne.0) goto 888   ! ERROR DETECTED

      if(ifield.eq.1) then

! CHECK GRID SPECIFICATIONS

        if(isec2(2).ne.nxfield) stop 'READ_INCT: NX NOT CONSISTENT'
        if(isec2(3).ne.ny) stop 'READ_INCT: NY NOT CONSISTENT'
        if(isec2(12)/2-1.ne.nlev_ec) &
          stop 'READ_INCT: VERTICAL DISCRETIZATION NOT CONSISTENT'
        xaux=float(isec2(5))/1000.
        yaux=float(isec2(7))/1000.
        xaux0=xlon0
        yaux0=ylat0
        if(xaux<0.) xaux=xaux+360.
        if(yaux<0.) yaux=yaux+360.
        if(xaux0<0.) xaux0=xaux0+360.
        if(yaux0<0.) yaux0=yaux0+360.
        if(abs(xaux-xaux0)>1.e-3) then
          print *, xaux, xaux0
          stop 'READ_INCT: LOWER LEFT LONGITUDE NOT CONSISTENT'
        endif
        if(abs(yaux-yaux0)>1.e-3) then
          print *,yaux,yaux0
          stop 'READ_INCT: LOWER LEFT LATITUDE NOT CONSISTENT'
        endif
      endif

!     ACTHTUNG: The temperature increments are stored like horizontal wind
!     leaving here the bottom level defined to zero in absence of other
!     prescription. It is probably not a good idea, anyway, to use
!     diabatic velocities in the boundary layer.
!     It is important for DIAB runs to set properly upper and lower
!     theta levels in COMMAND file
      k=isec1(8)
      field_identifier: select case (isec1(6))
        case (130)              !! Tendency of clear sky SWR
          do j=0,ny-1 ; do i=0,nxfield-1
            w_inct(i,j,nlev_ec-k+2,n) = zsec4(nxfield*(ny-j-1)+i+1)
          enddo       ; enddo
          count130=count130+1  
        case  default
 !      stop 'READ_DIAB DOES NOT FIND T TENDENCY IN FILE READ'
        end select field_identifier     

      goto 10              !! READ NEXT LEVEL OR PARAMETER
!
! CLOSING OF INPUT DATA FILE
!
50    call pbclose(lunit,ierr)     !! FINISHED READING / CLOSING GRIB FILE
      

! For global fields, assign rightmost grid point the value of the
! leftmost point
!****************************************************************

      if (xglobal) then
        w_inct(nx-1,:,:,n)=w_inct(0,:,:,n)
      endif

      deallocate (zsec4,inbuff,isec2,zsec2)
     
      if(count130==0) go to 889
!      print *,'read_inct> file ',trim(wfname_inct(indj)), count130

      return    
888   write(*,*) ' #### FLEXPART MODEL ERROR! INCT FIELD        #### ' 
      write(*,*) ' #### ',wfname_inct(indj),'                   #### '
      write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
      stop 'Execution terminated'
889   write(*,*) ' #### FLEXPART MODEL ERROR! INCT FIELD        #### '
      write(*,*) ' #### ',wfname_inct(indj),'                   #### '
      write(*,*) ' #### does not contain required fields        #### '
      print *,count130
      stop 'Execution terminated'
999   write(*,*) ' #### FLEXPART MODEL ERROR! INCT FIELD        #### ' 
      write(*,*) ' #### ',wfname_inct(indj),'                   #### '
      write(*,*) ' #### CANNOT BE OPENED !!!                    #### '
      stop 'Execution terminated'

      return
      end subroutine read_inct

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine readpaths_inct(pathfile,error)
!                            o
!*******************************************************************************
!                                                                              *
!     Reads the pathnames, where input/output files are expected to be.        *
!     The file pathnames must be available in the current working directory.   *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     1 February 1994                                                          *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! error              .true., if file pathnames does not exist                  *
! len(numpath)       lengths of the path names                                 *
! path(numpath)      pathnames of input/output files                           *
!                                                                              *
! Constants:                                                                   *
! numpath            number of pathnames to be read in                         *
!                                                                              *
!*******************************************************************************

      integer i
      logical, intent(out)::error

      character(len=*), intent(in):: pathfile

      error=.false.

! Read the pathname information stored in unitpath
!*************************************************

      open(unitpath,file=pathfile,status='old',err=999)
      print *,'opening pathfile ',trim(pathfile)

!     skip previously read data
      do i=1,numpath+2
        read(unitpath,*)
      enddo
      do i=1,2
        read(unitpath,'(a)',err=998) path_inct(i) 
        len_inct(i)=index(path_inct(i),' ')-1
      enddo

      close(unitpath)
      return    

998   write(*,*) ' #### TRAJECTORY MODEL ERROR! ERROR WHILE     #### ' 
      write(*,*) ' #### READING INCT     FILE PATHNAMES.        #### ' 

999   write(*,*) ' #### TRAJECTORY MODEL ERROR! PATH FILE       #### ' 
      write(*,*) ' #### CANNOT BE OPENED IN THE CURRENT WORKING #### '
      write(*,*) ' #### DIRECTORY.                              #### '
      error=.true.

      return
      end subroutine readpaths_inct

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine readavailable_inct(error)
!                                o
!*******************************************************************************
!                                                                              *
!   This routine reads the dates and times for which windfields are available. *
!                                                                              *
!     Authors: A. Stohl                                                        *
!                                                                              *
!     6 February 1994                                                          *
!     8 February 1999, Use of nested fields, A. Stohl    
!     7 April    2002, Suppression of code related to nested fields, B. Legras
!       August   2005, Reads available fields for temperatue tendencies
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! bdate                beginning date as Julian date                           *
! beg                  beginning date for windfields                           *
! end                  ending date for windfields                              *
! error                .true., if error ocurred in subprogram, else .false.    *
! fname                filename of wind field, help variable                   *
! ideltas [s]          duration of modelling period                            *
! idiff                time difference between 2 T tendency fields             *
! idiffnorm            normal time difference between 2 T tendency fields      *
! idiffmax_inct [s]         maximum allowable time between 2 T tendency fields      *
! jul                  julian date, help variable                              *
! numbwf_inct          actual number of T tendency fields                      *
! wfname_inct(maxwf)   file names of needed T tendency fields                  *
! wfspec_inct(maxwf)   file specifications of wind fields (e.g., if on disc)   *
! wftime_inct(maxwf) [s] times of T tendency fields relative to beginning time *
! wfname1,wfspec1,wftime1 = same as above, but only local (help variables)     *
!                                                                              *
! Constants:                                                                   *
! maxwf                maximum number of wind fields                           *
! unitavailab          unit connected to file AVAILABLE                        *
!                                                                              *
!*******************************************************************************
   
      use date
      logical error, toobig
      integer i,j,idiff,ldat,ltim,wftime1(maxwf)
      integer year,month,day
      double precision jul,beg,end
      character(len=32):: fname,spec,wfname1(maxwf),wfspec1(maxwf)

      error=.false.

! Windfields are only used, if they are within the modelling period.
! However, 1 additional day at the beginning and at the end is used for
! interpolation. -> Compute beginning and ending date for the windfields.
!************************************************************************

      if (ideltas.gt.0) then         ! forward trajectories
        beg=bdate-1.                 ! idiffmax_inct should be used here too               
        end=bdate+dble(float(ideltas)/86400.) &
                 +dble(float(idiffmax_inct)/86400.)
      else                           ! backward trajectories
        beg=bdate+dble(float(ideltas)/86400.) &
                 -dble(float(idiffmax_inct)/86400.)
        end=bdate+1.                 ! idiffmax_inct shoukd be used here too
      endif

! Open the wind field availability file and read available wind fields
! within the modelling period.
!*********************************************************************

      open(unitavailab,file=path_inct(2)(1:len_inct(2)),status='old',err=999)

      do i=1,3
        read(unitavailab,*)
      enddo
      
      numbwf_inct=0
100   read(unitavailab,'(i8,1x,i6,6x,a16,6x,a10)',end=99) &
           ldat,ltim,fname,spec
      jul=juldate(ldat,ltim)
      if ((jul.ge.beg).and.(jul.le.end)) then
        numbwf_inct=numbwf_inct+1
        if (numbwf_inct.gt.maxwf) then      ! check exceedance of dimension
           write(*,*) 'Number of needed inct fields is too great.'
           write(*,*) 'actual, allowed ',numbwf_inct, maxwf
           write(*,*) 'Reduce modelling period (file "COMMAND") or'
           write(*,*) 'reduce number of inct fields (file "AVAILABLE").'
           goto 1000
        endif

!  set list of filenames and times
!  offset time by +6h as indicated in AVAILABLE takes 
!  into account the fact that assimilation is performed 
!  over an interval of 12h
        wfname1(numbwf_inct)=fname
        wfspec1(numbwf_inct)=spec
        wftime1(numbwf_inct)=nint((jul-bdate)*86400.)
      endif
      goto 100       ! next wind field

99    continue

      close(unitavailab)

! Check wind field times of file AVAILABLE (expected to be in temporal order)
!****************************************************************************

      if (numbwf_inct.eq.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! NO INCT FIELDS    #### '
        write(*,*) ' #### AVAILABLE FOR SELECTED TIME PERIOD.     #### '
        error=.TRUE.
        return
      endif

      do i=2,numbwf_inct
        if (wftime1(i).le.wftime1(i-1)) then
          write(*,*) 'FLEXPART ERROR: FILE AVAILABLE-INCT IS CORRUPT.'
          write(*,*) 'THE TEMPERATURE TENDENCIES ARE NOT IN TEMPORAL ORDER.'
          write(*,*) 'PLEASE CHECK FIELD ',wfname1(i)
          error=.TRUE.
          return
        endif
      enddo

! For backward trajectories, reverse the order of the windfields
!***************************************************************

      if (ideltas.ge.0) then
        wfname_inct(1:numbwf_inct)=wfname1(1:numbwf_inct)
        wfspec_inct(1:numbwf_inct)=wfspec1(1:numbwf_inct)
        wftime_inct(1:numbwf_inct)=wftime1(1:numbwf_inct)
      else
        do i=1,numbwf_inct
          wfname_inct(numbwf_inct-i+1)=wfname1(i)
          wfspec_inct(numbwf_inct-i+1)=wfspec1(i)
          wftime_inct(numbwf_inct-i+1)=wftime1(i)
        enddo
      endif

! Check the time difference between the wind fields. If it is big, 
! write a warning message. If it is too big, terminate the trajectory. 
!*********************************************************************

      wftime1(:)=wftime_inct(:)
      do i=2,numbwf_inct
        idiff=abs(wftime_inct(i)-wftime_inct(i-1))
        toobig=.false.
        if (idiff.gt.idiffmax_inct) then
!         detect end of february for bissextil year in perpetual run
!         shift all remaining times in such case to satisfy check done in
!         getfield
          if(perpetual) then
            ! check done on unshifted times to process several bissextil years
            if(ideltas > 0) call caldate(wftime1(i-1)/86400.+bdate,ldat,ltim)
            if(ideltas < 0) call caldate(wftime1(i  )/86400.+bdate,ldat,ltim)
            year=int(ldat/10000)
            month=int((ldat-10000*year)/100)
            day=int(ldat-10000*year-100*month)
            if(month==2.and.day==28.and.mod(year,4)==0) then
              print *,'readavailable_inct > jump due to bissextil year in perpetual run'
              print *,'readavailable_inct > apply one-day shift to all remaining times'
              if(ideltas > 0) then
                do j=i,numbwf_inct
                  wftime_inct(j)=wftime_inct(j)-86400
                enddo
              else
                do j=i,numbwf_inct
                  wftime_inct(j)=wftime_inct(j)+86400
                enddo
              endif
            else
              toobig=.true.
            endif
          endif
          if(.not.perpetual.or.toobig) then
            write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
            write(*,*) 'INCT FIELDS IS TOO BIG FOR TRANSPORT CALCULATION.'
            write(*,*) 'THEREFORE, TRAJECTORIES HAVE TO BE SKIPPED.'
            print *, wfname_inct(i),wfname_inct(i-1),idiffmax_inct
          endif
        else if (idiff.gt.idiffnorm_inct) then
          write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
          write(*,*) 'INCT FIELDS IS BIG. THIS MAY CAUSE A DEGRADATION'
          write(*,*) 'OF SIMULATION QUALITY.'
          print *,  wfname_inct(i),wfname_inct(i-1),idiff,idiffnorm
        endif
      enddo


! Reset the times of the wind fields that are kept in memory to no time
!**********************************************************************

      do i=1,2
        memind_inct(i)=i
        memtime_inct(i)=999999999
      enddo

      return    

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE #### '
      write(*,'(a)') '     '//path_inct(2)(1:len_inct(2)) 
      write(*,*) ' #### CANNOT BE OPENED           #### '
1000  error=.true.

      return
      end subroutine readavailable_inct

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine verttransform_inct(n,m)

!*******************************************************************************
!                                                                         *
!     This subroutine calculates the potential temperature tendencies from     *
!     the ECMWF temperature tendencies                                         *
!                                                                              *
!     Author: B. Legras                                                        *
!     date:   August 2005                                                      *
!                                                                              *
!     from a previous version of verttransform.f                               *
!                                                                              *
!                                                                              *
! Variables:                                                                   *
! nx,ny,nz                        field dimensions in x,y and z direction      *
! w_inct(0:nxmax,0:nymax,nuvz,2)  potential temperature tendencies [K/s]       *
! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                        *
!                                                                              *
!*******************************************************************************

  integer, intent(in) :: n,m
  integer ix, jy, iz
  real pp

! Loop over the whole grid
!*************************
! Rescale increment over 12h to increment per s

  do jy=0,ny-1
    do ix=0,nx-1
      do iz=2, nuvz
        pp=akz(iz)+bkz(iz)*ps(ix,jy,1,n)
        w_inct(ix,jy,iz,n)=w_inct(ix,jy,iz,n)*((p0/pp)**kappa)/43200.
      enddo
    enddo
  enddo
! Note it is not correct to use the index n for pp as it is not sure that it
! correponds to the inc centered time or to 3h before or after

! Set the tendency at 10m to be that at next level
  w_inct(:,:,1,n)=w_inct(:,:,2,n)

!  Perform mass correction if required
  if(mass_correction) call diab_mass_inct(n,m)

  return
  end subroutine verttransform_inct

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass_inct(n,m)

!*******************************************************************************
! This subroutine calculates the mass flux across pure pressure levels in the 
! stratosphere from heating rates, pressure and temperature.
! The isentropic density sigma is calculated and the mass flux is calculated
! by taking its product with the heating rate (in d theta/dt) and integrated 
! over the sphere.
! The integral is then divided by the integral of sigma to get the the mass 
! averaged heating rate to be removed on each vertical in order to achieve the 
! mass conservation acrosss the pressure surface.
! Vertical derivate involved in sigma is calculated by centered differences
! except on the last upper level.
! (note that the derivative on the first pure pressure level neglects the bkz
! coefficient on the layer just below)
!
! B. Legras
! 01/03/2013
!
!*******************************************************************************

  integer, intent(in) :: n,m
  real, allocatable :: theta(:,:,:), sigma(:,:,:), flux(:,:), flux_lat(:)
  real, allocatable :: mass_flux(:), mean_sigma(:), mean_sigma_lat(:)
  real, allocatable :: mean_w(:)
  integer :: k

  allocate(theta(0:nx-1,0:ny-1,nuvz))
! Calculation of the potential temperature
  do k=nuvz-NPureP,nuvz
     theta(:,:,k)=0.5*(tth(:,:,k,1)+tth(:,:,k,2))*pif(k)
  enddo
! Calculation of the density sigma
  allocate(sigma(0:nx-1,0:ny-1,nuvz))
  do k=nuvz-NPureP+1,nuvz-1
     sigma(:,:,k)=pmc(k)/(theta(:,:,k+1)-theta(:,:,k-1))
  enddo
  sigma(:,:,nuvz)=pmc(nuvz)/(theta(:,:,nuvz)-theta(:,:,nuvz-1))
  allocate(mean_sigma(nuvz))
  allocate(mean_sigma_lat(0:ny-1))
  
! Calculation of the mass flux across the surface
  allocate(flux(0:nx-1,0:ny-1))
  allocate(flux_lat(0:ny-1))
  allocate(mass_flux(nuvz),mean_w(nuvz))
  do k=nuvz-NPureP+1,nuvz
     flux(:,:)=w_inct(:,:,k,n)*sigma(:,:,k)
     flux_lat(:)=sum(flux(0:nx-2,:),DIM=1)
     mean_sigma_lat(:)=sum(sigma(0:nx-2,:,k),DIM=1)
     mass_flux(k)=0.5*dot_product(flux_lat(:),area_coefft)
     mean_sigma(k)=0.5*dot_product(mean_sigma_lat(:),area_coefft)
     mean_w(k)=mass_flux(k)/mean_sigma(k)
  enddo

! Output of the ass averaged diab heating
  if(mean_diab_output) then
    write(unitflux2) memtime_inct(m),mean_w(nuvz-NPureP+1:nuvz)
    flush(unitflux2)
  endif

! Correction step

  do k=nuvz-NPureP+1,nuvz
    w_inct(:,:,k,n)=w_inct(:,:,k,n)-mean_w(k)
  enddo

  deallocate(theta,sigma,mean_sigma,mean_sigma_lat)
  deallocate(flux,flux_lat,mass_flux,mean_w)

  return
  end subroutine diab_mass_inct

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass_inct_init

!*************************************************************************************
! This routine initializes the output file
! Coefficients are initialized in diab_mass_init belonging to the module
! ecmwf_diab which must be used together with this module
!*************************************************************************************

  NPureP=23

  if(.not.(xglobal.and.nglobal.and.sglobal)) then
      write(*,*) ' #### TRACZILLA MODEL ERROR! DIAB_MASS        #### ' 
      write(*,*) ' #### mass flux cannot be balanced            #### '
      write(*,*) ' #### on not global grid                      #### '
      stop  
  endif

! Create the output file of the mas averaged heating flux
  open(unitflux2,file=trim(path(2))//'MassMeanInct2.dat',form='unformatted',position='append')

  return
  end subroutine diab_mass_inct_init


!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine interpol_wind_theta_inct   &
         (itime,xt,yt, ztheta, dxdt,dydt,dzdt, ngrid, theta_inf, theta_sup, & 
          z_factor,tint,nstop)

!*******************************************************************************
!                                                                              *
!  This subroutine interpolates the wind data to current trajectory position.  *
!                                                                              *
!    Author: A. Stohl                                                          *
!                                                                              *
!    16 December 1997  
! 
!    Changes: B. Legras, April 2002
!             interpolation from eta winds
!             variance calculation cancelled
!             B. Legras, June 2002
!             optimisation of log calculations
!             new calculation of z_factor for z and theta diffusion
!             B. Legras, June 2005
!             isentropc version
!             B. Legras, August 2005
!             diabatic version
!             B. Legras, March 2006
!             test on theta boundaries to stop the parcel integration
!             B. Legras, June 2009
!             add temperature increment
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! dxdt,dydt          horizontal wind components in grid units per second       *
! itime [s]          current temporal position                                 *
! memtime(3) [s]     times of the wind fields in memory                        *
! xt,yt,theta        coordinates position for which wind data shall be calculat*
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!*******************************************************************************

      integer, intent(in) :: itime,ngrid
      integer, intent(inout):: nstop
      real, intent(in) :: xt,yt,ztheta
      real, intent(out) :: dxdt,dydt,dzdt,z_factor,tint, theta_inf, theta_sup

! Auxiliary variables needed for interpolation
      real theta
      real u1(2),v1(2),w1d(2),w1i(2),dt1,dt2,dtt,tp1(2)
      real dt1_diab,dt2_diab,dtt_diab
      real dt1_inct,dt2_inct,dtt_inct
      real tr(4,2),trp(4,2),u(4,2),v(4,2),wd(4,2),wi(4,2),tp(4,2)
      integer m,indexh,indz(4,2)
      integer ix,jy,ixp,jyp,i0,j0,idxy
      real ddx,ddy,rddx,rddy,p1,p2,p3,p4
      real psl0,psup0,pinf0,pisup0,piinf0

!********************************************
! Multilinear interpolation in time and space
!********************************************

! Determine the lower left corner and its distance to the current position
!*************************************************************************

      ! min and max required for the points just falling on the boundary
      ! as it may happen typically as a result of initialization
      ix=min(floor(xt),nx-2)       ;  jy=min(floor(yt),ny-2)
      ixp=ix+1         ;  jyp=jy+1
      ddx=modulo(xt-float(ix),1.) ;  ddy=yt-float(jy)
      rddx=1.-ddx      ;  rddy=1.-ddy
      p1=rddx*rddy     ;  p2=ddx*rddy
      p3=rddx*ddy      ;  p4=ddx*ddy
      
! Calculate coefficients for temporal interpolation
!**************************************************

      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)
      dt1_diab=float(itime-memtime_diab(1))
      dt2_diab=float(memtime_diab(2)-itime)
      dtt_diab=1./(dt1_diab+dt2_diab)
      dt1_inct=float(itime-memtime_inct(1))
      dt2_inct=float(memtime_inct(2)-itime)
      dtt_inct=1./(dt1_inct+dt2_inct)

! Calculates the theta values on the four adjacent columns if required
!*********************************************************************

      if(.not.theta_col(ix,jy,memind(1)))   call calc_col_theta(ix,jy,memind(1))
      if(.not.theta_col(ix,jy,memind(2)))   call calc_col_theta(ix,jy,memind(2))
      if(.not.theta_col(ix,jyp,memind(1)))  call calc_col_theta(ix,jyp,memind(1))
      if(.not.theta_col(ix,jyp,memind(2)))  call calc_col_theta(ix,jyp,memind(2))
      if(.not.theta_col(ixp,jy,memind(1)))  call calc_col_theta(ixp,jy,memind(1))
      if(.not.theta_col(ixp,jy,memind(2)))  call calc_col_theta(ixp,jy,memind(2))
      if(.not.theta_col(ixp,jyp,memind(1))) call calc_col_theta(ixp,jyp,memind(1))
      if(.not.theta_col(ixp,jyp,memind(2))) call calc_col_theta(ixp,jyp,memind(2))

! Determine the level below the current position for u,v
!*******************************************************

!  Locates lower left corner

      theta=ztheta
      indz(1,1) = locisent(theta,ix,jy,lower_theta_level,upper_theta_level,memind(1))
      
!  Locates other points by assuming they are close to the first one
            
      indz(2,1) = locisent2(theta,ix,jyp,indz(1,1),memind(1))
      indz(3,1) = locisent2(theta,ixp,jy,indz(1,1),memind(1))
      indz(4,1) = locisent2(theta,ixp,jyp,indz(1,1),memind(1))
      indz(1,2) = locisent2(theta,ix,jy,indz(1,1),memind(2))
      indz(2,2) = locisent2(theta,ix,jyp,indz(1,1),memind(2))
      indz(3,2) = locisent2(theta,ixp,jy,indz(1,1),memind(2))
      indz(4,2) = locisent2(theta,ixp,jyp,indz(1,1),memind(2))

!  Defines potential temperature at the 8 nearby meshpoints for
!  the two times

      tr (1,1) = theta_g(indz(1,1),ix,jy,memind(1))
      tr (1,2) = theta_g(indz(1,2),ix,jy,memind(2))
      tr (2,1) = theta_g(indz(2,1),ix,jyp,memind(1))
      tr (2,2) = theta_g(indz(2,2),ix,jyp,memind(2))
      tr (3,1) = theta_g(indz(3,1),ixp,jy,memind(1))
      tr (3,2) = theta_g(indz(3,2),ixp,jy,memind(2))
      tr (4,1) = theta_g(indz(4,1),ixp,jyp,memind(1))
      tr (4,2) = theta_g(indz(4,2),ixp,jyp,memind(2))
      trp(1,1) = theta_g(indz(1,1)+1,ix,jy,memind(1))
      trp(1,2) = theta_g(indz(1,2)+1,ix,jy,memind(2))
      trp(2,1) = theta_g(indz(2,1)+1,ix,jyp,memind(1))
      trp(2,2) = theta_g(indz(2,2)+1,ix,jyp,memind(2))
      trp(3,1) = theta_g(indz(3,1)+1,ixp,jy,memind(1))
      trp(3,2) = theta_g(indz(3,2)+1,ixp,jy,memind(2))
      trp(4,1) = theta_g(indz(4,1)+1,ixp,jyp,memind(1))
      trp(4,2) = theta_g(indz(4,2)+1,ixp,jyp,memind(2))
!  Patch to avoid extrapolations at top and bottom of the domain     
      theta = max(theta,maxval(tr))
      theta = min(theta,minval(trp))

! Provides upper and lower theta bounds
!**************************************
      
      if(theta_bounds) then
        theta_inf = ((theta_g(lower_theta_level,ix,jy,memind(1))*p1 &
                    + theta_g(lower_theta_level,ix,jyp,memind(1))*p2 &
                    + theta_g(lower_theta_level,ixp,jy,memind(1))*p3 &
                    + theta_g(lower_theta_level,ixp,jyp,memind(1))*p4)*dt1 &
                   + (theta_g(lower_theta_level,ix,jy,memind(2))*p1 &
                    + theta_g(lower_theta_level,ix,jyp,memind(2))*p2 &
                    + theta_g(lower_theta_level,ixp,jy,memind(2))*p3 &
                    + theta_g(lower_theta_level,ixp,jyp,memind(2))*p4)*dt2) &
                   * dtt
        theta_sup = ((theta_g(upper_theta_level,ix,jy,memind(1))*p1 &
                    + theta_g(upper_theta_level,ix,jyp,memind(1))*p2 &
                    + theta_g(upper_theta_level,ixp,jy,memind(1))*p3 &
                    + theta_g(upper_theta_level,ixp,jyp,memind(1))*p4)*dt1 &
                   + (theta_g(upper_theta_level,ix,jy,memind(2))*p1 &
                    + theta_g(upper_theta_level,ix,jyp,memind(2))*p2 &
                    + theta_g(upper_theta_level,ixp,jy,memind(2))*p3 &
                    + theta_g(upper_theta_level,ixp,jyp,memind(2))*p4)*dt2) &
                   * dtt
       
      endif

! Halt trajectories which are too close from lower boundary
! and flag those which reach the upper boundary
!**********************************************************************

      if(minval(indz)==lower_theta_level) then
          nstop=2
          dxdt=0. ; dydt=0.; dzdt=0.; tint=275.; z_factor=1.
          return
      endif
      if(maxval(indz)==upper_theta_level-1) then
          nstop=-4
      endif


!**********************************************************************
! 1.) Bilinear horizontal interpolation
! This has to be done separately for 4 fields (Temporal(2)*Vertical(2))
!**********************************************************************

! Loop over 2 time steps and 2 levels
!************************************

      if (ngrid < 0) then  ! polar region
      
      select case(vert_interpol)
      
      case('log')
        do m=1,2
          indexh=memind(m)
            
            u(1,m)=(uupol(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + uupol(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            v(1,m)=(vvpol(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + vvpol(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            wd(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            wi(1,m)=(w_inct(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + w_inct(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            wd(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))   
            wi(2,m)=(w_inct(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + w_inct(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            wd(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            wi(3,m)=(w_inct(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + w_inct(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            wd(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            wi(4,m)=(w_inct(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + w_inct(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1d(m)=p1*wd(1,m)+p2*wd(2,m)+p3*wd(3,m)+p4*wd(4,m)
          w1i(m)=p1*wi(1,m)+p2*wi(2,m)+p3*wi(3,m)+p4*wi(4,m)
         
        enddo
        
      case('lin')
      
        do m=1,2
          indexh=memind(m)
            
            u(1,m)=(uupol(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + uupol(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            v(1,m)=(vvpol(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + vvpol(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            wd(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            wi(1,m)=(w_inct(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + w_inct(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            wd(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            wi(2,m)=(w_inct(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + w_inct(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            wd(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            wi(3,m)=(w_inct(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + w_inct(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            wd(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            wi(4,m)=(w_inct(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + w_inct(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1d(m)=p1*wd(1,m)+p2*wd(2,m)+p3*wd(3,m)+p4*wd(4,m)
          w1i(m)=p1*wi(1,m)+p2*wi(2,m)+p3*wi(3,m)+p4*wi(4,m)
          
        enddo
      
      end select 
                
      else     ! non polar region
      
      select case(vert_interpol)
      
      case('log')      
        do m=1,2
          indexh=memind(m)

            u(1,m)=(uuh(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) & 
                  + uuh(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            v(1,m)=(vvh(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + vvh(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            wd(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            wi(1,m)=(w_inct(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + w_inct(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            wd(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            wi(2,m)=(w_inct(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + w_inct(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            wd(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            wi(3,m)=(w_inct(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + w_inct(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            wd(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            wi(4,m)=(w_inct(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + w_inct(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1d(m)=p1*wd(1,m)+p2*wd(2,m)+p3*wd(3,m)+p4*wd(4,m)
          w1i(m)=p1*wi(1,m)+p2*wi(2,m)+p3*wi(3,m)+p4*wi(4,m)
  
        enddo

      case('lin')
      
        do m=1,2
          indexh=memind(m)

            u(1,m)=(uuh(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m)) &
                  + uuh(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))
            v(1,m)=(vvh(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m)) &
                  + vvh(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))
            wd(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m)) &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))  
            wi(1,m)=(w_inct(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m)) &
                  + w_inct(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))  
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            wd(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            wi(2,m)=(w_inct(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + w_inct(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            wd(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            wi(3,m)=(w_inct(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + w_inct(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            wd(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            wi(4,m)=(w_inct(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + w_inct(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))

          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1d(m)=p1*wd(1,m)+p2*wd(2,m)+p3*wd(3,m)+p4*wd(4,m)
          w1i(m)=p1*wi(1,m)+p2*wi(2,m)+p3*wi(3,m)+p4*wi(4,m)    
        
        enddo
      
      end select     

      endif

! Calculation of z_factor for vertical diffusion

      select case (diftype)
      case (1)          ! diffusion in z (cf p.46, book C, part2)
!       d theta / dz = - (g/theta) d theta / d Pi = -g d Log(theta) / d Pi
!       where Pi = (p/p0)**kappa
!       estimated from the data on the lower left corner at first time
!       of the interval, for a better estimate using the closest point
!       activate the first following line and deactivate the second one 
!       call sort_hor_distance       
        i0=ix; j0=jy; idxy=1
        psl0 = ps(ix ,jy ,1,memind(1))
        psup0 = akz(indz(idxy,1)+1)+bkz(indz(idxy,1)+1)*psl0
        pisup0 = cpa*(psup0/p0)**kappa
        pinf0 = akz(indz(idxy,1))+bkz(indz(idxy,1))*psl0
        piinf0 = cpa*(pinf0/p0)**kappa
        z_factor = -ga*(log(trp(idxy,1))-log(tr(idxy,1)))/(pisup0-piinf0)
      case (2)          ! diffusion in theta
        z_factor = 1.
      case default
        z_factor = 0.
      end select
 
! Accurate calculation of the temperature if needed
 
      if (AccurateTemp) then
      
      select case(vert_interpol)
      
      case('log')
        do m=1,2
          indexh=memind(m)
            
            tp(1,m)=(tth(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + tth(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m))) 
            tp(2,m)=(tth(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + tth(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            tp(3,m)=(tth(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + tth(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            tp(4,m)=(tth(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + tth(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
          
          tp1(m)=p1*tp(1,m)+p2*tp(2,m)+p3*tp(3,m)+p4*tp(4,m)
          
        enddo

      case('lin')
      
        do m=1,2
          indexh=memind(m)
            
            tp(1,m)=(tth(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + tth(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            tp(2,m)=(tth(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + tth(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m)) 
            tp(3,m)=(tth(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + tth(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m)) 
            tp(4,m)=(tth(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + tth(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))      
          
          tp1(m)=p1*tp(1,m)+p2*tp(2,m)+p3*tp(3,m)+p4*tp(4,m)
          
        enddo
      
      end select
      endif
      
!************************************
! 3.) Temporal interpolation (linear)
!************************************

      dxdt=(u1(1)*dt2+u1(2)*dt1)*dtt
      dydt=(v1(1)*dt2+v1(2)*dt1)*dtt
      dzdt=(w1d(1)*dt2_diab+w1d(2)*dt1_diab)*dtt_diab &
          +(w1i(1)*dt2_inct+w1i(2)*dt1_inct)*dtt_inct
      if(AccurateTemp) tint=(tp1(1)*dt2+tp1(2)*dt1)*dtt
      !if(debug_out) &
      ! print "('interpol>',i3,' P ',3f7.0,' T ',3f7.2,' TH ',3f7.2)", & ! test
      !   indz(1,1),p0*(theta/tint)**(-1/kappa),&
      !   p0*(theta_g(indz(1,1),ix,jy,memind(1))/tth(ix,jy,indz(1,1),memind(1)))**(-1/kappa),&
      !   p0*(theta_g(indz(1,1)+1,ix,jy,memind(1))/tth(ix,jy,indz(1,1)+1,memind(1)))**(-1/kappa),&
      !   tint,tth(ix,jy,indz(1,1),memind(1)),tth(ix,jy,indz(1,1)+1,memind(1)), &
      !   theta,theta_g(indz(1,1),ix,jy,memind(1)),theta_g(indz(1,1)+1,ix,jy,memind(1))

      return

      end subroutine interpol_wind_theta_inct

end module ecmwf_inct
!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log:
