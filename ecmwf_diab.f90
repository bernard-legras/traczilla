
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

module ecmwf_diab

!*******************************************************************************
! This modules allow to use temperature tendencies archived at ECMWF as
! vertical velocities.
! Must be combined with isentropic interpolation
! (in fact uses only lociso in the isentrop_m module)
!*******************************************************************************

use commons
use isentrop_m
use netcdf
implicit none

logical, save :: ecmwf_diabatic_w

! other variables now declared in commons which are also
! used by merra
!logical clear_sky, cloud_sky
!logical mass_correction, mean_diab_output

! ecmwf_diabatic_w        Specify that diabatic velocities are extracted
!                         from ECMWF temperature tendencies
! clear_sky		  Use of clear sky tendency (LW+SW)
! cloud_sky               Use of cloud sky radiative tendency (LW+SW) no latent
!                         heat
! mass_correction         apply a correction to heating on isobaric surfaces
!                         such that mean mass flux is zero across this surface
! mean_diab_output	  output of the mean mass flux before correction	

!integer numbwf_diab,wftime_diab(maxwf),lwindinterv_diab
!character(len=11):: wfname_diab(maxwf),wfspec_diab(maxwf)

! lwindinterv_diab [s]    Interval between wind fields currently in memory
! numbwf_diab             actual number of wind fields
! wftime_diab(maxwf) [s]  times relative to beginning time of wind fields
! wfname_diab(maxwf)      file names of wind fields
! wfspec_diab(maxwf)      specifications of wind field file, e.g. if on hard 
!                         disc or on tape

!integer memtime_diab(2),memind_diab(2)

! memtime_diab [s]        validation times of wind fields in memory
! memind_diab             pointer to wind field, in order to avoid sh!real, allocatable :: w_diab(:,:,:,:)

! w_diab [K/s]            vertical velocity as tendency in potential temperatureuffling
!                         of wind fields


 
!character(len=128):: path_diab(2)
!integer len_diab(2)

! path_diab               paths for diabatic directory and associated AVAILABLE file
! len_diab                lengths of the two previous strings

real(dp), save, allocatable :: area_coefft(:),pmc(:),pif(:)
integer, save :: NPureP

! coefficients for the correction of the mass conservation across pressure surfaces

contains

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine alloc_ecmwf_diab

print *,'alloc_ecmwf_diab w_diab'
allocate (w_diab(0:nx-1,0:ny-1,nuvz,2))
w_diab(:,:,:,:) = MISSING
end subroutine alloc_ecmwf_diab


!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine getfields_diab(itime,nstop)

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
! idiffmax             maximum allowable time difference between 2 wind fields *
!                                                                            *
!*******************************************************************************

      integer :: indj,indmin,itime,nstop,memaux
      save :: indmin

      data indmin/1/
 
! Check, if wind fields are available for the current time step
!**************************************************************

      nstop=0

      if ((ldirect*wftime_diab(1).gt.ldirect*itime).or.    &
        (ldirect*wftime_diab(numbwf_diab).lt.ldirect*itime)) then
        write(*,*) 'FLEXPART WARNING: NO DIABATIC VELOCITIES ARE AVAILABLE.'
        write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
        nstop=4
        return
      endif

      if ((ldirect*memtime_diab(1).le.ldirect*itime).and.   &
        (ldirect*memtime_diab(2).gt.ldirect*itime)) then

! The right wind fields are already in memory -> don't do anything
!*****************************************************************

        continue

      else if ((ldirect*memtime_diab(2).le.ldirect*itime).and. &
        (memtime_diab(2).ne.999999999)) then
 
! Current time is after 2nd wind field
! -> Resort wind field pointers, so that current time is between 1st and 2nd
!***************************************************************************

        memaux=memind_diab(1)
        memind_diab(1)=memind_diab(2)
        memind_diab(2)=memaux
        memtime_diab(1)=memtime_diab(2)

! Read a new wind field and store it on place memind(2)
!******************************************************

        do indj=indmin,numbwf_diab-1
           if (ldirect*wftime_diab(indj+1).gt.ldirect*itime) then
              call read_diab(indj+1,memind_diab(2))
              memtime_diab(2)=wftime_diab(indj+1)
              call verttransform_diab(memind_diab(2),2)
              write(*,'(a,a,a,i11,a,i11)') &
                 ' getfields_diab> ',trim(wfname_diab(indj+1)), &
                 '  memtime ',memtime_diab(2),'  time ',itime 
              nstop = 1
              goto 40
           endif
       enddo
 40    indmin=indj

      else

! No wind fields, which can be used, are currently in memory 
! -> read both wind fields
!***********************************************************

         do indj=indmin,numbwf_diab-1
            if ((ldirect*wftime_diab(indj).le.ldirect*itime).and.   &
                  (ldirect*wftime_diab(indj+1).gt.ldirect*itime)) then
               memind_diab(1)=1
               call read_diab(indj,memind_diab(1))
               memtime_diab(1)=wftime_diab(indj)
               write(*,'(a,a,a,i11,a,i11)') &
                 ' getfields_diab> ',trim(wfname_diab(indj)), &
                 '  memtime ',memtime_diab(1),'  time ',itime 
               memind_diab(2)=2
               call read_diab(indj+1,memind_diab(2))
               memtime_diab(2)=wftime_diab(indj+1)
               write(*,'(a,a,a,i11,a,i11)') &
                 ' getfields_diab> ',trim(wfname_diab(indj+1)), &
                 '  memtime ',memtime_diab(2),'  time ',itime 
               call verttransform_diab(memind_diab(1),1)
               call verttransform_diab(memind_diab(2),2)
               nstop = 1
               goto 60
            endif
         enddo
 60      indmin=indj

      endif

      lwindinterv_diab=abs(memtime_diab(2)-memtime_diab(1))
 
      if (lwindinterv_diab.gt.idiffmax) then
        print*,'getfields_diab> lwindinterv_diab, idiffmax ',lwindinterv_diab, idiffmax
        nstop=3
      endif

      return
      end subroutine getfields_diab

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine read_diab(indj,n)
      integer, intent(in):: indj,n
      
      if (ecmwf_netcdf) then 
        call read_diab_nc(indj,n)
      else 
        call read_diab_grb(indj,n)
      endif
      
      return
      end subroutine read_diab
      
      subroutine read_diab_grb(indj,n)
      
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READ_DIAB                   *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      B. LEGRAS                                  *
!             DATE:        2005   
!             Rewritten using grib_index and new ecccodes: Dec 2017   *
!                                                                     * 
!**********************************************************************
!  FROM READWIND
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
! IMPORTANT VARIABLES FROM COMMONS:                          *
!                                                                     *
! wfname_diab        File name of data to be read in                  *
! nx,ny,nuvz,nwz     expected field dimensions                        *
! nlev_ec            number of vertical levels ecmwf model            *
! w_diab             temperature tendency over x h                    *
!                                                                     *
!**********************************************************************
      use grib_api
      integer, intent(in):: indj,n
      
      integer :: count1, count2
      integer :: iret, idx, igrib, ifield
      integer :: nxf, nyf, nlevf
      integer :: i,j,k, off_bot
      integer :: parLev, paramSize, par1, par2
      integer, allocatable :: param(:)
      real(kind=4), allocatable:: zsec4(:)
      character (len=20):: gribFunction = 'read_diab_gr'
!     table of the paramID for ERA-Interim and ERA5
!     respectively sw, lw, swcs, lwcs, all
      integer :: table(5), tableEI(5), tableEA(5)       
      data tableEI /162100,162101,162102,162103,162110/
      data tableEA /235001,235002,235003,235004,235005/
      logical :: twopar
      allocate(zsec4(jpunp))
      
! INITIALIZING

      w_diab(0:nxfield-1,0:ny-1,:,n) = MISSING
      
      if (u_bot) then
        off_bot=1
      else
        off_bot=0
      endif
      
! GENERATE THE INDEX
      
      call grib_index_create(idx,path_diab(1)(1:len_diab(1))//wfname_diab(indj),'paramId',iret)
      if(iret.ne.GRIB_SUCCESS) go to 999
      call grib_index_get_size(idx,'paramId',paramSize)
      allocate(param(paramSize))
      call grib_index_get(idx,'paramId',param)
      
! DETERMINE WHICH TABLE TO BE USED

      if (any(param > 200000)) then
        table(:) = tableEA(:)
      else
        table(:) = tableEI(:)                   
      endif
      deallocate(param) 
      
! SELECT THE PARAMETERS TO BE READ

      if (clear_sky) then
        twopar = .true.
        par1 = table(3)
        par2 = table(4)
      else if (cloud_sky) then
        twopar = .true.
        par1 = table(1)
        par2 = table(2)
      else
        twopar = .false.
        par1 = table(5)
      endif
      
! GET FIRST FIELD
 
      call grib_index_select(idx,'paramId',par1)
      count1 = 0
      call grib_new_from_index(idx,igrib,iret)
      loop1: do while (iret /= GRIB_END_OF_INDEX)
        count1 = count1+1
        call grib_get_int(igrib,'level',parLev,iret)
      
! check grid specifications

        if (count1==1) then
          call grib_get_int(igrib,'numberOfPointsAlongAParallel',nxf,iret)
          call grib_get_int(igrib,'numberOfPointsAlongAMeridian',nyf,iret)
          call grib_get_int(igrib,'numberOfVerticalCoordinateValues',nlevf,iret)
          if (nxf/=nxfield) stop 'READDIABGRB: NX NOT CONSISTENT'
          if (nyf/=ny)      stop 'READDIABGRB: NY NOT CONSISTENT'  
          if (nlevf/2-1/=nlev_ec) stop 'READDIAB VERTICAL LEVELS NON CONSISTENT'
        endif
        
! read the data

        call grib_get_real4_array(igrib,'values',zsec4,iret)
        call grib_check(iret,gribFunction,'ERROR reading values')                 

!     ACTHTUNG: The diabatic velocities are stored like horizontal wind
!     leaving here the bottom level defined to zero in absence of other
!     prescription. It is probably not a good idea, anyway, to use
!     diabatic velocities in the boundary layer.
!     It is important for DIAB runs to set properly upper and lower
!     theta levels in COMMAND file 

        k = parLev
        do j=0,ny-1 ; do i=0,nxfield-1
              w_diab(i,j,nlev_ec-k+1+off_bot,n) = zsec4(nxfield*(ny-j-1)+i+1)
        enddo       ; enddo
        
        call grib_release(igrib)
        call grib_new_from_index(idx,igrib,iret)
        
      enddo loop1
      
! GET SECOND FIELD      
      
      if (twopar) then
      
        call grib_index_select(idx,'paramId',par2)
        count2 = 0
        call grib_new_from_index(idx,igrib,iret)
        loop2: do while (iret /= GRIB_END_OF_INDEX)
          count2 = count2+1
          call grib_get_int(igrib,'level',parLev,iret)
     
! read the data

          call grib_get_real4_array(igrib,'values',zsec4,iret)
          call grib_check(iret,gribFunction,'ERROR reading values')                 

          k = parLev
          do j=0,ny-1 ; do i=0,nxfield-1
              w_diab(i,j,nlev_ec-k+1+off_bot,n) = w_diab(i,j,nlev_ec-k+1+off_bot,n) + &
                                        zsec4(nxfield*(ny-j-1)+i+1)
          enddo       ; enddo
        
        call grib_release(igrib)
        call grib_new_from_index(idx,igrib,iret)
        
        enddo loop2     
      
      endif

! CLOSING OF INPUT DATA FILE

      call grib_index_release(idx)
      
! For global fields, assign rightmost grid point the value of the
! leftmost point
!****************************************************************

      if (xglobal) then
        w_diab(nx-1,:,:,n) = w_diab(0,:,:,n)
      endif

      deallocate (zsec4)

      if(clear_sky.or.cloud_sky) then
        if(count1/=count2) go to 889
!        print *,'read_diab> file ',wfname_diab(indj), count102, count103
      else
        if(count1==0) go to 889
!        print *,'read_diab> file ',wfname_diab(indj), count110
      endif
      !print *,count100,count101,count110

      return    

889   write(*,*) ' #### TRACZILLA MODEL ERROR! DIABFIELD        #### '
      write(*,*) ' #### ',wfname_diab(indj),'                   #### '
      write(*,*) ' #### does not contain required fields        #### '
      print *,clear_sky,cloud_sky
      print *,count1,count2
      stop 'Execution terminated'
999   write(*,*) ' #### TRACZILLA MODEL ERROR! DIABFIELD        #### ' 
      write(*,*) ' #### ',wfname_diab(indj),'                   #### '
      write(*,*) ' #### CANNOT BE OPENED !!!                    #### '
      stop 'Execution terminated'      

      return
      end subroutine read_diab_grb
 
 !**********************************************************************
      
      subroutine check(etat,code)
      integer, intent(in):: etat, code
      if(etat /= nf90_noerr) then
        print *,'Error in read_diab_nc'
        print *, trim(nf90_strerror(etat)),code
        stop 2
      endif
      end subroutine check
      
 !**********************************************************************     
      
      subroutine read_diab_nc(indj,n)
      
      integer, intent(in):: indj,n
      integer :: ncid
      integer :: NbLev, NbLon, NbLat, LevId, LonId, LatId, RVId
      integer :: i,j,k
      integer, allocatable :: Levlist(:)
      real(kind=4), allocatable :: buffer(:,:,:)
      
      ! Open file
      call check(NF90_OPEN(path_diab(1)(1:len_diab(1))//wfname_diab(indj)//'.nc4',&
          NF90_NOWRITE,ncid),1)
      ! get long Id & number of longitudes
      call check(NF90_INQ_DIMID(ncid,'longitude',LonId),0)
      call check(NF90_INQUIRE_DIMENSION(ncid,LonId,LEN=NbLon),1)
      if (NbLon/=nxfield) stop 'READDIABNC: NX NOT CONSISTENT'
      ! get lat Id & number of latitudes
      call check(NF90_INQ_DIMID(ncid,'latitude',LatId),0)
      call check(NF90_INQUIRE_DIMENSION(ncid,LatId,LEN=NbLat),1)
      if (NbLat/=ny)      stop 'READDIABNC: NY NOT CONSISTENT'
      ! TODO: some check that these quantities are in pace with nx, ny
      ! get level Id & number of levels
      call check(NF90_INQ_DIMID(ncid, 'level', LevId),2)
      call check(NF90_INQUIRE_DIMENSION(ncid, LevId, LEN=NbLev),3)
      ! read level list
      allocate(LevList(NbLev))
      call check(NF90_INQ_VARID(ncid, 'level', LevId),4)
      call check(NF90_GET_VAR(ncid,LevId,LevList),5)
      ! allocate the buffer
      allocate(buffer(0:nx-1,0:ny-1,NbLev))
      ! read the data
      if (clear_sky) then
        call check(NF90_INQ_VARID(ncid, 'mttswrcs', RVId),6)
        call check(NF90_GET_VAR(ncid,RVId,buffer,&
          start=(/1,1,1,1/),count=(/NbLon,NbLat,NbLev,1/)),7)
        do k=1,NbLev; do j=0,ny-1 ; do i=0,nxfield-1
          w_diab(i,ny-j-1,nlev_ec-LevList(k)+1,n) = buffer(i,j,k)
        enddo;        enddo       ; enddo
        call check(NF90_INQ_VARID(ncid, 'mttlwrcs', RVId),8)
        call check(NF90_GET_VAR(ncid,RVId,buffer,&
          start=(/1,1,1,1/),count=(/NbLon,NbLat,NbLev,1/)),9) 
        do k=1,NbLev; do j=0,ny-1 ; do i=0,nxfield-1
          w_diab(i,ny-j-1,nlev_ec-LevList(k)+1,n) = buffer(i,j,k) &
             + w_diab(i,ny-j-1,nlev_ec-LevList(k)+1,n)    
        enddo;        enddo       ; enddo
      else if (cloud_sky) then
        call check(NF90_INQ_VARID(ncid, 'mttswr',RVId),10)
        call check(NF90_GET_VAR(ncid,RVId,buffer,&
          start=(/1,1,1,1/),count=(/NbLon,NbLat,NbLev,1/)),11)
        do k=1,NbLev; do j=0,ny-1 ; do i=0,nxfield-1
          w_diab(i,ny-j-1,nlev_ec-LevList(k)+1,n) = buffer(i,j,k)
        enddo;        enddo       ; enddo
        call check(NF90_INQ_VARID(ncid, 'mttlwr', RVId),12)
        call check(NF90_GET_VAR(ncid,RVId,buffer,&
          start=(/1,1,1,1/),count=(/NbLon,NbLat,NbLev,1/)),13) 
        do k=1,NbLev; do j=0,ny-1 ; do i=0,nxfield-1
          w_diab(i,ny-j-1,nlev_ec-LevList(k)+1,n) = buffer(i,j,k) &
             + w_diab(i,ny-j-1,nlev_ec-LevList(k)+1,n)   
        enddo;        enddo       ; enddo
      else
        call check(NF90_INQ_VARID(ncid, 'mttpm', RVId),14)
        call check(NF90_GET_VAR(ncid,RVId,buffer,&
          start=(/1,1,1,1/),count=(/NbLon,NbLat,NbLev,1/)),15)
        do k=1,NbLev; do j=0,ny-1 ; do i=0,nxfield-1
          w_diab(i,ny-j-1,nlev_ec-LevList(k)+1,n) = buffer(i,j,k)
        enddo;        enddo       ; enddo
      endif    
      ! Close file and deallocate the buffer
      call check(NF90_CLOSE(ncid),16)
      deallocate(buffer,LevList)
      
      ! For global fields, assign rightmost grid point the value of the
      ! leftmost point

      if (xglobal) then
        w_diab(nx-1,:,:,n) = w_diab(0,:,:,n)
      endif

      return         
      end subroutine read_diab_nc    


!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine readpaths_diab(pathfile,error)
!                            o
!*******************************************************************************
!                                                                              *
!     Reads the pathnames, where input/output files are expected to be.        *
!     The file pathnames must be available in the current working directory.   *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     1 February 1994                                                          *
!
!     modified by B. Legras
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

      integer :: i
      logical, intent(out)::error

      character(len=*), intent(in):: pathfile

      error=.false.

! Read the pathname information stored in unitpath
!*************************************************

      open(unitpath,file=pathfile,status='old',err=999)
      print *,'opening pathfile ',trim(pathfile)

!     skip previously read data
      do i=1,numpath
        read(unitpath,*)
      enddo
      do i=1,2
        read(unitpath,'(a)',err=998) path_diab(i) 
        len_diab(i)=index(path_diab(i),' ')-1
      enddo

      close(unitpath)
      return    

998   write(*,*) ' #### TRAJECTORY MODEL ERROR! ERROR WHILE     #### ' 
      write(*,*) ' #### READING DIABATIC FILE PATHNAMES.        #### ' 

999   write(*,*) ' #### TRAJECTORY MODEL ERROR! PATH FILE       #### ' 
      write(*,*) ' #### CANNOT BE OPENED IN THE CURRENT WORKING #### '
      write(*,*) ' #### DIRECTORY.                              #### '
      error=.true.

      return
      end subroutine readpaths_diab

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine readavailable_diab(error)
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
! idiffmax [s]         maximum allowable time between 2 T tendency fields      *
! jul                  julian date, help variable                              *
! numbwf_diab          actual number of T tendency fields                      *
! wfname_diab(maxwf)   file names of needed T tendency fields                  *
! wfspec_diab(maxwf)   file specifications of wind fields (e.g., if on disc)   *
! wftime_diab(maxwf) [s] times of T tendency fields relative to beginning time *
! wfname1,wfspec1,wftime1 = same as above, but only local (help variables)     *
!                                                                              *
! Constants:                                                                   *
! maxwf                maximum number of wind fields                           *
! unitavailab          unit connected to file AVAILABLE                        *
!                                                                              *
!*******************************************************************************

      use date
      logical :: error, toobig
      integer :: i,j,idiff,ldat,ltim
      integer :: year,month,day
      real(dbl) :: jul,beg,end
      character(len=16):: fname,spec
      character(len=16), allocatable:: wfname1(:),wfspec1(:)
      integer, allocatable:: wftime1(:)

      allocate(wfname1(maxwf),wfspec1(maxwf),wftime1(maxwf))

      error=.false.

! Windfields are only used, if they are within the modelling period.
! However, 1 additional day at the beginning and at the end is used for
! interpolation. -> Compute beginning and ending date for the windfields.
!************************************************************************

      if (ideltas.gt.0) then         ! forward trajectories
        beg=bdate-1.                 ! idiffmax should be used here too               
        end=bdate+dble(float(ideltas)/86400._dp) &
                 +dble(float(idiffmax)/86400._dp)
      else                           ! backward trajectories
        beg=bdate+dble(float(ideltas)/86400._dp) &
                 -dble(float(idiffmax)/86400._dp)
        end=bdate+1.                 ! idiffmax shoukd be used here too
      endif

! Open the wind field availability file and read available wind fields
! within the modelling period.
!*********************************************************************

      open(unitavailab,file=path_diab(2)(1:len_diab(2)),status='old',err=999)

      do i=1,3
        read(unitavailab,*)
      enddo
      
      numbwf_diab=0
100   read(unitavailab,'(i8,1x,i6,3x,a16,3x,a10)',end=99) &
           ldat,ltim,fname,spec
      jul=juldate(ldat,ltim)
      if ((jul.ge.beg).and.(jul.le.end)) then
        numbwf_diab=numbwf_diab+1
        if (numbwf_diab.gt.maxwf) then      ! check exceedance of dimension
           write(*,*) 'Number of needed wind fields is too great.'
           write(*,*) 'actual, allowed ',numbwf_diab, maxwf
           write(*,*) 'Reduce modelling period (file "COMMAND") or'
           write(*,*) 'reduce number of wind fields (file "AVAILABLE").'
           goto 1000
        endif

!  set list of filenames and times
!  offset time by +1h30 as indicated in AVAILABLE takes 
!  into account the fact that temperature
!  tendencies are accumulated over a time interval of 3h
        wfname1(numbwf_diab)=adjustl(fname)
        wfspec1(numbwf_diab)=trim(spec)
        wftime1(numbwf_diab)=nint((jul-bdate)*86400._dp)
      endif
      goto 100       ! next wind field

99    continue

      close(unitavailab)

! Check wind field times of file AVAILABLE (expected to be in temporal order)
!****************************************************************************

      if (numbwf_diab.eq.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! NO WIND FIELDS    #### '
        write(*,*) ' #### AVAILABLE FOR SELECTED TIME PERIOD.     #### '
        error=.TRUE.
        return
      endif

      do i=2,numbwf_diab
        if (wftime1(i).le.wftime1(i-1)) then
          write(*,*) 'FLEXPART ERROR: FILE AVAILABLE-DIAB IS CORRUPT.'
          write(*,*) 'THE TEMPERATURE TENDENCIES ARE NOT IN TEMPORAL ORDER.'
          write(*,*) 'PLEASE CHECK FIELD ',wfname1(i)
          error=.TRUE.
          return
        endif
      enddo

! For backward trajectories, reverse the order of the windfields
!***************************************************************

      if (ideltas.ge.0) then
        wfname_diab(1:numbwf_diab)=wfname1(1:numbwf_diab)
        wfspec_diab(1:numbwf_diab)=wfspec1(1:numbwf_diab)
        wftime_diab(1:numbwf_diab)=wftime1(1:numbwf_diab)
      else
        do i=1,numbwf_diab
          wfname_diab(numbwf_diab-i+1)=wfname1(i)
          wfspec_diab(numbwf_diab-i+1)=wfspec1(i)
          wftime_diab(numbwf_diab-i+1)=wftime1(i)
        enddo
      endif

! Check the time difference between the wind fields. If it is big, 
! write a warning message. If it is too big, terminate the trajectory. 
!*********************************************************************

      wftime1(:)=wftime_diab(:)
      do i=2,numbwf_diab
        idiff=abs(wftime_diab(i)-wftime_diab(i-1))
        toobig=.false.
        if (idiff.gt.idiffmax) then
!         detect end of february for bissextil year in perpetual run
!         shift all remaining times in such case to satisfy check done in
!         getfield
          if(perpetual) then
            ! check done on unshifted times to process several bissextil years
            if(ideltas > 0) call caldate(wftime1(i-1)/86400._dp+bdate,ldat,ltim)
            if(ideltas < 0) call caldate(wftime1(i  )/86400._dp+bdate,ldat,ltim)
            year=int(ldat/10000)
            month=int((ldat-10000*year)/100)
            day=int(ldat-10000*year-100*month)
            if(month==2.and.day==28.and.mod(year,4)==0) then
              print *,'readavailable_diab > jump due to bissextil year in perpetual run'
              print *,'readavailable_diab > apply one-day shift to all remaining times'
              if(ideltas > 0) then
                do j=i,numbwf_diab
                  wftime_diab(j)=wftime_diab(j)-86400
                enddo
              else
                do j=i,numbwf_diab
                  wftime_diab(j)=wftime_diab(j)+86400
                enddo
              endif
            else
              toobig=.true.
            endif
          endif
          if(.not.perpetual.or.toobig) then
            write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
            write(*,*) 'WIND FIELDS IS TOO BIG FOR TRANSPORT CALCULATION.'
            write(*,*) 'THEREFORE, TRAJECTORIES HAVE TO BE SKIPPED.'
            print *, wfname_diab(i),wfname_diab(i-1),idiffmax
          endif
        else if (idiff.gt.idiffnorm) then
          write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
          write(*,*) 'DIAB FIELDS IS BIG. THIS MAY CAUSE A DEGRADATION'
          write(*,*) 'OF SIMULATION QUALITY.'
          print *,  wfname_diab(i),wfname_diab(i-1),idiff,idiffnorm
        endif
      enddo


! Reset the times of the wind fields that are kept in memory to no time
!**********************************************************************

      do i=1,2
        memind_diab(i)=i
        memtime_diab(i)=999999999
      enddo

      deallocate (wfname1,wfspec1,wftime1)

      print *,'readavailable_diab> done ',numbwf_diab

      return    

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE #### '
      write(*,'(a)') '     '//path_diab(2)(1:len_diab(2)) 
      write(*,*) ' #### CANNOT BE OPENED           #### '
1000  error=.true.

      return
      end subroutine readavailable_diab

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine verttransform_diab(n,m)

!*******************************************************************************
!                                                                              *
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
! w_diab(0:nxmax,0:nymax,nuvz,2)  potential temperature tendencies [K/s]       *
! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]    
! 
! Modif: 27/02/2013, B. Legras
! As the heating rate is read between two pressure fields, the surface pressure
! used in the calculation of w_diab should be the average of the two adjacent
! time. The previous pressure was taken with last index n which means 1:30
! before or after
!                                                                              *
!*******************************************************************************
  integer, intent(in) :: n,m
  integer :: ix, jy, iz
  real(dp) :: pp, ff

! Loop over the whole grid
!*************************

! This is not a good way to distinguish ERA40/ERA-Interim from ERA5
  if (hour_accu==3) then
    ff = 1./10800._dp
  else
    ff = 1.
  endif
  do jy=0,ny-1
    do ix=0,nx-1
      do iz=2, nuvz
        pp=akz(iz)+bkz(iz)*(ps(ix,jy,1,1)+(ps(ix,jy,1,2)))/2
        w_diab(ix,jy,iz,n) = ff*w_diab(ix,jy,iz,n)*((p0/pp)**kappa)
      enddo
    enddo
  enddo

! Set the tendency at 10m to be that at next level
  w_diab(:,:,1,n)=w_diab(:,:,2,n)

! Perform mass correction if required
  if(mass_correction) call diab_mass(n,m)

  return
  end subroutine verttransform_diab

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass(n,m)

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
! 28/02/2013
!
!*******************************************************************************

  integer, intent(in) :: n,m
  real(dp), allocatable :: theta(:,:,:), sigma(:,:,:), flux(:,:), flux_lat(:)
  real(dp), allocatable :: mass_flux(:), mean_sigma(:), mean_sigma_lat(:)
  real(dp), allocatable :: mean_w(:)
  integer :: k

  allocate(theta(0:nx-1,0:ny-1,nuvz))
  theta(:,:,:) = MISSING
  
! Calculation of the potential temperature
  do k=nuvz-NPureP,nuvz
     theta(:,:,k)=0.5*(tth(:,:,k,1)+tth(:,:,k,2))*pif(k)
!    theta(:,:,k)=0.5*(theta_g(:,:,k,1)+theta_g(:,:,k,2))
  enddo
! 
  allocate(sigma(0:nx-1,0:ny-1,nuvz))
  sigma(:,:,:) = MISSING
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
     flux(:,:)=w_diab(:,:,k,n)*sigma(:,:,k)
     flux_lat(:)=sum(flux(0:nx-2,:),DIM=1)
     mean_sigma_lat(:)=sum(sigma(0:nx-2,:,k),DIM=1)
     mass_flux(k)=0.5*dot_product(flux_lat(:),area_coefft)
     mean_sigma(k)=0.5*dot_product(mean_sigma_lat(:),area_coefft)
     mean_w(k)=mass_flux(k)/mean_sigma(k)
  enddo
  deallocate(flux,flux_lat)

! Output of the ass averaged diab heating
  if(mean_diab_output) then
    write(unitflux) memtime_diab(m),mean_w(nuvz-NPureP+1:nuvz)
    flush(unitflux)
  endif

! Correction step

!if(m==1) then ! old usage of m to print tests 
!    print *,'diab_mass min max : ',minval(mean_w(nuvz-NPureP+1:nuvz)),maxval(mean_w(nuvz-NPureP+1:nuvz))
!    print *,'diab_mass sigma   : ',minval(mean_sigma(nuvz-NPureP+1:nuvz)),maxval(mean_sigma(nuvz-NPureP+1:nuvz))
!    print *,'diab_mass w_diab  : ',minval(w_diab(:,:,nuvz-NPureP+1:nuvz,n)),maxval(w_diab(:,:,nuvz-NPureP+1:nuvz,n))
!    print *,'mean_w'
!    write(*,'(8g11.3)') mean_w(nuvz-NPureP+1:nuvz)
!    print *,'mass_flux'
!    write(*,'(8g11.3)') mass_flux(nuvz-NPureP+1:nuvz)
!    print *,'mean_sigma'
!    write(*,'(8g11.3)') mean_sigma(nuvz-NPureP+1:nuvz)
!endif

  do k=nuvz-NPureP+1,nuvz
    w_diab(:,:,k,n)=w_diab(:,:,k,n)-mean_w(k)
  enddo

  deallocate(theta,sigma,mean_sigma,mean_sigma_lat)
  deallocate(mass_flux,mean_w)

  return
  end subroutine diab_mass

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass_init

!*************************************************************************************
! This routine initializes coefficients used in the calculation of the masss correction 
! The coefficients are p delta log(p) for the calculation of the vertical derivative
! and sin(phi + dphi) - sin(phi - dphi) = 2 sin(dphi) cos(phi) for the surface integral
! where phi si a grid latitude and dphi=0.5*dy. This quantity is the diffference 
! (up to factor 2 pi which is not accounted) of area
! between to caps bouded by phi-dphi and phi+dphi latitudes. The first and last values 
! are the area of polar caps of angle dphi/2, that is 2 sin^2(dphi/4)
! Conversion to radian is applied
!*************************************************************************************
   
  integer :: i
  
  NPureP=23
 
  if((upper_theta_level < nuvz).OR. nuvz-NPureP+1 < lower_theta_level) then
      write(*,*) 'NPureP mismatch ',lower_theta_level,upper_theta_level,nuvz,&
                                    NPureP
  endif

  if(.not.(xglobal.and.nglobal.and.sglobal)) then
      write(*,*) ' #### TRACZILLA MODEL ERROR! DIAB_MASS        #### ' 
      write(*,*) ' #### mass flux cannot be balanced            #### '
      write(*,*) ' #### on not global grid                      #### '
      stop  
  endif

! calculate area coefficient assuming pole to pole regular grid
  allocate(area_coefft(0:ny-1))
! it is assumed that dy = 180/(ny-1)
  do i=2,ny-1
     area_coefft(i-1)=2*sin(pi/(ny-1)/2)*cos(pi*(2*i-ny-1)/(ny-1)/2)
  enddo
  area_coefft(0)=2*(sin(pi/(ny-1)/8))**2
  area_coefft(ny-1)=area_coefft(1)
! test : the sum should be equal to 2
  write(*,'(" diab_mass_init > check sum area  ",g10.3 )') sum(area_coefft)

! calculate weighting pressure factors in the vertical
  allocate(pmc(nuvz),pif(nuvz))
  pmc=0.
  pif=0.
  do i=nuvz-NPureP+1,nuvz-1
     pmc=(akz(i)*(Log(akz(i-1))-Log(akz(i+1))))
  enddo
  pmc(nuvz)=(akz(nuvz)*(Log(akz(nuvz-1))-Log(akz(nuvz))))

  do i=nuvz-NPureP+1,nuvz
     pif(i)=(p0/akz(i))**kappa
  enddo

! Create the output file of the mas averaged heating flux
  open(unitflux,file=trim(path(2))//'MassMeanDiab2.dat',form='unformatted',position='append')

  return
  end subroutine diab_mass_init

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine interpol_wind_theta_diab   &
         (itime,xt,yt, ztheta, dxdt,dydt,dzdt, ngrid, theta_inf, theta_sup, & 
          z_factor,tint,nstop,j)

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
      integer, intent(in):: j
      real(dp), intent(in) :: xt,yt,ztheta
      real(dp), intent(out) :: dxdt,dydt,dzdt,z_factor,tint, theta_inf, theta_sup

! Auxiliary variables needed for interpolation
      real(dp) :: theta
      real(dp) :: u1(2),v1(2),w1(2),dt1,dt2,dtt,tp1(2)
      real(dp) :: dt1_diab,dt2_diab,dtt_diab
      real(dp) :: tr(4,2),trp(4,2),u(4,2),v(4,2),w(4,2),tp(4,2)
      integer :: m,indexh,indexh_diab,indz(4,2)
      integer :: ix,jy,ixp,jyp,i0,j0,idxy
      real(dp) :: ddx,ddy,rddx,rddy,p1,p2,p3,p4
      real(dp) :: psl0,psup0,pinf0,pisup0,piinf0

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

      ! note that dt1 and dt2 are negative in backward time but the product 
      ! with dtt which is also negative is positive and provides the right weight
      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)
      dt1=dt1*dtt
      dt2=dt2*dtt
      dt1_diab=float(itime-memtime_diab(1))
      dt2_diab=float(memtime_diab(2)-itime)
      dtt_diab=1./(dt1_diab+dt2_diab)
      dt1_diab=dt1_diab*dtt_diab
      dt2_diab=dt2_diab*dtt_diab
 
      !if (debug_out) print *,p1,p2,p3,p4
      !if (debug_out) print *,dt1,dt2,dt1_diab,dt2_diab

! Calculates the theta values on the four adjacent columns if required
!*********************************************************************

      !if(ix < 0) then
      !  print *,ix,jy,xt
      !  flush(6)
      !endif

#if defined(PAR_RUN)
#else
      if(.not.theta_col(ix,jy,memind(1)))   call calc_col_theta(ix,jy,memind(1))
      if(.not.theta_col(ix,jy,memind(2)))   call calc_col_theta(ix,jy,memind(2))
      if(.not.theta_col(ix,jyp,memind(1)))  call calc_col_theta(ix,jyp,memind(1))
      if(.not.theta_col(ix,jyp,memind(2)))  call calc_col_theta(ix,jyp,memind(2))
      if(.not.theta_col(ixp,jy,memind(1)))  call calc_col_theta(ixp,jy,memind(1))
      if(.not.theta_col(ixp,jy,memind(2)))  call calc_col_theta(ixp,jy,memind(2))
      if(.not.theta_col(ixp,jyp,memind(1))) call calc_col_theta(ixp,jyp,memind(1))
      if(.not.theta_col(ixp,jyp,memind(2))) call calc_col_theta(ixp,jyp,memind(2))
#endif

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
      
!      if(debug_out) then
!        print *,'interpol'
!        write(*,'("j, lon, lat, theta ",I7,3 F8.2)')j,xlon0+dx*xt,ylat0+dy*yt,theta
!        write(*,'("indz  ",8I4)')indz
!        write(*,'("tr    ",8F8.2)')tr
!        write(*,'("theta_g(1) ",6F8.2)') theta_g(indz(1,1)-2,ix,jy,memind(1)),&
!         theta_g(indz(1,1)-1,ix,jy,memind(1)),theta_g(indz(1,1),ix,jy,memind(1)),& 
!         theta_g(indz(1,1)+1,ix,jy,memind(1)),theta_g(indz(1,1)+2,ix,jy,memind(1)),&
!         theta_g(indz(1,1)+3,ix,jy,memind(1))
!        write(*,'("tth(1)     ",6F8.2)') tth(ix,jy,indz(1,1)-2,memind(1)),&
!         tth(ix,jy,indz(1,1)-1,memind(1)),tth(ix,jy,indz(1,1),memind(1)),& 
!         tth(ix,jy,indz(1,1)+1,memind(1)),tth(ix,jy,indz(1,1)+2,memind(1)),&
!         tth(ix,jy,indz(1,1)+3,memind(1))
!        write(*,'("pres(1)    ",6F8.0)') akz(indz(1,1)-2)+bkz(indz(1,1)-2)*ps(ix,jy,1,memind(1)),&
!                                         akz(indz(1,1)-1)+bkz(indz(1,1)-1)*ps(ix,jy,1,memind(1)),&
!                                         akz(indz(1,1)  )+bkz(indz(1,1)  )*ps(ix,jy,1,memind(1)),&
!                                         akz(indz(1,1)+1)+bkz(indz(1,1)+1)*ps(ix,jy,1,memind(1)),&
!                                         akz(indz(1,1)+2)+bkz(indz(1,1)+2)*ps(ix,jy,1,memind(1)),&
!                                         akz(indz(1,1)+3)+bkz(indz(1,1)+3)*ps(ix,jy,1,memind(1)) 
!        write(*,'("pres(1)rec ",6F8.0)') p0*(tth(ix,jy,indz(1,1)-2,memind(1))/theta_g(indz(1,1)-2,ix,jy,memind(1)))**(1/kappa),&
!                                         p0*(tth(ix,jy,indz(1,1)-1,memind(1))/theta_g(indz(1,1)-1,ix,jy,memind(1)))**(1/kappa),&
!                                         p0*(tth(ix,jy,indz(1,1)  ,memind(1))/theta_g(indz(1,1)  ,ix,jy,memind(1)))**(1/kappa),&
!                                         p0*(tth(ix,jy,indz(1,1)+1,memind(1))/theta_g(indz(1,1)+1,ix,jy,memind(1)))**(1/kappa),&
!                                         p0*(tth(ix,jy,indz(1,1)+2,memind(1))/theta_g(indz(1,1)+2,ix,jy,memind(1)))**(1/kappa),&
!                                         p0*(tth(ix,jy,indz(1,1)+3,memind(1))/theta_g(indz(1,1)+3,ix,jy,memind(1)))**(1/kappa)             
!      endif

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
                    + theta_g(lower_theta_level,ixp,jyp,memind(2))*p4)*dt2)
                
        theta_sup = ((theta_g(upper_theta_level,ix,jy,memind(1))*p1 &
                    + theta_g(upper_theta_level,ix,jyp,memind(1))*p2 &
                    + theta_g(upper_theta_level,ixp,jy,memind(1))*p3 &
                    + theta_g(upper_theta_level,ixp,jyp,memind(1))*p4)*dt1 &
                   + (theta_g(upper_theta_level,ix,jy,memind(2))*p1 &
                    + theta_g(upper_theta_level,ix,jyp,memind(2))*p2 &
                    + theta_g(upper_theta_level,ixp,jy,memind(2))*p3 &
                    + theta_g(upper_theta_level,ixp,jyp,memind(2))*p4)*dt2)
                         
      endif

! Halt trajectories which are too close from lower boundary
! and flag those which reach the upper boundary
!**********************************************************************

      if(minval(indz)<=lower_theta_level) then
          nstop=2
          dxdt=0. ; dydt=0.; dzdt=0.; tint=275.; z_factor=1.
!         print *,minval(indz),ix,jy,theta
!         print *,indz(:,1)
!         print *,indz(:,2)
!         print *,theta_g(10,ix,jy,1),theta_g(10,ixp,jy,1),theta_g(10,ix,jyp,1),theta_g(10,ixp,jyp,1)
!         print *,theta_g(12,ix,jy,1),theta_g(12,ixp,jy,1),theta_g(12,ix,jyp,1),theta_g(12,ixp,jyp,1)
!         print *,theta_g(14,ix,jy,1),theta_g(14,ixp,jy,1),theta_g(14,ix,jyp,1),theta_g(14,ixp,jyp,1)
!         print *,theta_g(16,ix,jy,1),theta_g(16,ixp,jy,1),theta_g(16,ix,jyp,1),theta_g(16,ixp,jyp,1)
!         print *,theta_g(18,ix,jy,1),theta_g(18,ixp,jy,1),theta_g(18,ix,jyp,1),theta_g(18,ixp,jyp,1)
!         print *,theta_g(20,ix,jy,1),theta_g(20,ixp,jy,1),theta_g(20,ix,jyp,1),theta_g(20,ixp,jyp,1)
          return
      endif
      if(maxval(indz)>=upper_theta_level-1) then
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
          indexh_diab=memind_diab(m)
                     
            u(1,m)=(uupol(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + uupol(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            v(1,m)=(vvpol(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + vvpol(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            w(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh_diab)*(log(theta)-log(trp(1,m))) &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh_diab)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            w(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh_diab)*(log(theta)-log(trp(2,m))) &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh_diab)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            w(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh_diab)*(log(theta)-log(trp(3,m))) &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh_diab)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            w(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh_diab)*(log(theta)-log(trp(4,m))) &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh_diab)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)
         
        enddo
        
      case('lin')
      
        do m=1,2
          indexh=memind(m)
          indexh_diab=memind_diab(m)
            
            u(1,m)=(uupol(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + uupol(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            v(1,m)=(vvpol(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + vvpol(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            w(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh_diab)*(theta-trp(1,m))  &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh_diab)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            w(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh_diab)*(theta-trp(2,m))  &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh_diab)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            w(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh_diab)*(theta-trp(3,m))  &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh_diab)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            w(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh_diab)*(theta-trp(4,m))  &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh_diab)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)
          
        enddo
      
      end select 
                
      else     ! non polar region
      
      select case(vert_interpol)
      
      case('log')      
        do m=1,2
          indexh=memind(m)
          indexh_diab=memind_diab(m)

            u(1,m)=(uuh(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) & 
                  + uuh(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            v(1,m)=(vvh(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + vvh(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            w(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh_diab)*(log(theta)-log(trp(1,m))) &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh_diab)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            w(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh_diab)*(log(theta)-log(trp(2,m))) &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh_diab)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            w(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh_diab)*(log(theta)-log(trp(3,m))) &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh_diab)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            w(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh_diab)*(log(theta)-log(trp(4,m))) &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh_diab)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)
          
        enddo

      case('lin')
      
        do m=1,2
          indexh=memind(m)
          indexh_diab=memind_diab(m)

            u(1,m)=(uuh(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m)) &
                  + uuh(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))
            v(1,m)=(vvh(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m)) &
                  + vvh(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))
            w(1,m)=(w_diab(ix ,jy ,indz(1,m)  ,indexh_diab)*(theta-trp(1,m)) &
                  + w_diab(ix ,jy ,indz(1,m)+1,indexh_diab)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))  
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh )*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            w(2,m)=(w_diab(ixp,jy ,indz(2,m)  ,indexh_diab)*(theta-trp(2,m)) &
                  + w_diab(ixp,jy ,indz(2,m)+1,indexh_diab)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            w(3,m)=(w_diab(ix ,jyp,indz(3,m)  ,indexh_diab)*(theta-trp(3,m)) &
                  + w_diab(ix ,jyp,indz(3,m)+1,indexh_diab)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            w(4,m)=(w_diab(ixp,jyp,indz(4,m)  ,indexh_diab)*(theta-trp(4,m)) &
                  + w_diab(ixp,jyp,indz(4,m)+1,indexh_diab)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))

          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)    
        
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

      dxdt=u1(1)*dt2+u1(2)*dt1
      dydt=v1(1)*dt2+v1(2)*dt1
      dzdt=w1(1)*dt2_diab+w1(2)*dt1_diab
      if(AccurateTemp) tint=tp1(1)*dt2+tp1(2)*dt1
      !if(debug_out) write(*,"('j tint p ',I7,F8.2,F8.0)")j,tint,p0*(tint/theta)**(1/kappa)
      !if(debug_out) &
      ! print "('interpol>',i3,' P ',3f7.0,' T ',3f7.2,' TH ',3f7.2)", & ! test
      !   indz(1,1),p0*(theta/tint)**(-1/kappa),&
      !   p0*(theta_g(indz(1,1),ix,jy,memind(1))/tth(ix,jy,indz(1,1),memind(1)))**(-1/kappa),&
      !   p0*(theta_g(indz(1,1)+1,ix,jy,memind(1))/tth(ix,jy,indz(1,1)+1,memind(1)))**(-1/kappa),&
      !   tint,tth(ix,jy,indz(1,1),memind(1)),tth(ix,jy,indz(1,1)+1,memind(1)), &
      !   theta,theta_g(indz(1,1),ix,jy,memind(1)),theta_g(indz(1,1)+1,ix,jy,memind(1))

      return

  end subroutine interpol_wind_theta_diab

end module ecmwf_diab
!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log:
