
!**********************************************************************
! Copyright 2007, 2012, 2013                                          *
! Bernard Legras                                                      *
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

module mass_iso

!*******************************************************************************
! This modules allow to use temperature tendencies archived at ECMWF as
! vertical velocities.
! Must be combined with isentropic interpolation
!*******************************************************************************

use commons
use sphereharmSPE
use coord
implicit none
save

logical iso_mass, mass_diabat, mass_isent

! iso_mass                 Specify that we use data interpolated 
!                          on isentropic levels with mass equilibrated winds 
! mass_isent               Isentropic run (not implemented)
! mass_diabat              Run with diabatic velocities 

integer :: NbTheta
real, allocatable :: ThetaLev(:)
real, allocatable :: grid_lath(:), icos_lath(:), cosl(:), sinl(:)
real, allocatable :: w_iso(:,:,:,:)

! ThetaLev [K]           niveaux theta
! grid_lath (degree)     latitudes of the regular grid (in NH)
! icos_lath              inverse of the cosines of lat grids
! cosl, sinl             cos and sin of regular longitudes starting at 0
! w_iso [K/s]            vertical velocity as tendency in potential temperature
 
contains

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine alloc_iso

print *,'alloc_iso'
print *,nx,ny,nuvz
allocate (w_iso(0:nx-1,0:ny-1,nuvz,2))
print *,'alloc_iso w_iso'

end subroutine alloc_iso


!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine gridcheck_iso(error)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE GRIDCHECK_ISO               
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      B. LEGRAS      
!             from version of gridcheck by G. Wotawa
!             DATE:                                       *
!             LAST UPDATE:                               *
!                                                                     * 
!                                                                     * 
!**********************************************************************
!                                                                     *
! DESCRIPTION:        A MODIFIER                                      *
!                                                                     *
! THIS SUBROUTINE DETERMINES THE GRID SPECIFICATIONS (LOWER LEFT      *
! LONGITUDE, LOWER LEFT LATITUDE, NUMBER OF GRID POINTS, GRID DIST-   *
! ANCE AND VERTICAL DISCRETIZATION OF THE ECMWF MODEL) FROM THE       *
! GRIB HEADER OF THE FIRST INPUT FILE. THE CONSISTANCY (NO CHANGES    *
! WITHIN ONE FLEXPART RUN) IS CHECKED IN THE ROUTINE "READWIND" AT    *
! ANY CALL.                                                           *
!                                                                     *
! OUTPUT       error .true.   - can not read grid specifications      *
!              error .false.  - normal                                *
!
!                                                                     * 
! XLON0                geographical longitude of lower left gridpoint *
! YLAT0                geographical latitude of lower left gridpoint  *
! NX                   number of grid points x-direction              *
! NY                   number of grid points y-direction              *
! DX                   grid distance x-direction                      *
! DY                   grid distance y-direction                      *
! NUVZ                 number of grid points for horizontal wind      *
!                      components in z direction                      *
! NWZ                  number of grid points for vertical wind        *
!                      component in z direction                       *
! sizesouth, sizenorth give the map scale (i.e. number of virtual grid*
!                      points of the polar stereographic grid):       *
!                      used to check the CFL criterion                *
!                                                                     *
!**********************************************************************

  real sizesouth,sizenorth
  integer ifn, lf, Nt, Nlt, Nlg, jl
 
  logical error

  error=.false.

! Reads first or last field according to the direction of integration
  if(ideltas.gt.0) then
    ifn=1
  else
    ifn=numbwf
  endif

! Open the first file
!  print *,ifn,wfname(ifn)
  open(15,file=path(3)(1:len_path(3))//wfname(ifn),status='OLD', &
      form='UNFORMATTED',access='SEQUENTIAL')

! Check that N_trunc, N_long and N_lat fit the values in sphereharm
  read(15) lf
  if(lf.ne.201) then
    print *, 'gridcheck_iso > data file with wrong format ',lf
    error=.true.
    return
  endif
  rewind (15)
  read (15) lf, Nt, Nlg, Nlt, NbTheta
  if( (Nt.ne.N_trunc) ) then
    print *,'gridcheck_iso > data disagree with sphereharm parameters'
    print *, Nt
    error=.true.
    return
  endif
! Get theta levels
  allocate (ThetaLev(NbTheta))
  read(15) ThetaLev
  close(15)
! Fix grid parameters of the run
! assuming global grid at the moment
! and replication of the first longitude
! N_long and NH_lat in sphereharm do not need to be the same
! as when the spectral field has been generated  
  nx = N_long+1
  ny = 2*NH_lat-1
  nuvz=NbTheta
  nwz=NbTheta
  xlon0=0.
  ylat0=-90.
  xglobal=.true.
  dx=360./(nx-1)
  dy=180./(ny-1)
  dxconst=180./(dx*r_earth*pi)
  dyconst=180./(dy*r_earth*pi)
! Imposes south pole
  sglobal=.true.               ! field contains south pole
! Enhance the map scale by factor 3 (*2=6) compared to north-south
! map scale
  sizesouth=6.*(switchsouth+90.)/dy
  call stlmbr(southpolemap,-90.,0.)
  call stcm2p(southpolemap,0.,0.,switchsouth,0.,sizesouth, &
  sizesouth,switchsouth,180.)
  switchsouthg=(switchsouth-ylat0)/dy
! Imposes north pole
  nglobal=.true.               ! field contains north pole
! Enhance the map scale by factor 3 (*2=6) compared to north-south
! map scale
  sizenorth=6.*(90.-switchnorth)/dy
  call stlmbr(northpolemap,90.,0.)
  call stcm2p(northpolemap,0.,0.,switchnorth,0.,sizenorth, &
  sizenorth,switchnorth,180.)
  switchnorthg=(switchnorth-ylat0)/dy
  write(*,*)
  write(*,'(a,2i5,3L3)')' gribcheck_iso> nx,ny,xglobal,nglobal,sglobal ', &
       nx, ny, xglobal, nglobal, sglobal
  write(*,'(a,2f10.2)')' gribcheck> switchsouthg, switchnorthg ', &
       switchsouthg, switchnorthg
  write(*,'(a,2i7)') &
        ' gribcheck> # of vertical levels in ECMWF data: ',  &
         nuvz,nwz
  write(*,'(a)') ' Mother domain:'
  write(*,'(a,f10.1,a1,f10.1,a,f10.1,a,i4)') '  Longitude range:', &
         xlon0,' to ',xlon0+(nx-1)*dx,'   Grid distance: ',dx,' # ',nx
  write(*,'(a,f10.1,a1,f10.1,a,f10.1,a,i4)') '  Latitude range: ',  &
         ylat0,' to ',ylat0+(ny-1)*dy,'   Grid distance: ',dy,' # ',ny
  write(*,*)

! Initialize the spectral calculations
! Calculate addressing tables for the truncature, and other coeffts
! Initialize COMFFT
  call init_ctrl
  call init_trsf_short
! Define latitude grid and its inverse cosine
! Define cos and sin of longitude
  allocate (grid_lath(nh_lat),icos_lath(nh_lat))
  do jl=1,nh_lat-1
     grid_lath(jl) = (jl-1)*dy
     icos_lath(jl) = 1./cos(pi*grid_lath(jl)/180.)
  enddo
  grid_lath(nh_lat) = 90.
  icos_lath(Nh_lat) = huge(1.)
  allocate (sinl(n_long),cosl(n_long))
  do jl=1,n_long
     cosl(jl)=cos((jl-1)*dx*pi/180.)
     sinl(jl)=sin((jl-1)*dx*pi/180.)
  enddo
! Calculate associated Legendre functions on the regular grid
  do jl=1,nh_lat
     call bsharm(plm(:,jl),sin(grid_lath(jl)*pi/180.),N_trunc)
  enddo
  print *,'gribcheck_iso > Spectral initialization done'
  
  return        
  end subroutine gridcheck_iso

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

      subroutine getfields_iso(itime,nstop)

!******************************************************************************
!                                                                              *

!                                                                              *
!*******************************************************************************

!        B. Legras, June 2006 (from getfields)
!        iso version: read by read_iso and no call to verttransform
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
!                                                                              *
! Constants:                                                                   *
! idiffmax             maximum allowable time difference between 2 wind fields *
!                                                                            *
!*******************************************************************************

      integer indj,indmin,itime,nstop,memaux
      save indmin

      data indmin/1/
 
! Check, if wind fields are available for the current time step
!**************************************************************

      nstop=0

      if ((ldirect*wftime(1).gt.ldirect*itime).or.    &
        (ldirect*wftime(numbwf).lt.ldirect*itime)) then
        write(*,*) 'FLEXPART WARNING: NO ISO FIELD ARE AVAILABLE.'
        write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
        write(*,*) ldirect*wftime(1)
        write(*,*) ldirect*itime
        write(*,*) ldirect*wftime(numbwf)
        nstop=4
        return
      endif

      if ((ldirect*memtime(1).le.ldirect*itime).and.   &
        (ldirect*memtime(2).gt.ldirect*itime)) then

! The right wind fields are already in memory -> don't do anything
!*****************************************************************

        continue

      else if ((ldirect*memtime(2).le.ldirect*itime).and. &
        (memtime(2).ne.999999999)) then
 
! Current time is after 2nd wind field
! -> Resort wind field pointers, so that current time is between 1st and 2nd
!***************************************************************************

        memaux=memind(1)
        memind(1)=memind(2)
        memind(2)=memaux
        memtime(1)=memtime(2)

! Read a new wind field and store it on place memind(2)
!******************************************************

        do indj=indmin,numbwf-1
           if (ldirect*wftime(indj+1).gt.ldirect*itime) then
              call read_iso(indj+1,memind(2))
              call verttransform_iso(memind(2))
              memtime(2)=wftime(indj+1)
              print *, 'getfields> reads new field at time ',memtime(2)           
              nstop = 1
              goto 40
           endif
       enddo
 40    indmin=indj

      else

! No wind fields, which can be used, are currently in memory 
! -> read both wind fields
!***********************************************************

         do indj=indmin,numbwf-1
            if ((ldirect*wftime(indj).le.ldirect*itime).and.   &
                  (ldirect*wftime(indj+1).gt.ldirect*itime)) then
               memind(1)=1
               call read_iso(indj,memind(1))
               call verttransform_iso(memind(1))
               memtime(1)=wftime(indj)
               print*,'getfields_iso> reads first field:  time ',memtime(1)
               memind(2)=2
               call read_iso(indj+1,memind(2))
               call verttransform_iso(memind(2))
               memtime(2)=wftime(indj+1)
               print*,'getfields_iso> reads second field: time ',memtime(2)
               nstop = 1
               goto 60
            endif
         enddo
 60      indmin=indj

      endif

      lwindinterv=abs(memtime(2)-memtime(1))
 
      if (lwindinterv.gt.idiffmax) nstop=3

      return
      end subroutine getfields_iso

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine read_iso(indj,n)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READWIND                    *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      B. Legras
!             from FLEXPART version by G. Wotawa                      *
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS INTERPOLATED ON THETA LEVELS *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN SPECTRAL FORMAT    *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!                                                                     *
! IMPORTANT VARIABLES FROM SHARED VARIABLES:                          *
!                                                                     *
! wfname           File name of data to be read in                  *
! nx,ny,nuvz,nwz     expected field dimensions                        *
! nlev_ec            number of vertical levels ecmwf model            *
! w_iso             temperature tendency over 3h                     *
!                                                                     *
!**********************************************************************

  integer, intent(in):: indj,n
  integer :: k, jl
  real, allocatable :: buffer_g(:,:), buffer_s(:,:)
  integer:: infcount, nancount, l
  logical isnan, isinf

! Allocate fields
  allocate(buffer_s(2,nb_deglib),buffer_g(n2_long,n_lat))
! Open file
  open(15,file=path(3)(1:len_path(3))//wfname(indj), status='OLD', &
      form='UNFORMATTED',access='SEQUENTIAL')
  print *,'read_iso> file ',wfname(indj)
! Skip first records (perform check if paranoiac)
  read(15); read(15)
! Loop on levels
  do k = 1, NbTheta
! Read the fields, transform them to grid data
!    Zonal wind
     read(15) buffer_s
!  TEST TEST
     infcount=0 ; nancount=0
     do l=1,nb_deglib
        if(isnan(buffer_s(1,l))) nancount=nancount+1
        if(isinf(buffer_s(1,l))) infcount=infcount+1
        if(isnan(buffer_s(2,l))) nancount=nancount+1
        if(isinf(buffer_s(2,l))) infcount=infcount+1
     enddo
     if(nancount>0) print *,'Level ',k,nancount,'NaN'
     if(infcount>0) print *,'Level ',k,infcount,'Inf'  
! end TEST   

     call spectogrid(buffer_g,buffer_s)

!    set off-equator values
     do jl = 2, nh_lat-1
       uuh(0:nx-2,jl-1,k,n) = buffer_g(1:n_long,n_lat-jl+1) * icos_lath(nh_lat-jl+1)
       uuh(0:nx-2,nh_lat+jl-2,k,n) = buffer_g(1:n_long,jl)  * icos_lath(jl)
     enddo
!    equator
     uuh(0:nx-2,nh_lat-1,k,n) =  buffer_g(1:n_long,1)
!    quick n' dirty way to set polar boundary
     uuh(0:nx-2,0,k,n) = uuh(0:nx-2,1,k,n)
     uuh(0:nx-2,ny-1,k,n) = uuh(0:nx-2,ny-2,k,n)
!    duplicate last longitude
     uuh(nx-1,0:ny-1,k,n) = uuh(0,0:ny-1,k,n)
!    Meridional wind
     read(15) buffer_s
     call spectogrid(buffer_g,buffer_s)
!    set off-equator values
     do jl = 2, nh_lat-1
       vvh(0:nx-2,jl-1,k,n) = buffer_g(1:n_long,n_lat-jl+1) * icos_lath(nh_lat-jl+1)
       vvh(0:nx-2,nh_lat+jl-2,k,n) = buffer_g(1:n_long,jl)  * icos_lath(jl)
     enddo
!    equator
     vvh(0:nx-2,nh_lat-1,k,n) =  buffer_g(1:n_long,1)
!    quick n' dirty way to set polar boundary
     vvh(0:nx-2,0,k,n) = vvh(0:nx-2,1,k,n)
     vvh(0:nx-2,ny-1,k,n) = vvh(0:nx-2,ny-2,k,n)
!    duplicate last longitude
     vvh(nx-1,0:ny-1,k,n) = vvh(0,0:ny-1,k,n)
!    Heating
     read(15) buffer_s
     call spectogrid(buffer_g,buffer_s)
!    set off-equator values
     do jl = 2, nh_lat-1
       w_iso(0:nx-2,jl-1,k,n) = buffer_g(1:n_long,n_lat-jl+1)
       w_iso(0:nx-2,nh_lat+jl-2,k,n) = buffer_g(1:n_long,jl)
     enddo
!    equator
     w_iso(0:nx-2,nh_lat-1,k,n) =  buffer_g(1:n_long,1)
!    polar boundary
     w_iso(0:nx-2,0,k,n) = buffer_g(1:n_long,n_lat)
     w_iso(0:nx-2,ny-1,k,n) =  buffer_g(1:n_long,nh_lat)
!    duplicate last longitude
     w_iso(nx-1,0:ny-1,k,n) = w_iso(0,0:ny-1,k,n)
!    Temperature
     read(15) buffer_s
     call spectogrid(buffer_g,buffer_s)
!    set off-equator values
     do jl = 2, nh_lat-1
       tth(0:nx-2,jl-1,k,n) = buffer_g(1:n_long,n_lat-jl+1)
       tth(0:nx-2,nh_lat+jl-2,k,n) = buffer_g(1:n_long,jl)
     enddo
!    equator
     tth(0:nx-2,nh_lat-1,k,n) =  buffer_g(1:n_long,1)
!    polar boundary
     tth(0:nx-2,0,k,n) = buffer_g(1:n_long,n_lat)
     tth(0:nx-2,ny-1,k,n) =  buffer_g(1:n_long,nh_lat)
!    duplicate last longitude
     tth(nx-1,0:ny-1,k,n) = tth(0,0:ny-1,k,n)
  enddo
  print *,'read_iso max min u v T H' ! corrigé le 14/9/06
  print *,maxval(uuh(:,:,:,n)),minval(uuh(:,:,:,n))
  print *,maxval(vvh(:,:,:,n)),minval(vvh(:,:,:,n))
  print *,maxval(tth(:,:,:,n)),minval(tth(:,:,:,n))
  print *,maxval(w_iso(:,:,:,n)),minval(w_iso(:,:,:,n))
!  print *,uuh(201:202,153:154,2:3,n)

  deallocate (buffer_s,buffer_g)

  return

  end subroutine read_iso

  subroutine verttransform_iso(n)
  
  integer, intent(in):: n
  integer :: ix, jy, iz
  real :: xlon, ylat, uupolaux, vvpolaux

     if (nglobal) then
        do jy=int(switchnorthg)-2,ny-1
          ylat=ylat0+float(jy)*dy
          do ix=0,nx-1
            xlon=xlon0+float(ix)*dx
            do iz=1,nuvz
              call cc2gll(northpolemap,ylat,xlon,uuh(ix,jy,iz,n), &
                vvh(ix,jy,iz,n),uupol(ix,jy,iz,n),                &
                vvpol(ix,jy,iz,n))
            enddo
          enddo
        enddo
        do iz=1,nuvz
          xlon=xlon0+float(nx/2-1)*dx
          call cc2gll(northpolemap,90.,xlon,             &
            uuh(nx/2-1,ny-1,iz,n),vvh(nx/2-1,ny-1,iz,n), &
            uupolaux,vvpolaux)
          uupol(0:nx-1,ny-1,iz,n)=uupolaux
          vvpol(0:nx-1,ny-1,iz,n)=vvpolaux
        enddo
     endif

     if (sglobal) then
        do jy=0,int(switchsouthg)+3
          ylat=ylat0+float(jy)*dy
          do ix=0,nx-1
            xlon=xlon0+float(ix)*dx
            do iz=1,nuvz
              call cc2gll(southpolemap,ylat,xlon,uuh(ix,jy,iz,n), &
                vvh(ix,jy,iz,n),uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
            enddo
          enddo
        enddo
        do iz=1,nuvz
          xlon=xlon0+float(nx/2-1)*dx
          call cc2gll(southpolemap,90.,xlon,             &
            uuh(nx/2-1,0,iz,n),vvh(nx/2-1,0,iz,n), &
            uupolaux,vvpolaux)
            uupol(0:nx-1,0,iz,n)=uupolaux
            vvpol(0:nx-1,0,iz,n)=vvpolaux
        enddo
      endif
 
  return
  end subroutine verttransform_iso

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine interpol_wind_iso&
         (itime,xt,yt, theta, dxdt,dydt,dzdt, ngrid, & 
          z_factor,tint,nstop)

!*******************************************************************************
!                                                                              *
!  This subroutine interpolates the wind data to current trajectory position.  *
!                                                                              *
!    Author: B. Legras (derived from original interpol_wind)                
!                                                                              *
!
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
      real, intent(in) :: xt,yt,theta
      real, intent(out) :: dxdt,dydt,dzdt,z_factor,tint

! Auxiliary variables needed for interpolation
      real u1(2),v1(2),w1(2),dt1,dt2,dtt,tp1(2)
      real tr(4,2),trp(4,2),u(4,2),v(4,2),w(4,2),tp(4,2)
      integer m,indexh,indz(4,2)
      integer ix,jy,ixp,jyp,i0,j0,idxy
      real ddx,ddy,rddx,rddy,p1,p2,p3,p4
      real pisup0,piinf0

!********************************************
! Multilinear interpolation in time and space
!********************************************

! Determine the lower left corner and its distance to the current position
!*************************************************************************

      ! min and max required for the points just falling on the boundary
      ! as it may happen typically as a result of initialization
      ix=floor(xt)
      if(xglobal) ix=modulo(ix,nx-1)       
      jy=max(min(floor(yt),ny-2),0)
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

! Determine the level below the current position for u,v
!*******************************************************

!  Locates the potential temperature in the grid

      indz(1,1) = lociso(theta,1,NbTheta)
      indz(2,1) = indz(1,1)
      indz(3,1) = indz(1,1)
      indz(4,1) = indz(1,1)
      
!  Locates next time by assuming proximity
            

      indz(1,2) = lociso2(theta,indz(1,1))
      indz(2,2) = indz(1,2)
      indz(3,2) = indz(1,2)
      indz(4,2) = indz(1,2)

!  Defines potential temperature at the 8 nearby meshpoints for
!  the two times (quick n' dirty)

      tr (1,1) = ThetaLev(indz(1,1))
      tr (1,2) = ThetaLev(indz(1,2))
      tr (2,1) = tr(1,1)
      tr (2,2) = tr(1,2)
      tr (3,1) = tr(1,1)
      tr (3,2) = tr(1,2)
      tr (4,1) = tr(1,1)
      tr (4,2) = tr(1,2)
      trp(1,1) = ThetaLev(indz(1,1)+1)
      trp(1,2) = ThetaLev(indz(1,2)+1)
      trp(2,1) = trp(1,1)
      trp(2,2) = trp(1,2)
      trp(3,1) = trp(1,1)
      trp(3,2) = trp(1,2)
      trp(4,1) = trp(1,1)
      trp(4,2) = trp(1,2)
      
! Halt trajectories which are too close from lower boundaries
!************************************************************
! parcels cannot exceed the max altitude due to bound in 
! advanceB

      if((minval(indz)==1)) then
         nstop=2
         dxdt=0. ; dydt=0.; dzdt=0.; tint=275.;
         z_factor=1.
         return
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
            w(1,m)=(w_iso(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + w_iso(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            w(2,m)=(w_iso(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + w_iso(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            w(3,m)=(w_iso(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + w_iso(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            w(4,m)=(w_iso(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + w_iso(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)
         
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
            w(1,m)=(w_iso(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m))  &
                  + w_iso(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta))  &
                  / (tr(1,m)-trp(1,m))
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            w(2,m)=(w_iso(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + w_iso(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            w(3,m)=(w_iso(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + w_iso(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            w(4,m)=(w_iso(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + w_iso(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
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
          
            u(1,m)=(uuh(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) & 
                  + uuh(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            v(1,m)=(vvh(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + vvh(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            w(1,m)=(w_iso(ix ,jy ,indz(1,m)  ,indexh)*(log(theta)-log(trp(1,m))) &
                  + w_iso(ix ,jy ,indz(1,m)+1,indexh)*(log(tr(1,m))-log(theta))) &
                  / (log(tr(1,m))-log(trp(1,m)))
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            w(2,m)=(w_iso(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + w_iso(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            w(3,m)=(w_iso(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + w_iso(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            w(4,m)=(w_iso(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + w_iso(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))

          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)
  
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
            w(1,m)=(w_iso(ix ,jy ,indz(1,m)  ,indexh)*(theta-trp(1,m)) &
                  + w_iso(ix ,jy ,indz(1,m)+1,indexh)*(tr(1,m)-theta)) &
                  / (tr(1,m)-trp(1,m))  
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            w(2,m)=(w_iso(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + w_iso(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            w(3,m)=(w_iso(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + w_iso(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            w(4,m)=(w_iso(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + w_iso(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))

          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)    
        
        enddo
      
      end select

      endif

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
      dzdt=(w1(1)*dt2+w1(2)*dt1)*dtt
      if(AccurateTemp) tint=(tp1(1)*dt2+tp1(2)*dt1)*dtt
 
! TEST    
!      if(jp==27385 .and. itime<=-8611200) then
!        print *,'interpol'
!        print *,dxdt,dydt,dzdt
!        print *,u1,V1
!        print *,u(:,2),v(:,2)
!        print *,w1
!        print *,ix,jy,indz(1,1),indz(1,2)
!        print *,uuh(ix:ix+1,jy:jy+1,indz(1,1):indz(1,1)+1,2)
!      endif
! TEST

!***************************************************
! 4.) Calculation of z_factor for vertical diffusion
!***************************************************

      select case (diftype)
      case (1)          ! diffusion in z (cf p.46, book C, part2)
!       d theta / dz = - (g/theta) d theta / d Pi = -g d Log(theta) / d Pi
!       where Pi = Cp (p/p0)**kappa = Cp T/theta
!       estimated from the data on the lower left corner at first time
!       of the interval, for a better estimate using the closest point
!       activate the first following line and deactivate the second one 
!       call sort_hor_distance       
        i0=ix; j0=jy; idxy=1
        pisup0 = cpa * tth(i0,j0,indz(idxy,1)+1,indexh)/trp(idxy,1)
        piinf0 = cpa * tth(i0,j0,indz(idxy,1),indexh)/tr(idxy,1)
        if(debug_out) then 
                print *,idxy
                print *,pisup0,piinf0
                print *,tth(i0,j0,indz(idxy,1),indexh),tth(i0,j0,indz(idxy,1)+1,indexh)
                print *,tr(idxy,1),trp(idxy,1)
        endif
        z_factor = -ga*(log(trp(idxy,1))-log(tr(idxy,1)))/(pisup0-piinf0)
      case (2)          ! diffusion in theta
        z_factor = 1.
      case default
        z_factor = 0.
      end select
     
      return

      end subroutine interpol_wind_iso

      function lociso(theta,ib1,iu1)
      integer, intent(in) :: ib1,iu1
      real, intent(in) :: theta
      integer :: lociso,ib,iu,im
      ib=ib1 ; iu=iu1
      do while(iu-ib>1)
        im = (iu+ib)/2
        if( theta >= ThetaLev(im) ) then
          ib=im
        else
          iu=im
        endif
      enddo
      lociso = ib
      end function lociso
      
      function lociso2(theta,ib1)
      integer, intent(in) :: ib1
      real, intent(in) :: theta
      integer :: lociso2,ib
      ib=ib1
      if (theta < ThetaLev(ib)) then
        do while(ib>1)
          ib = ib-1
          if ( theta >= ThetaLev(ib) ) exit
        enddo
      else
        do while(ib<NbTheta-1)
          if (theta < ThetaLev(ib+1)) exit
          ib = ib+1
        enddo
      endif
      lociso2 = ib
      end function lociso2 

end module mass_iso
!
!=====|==1=========2=========3=========4=========5=========6=========7==
