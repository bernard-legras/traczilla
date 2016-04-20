!**********************************************************************
! Copyright 1996, 1997, 2001, 2002, 2006, 2007, 2012, 2013, 2015      *
! Andreas Stohl, G. Wotawa, Bernard Legras                            *
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

module jra55

!*******************************************************************************
! This modules allow to use JRA-55 data.
! Does not need isentrop;h as interpolation is made in eta
!*******************************************************************************

use commons
use grib_api
use isentrop_h ! needs upper_theta_layer
implicit none
private :: locpj
private :: NbLon, NbLat, NbHLat, NbLev
private :: idx_uvwt, idx_ps, idx_hr
private :: current_year,current_year_hr,current_ym,current_ym_hr
!private :: w1,w2,g1,g2

logical, save :: jra55_data, jra55_diab

integer, save :: NbLon, NbLat, NbHLat, NbLev

integer, save :: NPureP_jra55
real(dp), save, allocatable :: area_coefft_jra55(:), pmc_jra55(:), pif_jra55(:)
!integer :: PSId, UId, VId, TId,  OMEGAId, LWRId, SWRId, LWRCLRId, SWRCLRId
!real :: missing_value
!real, allocatable :: PressLev(:), facT(:), LogPressLev(:), pmc_jra55(:)
!real*8, allocatable :: grid_lat(:), grid_lon(:),grid_ver(:)


character(len=4), save:: current_year,current_year_hr
character(len=6), save:: current_ym,current_ym_hr
integer, save :: idx_uvwt(4), idx_ps, idx_hr(2)
real(dp), save, allocatable :: w1(:),w2(:)
integer, save, allocatable :: g1(:),g2(:)


! LogPressLev           log(p0/PressLev)
! pmc_jra55             delta p factor in the calculation of density sigma
! facT                  factor (p0/PressLev)**kappa
! NPMass                number of upper levels on which mass equilibration is performed
! area_coefft_jra55     area coefft for the surface integral of the mass flux   
! grid_lat (degree)     latitudes of the regular grid
! grid_lon (degree)     longitudes of the regular grid
! w_jra55 [X/s]         vertical velocity either in zlog coordinate (X=zlog)
!                       or potential temperature tendency (X=K) 

contains

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine alloc_jra55

print *,'alloc_jra55 uuh vvh uupol vvpol tth'
print *,nx,ny,nuvz
allocate (wwh(0:nx-1,-1:ny,nuvz,2))
allocate (uuh(0:nx-1,0:ny-1,nuvz,2),vvh(0:nx-1,0:ny-1,nuvz,2))
allocate (uupol(0:nx-1,-1:ny,nuvz,2),vvpol(0:nx-1,-1:ny,nuvz,2))
allocate (tth(0:nx-1,-1:ny,nuvz,2))
allocate (ps(0:nx-1,-1:ny,1,2))

end subroutine alloc_jra55


!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

!  
!  name: jra55.gridcheck_jra55
!  @param
!  @return
!  
  subroutine gridcheck_jra55(error)

!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE GRIDCHECK                   *
!             FLEXPART VERSION -> DO NOT USE IN FLEXTRA, PRAEPRO      *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      G. WOTAWA                                  *
!             DATE:        1997-08-06                                 *
!             LAST UPDATE: 1997-10-10                                 *
!                                                                     * 
!             Update:      1999-02-08, global fields allowed, A. Stohl*
!             
!             Rewritten    21-05-2013  copy new gridcheck 
!                          from flexpart9 with grib2, B. Legras           
!                                                                     * 
!**********************************************************************
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! THIS SUBROUTINE DETERMINES THE GRID SPECIFICATIONS (LOWER LEFT      *
! LONGITUDE, LOWER LEFT LATITUDE, NUMBER OF GRID POINTS, GRID DIST-   *
! ANCE AND VERTICAL DISCRETIZATION OF THE JRA-55 MODEL) FROM THE      *
! GRIB HEADER OF THE FIRST INPUT FILE. THE CONSISTANCY (NO CHANGES    *
! WITHIN ONE FLEXPART RUN) IS CHECKED IN THE ROUTINE "READWIND" AT    *
! ANY CALL.                                                           *
!                                                                     *
! OUTPUT       error .true.   - can not read grid specifications      *
!              error .false.  - normal                                *
!              oronew .true.  - Terrain heights given in grib files   *
!              oronew .false. - Terrain heights not specified in the  *
!                               grib files (old file standard)        *
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
! UVHEIGHT(1)-         heights of gridpoints where u and v are        *
! UVHEIGHT(NUVZ)       given                                          *
! WHEIGHT(1)-          heights of gridpoints where w is given         *
! WHEIGHT(NWZ)                                                        *
!                                                                     *
!**********************************************************************

! note !! is for statements related to variables which exist in flexpart9
! but not in traczilla   B. Legras
! these variables are nxshift, nxminl, nyminl

   use coord
   use grib_api
   logical, intent(out) :: error
   integer :: i,ifn,ll
   real(dbl) :: xaux1,xaux2,xaux3
   real(dp) :: yfirst,ylast,xlast
   real(dp) :: sizesouth,sizenorth
   real(dbl), allocatable :: y_gauss(:),w_gauss(:)
   real(dp), allocatable :: ylat_gauss(:)
   
   integer :: vvv, lat
   integer :: DateSize, TimeSize, ParamSize, LevelSize
   character(len=8), allocatable :: DateList(:),ParamList(:)
   character(len=4), allocatable :: TimeList(:)
   integer, allocatable:: Levelist(:)

   integer :: iret, igrib, nb_var
   integer :: gribVer
   character(len=24) :: gribErrorMsg = 'Error reading grib file'
   character(len=20) :: gribFunction = 'gridcheck'
   character(len=4) :: year
   character(len=2) :: hour,month
   character(len=8) :: fulldate

   error=.false.
!
!   if(ideltas.gt.0) then
!     ifn=1
!   else
!     ifn=numbwf
!   endif
!   print *,'ideltas ifn ',ideltas,ifn
   ifn=1

! OPENING INDEX FILE
!   
   print *,path(3)(1:len_path(3))//trim(wfname(ifn))
! Finds the first date from wfname
   ll=len_trim(wfname(ifn))
   hour=wfname(ifn)(ll-1:ll)
   year=wfname(ifn)(1:4)
   month=wfname(ifn)(5:6)
   current_year=year
   current_ym=year//month
   current_year_hr=year
   current_ym_hr=year//month
   fulldate=wfname(ifn)(1:8)
   nb_var=3
   if (z_motion) nb_var=4
   do i=1,nb_var
     call grib_index_read(idx_uvwt(i), &
     path(3)(1:len_path(3))//year//'/'//'uvwt-'//year//'-'//month//'.gribidx')
   enddo
   call grib_index_read(idx_ps,&
   path(3)(1:len_path(3))//year//'/'//'ps-'//year//'-'//month//'.gribidx')
   do i=1,2
     call grib_index_read(idx_hr(i),&
     path(3)(1:len_path(3))//year//'/'//'hr-'//year//'-'//month//'.gribidx')
   enddo
   print *,'gribcheck_jra55> all index opened'
   
! Test the uvwt index
  
  ! get the number of distinct values of Date in the index
  call grib_index_get_size(idx_uvwt(1),'mars.date',DateSize)
  ! allocate the array to contain the list of distinct Date
  allocate(DateList(DateSize))
  ! get the list of distinct Date from the index
  call grib_index_get(idx_uvwt(1),'mars.date',DateList)
  write(*,'(a,i3)') 'DateSize=',DateSize
  ! get the number of distinct values of Time in the index
  call grib_index_get_size(idx_uvwt(1),'mars.time',TimeSize)
  ! allocate the array to contain the list of distinct Time
  allocate(TimeList(TimeSize))
  ! get the list of distinct Time from the index
  call grib_index_get(idx_uvwt(1),'mars.time',TimeList)
  write(*,'(a,i3)') 'TimeSize=',TimeSize
  print *,TimeList
  ! get the number of distinct values of Param in the index
  call grib_index_get_size(idx_uvwt(1),'mars.param',ParamSize)
  ! allocate the array to contain the list of distinct Param
  allocate(ParamList(ParamSize))
  ! get the list of distinct Param from the index
  call grib_index_get(idx_uvwt(1),'mars.param',ParamList)
  write(*,'(a,i3)') 'ParamSize=',ParamSize
  print *,ParamList
  ! get the number of distinct values of Levelist in the index
  call grib_index_get_size(idx_uvwt(1),'mars.levelist',LevelSize)
  ! allocate the array to contain the list of distinct Levelist
  allocate(Levelist(LevelSize))
  ! get the list of distinct Levelist from the index
  call grib_index_get(idx_uvwt(1),'mars.levelist',Levelist)
  write(*,'(a,i3)') 'LevelSize=',LevelSize
  NbLev=LevelSize
  
! Test reading a time on uvwt index
  call grib_index_select(idx_uvwt(1),'mars.date',fulldate)
  call grib_index_select(idx_uvwt(1),'mars.time',TimeList(1))
  call grib_index_select(idx_uvwt(1),'mars.param',ParamList(1))
  call grib_index_select(idx_uvwt(1),'mars.levelist',1)
  call grib_new_from_index(idx_uvwt(1),igrib, iret)
  ! first check that we read GRIB1
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  if (gribVer/=1) then 
     print *,'gridcheck_jra55 expecting grib 1 file'
     stop 155
  endif
  call grib_get_int(igrib,'table2Version',vvv,iret)
  print *,'Table 2 Version ',vvv
  call grib_get_int(igrib,'indicatorOfParameter',vvv,iret)
  print *,'indicatorOfParameter ',vvv
  call grib_get_int(igrib,'paramId',vvv,iret)
  print *,'paramId ',vvv
  call grib_get_int(igrib,'Ni',NbLon,iret)
  call grib_get_int(igrib,'Nj',NbLat,iret)
  print *,'NbLon NbLat ',NbLon,NbLat
  call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees',xaux1,iret)
  call grib_get_real8(igrib,'longitudeOfLastGridPointInDegrees',xaux2,iret)
  call grib_get_real8(igrib,'iDirectionIncrementInDegrees',xaux3,iret)
  xlon0=xaux1
  xlast=xaux2
  dx=xaux3
  dx=0.5625_dp ! correcting inaccurate value read from file
  dxconst=180._dp/(dx*r_earth*pi)
  if (xlon0 > 180._dp) xlon0=xlon0-360._dp
  print *,'xlon0 xlast dx ',xlon0,xlast,dx
  call grib_get_real4(igrib,'latitudeOfFirstGridPointInDegrees',yfirst,iret)
  call grib_get_real4(igrib,'latitudeOfLastGridPointInDegrees',ylast,iret)
  print *,'yfirst ylast ',yfirst, ylast
  if (yfirst > ylast ) then
    print *,'Gaussian lats from North to South'
  else
    print *,'Gaussian lats from South to North'
  endif
  deallocate(DateList,TimeList,ParamList,Levelist)
  
! Set the vertical grid
  nuvz=NbLev
  nwz=NbLev
  allocate (akz(NbLev),bkz(NbLev))
  call setabl(akz,bkz)
  allocate (etakzlog(NbLev))
  etakzlog=log(akz+bkz)
  akz=p0*akz
  
! Check whether the grid is cyclic in longitudes
  if (xaux2-xaux1+xaux3-360.D0 < 0.001D0) then
    ! field is cyclic
    xglobal=.true.
    ! Add a last longitude that repeats the first for interpolations
    nx=NbLon+1
  else
    xglobal=.false.
    nx=NbLon
  endif
  
! Define the latitude grid (assuming that the grid goes pole to pole) 
  ny = NbLat
  NbHLat=NbLat/2
  ylat0 = min(yfirst,ylast)+0.001
  dy=2*abs(ylat0)/(NbLat-1)
  dyconst=180._dp/(dy*r_earth*pi)
  print *,'ylat0 dy ',ylat0,dy
  ! Calculate the Gaussian grid
  allocate(y_gauss(NbHLat),w_gauss(NbHLat),ylat_gauss(NbLat))
  call gauleg(y_gauss,w_gauss)
  do lat=1,NbHLat
    ylat_gauss(NbHLat+lat)=y_gauss(lat)
    ylat_gauss(NbHlat+1-lat)=-y_gauss(lat)
  enddo
! Calculation of the interpolation coefficients for the regular grid
! from the Gaussian grid 
  allocate (g1(1:Nblat),g2(1:Nblat),w1(1:Nblat),w2(1:Nblat))
  if(ylat0<ylat_gauss(1)) then
     print *,'ylat0 ylat_gauss(1) ',ylat0,ylat_gauss(1)
     print *,'gridcheck_jra55 > ERROR IN THE LAT GRID'
     stop 2815
  endif
  do lat=1,NbHLat
    g1(NbHLat+lat)=NbHLat+lat-1
    g2(NbHLat+lat)=NbHLat+lat
    w1(NbHLat+lat)=(ylat_gauss(NbHLat+lat)-(lat-0.5_dp)*dy)/ &
                   (ylat_gauss(NbHLat+lat)-ylat_gauss(NbHLat+lat-1))
    w2(NbHLat+lat)=((lat-0.5_dp)*dy-ylat_gauss(NbHLat+lat-1))/ &
                   (ylat_gauss(NbHLat+lat)-ylat_gauss(NbHLat+lat-1))
    g1(NbHLat+1-lat)=NbHLat+1-lat
    g2(NbHLat+1-lat)=NbHLat+2-lat
    w1(NbHLat+1-lat)=w2(NbHLat+lat)
    w2(NbHLat+1-lat)=w1(NbHLat+lat)
  enddo
  deallocate (y_gauss,w_gauss,ylat_gauss)
    
! Check whether the grid is global and whether the poles are included 
 
  if (xglobal.and.(min(yfirst,ylast)<-89._dp)) then
    ! field contains south pole
    sglobal=.true.
    ! Enhance the map scale by factor 3 (*2=6) compared to north-south
    ! map scale
    sizesouth=6._dp*(switchsouth+90._dp)/dy
    call stlmbr(southpolemap,-90._dp,0._dp)
    call stcm2p(southpolemap,0._dp,0._dp,switchsouth,0._dp,sizesouth, &
    sizesouth,switchsouth,180._dp)
    switchsouthg=(switchsouth-ylat0)/dy
  else
    sglobal=.false.
    switchsouthg=999999._dp
  endif
  
  if (xglobal.and.(max(yfirst,ylast)>89._dp)) then
    ! field contains south pole
    nglobal=.true.
    ! Enhance the map scale by factor 3 (*2=6) compared to north-south
    ! map scale
    sizenorth=6._dp*(90._dp-switchnorth)/dy
    call stlmbr(northpolemap,90._dp,0._dp)
    call stcm2p(northpolemap,0._dp,0._dp,switchnorth,0._dp,sizenorth, &
    sizenorth,switchnorth,180._dp)
    switchnorthg=(switchnorth-ylat0)/dy
  else
    nglobal=.false.
    switchnorthg=999999._dp
  endif
  write(*,'(a,2i5,2L3)')' gribcheck> nx,ny,nglobal,sglobal ', &
           nx, ny, nglobal, sglobal
  write(*,'(a,2f10.2)')' gribcheck> switchsouthg, switchnorthg ', &
           switchsouthg, switchnorthg
  
  return        
  end subroutine gridcheck_jra55
 
  subroutine gauleg(x,w)
! Calculation of the Gaussian abscissa and weights,
! according to Numerical Recipes.
! Calculation are made in double precision.
  real(dbl), intent(inout) :: x(:),w(:)
  integer :: m  
  integer :: n,j,i
  real(dbl) :: eps,pi,p1,p2,p3,pp,z,z1 
  m = size(x)
  n = 2*m
  eps = 1.D-14
  pi = 4.d0 * DACOS(1./dsqrt(2.d0))
  do i = 1, m
     z = dcos(pi*(dfloat(i)-.25D0)/(dfloat(n)+.5D0))
 1   continue
     p1 = 1.D0
     p2 = 0.D0
     do j = 1,n
        p3 = p2
        p2 = p1
        p1 = ((2.d0*dfloat(j)-1.d0)*z*p2 -(dfloat(j)-1.d0)*p3)/dfloat(j)
     enddo
     pp = n*(z*p1-p2)/(z*z-1.d0)
     z1 = z
     z = z1-p1/pp
     if(dabs(z-z1).gt.eps) go to 1
     x(m+1-i) = dasin(z)*180.D0/pi
     w(m+1-i) = 2.d0/((1.d0-z*z)*pp*pp)
  enddo
  return
  end subroutine gauleg
  
  subroutine setabl(akz,bkz)
! gives pressure as p=p0*akz+ps*bkz
! where p0 and ps must be in same units 
  real(dp), intent(inout) :: akz(:),bkz(:)
  real(dbl) :: aa(60),bb(60)
  data aa /0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
0.D0, 0.00066734297302D0, 0.00249386812530D0, 0.00500293194890D0, &
0.00797965319352D0, 0.01154802009928D0, 0.01570654044647D0, &
0.02044605865779D0, 0.02573226625069D0, 0.03160783983371D0, &
0.03819194447371D0, 0.04531113834805D0, 0.05286236151278D0, &
0.06079191356662D0, 0.06892602711506D0, 0.07715965104333D0, &
0.08524595100489D0, 0.09291889331031D0, 0.09999794335439D0, &
0.10628175673245D0, 0.11147896913969D0, 0.11539157373431D0, &
0.11784260347792D0, 0.118638D0, 0.117637D0, 0.114774D0, 0.110092D0, 0.103721D0, &
0.0958905D0, 0.0869323D0, 0.0770862D0, 0.0663751D0, 0.0550953D0, &
0.044665762162664D0, 0.036081493251726D0, 0.029145013769238D0, &
0.023529814041243D0, 0.018988988240733D0, 0.015320365348742D0, &
0.012351296182718D0, 0.009971208628290D0, 0.008049498872017D0, &
0.006492593496917D0, 0.005240019775690D0, 0.004221641249659D0, &
0.003382881250714D0, 0.002683635936493D0, 0.002088961054950D0, &
0.001584429547310D0, 0.0011595266408863D0, 0.0008046871791830D0, &
0.0005145078427194D0, 0.0002897808727181D0, 0.000100D0 /
  data bb/0.9984996244365, 0.9954996233047D0, 0.9914989494018D0, 0.9854979282879D0, &
0.9769957352197D0, 0.9659937887680D0, 0.9529914304700D0, 0.93631825158787D0, &
0.91548798405419D0, 0.89197459659135D0, 0.86549055061078D0, &
0.83541343807984D0, 0.80224765254850D0, 0.76649977884601D0, &
0.72870766385701D0, 0.68882141575136D0, 0.64622485946865D0, &
0.60210094054390D0, 0.55703923517824D0, 0.51110318922792D0, &
0.46496166119883D0, 0.41871948953482D0, 0.37313015634147D0, &
0.32895373629157D0, 0.28637066249204D0, 0.24558214564561D0, &
0.20738803908350D0, 0.17196995546443D0, 0.13952209373832D0, 0.11022D0, &
0.084224D0, 0.0620912D0, 0.0437779D0, 0.0291539D0, 0.0179913D0, 0.00995778D0, &
0.00455656D0, 0.0012626D0, 0.0000547673D0, 0.D0, 0.D0, 0.D0, &
0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/

  if(size(akz)/=60) then
    print *,'Bad number of vertical levels ',size(akz)
    stop
  endif
  
  akz=aa
  bkz=bb
  
  return
  end subroutine setabl

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine getfields_jra55(itime,nstop)

!*******************************************************************************

!        B. Legras, March 2013 (from getfields)
!        Unified version that process both velocity and diabatic files
!        get U, V, T and w if required
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

      integer indj,indmin,indmin_diab,itime,nstop,memaux
      save indmin,indmin_diab

      data indmin/1/,indmin_diab/1/

! 1st part
! Check, if wind fields are available for the current time step
!**************************************************************

      nstop=0

      if ((ldirect*wftime(1).gt.ldirect*itime).or.    &
        (ldirect*wftime(numbwf).lt.ldirect*itime)) then
        write(*,*) 'TRACZILLA WARNING: NO jra55 WINDS ARE AVAILABLE.'
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
              call read_jra55(indj+1,memind(2))
              call verttransform_jra55(memind(2))
              memtime(2)=wftime(indj+1)
              write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_jra55   > date ',trim(wfname(indj+1)),&
                      '  memtime ',memtime(2),'  time ',itime
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
               call read_jra55(indj,memind(1))
               call verttransform_jra55(memind(1))
               memtime(1)=wftime(indj)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_jra55   > date ',trim(wfname(indj)),&
                      '  memtime ',memtime(1),'  time ',itime
               memind(2)=2
               call read_jra55(indj+1,memind(2))
               call verttransform_jra55(memind(2))
               memtime(2)=wftime(indj+1)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_jra55   > date ',trim(wfname(indj+1)),&
                      '  memtime ',memtime(2),'  time ',itime
               nstop = 1
               goto 60
            endif
         enddo
 60      indmin=indj

      endif

      lwindinterv=abs(memtime(2)-memtime(1))
 
      if (lwindinterv.gt.idiffmax) nstop=3
      
! 2nd part
! Check, if heating rate fields are available for the current time step
!**************************************************************

      ifdiab: if (jra55_diab) then

        if ((ldirect*wftime_diab(1).gt.ldirect*itime).or.    &
            (ldirect*wftime_diab(numbwf_diab).lt.ldirect*itime)) then
          write(*,*) 'TRACZILLA WARNING: NO jra55 HEATINGS ARE AVAILABLE.'
          write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
          write(*,*) ldirect*wftime_diab(1)
          write(*,*) ldirect*itime
          write(*,*) ldirect*wftime_diab(numbwf_diab)
          nstop=4
          return
        endif

        ifread: if ((ldirect*memtime_diab(1).le.ldirect*itime).and.   &
                    (ldirect*memtime_diab(2).gt.ldirect*itime)) then

! The right wind fields are already in memory -> don't do anything
!*****************************************************************
          continue

        else if ((ldirect*memtime_diab(2).le.ldirect*itime).and. &
           (memtime_diab(2).ne.999999999)) then ifread
 
! Current time is after 2nd wind field
! -> Resort wind field pointers, so that current time is between 1st and 2nd
!***************************************************************************

          memaux=memind_diab(1)
          memind_diab(1)=memind_diab(2)
          memind_diab(2)=memaux
          memtime_diab(1)=memtime_diab(2)

! Read a new wind field and store it on place memind(2)
!******************************************************

          do indj=indmin_diab,numbwf_diab-1
             if (ldirect*wftime_diab(indj+1).gt.ldirect*itime) then    
                call read_jra55_diab(indj+1,memind_diab(2))
                call verttransform_jra55_diab(memind_diab(2),2)
                memtime_diab(2)=wftime_diab(indj+1)
                write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_jra55_hr> date ',trim(wfname_diab(indj+1)),&
                      '  memtime ',memtime_diab(2),'  time ',itime

                nstop = 1
                goto 45
             endif
          enddo
 45       indmin_diab=indj

        else ifread

! No wind fields, which can be used, are currently in memory 
! -> read both wind fields
!***********************************************************

          do indj=indmin_diab,numbwf_diab-1
             if ((ldirect*wftime_diab(indj).le.ldirect*itime).and.   &
                 (ldirect*wftime_diab(indj+1).gt.ldirect*itime)) then
               memind_diab(1)=1  
               print *,'first read a couple of fields'        
               call read_jra55_diab(indj,memind_diab(1))
               call verttransform_jra55_diab(memind_diab(1),1)
               memtime_diab(1)=wftime_diab(indj)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_jra55_hr> date ',trim(wfname_diab(indj)),&
                      '  memtime ',memtime_diab(1),'  time ',itime
               memind_diab(2)=2
               call read_jra55_diab(indj+1,memind_diab(2))
               call verttransform_jra55_diab(memind_diab(2),2)
               memtime_diab(2)=wftime_diab(indj+1)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_jra55_hr> date ',trim(wfname_diab(indj+1)),&
                      '  memtime ',memtime_diab(2),'  time ',itime
               nstop = 1

               goto 65
             endif
           enddo
 65        indmin_diab=indj

        endif ifread

        lwindinterv=abs(memtime_diab(2)-memtime_diab(1))
 
        if (lwindinterv.gt.idiffmax) nstop=3
 
 !    if z_motion, set memind_diab to memind for correct
 !    indexing of vertical velocities in the interpolation routine     
      else if (z_motion) then ifdiab
         memind_diab=memind         
      
      endif ifdiab

      return
      end subroutine getfields_jra55

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine read_jra55(indj,n)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READ_jra55                  *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      B. Legras                    *
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS INTERPOLATED ON THETA LEVELS *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN SPECTRAL FORMAT    *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!
! This routine uses indexes for the grib files containing t,u,v,w and ps  
! Wind and temperature fields are read into temporary arrays
! in the Gaussian latitude grid called val_gauss
! They are then interpolated to the regular grid by linear interpolation
! Using pre-computed interpolation coefficients calculated in the 
! gridcheck routine
! Reading is parallelized over the variables
! 
! IMPORTANT VARIABLES FROM SHARED VARIABLES:                          *
!                                                                     *
! wfname                date and hour of the data to be read          *
! Nblon, Nblat,nblev    expected field dimensions
! idx_uvwt, idx_ps      current index pointers of uvwt and ps files   *
! current_year          année courante de l'index                     *
! 
!**********************************************************************

  integer, intent(in) :: indj,n
  integer :: ll,nb_var,hh
  integer :: l,lat,i,err_loc 
  character (len=2) :: hour,month
  character (len=4) :: year
  character (len=6) :: yearmonth
  character (len=8) :: fulldate
!  integer :: OMP_GET_THREAD_NUM
  integer :: igrib_thread,iret
  !** real(dp), allocatable :: values(:),val_gauss(:,:)
  real(dp), allocatable :: values(:,:),val_gauss(:,:,:)
  !**** real(dp) :: values(204800,5),val_gauss(640,320,5)
  ! ACHTUNG: 4 should be changed into 5 if z_motion
  
! Finds the date from wfname
  ll=len_trim(wfname(indj))
  hour=wfname(indj)(ll-1:ll)
  year=wfname(indj)(1:4)
  month=wfname(indj)(5:6)
  fulldate=wfname(indj)(1:8)
  yearmonth=wfname(indj)(1:6)
  nb_var=3
  if (z_motion) nb_var=4
 
  select case(hour)
  case('00')
    hh=0
  case('06')
    hh=600
  case('12')
    hh=1200
  case('18')
    hh=1800
  case default
    hh=-1
    print *,'error in hh ',hh
  end select
  
! Loads new year when needed
! the index is read nb_var time as we do not know what contains index structure and
! how to copy it
  if (yearmonth /= current_ym) then
    do i=1,nb_var
      call grib_index_release(idx_uvwt(i))
      call grib_index_read(idx_uvwt(i), &
      path(3)(1:len_path(3))//year//'/'//'uvwt-'//year//'-'//month//'.gribidx');
    enddo
    call grib_index_release(idx_ps)
    call grib_index_read(idx_ps, &
    path(3)(1:len_path(3))//year//'/'//'ps-'//year//'-'//month//'.gribidx');
    print *,'New uvwt ps idx loaded, year month ',year,month
    current_year=year
    current_ym=year//month
  endif
  
! Allocate temp fields out of the parallel loop
  !*** if(z_motion) then
  !***  allocate(values(NbLon*NbLat,5),val_gauss(NbLon,NbLat,5))
  !*** else
  allocate(values(NbLon*NbLat,4),val_gauss(NbLon,NbLat,4))
  !*** endif
  
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) & 
!$OMP PRIVATE(l,lat,igrib_thread,iret)
! Read one variable (or file) in each section
!$OMP SECTION
! Temperature
  !** allocate (values(NbLon*NbLat),val_gauss(NbLon,NbLat))
  call grib_index_select(idx_uvwt(1),'mars.date',fulldate)
  call grib_index_select(idx_uvwt(1),'mars.time',hh)
  call grib_index_select(idx_uvwt(1),'mars.param','11.200')
  do l=1,NbLev
     call grib_index_select(idx_uvwt(1),'mars.levelist',l)
     call grib_new_from_index(idx_uvwt(1),igrib_thread, iret)
     !** call grib_get_real4_array(igrib_thread,'values',values)        
     !** val_gauss=reshape(values,(/NbLon,NbLat/))
     call grib_get_real4_array(igrib_thread,'values',values(:,1))
     call grib_release(igrib_thread)
     val_gauss(:,:,1)=reshape(values(:,1),(/NbLon,NbLat/))
! Interpolation to regular grid
! Non polar latitudes, poles are treated in verttransform 
! This loop also invert the latitude order    
     do lat=1,NbLat
        !** tth(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat)) &
        !**                               +w2(lat)*val_gauss(1:NbLon,g2(lat))
        tth(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat),1) &
                                      +w2(lat)*val_gauss(1:NbLon,g2(lat),1)
     enddo
     if(xglobal) tth(NbLon,0:NbLat-1,l,n) = tth(0,0:NbLat-1,l,n)
  enddo
  !** deallocate(values,val_gauss)
!$OMP SECTION
! Zonal wind
  !** allocate (values(NbLon*NbLat),val_gauss(NbLon,NbLat))
  call grib_index_select(idx_uvwt(2),'mars.date',fulldate)
  call grib_index_select(idx_uvwt(2),'mars.time',hh)
  call grib_index_select(idx_uvwt(2),'mars.param','33.200')
  do l=1,NbLev
     call grib_index_select(idx_uvwt(2),'mars.levelist',l)
     call grib_new_from_index(idx_uvwt(2),igrib_thread, iret)
     !** call grib_get_real4_array(igrib_thread,'values',values)        
     !** val_gauss=reshape(values,(/NbLon,NbLat/))
     call grib_get_real4_array(igrib_thread,'values',values(:,2))
     call grib_release(igrib_thread)
     val_gauss(:,:,2)=reshape(values(:,2),(/NbLon,NbLat/))  
     do lat=1,NbLat
        !** uuh(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat),2) &
        !**                               +w2(lat)*val_gauss(1:NbLon,g2(lat),2)
        uuh(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat),2) &
                                      +w2(lat)*val_gauss(1:NbLon,g2(lat),2)                              
     enddo
     if(xglobal) uuh(NbLon,0:NbLat-1,l,n) = uuh(0,0:NbLat-1,l,n)
  enddo
  !** deallocate(values,val_gauss)
!$OMP SECTION
! Meridional wind
  !** allocate (values(NbLon*NbLat),val_gauss(NbLon,NbLat))
  call grib_index_select(idx_uvwt(3),'mars.date',fulldate)
  call grib_index_select(idx_uvwt(3),'mars.time',hh)
  call grib_index_select(idx_uvwt(3),'mars.param','34.200')
  do l=1,NbLev
     call grib_index_select(idx_uvwt(3),'mars.levelist',l)
     call grib_new_from_index(idx_uvwt(3),igrib_thread, iret)
     !** call grib_get_real4_array(igrib_thread,'values',values)        
     !** val_gauss=reshape(values,(/NbLon,NbLat/))
     call grib_get_real4_array(igrib_thread,'values',values(:,3))
     call grib_release(igrib_thread)
     val_gauss(:,:,3)=reshape(values(:,3),(/NbLon,NbLat/))   
     do lat=1,NbLat
        !** vvh(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat)) &
        !**                               +w2(lat)*val_gauss(1:NbLon,g2(lat))
        vvh(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat),3) &
                                      +w2(lat)*val_gauss(1:NbLon,g2(lat),3)                              
     enddo
     if(xglobal) uuh(NbLon,0:NbLat-1,l,n) = vvh(0,0:NbLat-1,l,n)
  enddo
  !** deallocate(values,val_gauss)
!$OMP SECTION
! Surface pressure
  !** allocate (values(Nblon*Nblat),val_gauss(Nblon,Nblat))
  call grib_index_select(idx_ps,'mars.date',fulldate)
  call grib_index_select(idx_ps,'mars.time',hh)
  call grib_new_from_index(idx_ps,igrib_thread, iret)
  !** call grib_get_real4_array(igrib_thread,'values',values)        
  !** val_gauss=reshape(values,(/Nblon,Nblat/))
  call grib_get_real4_array(igrib_thread,'values',values(:,4))
  call grib_release(igrib_thread)       
  val_gauss(:,:,4)=reshape(values(:,4),(/Nblon,Nblat/))
  do lat=1,Nblat
      !** ps(0:NbLon-1,NbLat-lat,1,n) = w1(lat)*val_gauss(1:NbLon,g1(lat)) &
      !**                              +w2(lat)*val_gauss(1:NbLon,g2(lat))
      ps(0:NbLon-1,NbLat-lat,1,n) = w1(lat)*val_gauss(1:NbLon,g1(lat),4) &
                                   +w2(lat)*val_gauss(1:NbLon,g2(lat),4)                        
  enddo
  if(xglobal) ps(NbLon,0:NbLat-1,1,n) = ps(0,0:NbLat-1,1,n)
  !** deallocate(values,val_gauss)	
!$OMP SECTION
! Vertical Dp/Dt
  if (z_motion) then
    !** allocate (values(Nblon*Nblat),val_gauss(Nblon,Nblat))
    call grib_index_select(idx_uvwt(4),'mars.date',fulldate)
    call grib_index_select(idx_uvwt(4),'mars.time',hh)
    call grib_index_select(idx_uvwt(4),'mars.param','39.200')
    do l=1,NbLev
       call grib_index_select(idx_uvwt(4),'mars.levelist',l)
       call grib_new_from_index(idx_uvwt(4),igrib_thread, iret)
       !** call grib_get_real4_array(igrib_thread,'values',values)        
       !** val_gauss=reshape(values,(/Nblon,Nblat/))
       call grib_get_real4_array(igrib_thread,'values',values(:,5))
       call grib_release(igrib_thread)
       val_gauss(:,:,5)=reshape(values(:,5),(/NbLon,NbLat/))   
       do lat=1,Nblat
          !** wwh(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat)) &
          !**                               +w2(lat)*val_gauss(1:NbLon,g2(lat))
          wwh(0:NbLon-1,NbLat-lat,l,n) = w1(lat)*val_gauss(1:NbLon,g1(lat),5) &
                                        +w2(lat)*val_gauss(1:NbLon,g2(lat),5)
       enddo
       if(xglobal) wwh(NbLon,0:NbLat-1,l,n) = wwh(0,0:NbLat-1,l,n)
    enddo
    deallocate(values,val_gauss) 
  endif	
				
!$OMP END PARALLEL SECTIONS

  deallocate(values,val_gauss,STAT=err_loc)
  
  return

  end subroutine read_jra55
  
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine read_jra55_diab(indj,n)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READ_jra55                  *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      B. Legras                    *
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS INTERPOLATED ON THETA LEVELS *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN SPECTRAL FORMAT    *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!
! This routine uses indexes for the grib files containing heating rates  
! Heating rate fields are read into temporary arrays
! in the Gaussian latitude grid called val_gauss
! They are then interpolated to the regular grid by linear interpolation
! Using pre-computed interpolation coefficients calculated in the 
! gridcheck routine
! Reading is parallelized over the variables
! 
! IMPORTANT VARIABLES FROM SHARED VARIABLES:                          *
!                                                                     *
! wfname                date and hour of the data to be read          *
! Nblon, Nblat,nblev    expected field dimensions
! idxhr                 current index pointers of heating rate files  *
! current_year_hr       année courante de l'index                     *
! 
!**********************************************************************

  integer, intent(in) :: indj,n
  integer :: ll,nb_var,hh
  integer :: l,lat,i,err_loc  
  character (len=2) :: hour,month
  character (len=4) :: year
  character (len=6) :: yearmonth
  character (len=8) :: fulldate
!  integer :: OMP_GET_THREAD_NUM
  real(dp), allocatable :: values(:,:),swr_gauss(:,:,:),lwr_gauss(:,:,:)
  !*** real(dp) :: values(204800,2),swr_gauss(640,320,60),lwr_gauss(640,320,60)
  integer :: igrib_thread, iret
! Finds the date and time from wfname
  !flush(6)
  ll=len_trim(wfname_diab(indj))
  hour=wfname_diab(indj)(ll-1:ll)
  year=wfname_diab(indj)(1:4)
  month=wfname_diab(indj)(5:6)
  yearmonth=wfname_diab(indj)(1:6)
  fulldate=wfname_diab(indj)(1:8)
  nb_var=2
 
  select case(hour)
  case('00')
    hh=0
  case('06')
    hh=600
  case('12')
    hh=1200
  case('18')
    hh=1800
  case default
    hh=-1
    print *,'error in hh ',hh
  end select
  
! Loads new year when needed
! reading is performed twice as we do not know how to copy index structure
  if (yearmonth /= current_ym_hr) then
    do i=1,nb_var
      call grib_index_release(idx_hr(i))
      call grib_index_read(idx_hr(i), &
      path(3)(1:len_path(3))//year//'/'//'hr-'//year//'-'//month//'.gribidx')
    enddo
    print *,'New hr idx loaded, year month',year, month
    current_ym_hr=year//month
  endif
  
  allocate (lwr_gauss(NbLon,NbLat,NbLev),swr_gauss(NbLon,NbLat,NbLev))
  allocate (values(NbLon*NbLat,2))
  
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(l,igrib_thread,iret)
!$OMP SECTION
! Short wave
  !** allocate (values(Nblon*Nblat))
  call grib_index_select(idx_hr(1),'mars.date',fulldate)
  call grib_index_select(idx_hr(1),'mars.time',hh)
  call grib_index_select(idx_hr(1),'mars.param','250.200')
  do l=1,NbLev
     call grib_index_select(idx_hr(1),'mars.levelist',l)
     call grib_new_from_index(idx_hr(1),igrib_thread, iret)
     !** call grib_get_real4_array(igrib_thread,'values',values)        
     !** swr_gauss(:,:,l)=reshape(values,(/Nblon,Nblat/))
     call grib_get_real4_array(igrib_thread,'values',values(:,1))
     call grib_release(igrib_thread)       
     swr_gauss(:,:,l)=reshape(values(:,1),(/Nblon,Nblat/))
  enddo
  !** deallocate(values)
!$OMP SECTION
! Long wave
  !** allocate (values(Nblon*Nblat))
  call grib_index_select(idx_hr(2),'mars.date',fulldate)
  call grib_index_select(idx_hr(2),'mars.time',hh)
  call grib_index_select(idx_hr(2),'mars.param','251.200')
  do l=1,NbLev
     call grib_index_select(idx_hr(2),'mars.levelist',l)
     call grib_new_from_index(idx_hr(2),igrib_thread, iret)
     !** call grib_get_real4_array(igrib_thread,'values',values)        
     !** lwr_gauss(:,:,l)=reshape(values,(/Nblon,Nblat/))
     call grib_get_real4_array(igrib_thread,'values',values(:,2))
     call grib_release(igrib_thread)      
     lwr_gauss(:,:,l)=reshape(values(:,2),(/Nblon,Nblat/))
  enddo
  !** deallocate(values)
!$OMP END PARALLEL SECTIONS

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(lat)
  do lat=1,NbLat
     wwh(0:NbLon-1,NbLat-lat,1:NbLev,n) = &
                   w1(lat)*(swr_gauss(1:NbLon,g1(lat),1:NbLev) &
                           +lwr_gauss(1:NbLon,g1(lat),1:NbLev)) &
                  +w2(lat)*(swr_gauss(1:NbLon,g2(lat),1:NbLev) &
                           +lwr_gauss(1:NbLon,g2(lat),1:NbLev))
  enddo
!$OMP END PARALLEL DO

  if(xglobal) wwh(NbLon,0:NbLat-1,1:NbLev,n) = wwh(0,0:NbLat-1,1:NbLev,n)

  deallocate(values,STAT=err_loc)
  if(err_loc/=0) then
    print *,'error deallocating values ',err_loc
  endif
  deallocate(lwr_gauss,swr_gauss,STAT=err_loc)
  if(err_loc/=0) then
    print *,'error deallocating lwr and swr_gauss ',err_loc
  endif
 
  return

  end subroutine read_jra55_diab
  
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine verttransform_jra55(n)

  use coord
  integer, intent(in):: n
  integer :: ix, jy, iz
  real(dp) :: xlon, ylat

! north pole region
  if (nglobal) then
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(ix,jy,iz,xlon,ylat)
        do ix=0,nx-1
          xlon=xlon0+float(ix)*dx
          do jy=int(switchnorthg)-2,ny-1
            ylat=ylat0+float(jy)*dy
            do iz=1,nuvz
              call cc2gll(northpolemap,ylat,xlon,  &
                uuh(ix,jy,iz,n), vvh(ix,jy,iz,n),  &
                uupol(ix,jy,iz,n), vvpol(ix,jy,iz,n))
            enddo
          enddo
        enddo
!$OMP END DO
!      calculation at the north pole from an average over the closest circle
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(iz)
        do iz=1,nuvz
          uupol(0:nx-1,ny,iz,n)=sum(uupol(0:nx-2,ny-1,iz,n))/(nx-1)
          vvpol(0:nx-1,ny,iz,n)=sum(vvpol(0:nx-2,ny-1,iz,n))/(nx-1)
          tth(0:nx-1,ny,iz,n)=sum(tth(0:nx-2,ny-1,iz,n))/(nx-1)
          if(z_motion) wwh(0:nx-1,ny,iz,n)=sum(wwh(0:nx-2,ny-1,iz,n))/(nx-1)
        enddo
!$OMP ENDDO
!$OMP END PARALLEL
        ps(0:nx-1,ny,1,n)=sum(ps(0:nx-2,ny-1,1,n))/(nx-1)
  endif

!    south pole region
  if (sglobal) then
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(ix,jy,iz,xlon,ylat)
        do ix=0,nx-1
          xlon=xlon0+float(ix)*dx
          do jy=0,int(switchsouthg)+3
            ylat=ylat0+float(jy)*dy           
            do iz=1,nuvz
              call cc2gll(southpolemap,ylat,xlon, &
                uuh(ix,jy,iz,n), vvh(ix,jy,iz,n), &
                uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
            enddo
          enddo
        enddo
!$OMP ENDDO
!       calculation at the south pole from the average of the closest circle
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(iz)
        do iz=1,nuvz
          uupol(0:nx-1,-1,iz,n)=sum(uupol(0:nx-2,0,iz,n))/(nx-1)
          vvpol(0:nx-1,-1,iz,n)=sum(vvpol(0:nx-2,0,iz,n))/(nx-1)
          tth(0:nx-1,-1,iz,n)=sum(tth(0:nx-2,0,iz,n))/(nx-1)
          if (z_motion) wwh(0:nx-1,-1,iz,n)=sum(wwh(0:nx-2,0,iz,n))/(nx-1)
        enddo
!$OMP ENDDO
!$OMP END PARALLEL
        ps(0:nx-1,-1,1,n)=sum(ps(0:nx-2,0,1,n))/(nx-1)
  endif
  
  return
  end subroutine verttransform_jra55
 
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine verttransform_jra55_diab(n,m)
  
  integer, intent(in):: n,m
  integer :: iz,ix,jy
  real(dp) :: pp 
 
!     Conversion in degree of theta per second
!     Initially in temperature and K/day
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC) PRIVATE(iz,jy,ix,pp)      
  do iz=1,nuvz
     do jy=0,ny-1
        do ix=0,nx-1
           pp=akz(iz)+bkz(iz)*(ps(ix,jy,1,1)+(ps(ix,jy,1,2)))/2
           wwh(ix,jy,iz,n)=wwh(ix,jy,iz,n)*((p0/pp)**kappa)/86400.
        enddo
     enddo
     if(nglobal) wwh(0:nx-1,ny,iz,n)=sum(wwh(0:nx-2,ny-1,iz,n))/(nx-1)
     if(sglobal) wwh(0:nx-1,-1,iz,n)=sum(wwh(0:nx-2,0,iz,n))/(nx-1)
  enddo
!$OMP END PARALLEL DO
  !print *,'wwh just read, max, min ',maxval(maxval(wwh(:,:,35,n),DIM=2),DIM=1),&
  !                                   minval(minval(wwh(:,:,35,n),DIM=2),DIM=1)     
        
! call mass correction when required
  if(mass_correction) call diab_mass_jra55(n,m)
 
  return
  end subroutine verttransform_jra55_diab
  
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass_jra55(n,m)

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
  real(dp), allocatable :: theta(:,:,:), sigma(:,:,:), flux(:,:,:), flux_lat(:,:)
  real(dp), allocatable :: mass_flux(:), mean_sigma(:), mean_sigma_lat(:,:)
  real(dp), allocatable :: mean_w(:)
  integer :: k

  allocate(theta(0:nx-1,0:ny-1,nuvz),sigma(0:nx-1,0:ny-1,nuvz))
  allocate(mean_sigma(nuvz),mass_flux(nuvz),mean_w(nuvz))
  allocate(flux(0:nx-1,0:ny-1,nuvz))
  allocate(flux_lat(0:ny-1,nuvz))
  allocate(mean_sigma_lat(0:ny-1,nuvz))
  
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(k)
! Calculation of the potential temperature
  do k=nuvz-NPureP_jra55,nuvz
     theta(:,:,k)=0.5*(tth(:,:,k,1)+tth(:,:,k,2))*pif_jra55(k)
!    theta(:,:,k)=0.5*(theta_g(:,:,k,1)+theta_g(:,:,k,2))
  enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(k)
  do k=nuvz-NPureP_jra55+1,nuvz-1
     sigma(:,:,k)=pmc_jra55(k)/(theta(:,:,k+1)-theta(:,:,k-1))
  enddo
!$OMP END DO
  sigma(:,:,nuvz)=pmc_jra55(nuvz)/(theta(:,:,nuvz)-theta(:,:,nuvz-1))
  
! Calculation of the mass flux across the surface
! Note that mean_sigma is not a mean of sigma but the spherical
! integral of sigma (that divides the flux to get the correction)
  !allocate(flux(0:nx-1,0:ny-1))
  !allocate(flux_lat(0:ny-1))
  !allocate(mean_sigma_lat(0:ny-1))
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(k)
  do k=nuvz-NPureP_jra55+1,nuvz
     flux(:,:,k)=wwh(:,:,k,n)*sigma(:,:,k)
     flux_lat(:,k)=sum(flux(0:nx-2,:,k),DIM=1)
     mean_sigma_lat(:,k)=sum(sigma(0:nx-2,:,k),DIM=1)
     mass_flux(k)=0.5*dot_product(flux_lat(:,k),area_coefft_jra55)
     mean_sigma(k)=0.5*dot_product(mean_sigma_lat(:,k),area_coefft_jra55)
     mean_w(k)=mass_flux(k)/mean_sigma(k)
  enddo
!$OMP END DO
  !deallocate(flux,flux_lat,mean_sigma_lat)
!$OMP END PARALLEL
  deallocate(flux,flux_lat,mean_sigma_lat) 

! Output of the averaged diab heating
  if(mean_diab_output) then
    write(unitflux) memtime_diab(m),mean_w(nuvz-NPureP_jra55+1:nuvz), &
                                    mean_sigma(nuvz-NPureP_jra55+1:nuvz)
    flush(unitflux)
  endif

! Correction step
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(k)
  do k=nuvz-NPureP_jra55+1,nuvz
    wwh(:,:,k,n)=wwh(:,:,k,n)-mean_w(k)
  enddo
!$OMP END PARALLEL DO

  deallocate(theta,sigma,mean_sigma)
  deallocate(mass_flux,mean_w)

  return
  end subroutine diab_mass_jra55

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass_jra55_init

!*************************************************************************************
! This routine initializes coefficients used in the calculation of the masss correction 
! The coefficients are p delta log(p) for the calculation of the vertical derivative
! and sin(phi + dphi) - sin(phi - dphi) = 2 sin(dphi) cos(phi) for the surface integral
! where phi si a grid latitude and dphi=0.5*dy. This quantity is the diffference 
! (up to factor 2 pi which is not accounted) of area
! between two caps bounded by phi-dphi and phi+dphi latitudes. The first and last values 
! are the area of polar caps of angle dphi/2, that is 2 sin^2(dphi/4)
! Conversion to radian is applied
!*************************************************************************************
   
  integer :: i
  
  NPureP_jra55=21
 
  if(.not.(xglobal.and.nglobal.and.sglobal)) then
      write(*,*) ' #### TRACZILLA MODEL ERROR! DIAB_MASS        #### ' 
      write(*,*) ' #### mass flux cannot be balanced            #### '
      write(*,*) ' #### on not global grid                      #### '
      stop  
  endif

! calculate area coefficient assuming pole to pole regular grid
  allocate(area_coefft_jra55(0:ny-1))
! it is assumed that dy = 180/(ny-1)
  do i=1,ny/2-1
     area_coefft_jra55(ny/2-1+i)=2*sin(pi*dy/180._dp)*cos((i-1/2)*dy*pi/180._dp)
     area_coefft_jra55(ny/2-i)=area_coefft_jra55(ny/2-1+i)
  enddo
  area_coefft_jra55(0)=2*cos(ylat0*pi/180._dp)*sin((dy/2+90._dp+ylat0)*pi/180._dp)
  area_coefft_jra55(ny-1)=area_coefft_jra55(0)
! test : the sum should be equal to 2g10.3 )') sum(area_coefft_jra55)

! calculate weighting pressure factors in the vertical
  allocate(pmc_jra55(nuvz),pif_jra55(nuvz))
  pmc_jra55=0.
  print *,' diab_mass_init > check sum area  '
  pif_jra55=0.
  do i=nuvz-NPureP_jra55+1,nuvz-1
     pmc_jra55(i)=akz(i)*(Log(akz(i-1))-Log(akz(i+1)))
  enddo
  pmc_jra55(nuvz)=akz(nuvz)*(Log(akz(nuvz-1))-Log(akz(nuvz)))

  do i=nuvz-NPureP_jra55+1,nuvz
     pif_jra55(i)=(p0/akz(i))**kappa
  enddo

! Create the output file of the mas averaged heating flux
! or open it in the append mode
  open(unitflux,file=trim(path(2))//'MassMeanDiab2.dat',form='unformatted',position='append')

  return
  end subroutine diab_mass_jra55_init

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine interpol_wind_jra55 &
         (itime,xt,yt,zt, dxdt,dydt,dzdt, ngrid, & 
          theta_inf,theta_sup,psaver,z_factor,tint,nstop)

!*******************************************************************************
!                                                                              *
!  This subroutine interpolates the wind data to current trajectory position.  *
!                                                                              *
!    Author: B. Legras (derived from original interpol_wind)                   *
!                                                                              *
!    The interpolation is performed in eta coordinates which are linked to     *
!    model data
!    This routine is intended to be used as a universal diabatic interpolator  *
!    to be used for ERA-Interim and JRA-55 diabatic cases (no INC)
!    and for JRA-55 z cases
!    It cannot be used for ERA-Interim z case because the w grid is shifted 
!    in the vertical, but in this case the interpolation is anyway done
!    in log pressure
!    It cannot be used for MERRA because the vertical grid is in pressure 
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! dxdt,dydt          horizontal wind components in grid units per second       *
! itime [s]          current temporal position                                 *
! memtime(2) [s]     times of the wind fields in memory
! memind(2)          indexes of the wind fields in memory                      *
! memtime_diab(2)[s] times of the diabatic fields in memory
! memind(2)          indexes of the diabatic fields in memory      
! xt,yt,zt           coordinates position for which wind data shall be         *
!                    interpolated. zt can be theta or -log(p/p0)
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!*******************************************************************************

      use ecmwf_diab, only : ecmwf_diabatic_w
      integer, intent(in) :: itime,ngrid
      integer, intent(inout):: nstop
      real(dp), intent(in) :: xt,yt,zt
      real(dp), intent(inout) :: tint
      real(dp), intent(out) :: dxdt,dydt,dzdt,z_factor
      real(dp), intent(out) :: psaver
      real(dp), intent(out) :: theta_inf,theta_sup

! Auxiliary variables needed for interpolation
      real(dp) :: press, siglog, etalog, siglogr, siglogrp 
      real(dp) :: u1(2),v1(2),w1(2),dt1,dt2,dtt,tp1(2),dt1_w,dt2_w,dtt_w
      real(dp) :: fup,fbot,u(4,2),v(4,2),w(4,2),tp(4,2)
      integer :: m,indexh,indexh_diab,indz
      integer :: ix,jy
      real(dp) :: ddx,ddy,rddx,rddy,p1,p2,p3,p4

!********************************************
! Multilinear interpolation in time and space
!********************************************

! Determine the lower left corner and its distance to the current position
!*************************************************************************

      ! This section depends on the reanalysis because of the variations
      ! in the latitude grid
      ! min and max required for the points just falling on the boundary
      ! as it may happen typically as a result of initialization
      
      if(xglobal) ix=modulo(floor(xt),nx-1)       
      jy=max(min(floor(yt),ny-1),-1)
      ddx=modulo(xt-float(ix),1.)
      ! accounts for the two polar regions in JRA55
      ! the case should not occur for EI since 0 <= yt <= ny-1
      if (yt<0) then
         ddy=yt/(90._dp/dy-ny/2)+1._dp
      else if (yt>ny-1) then
         ddy=(yt+1._dp-ny)/(90._dp/dy-ny/2)
      else
         ddy=yt-float(jy)
      endif
      rddx=1.-ddx      ;  rddy=1.-ddy
      p1=rddx*rddy     ;  p2=ddx*rddy
      p3=rddx*ddy      ;  p4=ddx*ddy

! Calculate coefficients for temporal interpolation
!**************************************************

      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)
      dt1=dt1*dtt
      dt2=dt2*dtt
      if (jra55_diab.or.ecmwf_diabatic_w) then
        dt1_w=float(itime-memtime_diab(1))
        dt2_w=float(memtime_diab(2)-itime) 
        dtt_w=1/(dt1_w+dt2_w)
        dt1_w=dt1_w*dtt_w
        dt2_w=dt2_w*dtt_w
      else
        dt1_w=dt1
        dt2_w=dt2
!       memind_diab should be set to memind in getfields
      endif
 
! Determine the pressure and localize the parcel in the vertical grid
!********************************************************************

!  z pressure coordinates
      if (z_motion) then       
        press=p0*exp(-zt)
!  theta coordinates       
      else if (jra55_diab.or.ecmwf_diabatic_w) then
        !if (tint <=0) then
        !  print *,'tint < 0 ',tint,ix,jy,zt
        !endif
        press=p0*(tint/zt)**(1/kappa)
      else
        stop 999
      endif
      !if (press<=0) then 
      !   print *,'press <=0 ',press,tint
      !   print *,xt,yt,zt
      !endif
      

! surface pressure under the parcel at the required time (in Pascal)
      psaver=dt2*(p1*ps(ix,jy,1,memind(1))+p2*ps(ix+1,jy,1,memind(1))    &
                   +p3*ps(ix,jy+1,1,memind(1))+p4*ps(ix+1,jy+1,1,memind(1))) &
            +dt1*(p1*ps(ix,jy,1,memind(2))+p2*ps(ix+1,jy,1,memind(2))    &
                   +p3*ps(ix,jy+1,1,memind(2))+p4*ps(ix+1,jy+1,1,memind(2)))
                   
      if(theta_bounds) then
        theta_inf=280._dp  ! test value always less than minimum ztra1 
                           ! to avoid unecessary calculation
                 
        theta_sup = ((tth(ix, jy, upper_theta_level,memind(1))*p1 &
                    + tth(ix, jy+1,upper_theta_level,memind(1))*p2 &
                    + tth(ix+1,jy, upper_theta_level,memind(1))*p3 &
                    + tth(ix+1,jy+1,upper_theta_level,memind(1))*p4)*dt2 &
                   + (tth(ix, jy, upper_theta_level,memind(2))*p1 &
                    + tth(ix, jy+1,upper_theta_level,memind(2))*p2 &
                    + tth(ix+1,jy, upper_theta_level,memind(2))*p3 &
                    + tth(ix+1,jy+1,upper_theta_level,memind(2))*p4)*dt1) &
                   * (p0/akz(upper_theta_level))**kappa
      endif

! Locates the pressure in the grid
! Equal indz since it is a pure eta grid 

      indz = locpj(1,nuvz,psaver,press)
      siglog=log(press/psaver)

! Halt trajectories which are too close to lower boundary     
      if ((press > psaver) .or. (indz==1)) then
        nstop=2
        dxdt=0. ; dydt=0.; dzdt=0.; tint=275.;
         z_factor=1.
        return
      endif
! Mark trajectories which are too close to upper boundary
      if (indz>=upper_theta_level-1) nstop=-4   
           
! Calculation of log_eta according to the inversion formula in log-log
  
      siglogr=log(bkz(indz)+akz(indz)/psaver)
      siglogrp=log(bkz(indz+1)+akz(indz+1)/psaver)
      etalog=0.5_dp*(etakzlog(indz)+etakzlog(indz+1))& 
            +0.5_dp*(2*siglog-siglogr-siglogrp)/(siglogr-siglogrp) &
                   *(etakzlog(indz)-etakzlog(indz+1))
                   
! Calculate factors of bilinear interpolation as they are common 
! to all elements

      fbot=(etalog-etakzlog(indz+1))/(etakzlog(indz)-etakzlog(indz+1))
      fup=(etakzlog(indz)-etalog)/(etakzlog(indz)-etakzlog(indz+1)) 
                     
!**********************************************************************
! 1.) Bilinear horizontal interpolation
! This has to be done separately for 4 fields (Temporal(2)*Vertical(2))
!**********************************************************************

! Loop over 2 time steps and 2 levels
!************************************

      if (ngrid < 0) then  ! polar region
      
        do m=1,2
          indexh=memind(m)
          indexh_diab=memind_diab(m)
            
          u(1,m)=uupol(ix ,jy ,indz  ,indexh)*fbot  &
                + uupol(ix ,jy ,indz+1,indexh)*fup  
          v(1,m)=vvpol(ix ,jy ,indz  ,indexh)*fbot  &
                + vvpol(ix ,jy ,indz+1,indexh)*fup  
          w(1,m)=wwh(ix ,jy ,indz  ,indexh_diab)*fbot  &
                + wwh(ix ,jy ,indz+1,indexh_diab)*fup 
          u(2,m)=uupol(ix+1,jy ,indz  ,indexh)*fbot  &
                + uupol(ix+1,jy ,indz+1,indexh)*fup 
          v(2,m)=vvpol(ix+1,jy ,indz  ,indexh)*fbot  &
                + vvpol(ix+1,jy ,indz+1,indexh)*fup
          w(2,m)=wwh(ix+1,jy ,indz  ,indexh_diab)*fbot &
                + wwh(ix+1,jy ,indz+1,indexh_diab)*fup
          u(3,m)=uupol(ix ,jy+1,indz  ,indexh)*fbot &
                + uupol(ix ,jy+1,indz+1,indexh)*fup
          v(3,m)=vvpol(ix ,jy+1,indz  ,indexh)*fbot &
                + vvpol(ix ,jy+1,indz+1,indexh)*fup
          w(3,m)=wwh(ix ,jy+1,indz  ,indexh_diab)*fbot  &
                + wwh(ix ,jy+1,indz+1,indexh_diab)*fup
          u(4,m)=uupol(ix+1,jy+1,indz  ,indexh)*fbot &
                + uupol(ix+1,jy+1,indz+1,indexh)*fup
          v(4,m)=vvpol(ix+1,jy+1,indz  ,indexh)*fbot &
                + vvpol(ix+1,jy+1,indz+1,indexh)*fup
          w(4,m)=wwh(ix+1,jy+1,indz  ,indexh_diab)*fbot  &
                + wwh(ix+1,jy+1,indz+1,indexh_diab)*fup
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)
          
        enddo 
                
      else     ! non polar region

        do m=1,2
          indexh=memind(m)
          indexh_diab=memind_diab(m)

          u(1,m)=uuh(ix ,jy ,indz  ,indexh)*fbot &
                + uuh(ix ,jy ,indz+1,indexh)*fup
          v(1,m)=vvh(ix ,jy ,indz  ,indexh)*fbot &
                + vvh(ix ,jy ,indz+1,indexh)*fup
          w(1,m)=wwh(ix ,jy ,indz  ,indexh_diab)*fbot &
                + wwh(ix ,jy ,indz+1,indexh_diab)*fup
          u(2,m)=uuh(ix+1,jy ,indz  ,indexh)*fbot &
                + uuh(ix+1,jy ,indz+1,indexh)*fup
          v(2,m)=vvh(ix+1,jy ,indz  ,indexh)*fbot &
                + vvh(ix+1,jy ,indz+1,indexh)*fup
          w(2,m)=wwh(ix+1,jy ,indz  ,indexh_diab)*fbot&
                + wwh(ix+1,jy ,indz+1,indexh_diab)*fup
          u(3,m)=uuh(ix ,jy+1,indz  ,indexh)*fbot &
                + uuh(ix ,jy+1,indz+1,indexh)*fup
          v(3,m)=vvh(ix ,jy+1,indz  ,indexh)*fbot &
                + vvh(ix ,jy+1,indz+1,indexh)*fup
          w(3,m)=wwh(ix ,jy+1,indz  ,indexh_diab)*fbot &
                + wwh(ix ,jy+1,indz+1,indexh_diab)*fup
          u(4,m)=uuh(ix+1,jy+1,indz  ,indexh)*fbot &
                + uuh(ix+1,jy+1,indz+1,indexh)*fup
          v(4,m)=vvh(ix+1,jy+1,indz  ,indexh)*fbot &
                + vvh(ix+1,jy+1,indz+1,indexh)*fup
          w(4,m)=wwh(ix+1,jy+1,indz  ,indexh_diab)*fbot &
                + wwh(ix+1,jy+1,indz+1,indexh_diab)*fup

          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          w1(m)=p1*w(1,m)+p2*w(2,m)+p3*w(3,m)+p4*w(4,m)    
        
        enddo
          
      endif

 ! Accurate calculation of the temperature
      
      do m=1,2
          indexh=memind(m)          
          tp(1,m)=tth(ix ,jy ,indz  ,indexh)*fbot  &
                + tth(ix ,jy ,indz+1,indexh)*fup
          tp(2,m)=tth(ix+1,jy ,indz  ,indexh)*fbot &
                + tth(ix+1,jy ,indz+1,indexh)*fup
          tp(3,m)=tth(ix ,jy+1,indz  ,indexh)*fbot  &
                + tth(ix ,jy+1,indz+1,indexh)*fup
          tp(4,m)=tth(ix+1,jy+1,indz  ,indexh)*fbot  &
                + tth(ix+1,jy+1,indz+1,indexh)*fup
          tp1(m)=p1*tp(1,m)+p2*tp(2,m)+p3*tp(3,m)+p4*tp(4,m)
          
      enddo
      
!************************************
! 3.) Temporal interpolation (linear)
!************************************

      dxdt=u1(1)*dt2+u1(2)*dt1
      dydt=v1(1)*dt2+v1(2)*dt1
      dzdt=w1(1)*dt2_w+w1(2)*dt1_w
      tint=tp1(1)*dt2+tp1(2)*dt1
      
      !if (tint <=0. ) then
      !   print *,'tint <=0', tint
      !   print *,tp1
      !   print *,ix+1,jy+1,indz
      !endif


!***************************************************
! 4.) Calculation of z_factor for vertical diffusion
!***************************************************
!      deactivated sequence

!      select case (diftype)
!      case (1)          ! diffusion in z (cf p.46, book C, part2)
!       d theta / dz = - (g/theta) d theta / d Pi = -g d Log(theta) / d Pi
!       where Pi = Cp (p/p0)**kappa = Cp T/theta
!       estimated from the data on the lower left corner at first time
!       of the interval, for a better estimate using the closest point
!       activate the first following line and deactivate the second one 
!       call sort_hor_distance       
!        i0=ix; j0=jy; idxy=1
!        pisup0 = cpa * tth(i0,j0,indz(idxy,1)+1)/trp(idxy,1)
!        piinf0 = cpa * tth(i0,j0,indz(idxy,1))/tr(idxy,1)
!        z_factor = -ga*(trp(idxy,1)-tr(idxy,1))/(pisup0-piinf0)
!        TO BE REACTIVATED
!        z_factor=0._dp
!      case (2)          ! diffusion in theta
!        z_factor = 1._dp
!      case default
!        z_factor = 0._dp
!      end select
     
      return

      end subroutine interpol_wind_jra55

!*******************************************************************************
! Identical to locuv in readinterp.f90
      function locpj(ib1,iu1,psu,pint)
      integer, intent(in) :: ib1,iu1
      real(dp), intent(in) :: psu,pint
      integer :: locpj,ib,iu,im
      real(dp) :: pm
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
      locpj = ib
      end function locpj
!
!=====|==1=========2=========3=========4=========5=========6=========7==

end module jra55
