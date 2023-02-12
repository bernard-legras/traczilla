!**********************************************************************
! Copyright 1996, 1997, 2001, 2002, 2006, 2007, 2012, 2013, 2015, 2020*
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

module era5

!*******************************************************************************
! This modules allow to use ERA5 data with daily files
! Does not need isentrop.h except for upper_theta_layer as interpolation is made in eta
! In this version, the read are not done in parallel as all the fields
! are in the same files. Parallel reads must be done in separate files
! like in jra55.f90
!*******************************************************************************

use commons
use isentrop_h ! needs upper_theta_layer
implicit none
private :: locpj
private :: NbLon, NbLat, NbHLat, NbLev
private :: idx_uvwt, idx_hr
private :: current_year,current_year_hr,current_ym,current_ym_hr

logical, save :: era5_data, era5_diab

integer, save :: NbLon, NbLat, NbHLat, NbLev

integer, save :: NPureP_era5
real(dp), save, allocatable :: area_coefft_era5(:), pmc_era5(:), pif_era5(:)

character(len=4), save:: current_year,current_year_hr
character(len=6), save:: current_ym,current_ym_hr
integer, save :: idx_uvwt, idx_ps, idx_hr

! LogPressLev           log(p0/PressLev)
! pmc_era5              delta p factor in the calculation of density sigma
! facT                  factor (p0/PressLev)**kappa
! NPMass                number of upper levels on which mass equilibration is performed
! area_coefft_era5      area coefft for the surface integral of the mass flux
! w_era5 [X/s]          vertical velocity either in zlog coordinate (X=zlog)
!                       or potential temperature tendency (X=K) 

contains

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine alloc_era5

print *,'alloc_era5 uuh vvh uupol vvpol tth'
print *,nx,ny,nuvz
allocate (wwh(0:nx-1,0:ny-1,nuvz,2))
allocate (uuh(0:nx-1,0:ny-1,nuvz,2),vvh(0:nx-1,0:ny-1,nuvz,2))
allocate (uupol(0:nx-1,0:ny-1,nuvz,2),vvpol(0:nx-1,0:ny-1,nuvz,2))
allocate (tth(0:nx-1,0:ny-1,nuvz,2))
allocate (ps(0:nx-1,0:ny-1,1,2))
wwh(:,:,:,:)   = MISSING
uuh(:,:,:,:)   = MISSING
vvh(:,:,:,:)   = MISSING
uupol(:,:,:,:) = MISSING
vvpol(:,:,:,:) = MISSING
tth(:,:,:,:)   = MISSING
ps(:,:,:,:)    = MISSING

end subroutine alloc_era5

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#


!  name: era5.gridcheck_era5
!  @param
!  @return
!  
  subroutine gridcheck_era5(error)

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
!                          20-05-2020       
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
   use eccodes
   logical, intent(out) :: error
   integer :: ifn,ll
   real(dbl) :: xaux1,xaux2,xaux3
   real(dp) :: yfirst,ylast,xlast,ylat1
   real(dbl) :: sizesouth,sizenorth
   
   integer :: vvv
   integer :: DateSize, TimeSize, ParamSize, LevelSize
   character(len=8), allocatable :: DateList(:),ParamList(:)
   character(len=4), allocatable :: TimeList(:)
   integer, allocatable:: Levelist(:)

   integer :: iret, igrib
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
   print *,'gridcheck > ',path(3)(1:len_path(3))//trim(wfname(ifn))
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
   
   print *,'gridcheck > ',path(3)(1:len_path(3))//year//'/'//'uvwt-'//year//'-'//month//'.gribidx'
   call codes_index_read(idx_uvwt, &
     path(3)(1:len_path(3))//year//'/'//'uvwt-'//year//'-'//month//'.gribidx')
   !idx_ps = idx_uvwt
   !call codes_index_read(idx_hr,&
   !  path(3)(1:len_path(3))//year//'/'//'hr-'//year//'-'//month//'.gribidx')
   print *,'gribcheck_era5> index opened'
   
! Test the uvwt index
  
  ! get the number of distinct values of Date in the index
  call codes_index_get_size(idx_uvwt,'mars.date',DateSize)
  ! allocate the array to contain the list of distinct Date
  allocate(DateList(DateSize))
  ! get the list of distinct Date from the index
  call codes_index_get(idx_uvwt,'mars.date',DateList)
  write(*,'(a,i3)') 'DateSize=',DateSize
  ! get the number of distinct values of Time in the index
  call codes_index_get_size(idx_uvwt,'mars.time',TimeSize)
  ! allocate the array to contain the list of distinct Time
  allocate(TimeList(TimeSize))
  ! get the list of distinct Time from the index
  call codes_index_get(idx_uvwt,'mars.time',TimeList)
  write(*,'(a,i3)') 'TimeSize=',TimeSize
  print *,TimeList
  ! get the number of distinct values of Param in the index
  call codes_index_get_size(idx_uvwt,'mars.param',ParamSize)
  ! allocate the array to contain the list of distinct Param
  allocate(ParamList(ParamSize))
  ! get the list of distinct Param from the index
  call codes_index_get(idx_uvwt,'mars.param',ParamList)
  write(*,'(a,i3)') 'ParamSize=',ParamSize
  print *,ParamList
  ! get the number of distinct values of Levelist in the index
  call codes_index_get_size(idx_uvwt,'mars.levelist',LevelSize)
  ! allocate the array to contain the list of distinct Levelist
  allocate(Levelist(LevelSize))
  ! get the list of distinct Levelist from the index
  call codes_index_get(idx_uvwt,'mars.levelist',Levelist)
  write(*,'(a,i3)') 'LevelSize=',LevelSize
  NbLev=LevelSize
  
! Test reading a time on uvwt index
  call codes_index_select(idx_uvwt,'mars.date',fulldate)
  call codes_index_select(idx_uvwt,'mars.time',TimeList(1))
  call codes_index_select(idx_uvwt,'mars.param',ParamList(1))
  call codes_index_select(idx_uvwt,'mars.levelist',1)
  ! get the handle from the selection
  call codes_new_from_index(idx_uvwt,igrib, iret)
  ! first check that we read GRIB2
  call codes_get_int(igrib,'editionNumber',gribVer,iret)
  call codes_check(iret,gribFunction,gribErrorMsg)
  if (gribVer/=2) then 
     print *,'gridcheck_era5 expecting grib 2 file'
     stop 155
  endif
  call codes_get_int(igrib,'table2Version',vvv,iret)
  print *,'Table 2 Version ',vvv
  call codes_get_int(igrib,'indicatorOfParameter',vvv,iret)
  print *,'indicatorOfParameter ',vvv
  call codes_get_int(igrib,'paramId',vvv,iret)
  print *,'paramId ',vvv
  call codes_get_int(igrib,'Ni',NbLon,iret)
  call codes_get_int(igrib,'Nj',NbLat,iret)
  print *,'NbLon NbLat ',NbLon,NbLat
  call codes_get_real8(igrib,'longitudeOfFirstGridPointInDegrees',xaux1,iret)
  call codes_get_real8(igrib,'longitudeOfLastGridPointInDegrees',xaux2,iret)
  call codes_get_real8(igrib,'iDirectionIncrementInDegrees',xaux3,iret)
  xlon0=real(xaux1,kind=dp)
  xlast=real(xaux2,kind=dp)
  dx=real(xaux3,kind=dp)
  dx=1._dp ! correcting inaccurate value read from file
  dxconst=180._dp/(dx*r_earth*pi)
  !if (xlon0 > 180._dp) xlon0=xlon0-360._dp
  print *,'xlon0 xlast dx ',xlon0,xlast,dx
  call codes_get_real4(igrib,'latitudeOfFirstGridPointInDegrees',yfirst,iret)
  call codes_get_real4(igrib,'latitudeOfLastGridPointInDegrees',ylast,iret)
  print *,'yfirst ylast ',yfirst, ylast
  if (yfirst > ylast ) then
    print *,'lats from North to South'
  else
    print *,'lats from South to North'
  endif
  deallocate(DateList,TimeList,ParamList,Levelist)
  ! release the handle and the index
  call codes_release(igrib)
  !call codes_index_release(idx_uvwt)
  
! Set the vertical grid
! No need for a shifted grid for w as we do not intend for now to make
! kinematic calculations
  nuvz=NbLev
  nwz=NbLev
  allocate (akz(NbLev),bkz(NbLev))
  call setabl(akz,bkz)
  allocate (etakzlog(NbLev))
  etakzlog=log(akz+bkz)
! ACHTUNG, here akz becomes the pressure level  
  akz=p0*akz
  
! Check whether the grid is cyclic in longitudes
  if (xaux2-xaux1+xaux3-360.D0 < 0.001D0) then
    ! field is cyclic
    xglobal=.true.
    ! Add a last longitude that repeats the first for interpolations
    nx=NbLon+1
    print *,'Global longitude grid'
  else
    xglobal=.false.
    nx=NbLon
    print *,'Non global longitude grid'
  endif
  
! Define the latitude grid (assuming that the grid goes pole to pole) 
  ny = NbLat
  NbHLat=NbLat/2
  ylat0 = min(yfirst,ylast)
  ylat1 = max(yfirst,ylast)
  dy = (ylat1-ylat0)/(NbLat-1)
  dyconst = 180._dp/(dy*r_earth*pi)
  print *,'ylat0 dy ',ylat0,dy
  
! Set zmaxzmax=-log(akm(nwz)/p0) 
  zmax = -log(akz(nwz)/p0)
    
! Check whether the grid is global and whether the poles are included 
 
  if (xglobal.and.(min(yfirst,ylast)<-89._dp)) then
    ! field contains south pole
    sglobal=.true.
    ! Enhance the map scale by factor 3 (*2=6) compared to north-south
    ! map scale
    ! all polar transforms done in real*8
    sizesouth=6.d0*(switchsouth+90.d0)/dy
    if (oldpole) then
      call stlmbr(southpolemap,-90.d0, 0.d0)
      call stcm2p(southpolemap,0.d0,0.d0,switchsouth,0.d0,sizesouth, &
      sizesouth,switchsouth,180.d0)
    endif
    switchsouthg=(real(switchsouth,kind=dp)-ylat0)/dy
    print *,'South pole in the grid'
  else
    sglobal=.false.
    switchsouthg=999999._dp
  endif
  
  if (xglobal.and.(max(yfirst,ylast)>89._dp)) then
    ! field contains south pole
    nglobal=.true.
    ! Enhance the map scale by factor 3 (*2=6) compared to north-south
    ! map scale
    ! all polar transforms done in real*8
    sizenorth=6.d0*(90.d0-switchnorth)/dy
    if (oldpole) then
      call stlmbr(northpolemap,90.d0, 0.d0)
      call stcm2p(northpolemap,0.d0,0.d0,switchnorth,0.d0,sizenorth, &
      sizenorth,switchnorth,180.d0)
    endif
    switchnorthg=(real(switchnorth,kind=dp)-ylat0)/dy
    print *,'North pole in the grid'
  else
    nglobal=.false.
    switchnorthg=999999._dp
  endif
  write(*,'(a,2i5,2L3)')' gribcheck> nx,ny,nglobal,sglobal ', &
           nx, ny, nglobal, sglobal
  write(*,'(a,2f10.2)')' gribcheck> switchsouthg, switchnorthg ', &
           switchsouthg, switchnorthg
  
  return        
  end subroutine gridcheck_era5
  
  subroutine setabl(akz,bkz)
! gives pressure as p = akz+ps*bkz
! where pressure and ps must be in Pascal
! akz is here defined as a coefficient to be multiplied by p0 = 100000 Pa
! The multiplication is done in gridcheck 
  real(dp), intent(inout) :: akz(:),bkz(:)
  real(dbl) :: aa(137),bb(137)
  
data aa  / &
 0.00000000000000D0, 0.00001878906500D0, 0.00013296875500D0, &
 0.00042808594000D0, 0.00092441406000D0, 0.00162124999750D0, &
 0.00252312499750D0, 0.00363445312000D0, 0.00496238281000D0, &
 0.00651527344000D0, 0.00830750000000D0, 0.01034878906000D0, &
 0.01265398437500D0, 0.01523511718500D0, 0.01810488281000D0, &
 0.02127871094000D0, 0.02476691406500D0, 0.02858203125000D0, &
 0.03273250000000D0, 0.03722597656500D0, 0.04206667969000D0, &
 0.04725585937500D0, 0.05279089843500D0, 0.05866457031000D0, &
 0.06486476562500D0, 0.07137382812500D0, 0.07816859375000D0, &
 0.08521914062500D0, 0.09248984375000D0, 0.09993845703000D0, &
 0.10751740234500D0, 0.11517322266000D0, 0.12284798828000D0, &
 0.13048013671500D0, 0.13800546875000D0, 0.14535888672000D0, &
 0.15247574218500D0, 0.15929371093500D0, 0.16575367187500D0, &
 0.17180263672000D0, 0.17739382812500D0, 0.18248316406000D0, &
 0.18703583984500D0, 0.19103843750000D0, 0.19449399414000D0, &
 0.19741298828500D0, 0.19980555664500D0, 0.20168298828500D0, &
 0.20305664062500D0, 0.20393767578000D0, 0.20433898437500D0, &
 0.20427193359500D0, 0.20375085937500D0, 0.20278763671500D0, &
 0.20139797851500D0, 0.19959661133000D0, 0.19739716797000D0, &
 0.19481777344000D0, 0.19187400390500D0, 0.18858503906000D0, &
 0.18497076172000D0, 0.18105027344000D0, 0.17684617187500D0, &
 0.17238201171500D0, 0.16768055664000D0, 0.16276718750000D0, &
 0.15767186035000D0, 0.15241936035000D0, 0.14703877441000D0, &
 0.14156735351500D0, 0.13603000000000D0, 0.13045770996000D0, &
 0.12489210449500D0, 0.11935807617500D0, 0.11388364746000D0, &
 0.10850646972500D0, 0.10325305176000D0, 0.09814330566500D0, &
 0.09319541504000D0, 0.08842462891000D0, 0.08383939697500D0, &
 0.07943383056500D0, 0.07519640625000D0, 0.07111869873000D0, &
 0.06719408691000D0, 0.06341510498000D0, 0.05977209473000D0, &
 0.05625667725000D0, 0.05286442871500D0, 0.04959522217000D0, &
 0.04644983398500D0, 0.04342874023500D0, 0.04053210571500D0, &
 0.03775979370500D0, 0.03511105957000D0, 0.03258431518500D0, &
 0.03017755371000D0, 0.02788869751000D0, 0.02571559326500D0, &
 0.02365601074500D0, 0.02170763794000D0, 0.01986807617000D0, &
 0.01813484131000D0, 0.01650535950000D0, 0.01497696533500D0, &
 0.01354690247000D0, 0.01221231995000D0, 0.01097027405000D0, &
 0.00981773041000D0, 0.00875156372000D0, 0.00776855987500D0, &
 0.00686542023000D0, 0.00603876678500D0, 0.00528514740000D0, &
 0.00460104263500D0, 0.00398287430000D0, 0.00342701553500D0, &
 0.00292980163500D0, 0.00248754264500D0, 0.00209653755000D0, &
 0.00175308983000D0, 0.00145352455500D0, 0.00119420624000D0, &
 0.00097155811500D0, 0.00078208229000D0, 0.00062238019500D0, &
 0.00048917352500D0, 0.00037932476000D0, 0.00028985714000D0, &
 0.00021797324500D0, 0.00016107177500D0, 0.00011676195000D0, &
 0.00008287471500D0, 0.00005747030500D0, 0.00003884162500D0, &
 0.00002551303000D0, 0.00001000182500D0 /
data bb  / &
 0.9988150D0, 0.9963165D0, 0.9934935D0, 0.9902420D0, 0.9865210D0, &
 0.9823070D0, 0.9775750D0, 0.9722955D0, 0.9664325D0, 0.9599510D0, &
 0.9528070D0, 0.9449620D0, 0.9363705D0, 0.9269885D0, 0.9167720D0, &
 0.9056740D0, 0.8936540D0, 0.8806685D0, 0.8666805D0, 0.8516565D0, &
 0.8355685D0, 0.8183960D0, 0.8001265D0, 0.7807575D0, 0.7602975D0, &
 0.7387680D0, 0.7162040D0, 0.6926560D0, 0.6681895D0, 0.6428860D0, &
 0.6168420D0, 0.5901700D0, 0.5629965D0, 0.5354600D0, 0.5077095D0, &
 0.4799015D0, 0.4521970D0, 0.4247580D0, 0.3977440D0, 0.3713085D0, &
 0.3455965D0, 0.3207685D0, 0.2969760D0, 0.2742980D0, 0.2527430D0, &
 0.2322885D0, 0.2129120D0, 0.1945900D0, 0.1772995D0, 0.1610175D0, &
 0.1457190D0, 0.1313805D0, 0.1179765D0, 0.1054835D0, 0.0938740D0, &
 0.0831220D0, 0.0732030D0, 0.0640880D0, 0.0557505D0, 0.0481605D0, &
 0.0412870D0, 0.0351010D0, 0.0295700D0, 0.0246595D0, 0.0203365D0, &
 0.0165670D0, 0.0133110D0, 0.0105335D0, 0.0081970D0, 0.0062555D0, &
 0.0046745D0, 0.0034140D0, 0.0024245D0, 0.0016725D0, 0.0011215D0, &
 0.0007260D0, 0.0004510D0, 0.0002695D0, 0.0001555D0, 0.0000855D0, &
 0.0000415D0, 0.0000155D0, 0.0000035D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, 0.0000000D0, &
 0.0000000D0, 0.0000000D0  /

  if(size(akz)/=137) then
    print *,'Bad number of vertical levels ',size(akz)
    stop
  endif
  
  ! beware: akz is multiplied per p0 in gridcheck
  akz=real(aa,kind=dp)
  bkz=real(bb,kind=dp)
  
  
  
  return
  end subroutine setabl

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine getfields_era5(itime,nstop)

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
        write(*,*) 'TRACZILLA WARNING: NO era5 WINDS ARE AVAILABLE.'
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
              call read_era5(indj+1,memind(2))
              call verttransform_era5(memind(2))
              memtime(2)=wftime(indj+1)
              write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_era5   > date ',trim(wfname(indj+1)),&
                      '     memtime ',memtime(2),'  time ',itime
              nstop = 1
              goto 40
           endif
        enddo
 40     indmin=indj

      else

! No wind fields, which can be used, are currently in memory 
! -> read both wind fields
!***********************************************************
         do indj=indmin,numbwf-1
            if ((ldirect*wftime(indj).le.ldirect*itime).and.   &
                  (ldirect*wftime(indj+1).gt.ldirect*itime)) then      
               memind(1)=1
               call read_era5(indj,memind(1))
               call verttransform_era5(memind(1))
               memtime(1)=wftime(indj)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_era5   > date ',trim(wfname(indj)),&
                      '     memtime ',memtime(1),'  time ',itime
               memind(2)=2
               call read_era5(indj+1,memind(2))
               call verttransform_era5(memind(2))
               memtime(2)=wftime(indj+1)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_era5   > date ',trim(wfname(indj+1)),&
                      '     memtime ',memtime(2),'  time ',itime
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
      ifdiab: if (era5_diab) then

        if ((ldirect*wftime_diab(1).gt.ldirect*itime).or.    &
            (ldirect*wftime_diab(numbwf_diab).lt.ldirect*itime)) then
          write(*,*) 'TRACZILLA WARNING: NO era5 HEATINGS ARE AVAILABLE.'
          write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
          write(*,*) ldirect*wftime_diab(1)
          write(*,*) ldirect*itime
          write(*,*) ldirect*wftime_diab(numbwf_diab)
          nstop=4
          return
        endif

        ifread: if ((ldirect*memtime_diab(1).le.ldirect*itime).and.   &
                    (ldirect*memtime_diab(2).gt.ldirect*itime)) then

! The right heating rate fields are already in memory -> don't do anything
!*****************************************************************
          continue

        else if ((ldirect*memtime_diab(2).le.ldirect*itime).and. &
           (memtime_diab(2).ne.999999999)) then ifread
 
! Current time is after 2nd wind field
! -> Resort heating rate field pointers, so that current time is between 1st and 2nd
!***************************************************************************

          memaux=memind_diab(1)
          memind_diab(1)=memind_diab(2)
          memind_diab(2)=memaux
          memtime_diab(1)=memtime_diab(2)

! Read a new heating rate field and store it on place memind(2)
!******************************************************

          do indj=indmin_diab,numbwf_diab-1
             if (ldirect*wftime_diab(indj+1).gt.ldirect*itime) then
                call read_era5_diab(indj+1,memind_diab(2))
                call verttransform_era5_diab(memind_diab(2),2)
                memtime_diab(2)=wftime_diab(indj+1)
                write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_era5_hr> date ',trim(wfname_diab(indj+1)),&
                      '  memtime ',memtime_diab(2),'  time ',itime

                nstop = 1
                goto 45
             endif
          enddo
 45       indmin_diab=indj

        else ifread

! No wind fields, which can be used, are currently in memory 
! -> read both heating rate fields
!***********************************************************
        ! Force read of the index
        current_ym_hr = '190001'  
        do indj=indmin_diab,numbwf_diab-1
             if ((ldirect*wftime_diab(indj).le.ldirect*itime).and.   &
                 (ldirect*wftime_diab(indj+1).gt.ldirect*itime)) then
               memind_diab(1)=1  
               print *,'first read a couple of fields'        
               call read_era5_diab(indj,memind_diab(1))
               call verttransform_era5_diab(memind_diab(1),1)
               memtime_diab(1)=wftime_diab(indj)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_era5_hr> date ',trim(wfname_diab(indj)),&
                      '  memtime ',memtime_diab(1),'  time ',itime
               memind_diab(2)=2
               call read_era5_diab(indj+1,memind_diab(2))
               call verttransform_era5_diab(memind_diab(2),2)
               memtime_diab(2)=wftime_diab(indj+1)
               write(*,'(a,a,a,i11,a,i11)') &
                      ' getfields_era5_hr> date ',trim(wfname_diab(indj+1)),&
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
      end subroutine getfields_era5

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine read_era5(indj,n)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READ_era5                  *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      B. Legras                    *
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS ON MODEL LEVELS              *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!
! This routine uses indexes for the grib files containing t,u,v,w and ps  
! Wind and temperature fields are read into temporary arrays
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

  use eccodes
  integer, intent(in) :: indj,n
  integer :: ll,hh
  integer :: l,lat
  character (len=2) :: month
  character (len=4) :: year
  character (len=6) :: yearmonth
  character (len=8) :: fulldate
  integer :: igrib,iret
  real*4, allocatable :: values(:),vals(:,:)
   
! Finds the date from wfname
  !print *,'read_era5> wfname ',wfname(indj)
  ll=len_trim(wfname(indj))
  !hour=wfname(indj)(ll-1:ll)
  year=wfname(indj)(1:4)
  month=wfname(indj)(5:6)
  fulldate=wfname(indj)(1:8)
  yearmonth=wfname(indj)(1:6)
  ! decode hh
  read(wfname(indj)(ll-3:ll),'(i4)') hh
  !print *,'hh fulldate',hh,' ',fulldate
  
! Loads new month when needed
! the index is read nb_var time as we do not know what contains index structure and
! how to copy it
  !print *,'read_era5> ',current_ym,' ',yearmonth,' ',fulldate
  if (yearmonth /= current_ym) then
    ! release previous index
    call codes_index_release(idx_uvwt)
    ! get new one
    call codes_index_read(idx_uvwt, &        
      path(3)(1:len_path(3))//year//'/'//'uvwt-'//year//'-'//month//'.gribidx')
    !idx_ps = idx_uvwt
    print *,'New uvwt ps idx loaded, year month ',year,month
    current_year=year
    current_ym=year//month
  endif
  
! Allocate temp fields out of the parallel loop
  !*** if(z_motion) then
  !***  allocate(values(NbLon*NbLat,5),val_gauss(NbLon,NbLat,5))
  !*** else
  allocate(vals(NbLon,NbLat))
  !*** endif
  
! Read one variable (or file) in each section
! Temperature
  call codes_index_select(idx_uvwt,'mars.date',fulldate)
  call codes_index_select(idx_uvwt,'mars.time',hh)
  call codes_index_select(idx_uvwt,'mars.param',130)
  do l=1,NbLev
     call codes_index_select(idx_uvwt,'mars.levelist',l)
     call codes_new_from_index(idx_uvwt,igrib, iret)
     call codes_get_real4_array(igrib,'values',values)
     call codes_release(igrib)
     vals = reshape(values,(/NbLon,NbLat/))
! This loop inverts the latitude and the altitude orders
! npw from bottom to top and south to north
     do lat=1,NbLat
        tth(0:NbLon-1,NbLat-lat,Nblev-l+1,n) = vals(:,lat)
     enddo
     ! Replication of the first column as the last in the case of a global grid
     if(xglobal) tth(NbLon,0:NbLat-1,Nblev-l+1,n) = tth(0,0:NbLat-1,Nblev-l+1,n)
     deallocate(values)
  enddo
! Zonal wind
  !call codes_index_select(idx_uvwt,'mars.date',fulldate)
  !call codes_index_select(idx_uvwt,'mars.time',hh)
  call codes_index_select(idx_uvwt,'mars.param',131)
  do l=1,NbLev
     call codes_index_select(idx_uvwt,'mars.levelist',l)
     call codes_new_from_index(idx_uvwt,igrib,iret)
     call codes_get_real4_array(igrib,'values',values)
     call codes_release(igrib)
     vals = reshape(values,(/NbLon,NbLat/))  
     do lat=1,NbLat
        uuh(0:NbLon-1,NbLat-lat,Nblev-l+1,n) = vals(:,lat)                         
     enddo
     if(xglobal) uuh(NbLon,0:NbLat-1,Nblev-l+1,n) = uuh(0,0:NbLat-1,Nblev-l+1,n)
     deallocate(values)
  enddo
! Meridional wind
  !call codes_index_select(idx_uvwt,'mars.date',fulldate)
  !call codes_index_select(idx_uvwt,'mars.time',hh)
  call codes_index_select(idx_uvwt,'mars.param',132)
  do l=1,NbLev
     call codes_index_select(idx_uvwt,'mars.levelist',l)
     call codes_new_from_index(idx_uvwt,igrib, iret)
     call codes_get_real4_array(igrib,'values',values)
     call codes_release(igrib)
     vals = reshape(values,(/NbLon,NbLat/))
     do lat=1,NbLat
        vvh(0:NbLon-1,NbLat-lat,Nblev-l+1,n) = vals(:,lat)        
     enddo
     if(xglobal) vvh(NbLon,0:NbLat-1,Nblev-l+1,n) = vvh(0,0:NbLat-1,Nblev-l+1,n)
     deallocate(values)
  enddo
! Surface pressure (what is read is the log of pressure)
  !call codes_index_select(idx_uvwt,'mars.date',fulldate)
  !call codes_index_select(idx_uvwt,'mars.time',hh)
  call codes_index_select(idx_uvwt,'mars.param',152)
  call codes_index_select(idx_uvwt,'mars.levelist',1)
  call codes_new_from_index(idx_uvwt,igrib,iret)
  call codes_get_real4_array(igrib,'values',values)
  call codes_release(igrib)       
  vals = reshape(values,(/Nblon,Nblat/))
  do lat=1,Nblat
      ps(0:NbLon-1,NbLat-lat,1,n) = exp(vals(:,lat))                    
  enddo
  if(xglobal) ps(NbLon,0:NbLat-1,1,n) = ps(0,0:NbLat-1,1,n)
  deallocate(values)
! Vertical Dp/Dt
  if (z_motion) then
!    call codes_index_select(idx_uvwt(4),'mars.date',fulldate)
!    call codes_index_select(idx_uvwt(4),'mars.time',hh)
    call codes_index_select(idx_uvwt,'mars.param',135)
    do l=1,NbLev
       call codes_index_select(idx_uvwt,'mars.levelist',l)
       call codes_new_from_index(idx_uvwt,igrib, iret)
       call codes_get_real4_array(igrib,'values',values)        
       vals=reshape(values,(/Nblon,Nblat/))
       call codes_release(igrib)
       do lat=1,NbLat
         wwh(0:NbLon-1,NbLat-lat,Nblev-l+1,n) = vals(:,lat)        
       enddo
       if(xglobal) wwh(NbLon,0:NbLat-1,Nblev-l+1,n) = wwh(0,0:NbLat-1,Nblev-l+1,n)
       deallocate(values)
    enddo
  endif

  deallocate(vals) 
  
  return

  end subroutine read_era5
  
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine read_era5_diab(indj,n)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READ_era5                  *
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

  use eccodes
  integer, intent(in) :: indj,n
  integer :: ll,hh,step
  integer :: l,lat
  character (len=2) :: month
  character (len=4) :: year
  character (len=6) :: yearmonth
  character (len=8) :: fulldate
!  integer :: OMP_GET_THREAD_NUM
  real*4, allocatable :: values(:),swr(:,:,:),lwr(:,:,:)
  integer :: igrib, iret
! Finds the date and time from wfname
  ll=len_trim(wfname_diab(indj))
  !hour=wfname_diab(indj)(ll-1:ll)
  year=wfname_diab(indj)(1:4)
  month=wfname_diab(indj)(5:6)
  yearmonth=wfname_diab(indj)(1:6)
  fulldate=wfname_diab(indj)(1:8)
  
  read(wfname_diab(indj)(ll-6:ll),'(i4,1x,i2)')hh,step
  !print*,'hh step fulldate wfname',hh,step,' ',fulldate,' ',wfname_diab(indj)
  
! Loads new year when needed
! reading is performed twice as we do not know how to copy index structure
  if (yearmonth /= current_ym_hr) then
    ! release previous index
    call codes_index_release(idx_hr)
    call codes_index_read(idx_hr, &
      path(3)(1:len_path(3))//year//'/'//'hr-'//year//'-'//month//'.gribidx')
    print *,'New hr idx loaded, year month ',year, month
    current_ym_hr=year//month
  endif
 
  allocate (lwr(NbLon,NbLat,NbLev),swr(NbLon,NbLat,NbLev))
  
! Short wave
  call codes_index_select(idx_hr,'mars.date',fulldate)
  call codes_index_select(idx_hr,'mars.time',hh)
  call codes_index_select(idx_hr,'mars.step',step)
  call codes_index_select(idx_hr,'mars.param',235001)
  do l=1,NbLev
     call codes_index_select(idx_hr,'mars.levelist',l)
     call codes_new_from_index(idx_hr,igrib, iret)   
     call codes_get_real4_array(igrib,'values',values)
     call codes_release(igrib)       
     swr(:,:,Nblev-l+1)=reshape(values,(/Nblon,Nblat/))
     deallocate(values)
  enddo
! Long wave
  !call codes_index_select(idx_hr,'mars.date',fulldate)
  !call codes_index_select(idx_hr,'mars.time',hh)
  !call codes_index_select(idx_hr,'mars.step',step)
  call codes_index_select(idx_hr,'mars.param',235002)
  do l=1,NbLev
     call codes_index_select(idx_hr,'mars.levelist',l)
     call codes_new_from_index(idx_hr,igrib, iret)
     call codes_get_real4_array(igrib,'values',values)
     call codes_release(igrib)      
     lwr(:,:,Nblev-l+1)=reshape(values,(/Nblon,Nblat/))
     deallocate(values)
  enddo

  ! Add the two components and revert lat
  do lat=1,NbLat
     wwh(0:NbLon-1,NbLat-lat,1:NbLev,n) = swr(1:NbLon,lat,1:NbLev) &
                                        + lwr(1:NbLon,lat,1:NbLev)
  enddo
  if(xglobal) wwh(NbLon,0:NbLat-1,1:NbLev,n) = wwh(0,0:NbLat-1,1:NbLev,n)

  deallocate(lwr,swr)
 
  return

  end subroutine read_era5_diab
  
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine verttransform_era5(n)

  use coord
  use polarproj
  integer, intent(in):: n
  integer :: ix, jy, iz, magic
  !pb real(dbl) :: xlon, ylat
  real(dbl) :: xlon, ylat

! north pole region
  if (nglobal) then
  if (oldpole) then
#if defined(PAR_RUN)
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(ix,jy,iz,xlon,ylat)
#endif
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
#if defined(PAR_RUN)
!$OMP END DO
!$OMP END PARALLEL
#endif
! As a temporary fix, the velocity at the pole in polar coordinates is
! replaced by the average over the closest latitude circle.
! Even id only the horizontal component is at trouble, we include the 
! vertical component as well.
    do iz=1,nuvz
      uupol(0:nx-1,ny-1,iz,n) = sum(uupol(0:nx-1,ny-2,iz,n))/nx
      vvpol(0:nx-1,ny-1,iz,n) = sum(vvpol(0:nx-1,ny-2,iz,n))/nx
      wwh(0:nx-1,ny-1,iz,n) = sum(wwh(0:nx-1,ny-2,iz,n))/nx
    enddo
! New interpolation method with polarproj
  else
    if (magicpole) then
      if (uuh(0,ny-1,20,n) == uuh(int(nx/4),ny-1,20,n)) then
        magic = 1
      else
        magic = 0
      endif
    else
      magic = 2
    endif
    call uvinNpol(n,magic)
  endif
   
  endif

!    south pole region
  if (sglobal) then
  if (oldpole) then
#if defined(PAR_RUN)
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(ix,jy,iz,xlon,ylat)
#endif
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
#if defined(PAR_RUN)
!$OMP ENDDO
!$OMP END PARALLEL
#endif

   ! As a temporary fix, the velocity at the pole in polar coordinates is
   ! replaced by the average over the closest latitude circle.
   ! Even id only the horizontal component is at trouble, we include the 
   ! vertical component as well.
    do iz=1,nuvz
      uupol(0:nx-1,0,iz,n) = sum(uupol(0:nx-1,1,iz,n))/nx
      vvpol(0:nx-1,0,iz,n) = sum(vvpol(0:nx-1,1,iz,n))/nx
      wwh(0:nx-1,0,iz,n) = sum(wwh(0:nx-1,1,iz,n))/nx
    enddo
  ! New interpolation method with polarproj
  else
    if (magicpole) then
      if (uuh(0,0,20,n) == uuh(int(nx/4),0,20,n)) then
        magic = 1
      else
        magic = 0
      endif
    else
      magic = 2
    endif
    call uvinSpol(n,magic)
  endif

  endif

  return
  end subroutine verttransform_era5
 
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine verttransform_era5_diab(n,m)
  
  integer, intent(in):: n,m
  integer :: iz,ix,jy
  real(dp) :: pp 

! No conversion needed as data are in K/s
! Conversion in D theta / D t
#if defined(PAR_RUN)
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC) PRIVATE(iz,jy,ix,pp)
#endif      
  do iz=1,nuvz
     do jy=0,ny-1
        do ix=0,nx-1
           pp=akz(iz)+bkz(iz)*(ps(ix,jy,1,1)+(ps(ix,jy,1,2)))/2
           wwh(ix,jy,iz,n)=wwh(ix,jy,iz,n)*((p0/pp)**kappa)
        enddo
     enddo
     ! The nex statement is missing in jra_55 file
     if (xglobal) wwh(NbLon,0:NbLat-1,1:NbLev,n) = wwh(0,0:NbLat-1,1:NbLev,n)
  enddo
#if defined(PAR_RUN)
!$OMP END PARALLEL DO
#endif
  !print *,'wwh just read, max, min ',maxval(maxval(wwh(:,:,35,n),DIM=2),DIM=1),&
  !                                   minval(minval(wwh(:,:,35,n),DIM=2),DIM=1)     
        
! call mass correction when required
  if(mass_correction) call diab_mass_era5(n,m)
 
  return
  end subroutine verttransform_era5_diab
  
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass_era5(n,m)

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

#if defined(PAR_RUN)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(k)
#endif
! Calculation of the potential temperature
  do k=nuvz-NPureP_era5,nuvz
     theta(:,:,k)=0.5*(tth(:,:,k,1)+tth(:,:,k,2))*pif_era5(k)
!    theta(:,:,k)=0.5*(theta_g(:,:,k,1)+theta_g(:,:,k,2))
  enddo
#if defined(PAR_RUN)
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(k)
#endif
  do k=nuvz-NPureP_era5+1,nuvz-1
     sigma(:,:,k)=pmc_era5(k)/(theta(:,:,k+1)-theta(:,:,k-1))
  enddo
#if defined(PAR_RUN)
!$OMP END DO
#endif
  sigma(:,:,nuvz)=pmc_era5(nuvz)/(theta(:,:,nuvz)-theta(:,:,nuvz-1))
  
! Calculation of the mass flux across the surface
! Note that mean_sigma is not a mean of sigma but the spherical
! integral of sigma (that divides the flux to get the correction)
  !allocate(flux(0:nx-1,0:ny-1))
  !allocate(flux_lat(0:ny-1))
  !allocate(mean_sigma_lat(0:ny-1))
#if defined(PAR_RUN)
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(k)
#endif
  do k=nuvz-NPureP_era5+1,nuvz
     flux(:,:,k)=wwh(:,:,k,n)*sigma(:,:,k)
     flux_lat(:,k)=sum(flux(0:nx-2,:,k),DIM=1)
     mean_sigma_lat(:,k)=sum(sigma(0:nx-2,:,k),DIM=1)
     mass_flux(k)=0.5*dot_product(flux_lat(:,k),area_coefft_era5)
     mean_sigma(k)=0.5*dot_product(mean_sigma_lat(:,k),area_coefft_era5)
     mean_w(k)=mass_flux(k)/mean_sigma(k)
  enddo
#if defined(PAR_RUN)
!$OMP END DO
#endif
  !deallocate(flux,flux_lat,mean_sigma_lat)
#if defined(PAR_RUN)
!$OMP END PARALLEL
#endif
  deallocate(flux,flux_lat,mean_sigma_lat) 

! Output of the averaged diab heating
  if(mean_diab_output) then
    write(unitflux) memtime_diab(m),mean_w(nuvz-NPureP_era5+1:nuvz), &
                                    mean_sigma(nuvz-NPureP_era5+1:nuvz)
    flush(unitflux)
  endif

! Correction step
#if defined(PAR_RUN)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(k)
#endif
  do k=nuvz-NPureP_era5+1,nuvz
    wwh(:,:,k,n)=wwh(:,:,k,n)-mean_w(k)
  enddo
#if defined(PAR_RUN)
!$OMP END PARALLEL DO
#endif
  deallocate(theta,sigma,mean_sigma)
  deallocate(mass_flux,mean_w)

  return
  end subroutine diab_mass_era5

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

  subroutine diab_mass_era5_init

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
  
! Numbers of pure pressure levels on which the correction is applied  
  NPureP_era5 = 54
 
  if(.not.(xglobal.and.nglobal.and.sglobal)) then
      write(*,*) ' #### TRACZILLA MODEL ERROR! DIAB_MASS        #### ' 
      write(*,*) ' #### mass flux cannot be balanced            #### '
      write(*,*) ' #### on not global grid                      #### '
      stop  
  endif

! calculate area coefficient assuming pole to pole regular grid
  allocate(area_coefft_era5(0:ny-1))
! it is assumed that dy = 180/(ny-1)
  do i=1,ny/2-1
     area_coefft_era5(ny/2-1+i)=2*sin(pi*dy/180._dp)*cos((i-0.5)*dy*pi/180._dp)
     area_coefft_era5(ny/2-i)=area_coefft_era5(ny/2-1+i)
  enddo
  area_coefft_era5(0)=2*cos(ylat0*pi/180._dp)*sin((dy/2+90._dp+ylat0)*pi/180._dp)
  area_coefft_era5(ny-1)=area_coefft_era5(0)
! test : the sum should be equal to 2g10.3 )') sum(area_coefft_era5)

! calculate weighting pressure factors in the vertical
  allocate(pmc_era5(nuvz),pif_era5(nuvz))
  pmc_era5=0.
  print *,' diab_mass_init > check sum area  '
  pif_era5=0.
  do i=nuvz-NPureP_era5+1,nuvz-1
     pmc_era5(i)=akz(i)*(Log(akz(i-1))-Log(akz(i+1)))
  enddo
  pmc_era5(nuvz)=akz(nuvz)*(Log(akz(nuvz-1))-Log(akz(nuvz)))

  do i=nuvz-NPureP_era5+1,nuvz
     pif_era5(i)=(p0/akz(i))**kappa
  enddo

! Create the output file of the mas averaged heating flux
! or open it in the append mode
  open(unitflux,file=trim(path(2))//'MassMeanDiab2.dat',form='unformatted',position='append')

  return
  end subroutine diab_mass_era5_init

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine interpol_wind_era5 &
         (jp,itime,xt,yt,zt, dxdt,dydt,dzdt, ngrid, & 
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
!    to be used for ERA-Interim, ERA5 and JRA-55 
!    If used for ERA-Interim or ERA5 kinematic, beware the vertical shift
!    in w is neglected. A cheap way would be to interpolate w. An even cheapest 
!    way would be to ignore the problem. 
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

      integer, intent(in) :: itime,ngrid,jp
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
      
      
      if(xglobal) then 
        ix=modulo(max(floor(xt),0),nx-1)
      else
        ix=min(max(floor(xt),0),nx-2)
      endif
      jy=min(max(floor(yt),0),ny-2)
      ddx=modulo(xt-float(ix),1.)
      ! accounts for the two polar regions in era5
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
      if (era5_diab.or.diabatic_w) then
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
      else if (era5_diab.or.diabatic_w) then
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
! Same indz is assumed on the 4 corners since we interpolate in eta 

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

      end subroutine interpol_wind_era5

!*******************************************************************************
! Identical to locuv in readinterp.f90
! This piece of code is private to the module
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

end module era5
