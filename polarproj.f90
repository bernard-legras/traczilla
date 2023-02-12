! Copyright or © or Copr. Bernard Legras (2023)

! bernard.legras@lmd.ipsl.fr

! This software is a computer program whose purpose is to perform Legendre
! transforms on the sphere.

! This software is governed by the CeCILL-C license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL-C
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".

! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.

! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and,  more generally, to use and operate it in the
! same conditions as regards security.

! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-C license and that you accept its terms.

module polarproj

use commons
real(dbl), parameter :: rdbl_earth = 6.371e6_dbl 

contains

subroutine ll2polxyNP(lon,lat,x,y)
  real(dbl), intent(in) :: lon, lat
  real(dbl), intent(out) :: x, y
  real(dbl) :: r
  !lat_rad = lat * pi /180
  r = 2*rdbl_earth * sqrt((1-sind(lat))/(1+sind(lat)))
  x = r * cosd(lon)
  y = r * sind(lon)
  return
end subroutine ll2polxyNP
  
subroutine ll2polxySP(lon,lat,x,y)
  real(dbl), intent(in) :: lon, lat
  real(dbl), intent(out) :: x, y
  real(dbl) :: r
  !lat_rad = lat * pi /180
  r = 2*rdbl_earth * sqrt((1+sind(lat))/(1-sind(lat)))
  x = r * cosd(lon)
  y = r * sind(lon)
  return
end subroutine ll2polxySP

subroutine polxy2llNP(lon,lat,x,y)
  real(dbl), intent(out) :: lon, lat
  real(dbl), intent(in) :: x, y
  real(dbl) :: r2
  lon = atan2d(y,x)
  r2 = (x**2 + y**2)/(4*rdbl_earth**2)
  lat = asind((1 - r2)/(1 + r2))
  return
end subroutine polxy2llNP
  
subroutine polxy2llSP(lon,lat,x,y)
  real(dbl), intent(out) :: lon, lat
  real(dbl), intent(in) :: x, y
  real(dbl) :: r2
  lon = atan2d(y,x)
  r2 = (x**2 + y**2)/(4*rdbl_earth**2)
  lat = - asind((1 - r2)/(1 + r2))
  return
end subroutine polxy2llSP
  
subroutine uvinNpol(n,magic)
  integer, intent(in) :: n, magic
  real(dp) :: xlon, ylat
  do jy=int(switchnorthg)-2,ny-2
    ylat=ylat0+float(jy)*dy     
    do ix=0,nx-1
      xlon=xlon0+float(ix)*dx     
      do iz=1,nuvz
        uupol(ix,jy,iz,n) = 2/(1+sind(ylat)) * (-cosd(xlon) * vvh(ix,jy,iz,n) - sind(xlon) * uuh(ix,jy,iz,n))
        vvpol(ix,jy,iz,n) = 2/(1+sind(ylat)) * (-sind(xlon) * vvh(ix,jy,iz,n) + cosd(xlon) * uuh(ix,jy,iz,n))
      enddo
    enddo
  enddo
  ! Processing of the line at the pole
  ! magic = 0: uuh, vvh are sinusoids in continuity with nearby lats
  ! magic = 1: single value for uuh and vvh, corresponding to expected value at 30E
  ! magic = 2: polar value defined as the average on the closest latitude circle
  if (magic == 0) then
    do iz=1,nuvz
      uupol(0:nx-1,ny-1,iz,n) = - vvh(0,ny-1,iz,n) 
      vvpol(0:nx-1,ny-1,iz,n) =   uuh(0,ny-1,iz,n)
    enddo
  else if (magic == 1) then
    do iz=1,nuvz
      uupol(0:nx-1,ny-1,iz,n) = -0.5_dp * (uuh(0,ny-1,iz,n) + vvh(0,ny-1,iz,n)*sqrt(3._dp))
      vvpol(0:nx-1,ny-1,iz,n) =  0.5_dp * (uuh(0,ny-1,iz,n)*sqrt(3._dp) - vvh(0,ny-1,iz,n))
    enddo
  else if (magic == 2) then
    do iz=1,nuvz
      uupol(0:nx-1,ny-1,iz,n) = sum(uupol(0:nx-1,ny-2,iz,n))/nx
      vvpol(0:nx-1,ny-1,iz,n) = sum(vvpol(0:nx-1,ny-2,iz,n))/nx
    enddo
  endif
  return
end subroutine uvinNpol

subroutine uvinSpol(n,magic)
  integer, intent(in) :: n, magic
  real(dp) :: xlon, ylat 
  do jy=1,int(switchsouthg)+3
    ylat=ylat0+float(jy)*dy
    do ix=0,nx-1
      xlon=xlon0+float(ix)*dx     
      do iz=1,nuvz
        uupol(ix,jy,iz,n) = 2/(1-sind(ylat)) * (cosd(xlon) * vvh(ix,jy,iz,n) - sind(xlon) * uuh(ix,jy,iz,n))
        vvpol(ix,jy,iz,n) = 2/(1-sind(ylat)) * (sind(xlon) * vvh(ix,jy,iz,n) + cosd(xlon) * uuh(ix,jy,iz,n))
      enddo
    enddo
  enddo 
! Processing of the line at the pole
  ! magic = 0: uuh, vvh are sinusoids in continuity with nearby lats
  ! magic = 1: single value for uuh and vvh, corresponding to expected value at 90W
  ! magic = 2: polar value defined as the average on the closest latitude circle
  if (magic == 0) then
    do iz=1,nuvz
      uupol(0:nx-1,0,iz,n) = vvh(0,0,iz,n) 
      vvpol(0:nx-1,0,iz,n) = uuh(0,0,iz,n)
    enddo
  else if (magic == 1) then
    do iz=1,nuvz
      uupol(0:nx-1,0,iz,n) =   uuh(0,0,iz,n)
      vvpol(0:nx-1,0,iz,n) = - vvh(0,0,iz,n) 
    enddo
  else if (magic == 2) then
    do iz=1,nuvz
      uupol(0:nx-1,0,iz,n) = sum(uupol(0:nx-1,1,iz,n))/nx
      vvpol(0:nx-1,0,iz,n) = sum(vvpol(0:nx-1,1,iz,n))/nx
    enddo
  endif
  return
end subroutine uvinSpol
  
end module polarproj
