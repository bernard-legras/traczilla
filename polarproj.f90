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
