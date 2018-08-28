
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
!=====================================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@[ TRACZILLA ]@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8======
!
module isentrop_m

!**************************************
! Isentropic interpolation module (methods)
!
! B. Legras, June 2005 
! 
!**************************************

use isentrop_h
implicit none

contains 
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine interpol_wind_theta   &
         (itime,xt,yt, theta, dxdt,dydt, ngrid, theta_inf, theta_sup, z_factor,tint)

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
      real(dp), intent(in) :: xt,yt,theta
      real(dp), intent(out) :: dxdt,dydt,z_factor,tint, theta_inf, theta_sup

! Auxiliary variables needed for interpolation
      real(dp) :: u1(2),v1(2),dt1,dt2,dtt,tp1(2)
      real(dp) :: tr(4,2),trp(4,2),u(4,2),v(4,2),tp(4,2)
      integer :: m,indexh,indz(4,2)
      integer :: ix,jy,ixp,jyp,i0,j0,idxy
      real(dp) :: ddx,ddy,rddx,rddy,p1,p2,p3,p4
      real(dp) :: psl0,psup0,pinf0,pisup0,piinf0

!********************************************
! Multilinear interpolation in time and space
!********************************************

! Determine the lower left corner and its distance to the current position
!*************************************************************************

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

! Calculates the theta values on the four adjacent columns if required
!*********************************************************************

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
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
         
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
            u(2,m)=(uupol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + uupol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvpol(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m))  &
                  + vvpol(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta))  &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uupol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + uupol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvpol(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m))  &
                  + vvpol(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta))  &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uupol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + uupol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvpol(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m))  &
                  + vvpol(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta))  &
                  / (tr(4,m)-trp(4,m))
          
          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
          
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
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(log(theta)-log(trp(2,m))) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh)*(log(tr(2,m))-log(theta))) &
                  / (log(tr(2,m))-log(trp(2,m)))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(log(theta)-log(trp(3,m))) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(log(tr(3,m))-log(theta))) &
                  / (log(tr(3,m))-log(trp(3,m)))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(log(theta)-log(trp(4,m))) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(log(tr(4,m))-log(theta))) &
                  / (log(tr(4,m))-log(trp(4,m)))

          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)
  
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
            u(2,m)=(uuh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + uuh(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            v(2,m)=(vvh(ixp,jy ,indz(2,m)  ,indexh)*(theta-trp(2,m)) &
                  + vvh(ixp,jy ,indz(2,m)+1,indexh)*(tr(2,m)-theta)) &
                  / (tr(2,m)-trp(2,m))
            u(3,m)=(uuh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + uuh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            v(3,m)=(vvh(ix ,jyp,indz(3,m)  ,indexh)*(theta-trp(3,m)) &
                  + vvh(ix ,jyp,indz(3,m)+1,indexh)*(tr(3,m)-theta)) &
                  / (tr(3,m)-trp(3,m))
            u(4,m)=(uuh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + uuh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))
            v(4,m)=(vvh(ixp,jyp,indz(4,m)  ,indexh)*(theta-trp(4,m)) &
                  + vvh(ixp,jyp,indz(4,m)+1,indexh)*(tr(4,m)-theta)) &
                  / (tr(4,m)-trp(4,m))

          u1(m)=p1*u(1,m)+p2*u(2,m)+p3*u(3,m)+p4*u(4,m)
          v1(m)=p1*v(1,m)+p2*v(2,m)+p3*v(3,m)+p4*v(4,m)     
        
        enddo
      
      end select     

      endif

! Calculation of z_factor for vertical diffusion

      select case (diftype)
      case (1)          ! diffusion in z
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
        z_factor = -ga*(log(trp(idxy,1))-log(tr(idxy,1)))/(piinf0-pisup0)
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
      if(AccurateTemp) tint=(tp1(1)*dt2+tp1(2)*dt1)*dtt
      
      if(debug_out) &
        print "('interpol>',i3,' P ',3f7.0,' T ',3f7.2,' TH ',3f7.2)", &
               indz(1,1),p0*(theta/tint)**(-1/kappa),&
               p0*(theta_g(indz(1,1),ix,jy,memind(1))/tth(ix,jy,indz(1,1),memind(1)))**(-1/kappa),&
               p0*(theta_g(indz(1,1)+1,ix,jy,memind(1))/tth(ix,jy,indz(1,1)+1,memind(1)))**(-1/kappa),&
               tint,tth(ix,jy,indz(1,1),memind(1)),tth(ix,jy,indz(1,1)+1,memind(1)), &
               theta,theta_g(indz(1,1),ix,jy,memind(1)),theta_g(indz(1,1)+1,ix,jy,memind(1)) 
      return

!************************************
      contains
!************************************      
  
      ! unused routine at the moment
      ! provides in i0, y0 the coordinates of the closest
      ! point
      subroutine sort_hor_dist()
      if (p2<p1) then
         if (p3<p1) then
            if (p4<p1) then
               i0=ix ; j0=jy; idxy=1
            else
               i0=ixp ; j0=jyp; idxy=4
            endif
         elseif (p4<p3) then
            i0=ixp; j0=jy; idxy=3
         else
            i0=ixp; j0=jyp; idxy=4
         endif
      elseif (p3<p2) then  
         if (p2<p4) then
            i0=ixp; j0=jyp; idxy=4
         else 
            i0=ix; j0=jyp; idxy=2 
         endif
      elseif (p4<p3) then
         i0=ixp ; j0=jy; idxy=3
      else
         i0=ixp ; j0=jyp; idxy=4
      endif
      end subroutine sort_hor_dist

   end subroutine interpol_wind_theta

!===============================================================================
! other routines not contained in interpol_wind_theta

      subroutine calc_col_theta(x,y,ind)
      integer, intent(in) :: x, y, ind
!     theta_diff  must be defined here to avoid // conflicts
      real, allocatable :: theta_diff(:) ! 
      real psl, tol
      integer k, id, id1(1), nbunmix     
!     Basic calculation of the potential temperature   
      psl = ps(x, y ,1, ind)
      do k= lower_theta_level,upper_theta_level
        theta_g(k,x,y,ind) = (p0/(akz(k)+bkz(k)*psl))**kappa * tth(x,y,k,ind)
      enddo
!     Tolerance on the difference of theta between two successive levels
      tol = 0.001
      nbunmix=0
      theta_col(x,y,ind) = .true.
      if(debug_out) then
        write(*,*) 'dans calc_col_theta'
        write(*,*) 'x=',x
        write(*,*) 'y=',y
        write(*,*) 'k=',k
        write(*,*) 'ind=',ind
        write(*,*) 'theta_g(lower_theta_level:upper_theta_level,x,y,ind)=',theta_g(lower_theta_level:upper_theta_level,x,y,ind)
      endif
      allocate (theta_diff(1:upper_theta_level-lower_theta_level))
      theta_diff = theta_g(lower_theta_level+1:upper_theta_level,x,y,ind) &
                 - theta_g(lower_theta_level:upper_theta_level-1,x,y,ind)
      if(minval(theta_diff) < 0.) then
        theta_inv_col(x,y,ind) = .true.
        call quicksort(theta_g(lower_theta_level,x,y,ind),0, &
          upper_theta_level-lower_theta_level)
          ! quicksort written in C, for which indices start at 0, not 1
        theta_diff = theta_g(lower_theta_level+1:upper_theta_level,x,y,ind) &
                   - theta_g(lower_theta_level:upper_theta_level-1,x,y,ind)
      endif
      ! Recalculation of the profile in case of mixing with equal theta over 
      ! two or more successive levels
      ! The algorithm is borrowed from unmix with the important difference that
      ! it does not aim at smoothing the profile but only at removing spurious
      ! singularities.
      ! Hence min are replaced by max and conversely, with respect to unmix.
      ! Interpolation is here only used to ensure monotonicity. (to be checked)
      do while (minval(theta_diff) < tol)
        nbunmix=nbunmix+1
        id1=minloc(theta_diff(:))
        id=id1(1)+1
        if(id==2) then
          theta_g(lower_theta_level,x,y,ind)= &
                  theta_g(lower_theta_level+1,x,y,ind)-3.2*tol
        else
          theta_g(lower_theta_level+id-1,x,y,ind)= &
              min((2*theta_g(lower_theta_level+id-2,x,y,ind) &
                    +theta_g(lower_theta_level+id,x,y,ind))/3., &
                  theta_g(lower_theta_level+id-1,x,y,ind)-2*tol)
!                  theta_g(lower_theta_level+id-1,x,y,ind)-1.1*tol
!              min((2*theta_g(lower_theta_level+id-3,x,y,ind) &
!                    +theta_g(lower_theta_level+id,x,y,ind))/3., &
!                  theta_g(lower_theta_level+id-1,x,y,ind)-2*tol)
          theta_g(lower_theta_level+id,x,y,ind)= &
              max((theta_g(lower_theta_level+id,x,y,ind) &
                +2*theta_g(lower_theta_level+id+1,x,y,ind))/3., &
                  theta_g(lower_theta_level+id,x,y,ind)+2*tol) 
!                  theta_g(lower_theta_level+id,x,y,ind)+1.1*tol
!              max((theta_g(lower_theta_level+id-3,x,y,ind) &
!                +2*theta_g(lower_theta_level+id,x,y,ind))/3., &
!                  theta_g(lower_theta_level+id-2,x,y,ind)+2*tol) 
          call quicksort(theta_g(lower_theta_level,x,y,ind),0, id)
        endif
        theta_diff = theta_g(lower_theta_level+1:upper_theta_level,x,y,ind) &
                   - theta_g(lower_theta_level:upper_theta_level-1,x,y,ind)
        !print *,id,x,y,ind,nbunmix
      enddo
      deallocate(theta_diff)
      return
      end subroutine calc_col_theta
!*******************************************************************************          
      function locisent(theta,x,y,ib1,iu1,ind)
      integer, intent(in) :: ib1,iu1,x,y,ind
      real, intent(in) :: theta
      integer :: locisent,ib,iu,im
      ib=ib1 ; iu=iu1
      do while(iu-ib>1)
        im = (iu+ib)/2
        if( theta >= theta_g(im,x,y,ind) ) then
          ib=im
        else
          iu=im
        endif
      enddo
      locisent = ib
      end function locisent
!*******************************************************************************      
      function locisent2(theta,x,y,ib1,ind)
      integer, intent(in) :: x,y,ib1,ind
      real, intent(in) :: theta
      integer :: locisent2,ib
      ib=ib1
      if (theta < theta_g(ib,x,y,ind)) then
        do while(ib>lower_theta_level)
          ib = ib-1
          if ( theta >= theta_g(ib,x,y,ind) ) exit
        enddo
      else
        do while(ib<upper_theta_level-1)
          if (theta < theta_g(ib+1,x,y,ind)) exit
          ib = ib+1
        enddo
      endif
      locisent2 = ib
      end function locisent2 

end module isentrop_m

!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log:

