!**********************************************************************
! Copyright 2004                                                      *
! Ignacio Pisso, Bernard Legras                                       *
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
!=====================================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@[ TRACZILLA ]@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8======
!
module lyapunov

!**************************************
! Calculation of Lyapunov exponent
!
! B. Legras & I. Pisso, April 2004 
! 
!**************************************

use commons
implicit none
save
     
      logical activ_Lyapunov
      integer numpart_np, tau_orthog
      real sheet_slope, delta_hor, delta_ver, height_scale
      real, dimension(:,:), allocatable   ::  lambdaa, lambdab
      real, dimension(:,:,:), allocatable ::  Mp
      !test---
      real, dimension(:), allocatable     ::  mass_triad_b
      !-------       
      
 ! activ_Lyapunov          logical switch to activate Lyapunov calculations     
 ! numpart_np              number of center particles for Lyapunov runs 
 !                         (should be nsample*n_loc)
 ! sheet_slope	           assumed slope of tracer sheets used to define norms
 ! delta_hor		   horizontal displacement 
 ! delta_ver		   vertical displacement (= delta_hor / sheet_slope)
 ! height_scale		   height scale
 ! tau_orthog		   interval between 2 orthogonalizations
 
 contains
 
subroutine alloc_Lyapunov
   allocate (lambdaa(3,numpart_np), lambdab(3,numpart_np))
   lambdaa=0. ; lambdab=0.
   allocate (Mp(3,3,numpart_np))
   allocate (mass_triad_b(numpart_np))
end subroutine alloc_Lyapunov
 
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine ortho(itime) 

!*******************************************************************************
!   	 Fixes perturbations used in the calculation of lyapunov charateristic *
!	 exponents							       *     
!	 Ignacio PISSO April 2004                                              *
! 	 						                       *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

 integer i,j,j1
 integer, intent (in) :: itime
 real, dimension(7,4) ::  p
 !perturbations with respect to the particle before and after gram
 real, dimension(3,4) ::  h, b 
 real, dimension(3,4) ::  c
 real, dimension(3) ::  lambdaa_loc,lambdab_loc
 !test-----
 real :: mass_triad_f1, mass_triad_f2
 !---------
 
 j1 = numpart_np 
 ! main loop
 do i=1,numpart_np 
 
    call grid2cart(p(1,1),p(1,2),p(1,3),xtra1(i),ytra1(i))
    p(1,4) = ztra1(i)
    do j=1,6  
       call grid2cart(p(j+1,1),p(j+1,2),p(j+1,3),xtra1(j1+j),ytra1(j1+j))
       p(j+1,4)= ztra1(j1+j)
    enddo
     
    do j=1,3
       c(j,:)=p(1,:) ! expand central particle
    enddo

!1st triad
!*********
 
    h = p(2:4,:) - c
    call gramlam(h,b,lambdaa_loc)
    p(2:4,:) = c + b
    lambdaa(:,i) = lambdaa(:,i) + lambdaa_loc
    
    !---test
    mass_triad_f1 = exp(-ztra1(i)+lambdaa_loc(1)+lambdaa_loc(2)+lambdaa_loc(3))
    !-------    
    
!2nd triad
!*********
  
    h = p(5:7,:) - c
    call gramlam(h,b,lambdab_loc)
    p(5:7,:) = c + b
    lambdab(:,i) = lambdab(:,i) + lambdab_loc

    !---test
    mass_triad_f2 = exp(-ztra1(i)+lambdab_loc(1)+lambdab_loc(2)+lambdab_loc(3))
    write(*,'(a,i10,i5,2G14.5)') ' ortho> ', itime, i, &
       mass_triad_f1/mass_triad_b(i),mass_triad_f2/mass_triad_b(i)
    !-------

!New perturbations in grid coordinates
!*************************************

    do j=1,6
       call sphere(p(j+1,:))
       call cart2grid(xtra1(j1+j),ytra1(j1+j),p(j+1,1),p(j+1,2),p(j+1,3))
       ztra1(j1+j) = p(j+1,4)
    enddo
   
    j1=j1+6 
 enddo  !end main loop

 return
 end subroutine ortho
 
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

subroutine gramlam(hh,b,lambdagram) 

!*************************************************************************************
!   Orthogonalize perturbations used in the calculation of Lyapunov characteristic
!        exponent
!
!	 Ignacio PISSO April 2004
!
!*************************************************************************************
!
! Variables:
!
!*************************************************************************************

real, dimension(3,4) , intent (in)      :: hh 

real                                    :: area1, vol
real, dimension(4)                      :: h3p , h2p 

real, dimension(3,4) , intent (out)     :: b
real, dimension(3)   , intent (out)     :: lambdagram

!*************************************************************************************

b(1,:) = hh(1,:)/norm4(hh(1,:))
lambdagram(1) = log(norm4(hh(1,:)))                    
h2p = hh(2,:) - prod4(b(1,:),hh(2,:))*b(1,:) 
b(2,:) = h2p/norm4(h2p)
area1 = cross4(hh(1,:),hh(2,:))
lambdagram(2)= log(area1) - lambdagram(1)            
h3p =  hh(3,:) - prod4(b(1,:),hh(3,:))*b(1,:) - prod4(b(2,:),hh(3,:))*b(2,:)
b(3,:) = h3p/norm4(h3p)
vol = abs(area1*norm4(h3p))
lambdagram(3) =  log(vol) - lambdagram(2) - lambdagram(1)

!test-----
!write(*,'(a,3g12.5)') ' gramlam> ',lambdagram(1),lambdagram(2),lambdagram(3)
!write(*,'(a,3g12.5)') ' gramlam> ',norm4(b(1,:)),norm4(b(2,:)),norm4(b(3,:))
!write(*,'(a,3g12.5)') ' gramlam> ',prod4(b(1,:),b(2,:)),prod4(b(2,:),b(3,:)), &
!                        prod4(b(3,:),b(1,:))
write(*,'(a,3g12.5)') ' gramlam> ', &
   180/pi * asin(prod4(hh(1,:),hh(2,:))/(norm4(hh(1,:))*norm4(hh(2,:)))), &
   180/pi * asin(prod4(hh(2,:),hh(3,:))/(norm4(hh(2,:))*norm4(hh(3,:)))), &
   180/pi * asin(prod4(hh(3,:),hh(1,:))/(norm4(hh(3,:))*norm4(hh(1,:))))
!---------

return
end subroutine gramlam
 
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
 
subroutine fix_perturb


!********************************************************************************
!   	 Fixes perturbations used in the calculation of lyapunov charateristic 	*
! 	 exponents								*
!									      	*
!	 Ignacio PISSO April 2004				      	       	*
!********************************************************************************
!                                                                              	*
! Variables:                                                                   	*
!                                                                              	*
!********************************************************************************
 
integer i,j, j1 

real, dimension(7,4) ::  p , b
real, dimension(3)   :: v,vx,vy,v1,v2
real deltaz

!main boucle 
!************************************************************************************* 
j1 = numpart_np 
do i=1,numpart_np
 	                          
   do j=1,7  
        call grid2cart(p(j,1),p(j,2),p(j,3),xtra1(i),ytra1(i))
        p(j,4) = ztra1(i)
        call sphere(p(j,:))
   enddo
	    
! Initialisation 2D base on the surface of the sphere 
!*********************************************************************************

    v  = p(1,1:3)
    vx = (/ 1, 0, 0/)
    vy = (/ 0, 1, 0/)
  
    v1 =  vx - v*prod3(v,vx)
    v1 =  v1/norm3(v1)
    v2 =  vy - v*prod3(v,vy) - v1*prod3(vy,v1) 
    v2 =  v2/norm3(v2)
    				
!Perturbations
!*********************************************************************************
		
    v1 = v1*(delta_hor/r_earth)
    v2 = v2*(delta_hor/r_earth)
    deltaz = delta_ver/height_scale
								
    b        = 0.
    b(2,1:3) = v1
    b(3,1:3) = v2
    b(4,4)   = deltaz
    b(5,1:3) = -v1
    b(6,1:3) = -v2
    b(7,4)   = -deltaz
    
    p = p + b
    
    do j=1,6
      call sphere(p(j+1,:))
      call cart2grid(xtra1(j1+j),ytra1(j1+j),p(j+1,1),p(j+1,2),p(j+1,3))
      ztra1(j1+j) = p(j+1,4)
      itra1(j1+j) = itra1(i)
      itramem(j1+j) = itra1(i)
    enddo
    j1=j1+6

enddo    !end main boucle

print * , "Perturbations initialized "

return
end subroutine fix_perturb

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

!***********************************************************************************
!Scalar product in 3 D, within the unit sphere world
!***********************************************************************************   
 
 real function prod3(u,v)
 real, dimension(3) :: u
 real, dimension(3) :: v
 prod3 = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
 return
 end function prod3
 
!************************************************************************************
!Norm in 3 D, within the unit sphere world
!************************************************************************************
 real function norm3(v)
 real, dimension(3) :: v
 norm3 =  sqrt((v(1))**2 + (v(2))**2 + (v(3))**2 ) 
 return
 end function norm3
 
!*************************************************************************************
!*************************************************************************************
 
 real function prod4(u,v)
 real, dimension(4) :: u
 real, dimension(4) :: v
 prod4 = ((r_earth/delta_hor)**2)*(u(1)*v(1)+u(2)*v(2)+u(3)*v(3)) &
        +((height_scale/delta_ver)**2)*u(4)*v(4)
 return
 end function prod4
  
!*************************************************************************************
!Norm in 4 D
!*************************************************************************************
 real function norm4(v)
 real, dimension(4) :: v
 norm4 =  sqrt(((r_earth/delta_hor)**2)*(v(1)**2 + v(2)**2 + v(3)**2) &
             + ((height_scale/delta_ver)*v(4))**2) 
 return
 end function norm4

!*************************************************************************************
!*************************************************************************************
 subroutine sphere(v)
 real, dimension(4) :: v
 real arrow
 arrow=1/sqrt(v(1)**2 + v(2)**2 + v(3)**2)
 v(1) = v(1)*arrow
 v(2) = v(2)*arrow
 v(3) = v(3)*arrow
 end subroutine sphere
 
!*************************************************************************************
!Cross product
!*************************************************************************************
 real function cross4(v,w)
 real, dimension(4) :: v,w
 real prod, norm_v2, norm_w2
 norm_v2 = ((r_earth/delta_hor)**2)*(v(1)**2 + v(2)**2 + v(3)**2) &
         + ((height_scale/delta_ver)*v(4))**2
 norm_w2 = ((r_earth/delta_hor)**2)*(w(1)**2 + w(2)**2 + w(3)**2) &
         + ((height_scale/delta_ver)*w(4))**2
 prod = ((r_earth/delta_hor)**2)*(w(1)*v(1)+w(2)*v(2)+w(3)*v(3)) &
      + ((height_scale/delta_ver)**2)*w(4)*v(4)
 cross4 = sqrt(norm_v2 * norm_w2 - prod**2)
 return
 end function cross4

!*************************************************************************************
!************************************************************************************* 
 subroutine grid2cart(x,y,z,xtra,ytra)
 
 real, intent(in) :: xtra,ytra
 real, intent(out) :: x,y,z
 real xlon,ylat
 
 xlon = (xtra*dx + xlon0)*pi/180
 ylat = (ytra*dx + ylat0)*pi/180
 x = cos(ylat)*cos(xlon)
 y = cos(ylat)*sin(xlon) 
 z = sin(ylat)
 return
 end subroutine grid2cart

!*************************************************************************************
!*************************************************************************************
 subroutine cart2grid(xtra,ytra,x,y,z)
 
 real, intent(in):: x,y,z
 real, intent(out):: xtra,ytra
 real xlon,ylat
 
 ylat = (asin(z))*180/pi
 xlon = (atan2(y,x))*180/pi
 
 ytra = (ylat - ylat0)/dy
 xtra = (xlon - xlon0)/dx
 return
 end subroutine cart2grid

 end module lyapunov
