
!**********************************************************************
! Copyright 2001, 2002, 2006, 2007, 2012, 2013           *
! Bernard Legras                                       *
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
module isentrop_h

!**************************************
! Isentropic interpolation modulea (head)
! 
! B. Legras, June 2005 
! 
!**************************************

use commons
implicit none
save 

     logical diabatic_w, isentropic_motion, theta_bounds

! diabatic_w              vertical velocity is from heating rates
! isentropic_motion       use theta as vertical coordinate
!                         and interpolate wind to theta level
! theta_bounds            bounds theta variations within local limits
     
      integer lower_theta_level, upper_theta_level
      real, dimension(:), allocatable :: theta_part
      real, dimension(:,:,:,:), allocatable:: theta_g
      logical, dimension(:,:,:), allocatable:: theta_col,theta_inv_col
      
 ! theta_part          theta for each particle (unused)
 ! theta_g             theta for the ECMWF data, on the grid
 ! lower_theta_level   lower ecmwf hybrid level used in isentropic
 !                     interpolation
 ! upper_theta_level   upper ecmwf hybrid level used in isentropic
 !                     interpolation
 ! theta_col           true when theta_g has been calculated for a
 !                     given point in the horizontal grid
 ! theta_inv_col       true when an inversion has been detected in the
 !                     theta column for the corresponding grid point
 !                     (diagnostic purpose)
 

contains
 
subroutine alloc_isentrop_perm
print *,'alloc_isentrop_perm theta_g theta_col theta_inv_col'
   if(upper_theta_level > nuvz) then
     print *,'WARNING: upper_theta_level is reduced to nuvz'
     upper_theta_level = nuvz
! previously cutoff at nuvz-1
   endif
   allocate (theta_g(lower_theta_level:upper_theta_level,0:nx-1,0:ny-1,2))
   allocate (theta_col(0:nx-1,0:ny-1,2),theta_inv_col(0:nx-1,0:ny-1,2))
!   allocate  of theta_part delegated to initialization routine
end subroutine alloc_isentrop_perm

end module isentrop_h
 
!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log:
