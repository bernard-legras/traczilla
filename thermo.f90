
!**********************************************************************
! Copyright 2004, 2006, 2007, 2012, 2013                              *
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
!
!###############################################################################
!----------------------------------- THERMO ------------------------------------
!###############################################################################

module thermo
implicit none

contains
	
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SATRATIO @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

real function satratio(press,temp)

!*******************************************************************************
!                                                                              *
!     Calculates the saturation mixing ratio of water vapor over ice
!     Using formula provided in K. Emanuel's book
!     Calibrated for temperatures in the TTL
!     Result is in mass mixing ratio (kg/kg)
!                                                                              *
!     Author: B. Legras
!                                                                              *
!     31 Jan 2004
!     23 Aug 2007: correction: 
!                  press assumed to be in Pa is divided by 100
!                  to fit in the formula in hPa
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:     
!
!
!*******************************************************************************

 real, intent(in):: press, temp
 real :: estar, satppmv
 ! saturation pressure
 estar=1.0008*exp(23.33086-(6111.72784/(temp))+0.15215*alog(temp))
 ! saturation in volume mixing ratio (ppmv)
 satppmv=1e6*estar/(0.01*press-estar)
 ! saturation in mass mixing ratio (kg/kg)
 satratio=0.622e-6 * satppmv    
 return
 end function satratio
!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log: 
!$Id:
!
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ F_QVSAT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################

      REAL FUNCTION F_QVSAT( p, t ) 

!     PURPOSE:
!
!     Calculate the saturation specific humidity using enhanced Teten's
!     formula.
!
!     AUTHOR: Yuhe Liu
!     01/08/1998
!
!     MODIFICATION HISTORY:
!
!     INPUT :
!       p        Pressure (Pascal)
!       t        Temperature (K)
!     OUTPUT:
!       f_qvsat  Saturation water vapor specific humidity (kg/kg).
!
!     Variable Declarations.
!

      real,intent(in) :: p,t         ! Pressure (Pascal)
                                     ! Temperature (K)
      real :: fespt
      
      real :: rd        ! Gas constant for dry air  (m**2/(s**2*K))
      parameter( rd     = 287.0 )
      real :: rv        ! Gas constant for water vapor  (m**2/(s**2*K)).
      parameter( rv     = 461.0 )
      real :: rddrv    
      parameter( rddrv  = rd/rv )


! Change by A. Stohl to save computation time:
      fespt=f_es(p,t)
      f_qvsat = rddrv * fespt / (p-(1.0-rddrv)*fespt)

      RETURN
      END FUNCTION F_QVSAT
      
      REAL FUNCTION F_ES( p, t ) 
!#######################################################################
!
!     PURPOSE:
!
!     Calculate the saturation specific humidity using enhanced Teten's
!     formula.
!
!     AUTHOR: Yuhe Liu
!     01/08/1998
!
!     MODIFICATION HISTORY:
!
!#######################################################################
!     INPUT :
!       p        Pressure (Pascal)
!       t        Temperature (K)
!     OUTPUT:
!       f_es     Saturation water vapor pressure (Pa)
!#######################################################################
!     Variable Declarations.

      real, intent(in) :: p,t         ! Pressure (Pascal)
                                      ! Temperature (K)
      real :: f_esl, f_esi

!      IF ( t.ge.273.15 ) THEN      ! for water
      IF ( t.ge.253.15 ) THEN      ! modification Petra Seibert
                                   ! (supercooled water may be present)
        f_es = f_esl( p,t )
      ELSE                            ! for ice
        f_es = f_esi( p,t )
      ENDIF

      RETURN
      END FUNCTION F_ES

      REAL FUNCTION F_ESL( p, t )  

      real, intent(in) :: p,t         ! Pressure (Pascal)
                                      ! Temperature (K)
      real :: f

!#######################################################################
!
!     Saturation specific humidity parameters used in enhanced Teten's
!     formula. (See A. Buck, JAM 1981)
!     Saturation water vapor pressure over liquid water
!
!#######################################################################

      real :: satfwa, satfwb
      parameter ( satfwa = 1.0007 )
      parameter ( satfwb = 3.46e-8 )  ! for p in Pa

      real :: satewa, satewb, satewc
      parameter ( satewa = 611.21 )   ! es in Pa
      parameter ( satewb = 17.502 )
      parameter ( satewc = 32.18 )

      f = satfwa + satfwb * p
      f_esl = f * satewa * exp( satewb*(t-273.15)/(t-satewc) )

      RETURN
      END  FUNCTION F_ESL

      REAL FUNCTION F_ESI( p, t )  

      real, intent(in) :: p,t         ! Pressure (Pascal)
                                      ! Temperature (K)
      real :: f

!#######################################################################
!
!     Saturation specific humidity parameters used in enhanced Teten's
!     formula. (See A. Buck, JAM 1981)
!     Saturation water vapor pressure over ice (Pa)
!
!#######################################################################
!
      real :: satfia, satfib
      parameter ( satfia = 1.0003 )
      parameter ( satfib = 4.18e-8 )  ! for p in Pa

      real :: sateia, sateib, sateic
      parameter ( sateia = 611.15 )   ! es in Pa
      parameter ( sateib = 22.452 )
      parameter ( sateic = 0.6 )

      f = satfia + satfib * p
      f_esi = f * sateia * exp( sateib*(t-273.15)/(t-sateic) )

      RETURN
      END FUNCTION F_ESI

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ EW @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

      REAL FUNCTION EW(X)
!     ****************************************************************
!     SAETTIGUNGSDAMPFDRUCK UEBER WASSER IN PA. X IN KELVIN.
!     NACH DER GOFF-GRATCH-FORMEL.
!     ****************************************************************
      REAL, INTENT(IN)::X
      REAL :: Y,A,C,D
      EW=0.
      IF(X.LE.0.) STOP 'SORRY: T NOT IN [K]'
      Y=373.16/X
      A=-7.90298*(Y-1.)
      A=A+(5.02808*0.43429*ALOG(Y))
      C=(1.-(1./Y))*11.344
      C=-1.+(10.**C)
      C=-1.3816*C/(10.**7)
      D=(1.-Y)*3.49149
      D=-1.+(10.**D)
      D=8.1328*D/(10.**3)
      Y=A+C+D
      EW=101324.6*(10.**Y)       ! Saettigungsdampfdruck in Pa
      RETURN
      END FUNCTION EW

 end module thermo
!
!=====|==1=========2=========3=========4=========5=========6=========7==
!
!$Log: 
