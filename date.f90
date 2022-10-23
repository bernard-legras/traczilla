!**********************************************************************
! Copyright 1993, 1994, 1997, 2002,  2013                             *
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
!###############################################################################
!------------------------------------ DATE -------------------------------------
!###############################################################################

module date
implicit none

contains

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CALDATE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
!
      SUBROUTINE caldate(JULDATE,YYYYMMDD,HHMISS)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!                                                                             *
!     Calculates the Gregorian date from the Julian date                      *
!                                                                             *
!     AUTHOR: Andreas Stohl (21 January 1994), adapted from Numerical Recipes *
!                                                                             *
!     Variables:                                                              *
!     DD             Day                                                      *
!     HH             Hour                                                     *
!     HHMISS         Hour, Minute, Second                                     *
!     JA,JB,JC,JD,JE help variables                                           *
!     JALPHA         help variable                                            *
!     JULDATE        Julian Date                                              *
!     JULDAY         help variable                                            *
!     MI             Minute                                                   *
!     MM             Month                                                    *
!     SS             Seconds                                                  *
!     YYYY           Year                                                     *
!     YYYYMMDD       Year, Month, Day                                         *
!                                                                             *
!     Constants:                                                              *
!     IGREG          help constant                                            *
!                                                                             *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      INTEGER YYYYMMDD,YYYY,MM,DD,HHMISS,HH,MI,SS
      INTEGER JULDAY,JA,JB,JC,JD,JE,IGREG,JALPHA
      REAL (kind=8):: JULDATE
      PARAMETER (IGREG=2299161)

      JULDAY=INT(JULDATE)
      IF(JULDAY.GE.IGREG)THEN
        JALPHA=INT(((JULDAY-1867216)-0.25)/36524.25)
        JA=JULDAY+1+JALPHA-INT(0.25*JALPHA)
      ELSE
        JA=JULDAY
      ENDIF
      JB=JA+1524
      JC=INT(6680.+((JB-2439870)-122.1)/365.25)
      JD=365*JC+INT(0.25*JC)
      JE=INT((JB-JD)/30.6001)
      DD=JB-JD-INT(30.6001*JE)
      MM=JE-1
      IF (MM.GT.12) MM=MM-12
      YYYY=JC-4715
      IF (MM.GT.2) YYYY=YYYY-1
      IF (YYYY.LE.0) YYYY=YYYY-1

      YYYYMMDD=10000*YYYY+100*MM+DD
      HH=INT(24.*(JULDATE-FLOAT(JULDAY)))
      MI=INT(1440.*(JULDATE-FLOAT(JULDAY))-60.*FLOAT(HH))
      SS=NINT(86400.*(JULDATE-FLOAT(JULDAY))-3600.*FLOAT(HH)) &
         -60.*FLOAT(MI)
      IF (SS.EQ.60) THEN  ! 60 seconds = 1 minute
        SS=0
        MI=MI+1
      ENDIF
      IF (MI.EQ.60) THEN
        MI=0
        HH=HH+1
      ENDIF
      HHMISS=10000*HH+100*MI+SS

      RETURN
      END SUBROUTINE caldate

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ JULDATE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
!
      FUNCTION juldate(YYYYMMDD,HHMISS)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!                                                                             *
!     Calculates the Julian date                                              *
!                                                                             *
!     AUTHOR: Andreas Stohl (15 October 1993)                                 *
!                                                                             *
!     Variables:                                                              *
!     DD             Day                                                      *
!     HH             Hour                                                     *
!     HHMISS         Hour, minute + second                                    *
!     JA,JM,JY       help variables                                           *
!     JULDATE        Julian Date                                              *
!     JULDAY         help variable                                            *
!     MI             Minute                                                   *
!     MM             Month                                                    *
!     SS             Second                                                   *
!     YYYY           Year                                                     *
!     YYYYMMDDHH     Date and Time                                            *
!                                                                             *
!     Constants:                                                              *
!     IGREG          help constant                                            *
!                                                                             *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      INTEGER YYYYMMDD,YYYY,MM,DD,HH,MI,SS,HHMISS
      INTEGER JULDAY,JY,JM,JA,IGREG
      REAL (kind=8):: JULDATE
      PARAMETER (IGREG=15+31*(10+12*1582))

      YYYY=YYYYMMDD/10000
      MM=(YYYYMMDD-10000*YYYY)/100
      DD=YYYYMMDD-10000*YYYY-100*MM
      HH=HHMISS/10000
      MI=(HHMISS-10000*HH)/100
      SS=HHMISS-10000*HH-100*MI

      IF (YYYY.EQ.0) STOP 'juldate> There is no Year Zero.'
      IF (YYYY.LT.0) YYYY=YYYY+1
      IF (MM.GT.2) THEN
        JY=YYYY
        JM=MM+1
      ELSE
        JY=YYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+DD+1720995
      IF (DD+31*(MM+12*YYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
       
      JULDATE=DBLE(FLOAT(JULDAY))+DBLE(FLOAT(HH)/24.)+ &
         DBLE(FLOAT(MI)/1440.)+DBLE(FLOAT(SS)/86400.)

      RETURN
      END FUNCTION juldate
      
end module date
!
!=====|==1=========2=========3=========4=========5=========6=========7=========8
!
!$Log: 
