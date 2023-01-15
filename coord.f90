!$Id: 
!
!###############################################################################
!----------------------------------- COORD -------------------------------------
!###############################################################################

! CHANGES TO THE ROUTINES BY A. STOHL
! XI,XI0,ETA,ETA0 ARE DOUBLE PRECISION VARIABLES TO AVOID PROBLEMS
! AT POLES
! Modified by B. Legras to eliminate weird WMO behaviour
! within the last degree near the poles
! Unused routines eliminated
! 3-1-2023: Passed integrally in real*8 except the CC2GLL wind arguments which are used
! to generate uupol and vvpol from uuh and vvh 

module coord
implicit none
integer, parameter:: dbl=kind(0.d0),dp=kind(0.0),sp=kind(0.0)
REAL(dbl), PARAMETER :: PI=3.14159265358979D0,RADPDG=PI/180,DGPRAD=180/PI
REAL(dbl), PARAMETER :: REARTH=6367.47D0
REAL(dbl), PARAMETER :: ALMST1=.9999999D0
contains

!*  GENERAL CONFORMAL MAP ROUTINES FOR METEOROLOGICAL MODELERS
!*  WRITTEN ON 3/31/94 BY

!* Dr. Albion Taylor
!* NOAA / OAR / ARL                  Phone: (301) 713-0295 x 132
!* Rm. 3151, 1315 East-West Highway  Fax:   (301) 713-0119
!* Silver Spring, MD 20910           E-mail: ADTaylor@arlrisc.ssmc.noaa.gov 

!*  SUBROUTINE STLMBR (STRCMP, TNGLAT, CLONG)
!*    THIS ROUTINE INITIALIZES THE MAP STRUCTURE ARRAY STRCMP TO
!*    THE FORM OF A SPECIFIC MAP PROJECTION
!*  INPUTS:
!*    TNGLAT - THE LATITUDE AT WHICH THE PROJECTION WILL BE TANGENT
!*             TO THE EARTH.  +90. FOR NORTH POLAR STEREOGRAPHIC,
!*             -90. FOR SOUTH POLAR STEREOGRAPHIC, 0. FOR MERCATOR,
!*             AND OTHER VALUES FOR LAMBERT CONFORMAL. 
!*             -90 <= TNGLAT <= 90.
!*    CLONG -  A LONGITUDE IN THE REGION UNDER CONSIDERATION.  LONGITUDES
!*             BETWEEN CLONG-180. AND CLONG+180.  WILL BE MAPPED IN ONE
!*             CONNECTED REGION
!*  OUTPUTS:
!*    STRCMP - A 9-VALUE MAP STRUCTURE ARRAY FOR USE WITH SUBSEQUENT
!*             CALLS TO THE COORDINATE TRANSFORM ROUTINES.
!*
!*  REAL FUNCTION EQVLAT (XLAT1,XLAT2)
!*    THIS FUNCTION IS PROVIDED TO ASSIST IN FINDING THE TANGENT LATITUDE
!*    EQUIVALENT TO THE 2-REFERENCE LATITUDE SPECIFICATION IN THE LEGEND
!*    OF MOST LAMBERT CONFORMAL MAPS.  IF THE MAP SPECIFIES "SCALE 
!*    1:XXXXX TRUE AT 40N AND 60N", THEN EQVLAT(40.,60.) WILL RETURN THE
!*    EQUIVALENT TANGENT LATITUDE.
!*  INPUTS:
!*    XLAT1,XLAT2:  THE TWO LATITUDES SPECIFIED IN THE MAP LEGEND
!*  RETURNS:
!*    THE EQUIVALENT TANGENT LATITUDE
!*  EXAMPLE:  CALL STLMBR(STRCMP, EQVLAT(40.,60.), 90.)

!*  SUBROUTINE STCM2P (STRCMP, X1,Y1, XLAT1,XLONG1,
!*          X2,Y2, XLAT2,XLONG2)
!*  SUBROUTINE STCM1P (STRCMP, X1,Y1, XLAT1,XLONG1,
!*          XLATG,XLONGG, GRIDSZ, ORIENT)
!*    THESE ROUTINES COMPLETE THE SPECIFICATION OF THE MAP STRUCTURE
!*    ARRAY BY CONFORMING THE MAP COORDINATES TO THE SPECIFICATIONS
!*    OF A PARTICULAR GRID.  EITHER STCM1P OR STCM2P MUST BE CALLED,
!*    BUT NOT BOTH
!*  INPUTS:
!*    STRCMP - A 9-VALUE MAP STRUCTURE ARRAY, SET TO A PARTICULAR MAP
!*             FORM BY A PREVIOUS CALL TO STLMBR
!*    FOR STCM2P:
!*      X1,Y1, X2,Y2 - THE MAP COORDINATES OF TWO POINTS ON THE GRID
!*      XLAT1,XLONG1, XLAT2,XLONG2 - THE GEOGRAPHIC COORDINATES OF THE
!*             SAME TWO POINTS
!*    FOR STCM1P:
!*      X1,Y1 - THE MAP COORDINATES OF ONE POINT ON THE GRID
!*      XLAT1,XLONG1 - THE GEOGRAPHIC COORDINATES OF THE SAME POINT
!*      XLATG,XLONGG - LATITUDE AND LONGITUDE OF REFERENCE POINT FOR
!*             GRIDSZ AND ORIENTATION SPECIFICATION.
!*      GRIDSZ - THE DESIRED GRID SIZE IN KILOMETERS, AT XLATG,XLONGG
!*      ORIENT - THE ANGLE, WITH RESPECT TO NORTH, OF A Y-GRID LINE, AT
!*             THE POINT XLATG,XLONGG
!*  OUTPUTS:
!*    STRCMP - A 9-VALUE MAP STRUCTURE ARRAY, FULLY SET FOR USE BY
!*             OTHER SUBROUTINES IN THIS SYSTEM

!*  SUBROUTINE CLL2XY (STRCMP, XLAT,XLONG, X,Y)
!*  SUBROUTINE CXY2LL (STRCMP, X,Y, XLAT,XLONG)
!*     THESE ROUTINES CONVERT BETWEEN MAP COORDINATES X,Y
!*     AND GEOGRAPHIC COORDINATES XLAT,XLONG
!*  INPUTS:
!*     STRCMP(9) - 9-VALUE MAP STRUCTURE ARRAY
!*     FOR CLL2XY:  XLAT,XLONG - GEOGRAPHIC COORDINATES
!*     FOR CXY2LL:  X,Y - MAP COORDINATES
!*  OUTPUTS:
!*     FOR CLL2XY:  X,Y - MAP COORDINATES
!*     FOR CXY2LL:  XLAT,XLONG - GEOGRAPHIC COORDINATES

!*  SUBROUTINE CC2GXY (STRCMP, X,Y, UE,VN, UG,VG)
!*  SUBROUTINE CG2CXY (STRCMP, X,Y, UG,VG, UE,VN)
!*  SUBROUTINE CC2GLL (STRCMP, XLAT,XLONG, UE,VN, UG,VG)
!*  SUBROUTINE CG2CLL (STRCMP, XLAT,XLONG, UG,VG, UE,VN)
!*     THESE SUBROUTINES CONVERT VECTOR WIND COMPONENTS FROM
!*     GEOGRAPHIC, OR COMPASS, COORDINATES, TO MAP OR
!*     GRID COORDINATES.  THE SITE OF THE WIND TO BE
!*     CONVERTED MAY BE GIVEN EITHER IN GEOGRAPHIC OR
!*     MAP COORDINATES.  WIND COMPONENTS ARE ALL IN KILOMETERS
!*     PER HOUR, WHETHER GEOGRAPHIC OR MAP COORDINATES.
!*  INPUTS:
!*    STRCMP(9) - 9-VALUE MAP STRUCTURE ARRAY
!*    FOR CC2GXY AND CG2CXY:  X,Y        -  MAP COORDINATES OF SITE
!*    FOR CC2GLL AND CG2CLL:  XLAT,XLONG -  GEOGRAPHIC COORDINATES OF SITE
!*    FOR CC2GXY AND CC2GLL:  UE,VN - EAST AND NORTH WIND COMPONENTS
!*    FOR CG2CXY AND CG2CLL:  UG,VG - X- AND Y- DIRECTION WIND COMPONENTS
!*  OUTPUTS:
!*    FOR CC2GXY AND CC2GLL:  UG,VG - X- AND Y- DIRECTION WIND COMPONENTS
!*    FOR CG2CXY AND CG2CLL:  UE,VN - EAST AND NORTH WIND COMPONENTS

!*  SUBROUTINE CCRVXY (STRCMP, X, Y,       GX,GY)
!*  SUBROUTINE CCRVLL (STRCMP, XLAT,XLONG, GX,GY)
!*    THESE SUBROUTINES RETURN THE CURVATURE VECTOR (GX,GY), AS REFERENCED
!*    TO MAP COORDINATES, INDUCED BY THE MAP TRANSFORMATION.  WHEN
!*    NON-LINEAR TERMS IN WIND SPEED ARE IMPORTANT, A "GEODESIC" FORCE
!*    SHOULD BE INCLUDED IN THE VECTOR FORM [ (U,U) G - (U,G) U ] WHERE THE
!*    INNER PRODUCT (U,G) IS DEFINED AS UX*GX + UY*GY.
!*  INPUTS:
!*    STRCMP(9) - 9-VALUE MAP STRUCTURE ARRAY
!*    FOR CCRVXY:  X,Y        -  MAP COORDINATES OF SITE
!*    FOR CCRVLL:  XLAT,XLONG -  GEOGRAPHIC COORDINATES OF SITE
!*  OUTPUTS:
!*    GX,GY       - VECTOR COEFFICIENTS OF CURVATURE, IN UNITS RADIANS
!*                  PER KILOMETER

!*  REAL FUNCTION CGSZLL (STRCMP, XLAT,XLONG)
!*  REAL FUNCTION CGSZXY (STRCMP, X,Y)
!*    THESE FUNCTIONS RETURN THE SIZE, IN KILOMETERS, OF EACH UNIT OF
!*    MOTION IN MAP COORDINATES (GRID SIZE).  THE GRID SIZE AT ANY
!*    LOCATION DEPENDS ON THAT LOCATION; THE POSITION MAY BE GIVEN IN
!*    EITHER MAP OR GEOGRAPHIC COORDINATES.
!*  INPUTS:
!*    STRCMP(9) - 9-VALUE MAP STRUCTURE ARRAY
!*    FOR CGSZXY:  X,Y        -  MAP COORDINATES OF SITE
!*    FOR CGSZLL:  XLAT,XLONG -  GEOGRAPHIC COORDINATES OF SITE
!*  RETURNS:
!*    GRIDSIZE IN KILOMETERS AT GIVEN SITE.

!*  SUBROUTINE CPOLXY (STRCMP, X,Y, ENX,ENY,ENZ)
!*  SUBROUTINE CPOLLL (STRCMP, XLAT,XLONG, ENX,ENY,ENZ)
!*    THESE SUBROUTINES PROVIDE 3-D VECTOR COMPONENTS OF A UNIT VECTOR
!*    IN THE DIRECTION OF THE NORTH POLAR AXIS.  WHEN MULTIPLIED
!*    BY TWICE THE ROTATION RATE OF THE EARTH (2 * PI/24 HR), THE
!*    VERTICAL COMPONENT YIELDS THE CORIOLIS FACTOR.
!*  INPUTS:
!*    STRCMP(9) - 9-VALUE MAP STRUCTURE ARRAY
!*    FOR CPOLXY:  X,Y        -  MAP COORDINATES OF SITE
!*    FOR CPOLLL:  XLAT,XLONG -  GEOGRAPHIC COORDINATES OF SITE
!*  RETURNS:
!*    ENX,ENY,ENZ THE DIRECTION COSINES OF A UNIT VECTOR IN THE
!*    DIRECTION OF THE ROTATION AXIS OF THE EARTH

!*  SUBROUTINE CNLLXY (STRCMP, XLAT,XLONG, XI,ETA)
!*  SUBROUTINE CNXYLL (STRCMP, XI,ETA, XLAT,XLONG)
!*    THESE SUBROUTINES PERFORM THE UNDERLYING TRANSFORMATIONS FROM
!*    GEOGRAPHIC COORDINATES TO AND FROM CANONICAL (EQUATOR CENTERED)
!*    COORDINATES.  THEY ARE CALLED BY CXY2LL AND CLL2XY, BUT ARE NOT
!*    INTENDED TO BE CALLED DIRECTLY

!*  REAL FUNCTION CSPANF (VALUE, BEGIN, END)
!*    THIS FUNCTION ASSISTS OTHER ROUTINES IN PROVIDING A LONGITUDE IN
!*    THE PROPER RANGE.  IT ADDS TO VALUE WHATEVER MULTIPLE OF 
!*    (END - BEGIN) IS NEEDED TO RETURN A NUMBER BEGIN < CSPANF <= END

       SUBROUTINE CC2GLL (STRCMP, XLAT,XLONG, UE,VN, UG,VG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
       REAL(dbl) XPOLG,YPOLG,ALONG,SLONG,CLONG,ROT
       ! As this is used to convert {uuh, vvh} into {vvpol, uupol}
       ! single precision is assumed
       ! To be changed when all the code is passed into double precision
       REAL(dp), INTENT(OUT) :: UG,VG
       REAL(dbl), INTENT(IN) :: STRCMP(9), XLAT, XLONG 
       REAL(dp), INTENT(IN) :: UE, VN

       ALONG = CSPANF( XLONG - STRCMP(2), -180.D0, 180.D0)
       !* The following paragraph assumes the pole line cotains two fixed values
       !* which are the pole wind given according to WMO rule
       ! if (xlat.gt.89.985D0) then
       !*  North polar meteorological orientation: "north" along prime meridian
       !  rot = - strcmp(1) * along + xlong - 180.D0
       !elseif (xlat.lt.-89.985D0) then
       !*  South polar meteorological orientation: "north" along prime meridian
       !  rot = - strcmp(1) * along - xlong
       !else
       !  rot = - strcmp(1) * along
       !endif
       !* When the pole line is a double sinusoid in quadrature in continuity with
       !* neighbour latitudes, we use the ansatz (It is not totally clear this
       !* is correct but at least it avoids parcels to stick at the pole)
       ROT = - STRCMP(1) * ALONG

       SLONG = SIN( RADPDG * ROT )
       CLONG = COS( RADPDG * ROT )
       XPOLG = SLONG * STRCMP(5) + CLONG * STRCMP(6)
       YPOLG = CLONG * STRCMP(5) - SLONG * STRCMP(6)
       UG = YPOLG * UE + XPOLG * VN
       VG = YPOLG * VN - XPOLG * UE
       RETURN        
       END SUBROUTINE CC2GLL

       REAL(dbl) FUNCTION CGSZLL (STRCMP, XLAT,XLONG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
       REAL(dbl), INTENT(IN) :: STRCMP(9), XLAT, XLONG
       REAL(dbl) SLAT,YMERC,EFACT
       IF (XLAT .GT. 89.985D0) THEN
!* CLOSE TO NORTH POLE
         IF (STRCMP(1) .GT. 0.9999D0) THEN
!* AND TO GAMMA == 1.
           CGSZLL = 2. * STRCMP(7)
           RETURN
         ENDIF
         EFACT = COS(RADPDG * XLAT)
         IF (EFACT .LE. 0.D0) THEN
           CGSZLL = 0.D0
           RETURN
         ELSE
           YMERC = - LOG( EFACT /(1.D0 + SIN(RADPDG * XLAT)))
         ENDIF
       ELSE IF (XLAT .LT. -89.985D0) THEN
!* CLOSE TO SOUTH POLE
         IF (STRCMP(1) .LT. -0.9999D0) THEN
!* AND TO GAMMA == -1.0
           CGSZLL = 2 * STRCMP(7)
           RETURN
         ENDIF
         EFACT = COS(RADPDG * XLAT)
         IF (EFACT .LE. 0.D0) THEN
           CGSZLL = 0.D0
           RETURN
         ELSE
           YMERC = LOG( EFACT /(1.D0 - SIN(RADPDG * XLAT)))
         ENDIF
       ELSE
       SLAT = SIN(RADPDG * XLAT)
       YMERC = LOG((1.D0 + SLAT) / (1.D0 - SLAT))/2.
!       EFACT = EXP(YMERC)
!       CGSZLL = 2. * STRCMP(7) * EXP (STRCMP(1) * YMERC)
!     C        / (EFACT + 1./EFACT)
       ENDIF
       CGSZLL = STRCMP(7) * COS(RADPDG * XLAT) * EXP(STRCMP(1) *YMERC)
       RETURN
       END FUNCTION CGSZLL
                            
       SUBROUTINE CLL2XY (STRCMP, XLAT,XLONG, X,Y)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL 
       REAL(dbl), INTENT(IN) ::  STRCMP(9), XLAT, XLONG
       REAL(dbl), INTENT(OUT) :: X, Y
       REAL(dbl) :: XI, ETA
       CALL CNLLXY(STRCMP, XLAT,XLONG, XI,ETA)
       X = STRCMP(3) + REARTH/STRCMP(7) *                  &
             (XI * STRCMP(5) + ETA * STRCMP(6) )
       Y = STRCMP(4) + REARTH/STRCMP(7) *                  &
             (ETA * STRCMP(5) - XI * STRCMP(6) )
       RETURN
       END SUBROUTINE CLL2XY

       SUBROUTINE CNLLXY (STRCMP, XLAT,XLONG, XI,ETA)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
!  MAIN TRANSFORMATION ROUTINE FROM LATITUDE-LONGITUDE TO
!  CANONICAL (EQUATOR-CENTERED, RADIAN UNIT) COORDINATES  
       REAL(dbl), INTENT(IN) :: STRCMP(9), XLAT, XLONG
       REAL(dbl), INTENT(OUT) :: XI, ETA
       REAL(dbl) :: GAMMAV, GDLONG, SNDGAM, CSDGAM, RHOG1
       REAL(dbl) :: DLONG,DLAT,SLAT,MERCY,GMERCY
       GAMMAV = STRCMP(1)
       DLAT = XLAT
       DLONG = CSPANF(XLONG - STRCMP(2), -180.D0, 180.D0)
       DLONG = DLONG * RADPDG
       GDLONG = GAMMAV * DLONG
       IF (ABS(GDLONG) .LT. .01D0) THEN
!  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
!  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
         GDLONG = GDLONG * GDLONG
         SNDGAM = DLONG * (1.D0 - 1.D0/6.D0 * GDLONG *        &
                          (1.D0 - 1.D0/20.D0 * GDLONG *       &
                          (1.D0 - 1.D0/42.D0 * GDLONG )))
         CSDGAM = DLONG * DLONG * .5 *                  &
                          (1.D0 - 1.D0/12.D0 * GDLONG *       &
                          (1.D0 - 1.D0/30.D0 * GDLONG *       &
                          (1.D0 - 1.D0/56.D0 * GDLONG )))
       ELSE
! CODE FOR MODERATE VALUES OF GAMMA
         SNDGAM = SIN (GDLONG) /GAMMAV
         CSDGAM = (1. - COS(GDLONG) )/GAMMAV /GAMMAV
       ENDIF
       SLAT = SIN(RADPDG * DLAT)
       IF ((SLAT .GE. ALMST1) .OR. (SLAT .LE. -ALMST1)) THEN
         ETA = 1.D0/STRCMP(1)
         XI = 0.D0
         RETURN
       ENDIF
       MERCY = .5D0 * LOG( (1.D0 + SLAT) / (1.D0 - SLAT) )
       GMERCY = GAMMAV * MERCY
       IF (ABS(GMERCY) .LT. .001D0) THEN
!  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
!  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
         RHOG1 = MERCY * (1.D0 - .5D0 * GMERCY *          &
                         (1.D0 - 1.D0/3.D0 * GMERCY *       &
                         (1.D0 - 1.D0/4.D0 * GMERCY ) ) )
       ELSE
! CODE FOR MODERATE VALUES OF GAMMA
         RHOG1 = (1.D0 - EXP(-GMERCY)) / GAMMAV
       ENDIF
       ETA = RHOG1 + (1.D0 - GAMMAV * RHOG1) * GAMMAV * CSDGAM
       XI = (1. - GAMMAV * RHOG1 ) * SNDGAM
       END SUBROUTINE CNLLXY

       SUBROUTINE CNXYLL (STRCMP, XI,ETA, XLAT,XLONG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
!  MAIN TRANSFORMATION ROUTINE FROM CANONICAL (EQUATOR-CENTERED,
!  RADIAN UNIT) COORDINATES
       REAL(dbl), INTENT(IN) :: STRCMP(9)
       REAL(dbl), INTENT(IN) ::  XI, ETA
       REAL(dbl), INTENT(OUT):: XLAT, XLONG
       REAL(dbl) :: GAMMAV,TEMP,ARG1,ARG2,YMERC,ALONG,GXI,CGETA,ODIST
       GAMMAV = STRCMP(1)
!  CALCULATE EQUIVALENT MERCATOR COORDINATE
       ODIST = XI*XI + ETA*ETA
       ARG2 = 2. * ETA - GAMMAV * (XI*XI + ETA*ETA)
       ARG1 = GAMMAV * ARG2
! Change by A. Stohl to avoid problems close to the poles
! IF (ARG1 .GE. ALMST1) THEN
!  DISTANCE TO NORTH (OR SOUTH) POLE IS ZERO (OR IMAGINARY ;) )
! XLAT = SIGN(90.,STRCMP(1))
! XLONG = STRCMP(2)
! RETURN
! ENDIF
       IF (ABS(ARG1) .LT. .01D0) THEN
!  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
!  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
         TEMP = (ARG1 / (2.D0 - ARG1) )**2
         YMERC = ARG2 / (2.D0 - ARG1) * (1.D0    + TEMP *   &
                                      (1.D0/3.D0 + TEMP *   &
                                      (1.D0/5.D0 + TEMP *   &
                                      (1.D0/7.D0 ))))
       ELSE
! CODE FOR MODERATE VALUES OF GAMMA
         YMERC = - LOG ( 1.D0 - ARG1 ) /2.D0 / GAMMAV
       ENDIF
!  CONVERT YMERC TO LATITUDE
       TEMP = EXP( - ABS(YMERC) )
       XLAT = SIGN(ATAN2((1.D0 - TEMP) * (1.D0 + TEMP), 2.D0 * TEMP), YMERC)
!  FIND LONGITUDES
       GXI = GAMMAV*XI
       CGETA = 1.D0 - GAMMAV * ETA
       IF ( ABS(GXI) .LT. .01D0*CGETA ) THEN
!  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
!  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
         TEMP = ( GXI /CGETA )**2
         ALONG = XI / CGETA * (1.D0    - TEMP *     &
                              (1.D0/3.D0 - TEMP *     &
                              (1.D0/5.D0 - TEMP *     &
                              (1.D0/7.D0   ))))
       ELSE
! CODE FOR MODERATE VALUES OF GAMMA
         ALONG = ATAN2( GXI, CGETA) / GAMMAV
       ENDIF
       XLONG = SNGL(STRCMP(2) + DGPRAD * ALONG)
       XLAT = XLAT * DGPRAD
       RETURN
       END SUBROUTINE CNXYLL

       REAL(dbl) FUNCTION CSPANF (VALU, BEGIN, ENDVAL)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
!* REAL FUNCTION CSPANF RETURNS A VALUE IN THE INTERVAL (BEGIN,END]
!* WHICH IS EQUIVALENT TO VALUE, MOD (END - BEGIN).  IT IS USED TO
!* REDUCE PERIODIC VARIABLES TO A STANDARD RANGE.  IT ADJUSTS FOR THE
!* BEHAVIOR OF THE MOD FUNCTION WHICH PROVIDES POSITIVE RESULTS FOR
!* POSITIVE INPUT, AND NEGATIVE RESULTS FOR NEGATIVE INPUT
!* INPUT:
!*       VALUE - REAL NUMBER TO BE REDUCED TO THE SPAN
!*       BEGIN - FIRST VALUE OF THE SPAN
!*       END   - LAST VALUE OF THE SPAN
!* RETURNS:
!*       THE REDUCED VALUE
!* EXAMPLES:
!*      ALONG = CSPANF(XLONG, -180., +180.)
!*      DIR  = CSPANF(ANGLE, 0., 360.)
       REAL(dbl) :: FIRST,LAST,VAL
       REAL(dbl), INTENT(IN) :: VALU, BEGIN, ENDVAL
       FIRST = MIN(BEGIN,ENDVAL)
       LAST = MAX(BEGIN,ENDVAL)
       VAL = MOD( VALU - FIRST , LAST - FIRST)
       IF ( VAL .LE. 0.D0) THEN
         CSPANF = VAL + LAST
       ELSE
         CSPANF = VAL + FIRST
       ENDIF
       RETURN
       END FUNCTION CSPANF

       SUBROUTINE CXY2LL (STRCMP, X,Y, XLAT,XLONG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
       IMPLICIT NONE
       REAL(dbl), PARAMETER :: REARTH=6367.47D0
       REAL(dbl) XI0,ETA0,XI,ETA
       REAL(dbl), INTENT(IN) :: STRCMP(9), X, Y
       REAL(dbl), INTENT(OUT) :: XLAT, XLONG 
       XI0 = ( X - STRCMP(3) ) * STRCMP(7) / REARTH
       ETA0 = ( Y - STRCMP(4) ) * STRCMP(7) /REARTH
       XI = XI0 * STRCMP(5) - ETA0 * STRCMP(6)
       ETA = ETA0 * STRCMP(5) + XI0 * STRCMP(6)
       CALL CNXYLL(STRCMP, XI,ETA, XLAT,XLONG)
       XLONG = CSPANF(XLONG, -180.D0, 180.D0)
       RETURN
       END SUBROUTINE CXY2LL

       SUBROUTINE STCM1P(STRCMP, X1,Y1, XLAT1,XLONG1,  &
                         XLATG,XLONGG, GRIDSZ, ORIENT)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
       REAL(dbl), INTENT(INOUT) :: STRCMP(9)
       REAL(dbl), INTENT(IN) :: ORIENT, GRIDSZ, XLATG, XLONGG, X1, Y1, XLAT1, XLONG1
       REAL(dbl) :: TURN, X1A, Y1A 
       INTEGER :: K
       DO K=3,4
         STRCMP (K) = 0.D0
       ENDDO
       TURN = RADPDG * (ORIENT - STRCMP(1) *      &
                  CSPANF(XLONGG - STRCMP(2), -180.D0, 180.D0) )
       STRCMP (5) = COS (TURN)
       STRCMP (6) = - SIN (TURN)
       STRCMP (7) = 1.D0
       STRCMP (7) = GRIDSZ * STRCMP(7)            &
                   / CGSZLL(STRCMP, XLATG, STRCMP(2))
       CALL CLL2XY (STRCMP, XLAT1,XLONG1, X1A,Y1A)
       STRCMP(3) = STRCMP(3) + X1 - X1A
       STRCMP(4) = STRCMP(4) + Y1 - Y1A
       RETURN
       END SUBROUTINE STCM1P

       SUBROUTINE STCM2P(STRCMP, X1,Y1, XLAT1,XLONG1, &
                                  X2,Y2, XLAT2,XLONG2)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
       REAL(dbl), INTENT(INOUT) :: STRCMP(9)
       REAL(dbl), INTENT(IN) :: X1, Y1, X2, Y2, XLAT1, XLONG1, XLAT2, XLONG2
       REAL(dbl) :: X1A, Y1A, X2A, Y2A, DEN, DENA
       INTEGER :: K
       DO K=3,6
         STRCMP (K) = 0.D0
       ENDDO
       STRCMP (5) = 1.D0
       STRCMP (7) = 1.D0
       CALL CLL2XY (STRCMP, XLAT1,XLONG1, X1A,Y1A)
       CALL CLL2XY (STRCMP, XLAT2,XLONG2, X2A,Y2A)
       DEN = SQRT( (X1 - X2)**2 + (Y1 - Y2)**2 )
       DENA = SQRT( (X1A - X2A)**2 + (Y1A - Y2A)**2 )
       STRCMP(5) = ((X1A - X2A)*(X1 - X2) + (Y1A - Y2A) * (Y1 - Y2)) &
          /DEN /DENA
       STRCMP(6) = ((Y1A - Y2A)*(X1 - X2) - (X1A - X2A) * (Y1 - Y2)) &
          /DEN /DENA
       STRCMP (7) = STRCMP(7) * DENA / DEN
       CALL CLL2XY (STRCMP, XLAT1,XLONG1, X1A,Y1A)
       STRCMP(3) = STRCMP(3) + X1 - X1A
       STRCMP(4) = STRCMP(4) + Y1 - Y1A
       RETURN
       END SUBROUTINE STCM2P

       SUBROUTINE STLMBR(STRCMP, TNGLAT, XLONG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
       REAL(dbl), INTENT(OUT) :: STRCMP(9)
       REAL(dbl), INTENT(IN) :: TNGLAT, XLONG
       REAL(dbl) :: XI, ETA
       STRCMP(1) = SIN(RADPDG * TNGLAT)
!*  GAMMA = SINE OF THE TANGENT LATITUDE
       STRCMP(2) = CSPANF( XLONG, -180.D0, +180.D0)
!* LAMBDA_0 = REFERENCE LONGITUDE
       STRCMP(3) = 0.D0
!* X_0 = X- GRID COORDINATE OF ORIGIN (XI,ETA) = (0.,0.)
       STRCMP(4) = 0.D0
!* y_0 = Y-GRID COORDINATE OF ORIGIN (XI,ETA) = (0.,0.)
       STRCMP(5) = 1.D0
!* COSINE OF ROTATION ANGLE FROM XI,ETA TO X,Y
       STRCMP(6) = 0.D0
!* SINE OF ROTATION ANGLE FROM XI,ETA TO X,Y
       STRCMP(7) = REARTH
!* GRIDSIZE IN KILOMETERS AT THE EQUATOR
       CALL CNLLXY(STRCMP, 89.D0,XLONG, XI,ETA)
       STRCMP(8) = 2 * ETA - STRCMP(1) * ETA * ETA
!* RADIAL COORDINATE FOR 1 DEGREE FROM NORTH POLE
       CALL CNLLXY(STRCMP, -89.D0,XLONG, XI,ETA)
          STRCMP(9) = 2 * ETA - STRCMP(1) * ETA * ETA
!* RADIAL COORDINATE FOR 1 DEGREE FROM SOUTH POLE
       RETURN
       END SUBROUTINE STLMBR

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ COORDTRAFO @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

      subroutine coordtrafo(error)

!**********************************************************************
!                                                                     * 
!             FLEXPART MODEL SUBROUTINE COORDTRAFO                    *
!                                                                     *
!**********************************************************************
!                                                                     * 
!             AUTHOR:      G. WOTAWA                                  *
!             DATE:        1994-02-07                                 *
!             LAST UPDATE: 1996-05-18   A. STOHL                      *
!                                                                     * 
!**********************************************************************
!                                                                     *
! DESCRIPTION: This subroutine transforms x and y coordinates of      *
! particle release points to grid coordinates.                        *
!                                                                     *
!**********************************************************************

      use commons
      implicit none

      integer i,j
      logical error

      error=.false.

      print *,'coordtrafo> xlon0, ylat0 ',xlon0,ylat0
      print *,'coordtrafo> nx, ny ',nx,ny

      if (numpoint.eq.0) goto 30

! TRANSFORM X- AND Y- COORDINATES OF STARTING POINTS TO GRID COORDINATES
!***********************************************************************

      do i=1,numpoint
        xpoint1(i)=(xpoint1(i)-xlon0)/dx
        xpoint2(i)=(xpoint2(i)-xlon0)/dx
        ypoint1(i)=(ypoint1(i)-ylat0)/dy
        ypoint2(i)=(ypoint2(i)-ylat0)/dy
      enddo

      print *,'coordtrafo> xpoint1(1), ypoint1(1) ', xpoint1(1),ypoint1(1) 

15    continue

! CHECK IF RELEASE POINTS ARE WITHIN DOMAIN
!******************************************

      do i=1,numpoint
      if ((ypoint1(i).lt.1.e-6).or.(ypoint1(i).ge.float(ny-1)-1.e-6).or. &
         (ypoint2(i).lt.1.e-6).or.(ypoint2(i).ge.float(ny-1)-1.e-6).or.  &
         (xpoint1(i).lt.1.e-6).or.(xpoint1(i).ge.float(nx-1)-1.e-6).or.  &
         (xpoint2(i).lt.1.e-6).or.(xpoint2(i).ge.float(nx-1)-1.e-6)) then 
          write(*,*) ' NOTICE: RELEASE POINT OUT OF DOMAIN DETECTED.'
          write(*,*) ' IT IS REMOVED NOW ... '
          write(*,*) ' COMMENT: ',compoint(i)

          if (i.lt.numpoint) then
            do j=i+1,numpoint
              xpoint1(j-1)=xpoint1(j)
              ypoint1(j-1)=ypoint1(j)
              xpoint2(j-1)=xpoint2(j)
              ypoint2(j-1)=ypoint2(j)
              zpoint1(j-1)=zpoint1(j)
              zpoint2(j-1)=zpoint2(j)
              npart(j-1)=npart(j)
              compoint(j-1)=compoint(j)
            enddo         
          endif

          numpoint=numpoint-1
          if (numpoint.gt.0) goto 15
        endif
      enddo

30    if (numpoint.eq.0) then
        error=.true.
        write(*,*) ' FLEXPART MODEL SUBROUTINE COORDTRAFO: ERROR ! '
        write(*,*) ' NO TRAJECTORY STARTING POINTS ARE GIVEN !!!'
      endif

      return
      end subroutine coordtrafo
      
end module coord
!
!=====|==1=========2=========3=========4=========5=========6=========7=========8
!
!$Log: 
