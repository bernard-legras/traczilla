module interpol

integer :: icode

contains

      SUBROUTINE SLOPES(XTAB, YTAB, MTAB, NUM)                          
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C.STATE UNIVERSITY
!
      INTEGER :: NUM
      REAL, INTENT(IN) :: XTAB(:), YTAB(:)
      REAL, INTENT(OUT):: MTAB(:)

      REAL :: M1, M2, XBAR, XHAT, YDIF1,YDIF2, YXMID, XMID, M1S, M2S
!
! SLOPES CALCULATES THE DERIVATIVE AT EACH OF THE DATA POINTS.  THE
! SLOPES PROVIDED WILL INSURE THAT AN OSCULATORY QUADRATIC SPLINE WILL
! HAVE ONE ADDITIONAL KNOT BETWEEN TWO ADJACENT POINTS OF INTERPOLATION.
! CONVEXITY AND MONOTONICITY ARE PRESERVED WHEREVER THESE CONDITIONS
! ARE COMPATIBLE WITH THE DATA.
!
! ON INPUT--
!
!   XTAB CONTAINS THE ABSCISSAS OF THE DATA POINTS.
!
!   YTAB CONTAINS THE ORDINATES OF THE DATA POINTS.
!
!   NUM IS THE NUMBER OF DATA POINTS (DIMENSION OF XTAB,YTAB).
!
!
! ON OUTPUT--
!
!   MTAB CONTAINS THE VALUE OF THE FIRST DERIVATIVE AT EACH DATA POINT.
!
! AND
!
!   SLOPES DOES NOT ALTER XTAB,YTAB,NUM.
!
!-----------------------------------------------------------------------
!
      NUM1 = NUM - 1
      IM1 = 1
      I = 2
      I1 = 3
!
! CALCULATE THE SLOPES OF THE TWO LINES JOINING THE FIRST THREE DATA
! POINTS.
      YDIF1 = YTAB(2) - YTAB(1)
      YDIF2 = YTAB(3) - YTAB(2)
      M1 = YDIF1/(XTAB(2)-XTAB(1))
      M1S = M1
      M2 = YDIF2/(XTAB(3)-XTAB(2))
      M2S = M2
!
! IF ONE OF THE PRECEDING SLOPES IS ZERO OR IF THEY HAVE OPPOSITE SIGN,
! ASSIGN THE VALUE ZERO TO THE DERIVATIVE AT THE MIDDLE POINT.
   10 IF (M1.EQ.0.E0.OR.M2.EQ.0.E0.OR.(M1*M2).LE.0.E0) GO TO 20
      IF (ABS(M1).GT.ABS(M2)) GO TO 30
      GO TO 40
   20 MTAB(I) = 0.E0
      GO TO 50
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M1.
   30 XBAR = (YDIF2/M1) + XTAB(I)
      XHAT = (XBAR+XTAB(I1))/2.E0
      MTAB(I) = YDIF2/(XHAT-XTAB(I))
      GO TO 50
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M2.
   40 XBAR = (-YDIF1/M2) + XTAB(I)
      XHAT = (XTAB(IM1)+XBAR)/2.E0
      MTAB(I) = YDIF1/(XTAB(I)-XHAT)
!
! INCREMENT COUNTERS.
   50 IM1 = I
      I = I1
      I1 = I1 + 1
      IF (I.GT.NUM1) GO TO 60
!
! CALCULATE THE SLOPES OF THE TWO LINES JOINING THREE CONSECUTIVE DATA
! POINTS.
      YDIF1 = YTAB(I) - YTAB(IM1)
      YDIF2 = YTAB(I1) - YTAB(I)
      M1 = YDIF1/(XTAB(I)-XTAB(IM1))
      M2 = YDIF2/(XTAB(I1)-XTAB(I))
      GO TO 10
!
! CALCULATE THE SLOPE AT THE LAST POINT, XTAB(NUM).
   60 IF ((M1*M2).LT.0.E0) GO TO 80
      XMID = (XTAB(NUM1)+XTAB(NUM))/2.E0
      YXMID = MTAB(NUM1)*(XMID-XTAB(NUM1)) + YTAB(NUM1)
      MTAB(NUM) = (YTAB(NUM)-YXMID)/(XTAB(NUM)-XMID)
      IF ((MTAB(NUM)*M2).LT.0.E0) GO TO 70
      GO TO 90
   70 MTAB(NUM) = 0.E0
      GO TO 90
   80 MTAB(NUM) = 2.E0*M2
!
! CALCULATE THE SLOPE AT THE FIRST POINT, XTAB(1).
   90 IF ((M1S*M2S).LT.0.E0) GO TO 110
      XMID = (XTAB(1)+XTAB(2))/2.E0
      YXMID = MTAB(2)*(XMID-XTAB(2)) + YTAB(2)
      MTAB(1) = (YXMID-YTAB(1))/(XMID-XTAB(1))
      IF ((MTAB(1)*M1S).LT.0.E0) GO TO 100
      RETURN
  100 MTAB(1) = 0.E0
      RETURN
  110 MTAB(1) = 2.E0*M1S
      RETURN
!
      END SUBROUTINE SLOPES

      REAL FUNCTION SPLINE(XVALS, Z1, Z2, XTABS, YTABS, XTABS1, YTABS1, &
                      Y1, Y2, E2, W2, V2, NCASE)
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
      REAL XVALS, Z1, Z2, XTABS, YTABS, XTABS1, YTABS1, &
           V2, W2, Y1, Y2, E2
!
!   SPLINE FINDS THE IMAGE OF A POINT IN XVAL.
!
! ON INPUT--
!
!   XVALS CONTAINS THE VALUE AT WHICH THE SPLINE IS EVALUATED.
!
!   (XTABS,YTABS) ARE THE COORDINATES OF THE LEFT-HAND DATA POINT
!   USED IN THE EVALUATION OF XVALS.
!
!   (XTABS1,YTABS1) ARE THE COORDINATES OF THE RIGHT-HAND DATA POINT
!   USED IN THE EVALUATION OF XVALS.
!
!   Z1,Z2,Y1,Y2,E2,W2,V2 ARE THE PARAMETERS OF THE SPLINE.
!
!   NCASE CONTROLS THE EVALUATION OF THE SPLINE BY INDICATING WHETHER
!   ONE OR TWO KNOTS WERE PLACED IN THE INTERVAL (XTABS,XTABS1).
!
!
! ON OUTPUT--
!
!   SPLINE IS THE IMAGE OF XVALS.
!
! AND
!
!   SPLINE DOES NOT ALTER ANY OF THE INPUT PARAMETERS.
!
!-----------------------------------------------------------------------
!
! IF NCASE .EQ. 4, MORE THAN ONE KNOT WAS PLACED IN THE INTERVAL.
      IF (NCASE.EQ.4) GO TO 40
!
! CASES 1,2, OR 3.
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE KNOT.
      IF (Z1-XVALS) 10, 20, 30

   10 SPLINE = (Z2*(XTABS1-XVALS)**2+W2*2.E0*(XVALS-Z1)*(XTABS1-XVALS) &
               +YTABS1*(XVALS-Z1)**2)/(XTABS1-Z1)**2
      RETURN

   20 SPLINE = Z2
      RETURN

   30 SPLINE = (YTABS*(Z1-XVALS)**2+V2*2.E0*(XVALS-XTABS)*(Z1-XVALS)+ &
               Z2*(XVALS-XTABS)**2)/(Z1-XTABS)**2
      RETURN

! CASE 4.

! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE FIRST KNOT.
   40 IF (Y1-XVALS) 70, 60, 50

   50 SPLINE = (YTABS*(Y1-XVALS)**2+V2*2.E0*(XVALS-XTABS)*(Y1-XVALS)+ &
                Y2*(XVALS-XTABS)**2)/(Y1-XTABS)**2
      RETURN

   60 SPLINE = Y2
      RETURN

! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE SECOND KNOT.
   70 IF (Z1-XVALS) 100, 90, 80

   80 SPLINE = (Y2*(Z1-XVALS)**2+E2*2.E0*(XVALS-Y1)*(Z1-XVALS) &
               +Z2*(XVALS-Y1)**2)/(Z1-Y1)**2
      RETURN

   90 SPLINE = Z2
      RETURN

  100 SPLINE = (Z2*(XTABS1-XVALS)**2+W2*2.E0*(XVALS-Z1)*(XTABS1-XVALS) &
               +YTABS1*(XVALS-Z1)**2)/(XTABS1-Z1)**2
      RETURN


      END FUNCTION SPLINE      

      SUBROUTINE MEVAL_SPE(XVAL, YVAL, XTAB, YTAB, MTAB, NUM, NUME, &
                           EPS, ERR)
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
      INTEGER, INTENT(IN) :: NUM, NUME
      REAL, INTENT(IN) :: XVAL(:), YTAB(:), MTAB(:), EPS
      REAL, INTENT(INOUT) :: XTAB(:)
      REAL, INTENT(OUT):: YVAL(:)
      INTEGER, INTENT(OUT) :: ERR
      REAL :: V1, V2, W1, W2, Z1, Z2, Y1, Y2, E1, E2, mute
      INTEGER :: START, START1, END, END1, FND    
!
! MEVAL CONTROLS THE EVALUATION OF AN OSCULATORY QUADRATIC SPLINE.  THE
! USER MAY PROVIDE HIS OWN SLOPES AT THE POINTS OF INTERPOLATION OR USE
! THE SUBROUTINE 'SLOPES' TO CALCULATE SLOPES WHICH ARE CONSISTENT WITH
! THE SHAPE OF THE DATA.
!
!
!
! ON INPUT--
!
!   XVAL MUST BE A NONDECREASING VECTOR OF POINTS AT WHICH THE SPLINE
!   WILL BE EVALUATED.
!
!   XTAB CONTAINS THE ABSCISSAS OF THE DATA POINTS TO BE INTERPOLATED.
!   XTAB MUST BE INCREASING.
!
!   YTAB CONTAINS THE ORDINATES OF THE DATA POINTS TO BE INTERPOLATED.
!
!   MTAB CONTAINS THE SLOPE OF THE SPLINE AT EACH POINT OF INTERPOLA-
!   TION.
!
!   NUM IS THE NUMBER OF DATA POINTS (DIMENSION OF XTAB AND YTAB).
!
!   NUME IS THE NUMBER OF POINTS OF EVALUATION (DIMENSION OF XVAL AND
!   YVAL).
!
!   EPS IS A RELATIVE ERROR TOLERANCE USED IN SUBROUTINE 'CHOOSE'
!   TO DISTINGUISH THE SITUATION MTAB(I) OR MTAB(I+1) IS RELATIVELY
!   CLOSE TO THE SLOPE OR TWICE THE SLOPE OF THE LINEAR SEGMENT
!   BETWEEN XTAB(I) AND XTAB(I+1).  IF THIS SITUATION OCCURS,
!   ROUNDOFF MAY CAUSE A CHANGE IN CONVEXITY OR MONOTONICITY OF THE
!   RESULTING SPLINE AND A CHANGE IN THE CASE NUMBER PROVIDED BY
!   CHOOSE.  IF EPS IS NOT EQUAL TO ZERO, THEN EPS SHOULD BE GREATER
!   THAN OR EQUAL TO MACHINE EPSILON.
!
!
! ON OUTPUT--
!
! YVAL CONTAINS THE IMAGES OF THE POINTS IN XVAL.
!
!   ERR IS AN ERROR CODE--
!      ERR=0 - MEVAL RAN NORMALLY.
!      ERR=1 - XVAL(I) IS LESS THAN XTAB(1) FOR AT LEAST ONE I OR
!              XVAL(I) IS GREATER THAN XTAB(NUM) FOR AT LEAST ONE I.
!              MEVAL WILL EXTRAPOLATE TO PROVIDE FUNCTION VALUES FOR
!              THESE ABSCISSAS.
!      ERR=2 - XVAL(I+1) .LT. XVAL(I) FOR SOME I.
!
! AND
!
!   MEVAL DOES NOT ALTER XVAL,XTAB,YTAB,MTAB,NUM,NUME.
!
!
!   MEVAL CALLS THE FOLLOWING SUBROUTINES OR FUNCTIONS:
!      SEARCH
!      CASES
!      CHOOSE
!      SPLINE
!
!-----------------------------------------------------------------------
!
      START = 1
      END = NUME
      ERR = 0
      IF (NUME.EQ.1) GO TO 20
!
! DETERMINE IF XVAL IS NONDECREASING.
      NUME1 = NUME - 1
      DO I=1,NUME1
        IF (XVAL(I+1).GE.XVAL(I)) CYCLE
        ERR = 2
        GO TO 230
      ENDDO
   
! DETERMINE IF XTAB IS INCREASING
   20 DO I=1,num-1
        IF (xtab(i) > xtab(i+1)) THEN
          write(*,'(a,2G12.5)') 'Inverted levels ',xtab(i),xtab(i+1)
          !special version that invert xtab in this uncommon case!
          ! err=2
          mute = xtab(i)
          xtab(i)=xtab(i+1)
          xtab(i+1) = mute
        ENDIF 
      ENDDO
!
! IF XVAL(I) .LT. XTAB(1), THEN XVAL(I)=YTAB(1).
! IF XVAL(I) .GT. XTAB(NUM), THEN XVAL(I)=YTAB(NUM).
!
! DETERMINE IF ANY OF THE POINTS IN XVAL ARE LESS THAN THE ABSCISSA OF
! THE FIRST DATA POINT.
      DO I=1,NUME
        IF (XVAL(I).GE.XTAB(1)) GO TO 40
        START = I + 1
      ENDDO
!
!
   40 NUME1 = NUME + 1
!
! DETERMINE IF ANY OF THE POINTS IN XVAL ARE GREATER THAN THE ABSCISSA
! OF THE LAST DATA POINT.
      DO I=1,NUME
        IND = NUME1 - I
        IF (XVAL(IND).LE.XTAB(NUM)) GO TO 60
        END = IND - 1
      ENDDO

! CALCULATE THE IMAGES OF POINTS OF EVALUATION WHOSE ABSCISSAS
! ARE LESS THAN THE ABSCISSA OF THE FIRST DATA POINT.
   60 IF (START.EQ.1) GO TO 80
! SET THE ERROR PARAMETER TO INDICATE THAT EXTRAPOLATION HAS OCCURRED.
      ERR = 1
      CALL CHOOSE(XTAB(1), YTAB(1), MTAB(1), MTAB(2), XTAB(2), YTAB(2), &
                  EPS, NCASE)
      CALL CASES(XTAB(1), YTAB(1), MTAB(1), MTAB(2), XTAB(2), YTAB(2), &
                 E1, E2, V1, V2, W1, W2, Z1, Z2, Y1, Y2, NCASE)
      START1 = START - 1
      DO 70 I=1,START1
        YVAL(I) = SPLINE(XVAL(I),Z1,Z2,XTAB(1),YTAB(1),XTAB(2),YTAB(2), &
                         Y1,Y2,E2,W2,V2,NCASE)
   70 CONTINUE
      IF (NUME.EQ.1) GO TO 230

! SEARCH LOCATES THE INTERVAL IN WHICH THE FIRST IN-RANGE POINT OF
! EVALUATION LIES.
   80 IF ((NUME.EQ.1) .AND. (END.NE.NUME)) GO TO 200
      CALL SEARCH(XTAB, NUM, XVAL(START), LCN, FND)

      LCN1 = LCN + 1


! IF THE FIRST IN-RANGE POINT OF EVALUATION IS EQUAL TO ONE OF THE DATA
! POINTS, ASSIGN THE APPROPRIATE VALUE FROM YTAB.  CONTINUE UNTIL A
! POINT OF EVALUATION IS FOUND WHICH IS NOT EQUAL TO A DATA POINT.
      IF (FND.EQ.0) GO TO 130
   90 YVAL(START) = YTAB(LCN)
      START1 = START
      START = START + 1
      IF (START.GT.NUME) GO TO 230
      IF (XVAL(START1).EQ.XVAL(START)) GO TO 90

  100 IF (XVAL(START)-XTAB(LCN1)) 130, 110, 120

  110 YVAL(START) = YTAB(LCN1)
      START1 = START
      START = START + 1
      IF (START.GT.NUME) GO TO 230
      IF (XVAL(START).NE.XVAL(START1)) GO TO 120
      GO TO 110

  120 LCN = LCN1
      LCN1 = LCN1 + 1
      GO TO 100

! CALCULATE THE IMAGES OF ALL THE POINTS WHICH LIE WITHIN RANGE OF THE
! DATA.

  130 IF ((LCN.EQ.1) .AND. (ERR.EQ.1)) GO TO 140
      CALL CHOOSE(XTAB(LCN), YTAB(LCN), MTAB(LCN), MTAB(LCN1), &
                  XTAB(LCN1), YTAB(LCN1), EPS, NCASE)
      CALL CASES(XTAB(LCN), YTAB(LCN), MTAB(LCN), MTAB(LCN1), &
                 XTAB(LCN1), YTAB(LCN1), &
                 E1, E2, V1, V2, W1, W2, Z1, Z2, Y1, Y2, NCASE)

  140 DO 190 I=START,END
!
! IF XVAL(I) -XTAB(LCN1) IS NEGATIVE, DO NOT RECALCULATE THE PARAMETERS
! FOR THIS SECTION OF THE SPLINE SINCE THEY ARE ALREADY KNOWN.
      IF (XVAL(I)-XTAB(LCN1)) 150, 160, 170
!
  150 YVAL(I) = SPLINE(XVAL(I),Z1,Z2,XTAB(LCN),YTAB(LCN), &
                         XTAB(LCN1),YTAB(LCN1),Y1,Y2,E2,W2,V2,NCASE)
      GO TO 190
!
!  IF XVAL(I) IS A DATA POINT, ITS IMAGE IS KNOWN.
  160 YVAL(I) = YTAB(LCN1)
      GO TO 190
!
! INCREMENT THE POINTERS WHICH GIVE THE LOCATION IN THE DATA VECTOR.
  170 LCN = LCN1
      LCN1 = LCN + 1
!
! DETERMINE THAT THE ROUTINE IS IN THE CORRECT PART OF THE SPLINE.
      IF (XVAL(I)-XTAB(LCN1)) 180, 160, 170
!
! CALL CHOOSE TO DETERMINE THE APPROPRIATE CASE AND THEN CALL
! CASES TO COMPUTE THE PARAMETERS OF THE SPLINE.
  180 CALL CHOOSE(XTAB(LCN), YTAB(LCN), MTAB(LCN), MTAB(LCN1), &
                    XTAB(LCN1), YTAB(LCN1), EPS, NCASE)
      CALL CASES(XTAB(LCN), YTAB(LCN), MTAB(LCN), MTAB(LCN1), &
                   XTAB(LCN1), YTAB(LCN1), &
                   E1, E2, V1, V2, W1, W2, Z1, Z2, Y1, Y2, NCASE)
      GO TO 150
  190 CONTINUE
!
!
! CALCULATE THE IMAGES OF THE POINTS OF EVALUATION WHOSE ABSCISSAS
! ARE GREATER THAN THE ABSCISSA OF THE LAST DATA POINT.
      IF (END.EQ.NUME) GO TO 230
      IF ((LCN1.EQ.NUM) .AND. (XVAL(END).NE.XTAB(NUM))) GO TO 210
! SET THE ERROR PRARMETER TO INDICATE THAT EXTRAPOLATION HAS OCCURRED.
  200 ERR = 1
      NUM1 = NUM - 1
      CALL CHOOSE(XTAB(NUM1), YTAB(NUM1), MTAB(NUM1), XTAB(NUM), &
                  XTAB(NUM), YTAB(NUM), EPS, NCASE)
      CALL CASES(XTAB(NUM1), YTAB(NUM1), MTAB(NUM1), MTAB(NUM), &
                 XTAB(NUM), YTAB(NUM), &
                 E1, E2, V1, V2, W1, W2, Z1, Z2, Y1, Y2, NCASE)
  210 END1 = END + 1
      DO 220 I=END1,NUME
        YVAL(I) = SPLINE(XVAL(I),Z1,Z2,XTAB(NUM1),YTAB(NUM1), &
                         XTAB(NUM),YTAB(NUM),Y1,Y2,E2,W2,V2,NCASE)
  220 CONTINUE
!
!
  230 RETURN
  
      END SUBROUTINE MEVAL_SPE

      SUBROUTINE SEARCH(XTAB, NUM, S, LCN, FND)                         
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
      INTEGER, INTENT(IN) :: NUM
      REAL, INTENT(IN) :: XTAB(:), S
      INTEGER, INTENT(OUT) :: LCN, FND
      INTEGER FIRST
!
! SEARCH CONDUCTS A BINARY SEARCH FOR S.  SEARCH IS CALLED ONLY IF S IS
! BETWEEN XTAB(1) AND XTAB(NUM).
!
! ON INPUT--
!
!   XTAB CONTAINS THE ABSCISSAS OF THE DATA POINTS OF INTERPOLATION.
!
!   NUM IS THE DIMENSION OF XTAB.
!
!   S IS THE VALUE WHOSE RELATIVE POSITION IN XTAB IS LOCATED BY SEARCH.
!
!
! ON OUTPUT--
!
!   FND IS SET EQUAL TO 1 IF S IS FOUND IN XTAB AND IS SET EQUAL TO 0
!   OTHERWISE.
!
!   LCN IS THE INDEX OF THE LARGEST VALUE IN XTAB FOR WHICH XTAB(I)
!   .LT. S.
!
! AND
!
!   SEARCH DOES NOT ALTER XTAB,NUM,S.
!
!-----------------------------------------------------------------------
!
      FIRST = 1
      LAST = NUM
      FND = 0
!
      IF (XTAB(1).NE.S) GO TO 10
      LCN = 1
      FND = 1
      RETURN
!
   10 IF (XTAB(NUM).NE.S) GO TO 20
      LCN = NUM
      FND = 1
      RETURN
!
! IF (LAST-FIRST) .EQ. 1, S IS NOT IN XTAB.  SET POSITION EQUAL TO
! FIRST.
   20 IF ((LAST-FIRST).EQ.1) GO TO 30
!
      MIDDLE = (FIRST+LAST)/2
!
! CHECK IF S .EQ. XTAB(MIDDLE).  IF NOT, CONTINUE THE SEARCH IN THE
! APPROPRIATE HALF OF THE VECTOR XTAB.
      IF (XTAB(MIDDLE)-S) 40, 50, 60
!
   30 LCN = FIRST
      RETURN
   40 FIRST = MIDDLE
      GO TO 20
   50 LCN = MIDDLE
      FND = 1
      RETURN
   60 LAST = MIDDLE
      GO TO 20
      END SUBROUTINE SEARCH

      SUBROUTINE CHOOSE(P1, P2, M1, M2, Q1, Q2, EPS, NCASE)             
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!
      REAL, INTENT(IN) :: P1, P2, M1, M2, Q1, Q2, EPS
      INTEGER, INTENT(OUT) :: NCASE
      REAL MREF, MREF1, MREF2, SPQ, PROD, PROD1, PROD2
!
! CHOOSE DETERMINES THE CASE NEEDED FOR THE COMPUTATION OF THE PARAME-
! TERS OF THE QUADRATIC SPLINE AND RETURNS THE VALUE IN THE VARIABLE
! NCASE.
!
! ON INPUT--
!
!   (P1,P2) GIVES THE COORDINATES OF ONE OF THE POINTS OF INTERPOLATION.
!
!   M1 SPECIFIES THE DERIVATIVE CONDITION AT (P1,P2).
!
!   (Q1,Q2) GIVES THE COORDINATES OF ONE OF THE POINTS OF INTERPOLATION.
!
!   M2 SPECIFIES THE DERIVATIVE CONDITION AT (Q1,Q2).
!
!   EPS IS AN ERROR TOLERANCE USED TO DISTINGUISH CASES WHEN M1 OR M2 IS
!   RELATIVELY CLOSE TO THE SLOPE OR TWICE THE SLOPE OF THE LINE
!   SEGMENT JOINING (P1,P2) AND (Q1,Q2).  IF EPS IS NOT EQUAL TO ZERO,
!   THEN EPS SHOULD BE GREATER THAN OR EQUAL TO MACHINE EPSILON.
!
!
! ON OUTPUT--
!
!   NCASE CONTAINS THE VALUE WHICH CONTROLS HOW THE PARAMETERS OF THE
!   QUADRATIC SPLINE ARE EVALUATED.
!
! AND
!
!   CHOOSE DOES NOT ALTER P1,P2,Q1,Q2,M1,M2,EPS.
!
!-----------------------------------------------------------------------
!
! CALCULATE THE SLOPE SPQ OF THE LINE JOINING (P1,P2),(Q1,Q2).
      SPQ = (Q2-P2)/(Q1-P1)
!
! CHECK WHETHER OR NOT SPQ IS 0.
      IF (SPQ.NE.0.E0) GO TO 20
      IF ((M1*M2).GE.0.E0) GO TO 10
      NCASE = 1
      RETURN
   10 NCASE = 2
      RETURN
!
   20 PROD1 = SPQ*M1
      PROD2 = SPQ*M2
!
! FIND THE ABSOLUTE VALUES OF THE SLOPES SPQ,M1,AND M2.
      MREF = ABS(SPQ)
      MREF1 = ABS(M1)
      MREF2 = ABS(M2)
!
! IF THE RELATIVE DEVIATION OF M1 OR M2 FROM SPQ IS LESS THAN EPS, THEN
! CHOOSE CASE 2 OR CASE 3.
      IF ((ABS(SPQ-M1).LE.EPS*MREF) .OR. (ABS(SPQ-M2).LE.EPS*MREF)) &
         GO TO 30
!
! COMPARE THE SIGNS OF THE SLOPES SPQ,M1, AND M2.
      IF ((PROD1.LT.0.E0).OR.(PROD2.LT.0.E0)) GO TO 80
      PROD = (MREF-MREF1)*(MREF-MREF2)
      IF (PROD.GE.0.E0) GO TO 40
!
! L1, THE LINE THROUGH (P1,P2) WITH SLOPE M1, AND L2, THE LINE THROUGH
! (Q1,Q2) WITH SLOPE M2, INTERSECT AT A POINT WHOSE ABSCISSA IS BETWEEN
! P1 AND Q1.  THE ABSCISSA BECOMES A KNOT OF THE SPLINE.
      NCASE = 1
      RETURN
!
   30 IF ((PROD1.LT.0.E0).OR.(PROD2.LT.0.E0)) GO TO 80
   40 IF (MREF1.GT.(2.E0*MREF)) GO TO 50
      IF (MREF2.GT.(2.E0*MREF)) GO TO 60
!
! BOTH L1 AND L2 CROSS THE LINE THROUGH (P1+Q1/2.,P2) AND (P1+Q1/2.,Q2),
! WHICH IS THE MIDLINE OF THE RECTANGLE FORMED BY (P1,P2),(Q1,P2),
! (Q1,Q2), AND (P1,Q2), OR BOTH M1 AND M2 HAVE SIGNS DIFFERENT THAN THE
! SIGN OF SPQ, OR ONE OF M1 AND M2 HAS OPPOSITE SIGN FROM SPQ AND L1
! AND L2 INTERSECT TO THE LEFT OF P1 OR TO THE RIGHT OF Q1.  THE POINT
! (P1+Q1)/2. IS A KNOT OF THE SPLINE.
      NCASE = 2
      RETURN
!
! CHOOSE CASE 4 IF MREF2 IS GREATER THAN (2.-EPS)*MREF; OTHERWISE,
! CHOOSE CASE 3.
   50 IF (MREF2.GT.(2.E0-EPS)*MREF) GO TO 70
      NCASE = 3
      RETURN
!
! IN CASES 3 AND 4, SIGN(M1)=SIGN(M2)=SIGN(SPQ).
!
! EITHER L1 OR L2 CROSSES THE MIDLINE, BUT NOT BOTH.
! CHOOSE CASE 4 IF MREF1 IS GREATER THAN (2.-EPS)*MREF; OTHERWISE,
! CHOOSE CASE 3.
   60 IF (MREF1.GT.(2.E0-EPS)*MREF) GO TO 70
      NCASE = 3
      RETURN
!
! IF NEITHER L1 NOR L2 CROSSES THE MIDLINE, THE SPLINE REQUIRES TWO
! KNOTS BETWEEN P1 AND Q1.
   70 NCASE = 4
      RETURN
!
! THE SIGN OF AT LEAST ONE OF THE SLOPES M1,M2 DOES NOT AGREE WITH THE
! SIGN OF THE SLOPE SPQ.
   80 IF ((PROD1.LT.0.E0).AND.(PROD2.LT.0.E0)) GO TO 130
!
      IF (PROD1.LT.0.E0) GO TO 90
      GO TO 110
!
   90 IF (MREF2.GT.((1.E0+EPS)*MREF)) GO TO 100
      NCASE = 2
      RETURN
!
  100 NCASE = 1
      RETURN
!
  110 IF (MREF1.GT.((1.E0+EPS)*MREF)) GO TO 120
      NCASE = 2
      RETURN
!
  120 NCASE = 1
      RETURN
!
  130 NCASE = 2
      RETURN
!
      END SUBROUTINE CHOOSE

      SUBROUTINE CASES(P1, P2, M1, M2, Q1, Q2, E1, E2, V1, V2, W1, W2, & 
                       Z1, Z2, Y1, Y2, NCASE)
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
      REAL, INTENT(IN) ::  P1, P2, M1, M2, Q1, Q2
      INTEGER, INTENT(IN) :: NCASE
      REAL, INTENT(OUT) :: V1, V2, Z1, Z2, W1, W2, E1, E2,Y1, Y2
      REAL MBAR1, MBAR2, MBAR3, C1, D1, H1, J1, K1, ZTWO
!
! CASES COMPUTES THE KNOTS AND OTHER PARAMETERS OF THE SPLINE ON THE
! INTERVAL (P1,Q1).
!
!
! ON INPUT--
!
!   (P1,P2) AND (Q1,Q2) ARE THE COORDINATES OF THE POINTS OF
!   INTERPOLATION.
!
!   M1 IS THE SLOPE AT (P1,P2).
!
!   M2 IS THE SLOPE AT (Q1,Q2)
!
!   NCASE CONTROLS THE NUMBER AND LOCATION OF THE KNOTS.
!
!
! ON OUTPUT--
!
!   (V1,V2),(W1,W2),(Z1,Z2), AND (E1,E2) ARE THE COORDINATES OF THE
!   KNOTS AND OTHER PARAMETERS OF THE SPLINE ON (P1,Q1).  (E1,E2)
!   AND (Y1,Y2) ARE USED ONLY IF NCASE=4.
!
! AND
!
!   CASES DOES NOT ALTER P1,P2,M1,M2,Q1,Q2.
!
!-----------------------------------------------------------------------
!
      IF ((NCASE.EQ.3) .OR. (NCASE.EQ.4)) GO TO 20
      IF (NCASE.EQ.2) GO TO 10
!
! CALCULATE THE PARAMETERS FOR CASE 1.
      Z1 = (P2-Q2+M2*Q1-M1*P1)/(M2-M1)
      ZTWO = P2 + M1*(Z1-P1)
      V1 = (P1+Z1)/2.E0
      V2 = (P2+ZTWO)/2.E0
      W1 = (Z1+Q1)/2.E0
      W2 = (ZTWO+Q2)/2.E0
      Z2 = V2 + ((W2-V2)/(W1-V1))*(Z1-V1)
      RETURN
!
! CALCULATE THE PARAMETERS FOR CASE 2.
   10 Z1 = (P1+Q1)/2.E0
      V1 = (P1+Z1)/2.E0
      V2 = P2 + M1*(V1-P1)
      W1 = (Z1+Q1)/2.E0
      W2 = Q2 + M2*(W1-Q1)
      Z2 = (V2+W2)/2.E0
      RETURN
!
! CALCULATE THE PARAMETERS USED IN BOTH CASES 3 AND 4.
   20 C1 = P1 + (Q2-P2)/M1
      D1 = Q1 + (P2-Q2)/M2
      H1 = 2.E0*C1 - P1
      J1 = 2.E0*D1 - Q1
      MBAR1 = (Q2-P2)/(H1-P1)
      MBAR2 = (P2-Q2)/(J1-Q1)
!
      IF (NCASE.EQ.4) GO TO 50
!
! CALCULATE THE PARAMETERS FOR CASE 3.
      K1 = (P2-Q2+Q1*MBAR2-P1*MBAR1)/(MBAR2-MBAR1)
      IF (ABS(M1).GT.ABS(M2)) GO TO 30
      Z1 = (K1+Q1)/2.E0
      GO TO 40
   30 Z1 = (K1+P1)/2.E0
   40 V1 = (P1+Z1)/2.E0
      V2 = P2 + M1*(V1-P1)
      W1 = (Q1+Z1)/2.E0
      W2 = Q2 + M2*(W1-Q1)
      Z2 = V2 + ((W2-V2)/(W1-V1))*(Z1-V1)
      RETURN
!
! CALCULATE THE PARAMETERS FOR CASE 4.
   50 Y1 = (P1+C1)/2.E0
      V1 = (P1+Y1)/2.E0
      V2 = M1*(V1-P1) + P2
      Z1 = (D1+Q1)/2.E0
      W1 = (Q1+Z1)/2.E0
      W2 = M2*(W1-Q1) + Q2
      MBAR3 = (W2-V2)/(W1-V1)
      Y2 = MBAR3*(Y1-V1) + V2
      Z2 = MBAR3*(Z1-V1) + V2
      E1 = (Y1+Z1)/2.E0
      E2 = MBAR3*(E1-V1) + V2
      RETURN
!
      END SUBROUTINE CASES


      SUBROUTINE  UVIP3P(NP,ND,XD,YD,NI,XI, YI)
!
! Univariate Interpolation (Improved Akima Method)
!
! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 89/07/04
!
! This subroutine performs univariate interpolation.  It is based
! on the improved A method developed by Hiroshi Akima, 'A method
! of univariate interpolation that has the accuracy of a third-
! degree polynomial,' ACM TOMS, vol. xx, pp. xxx-xxx, 19xx.  (The
! equation numbers referred to in the comments below are those in
! the paper.)
!
! In this method, the interpolating function is a piecewise
! function composed of a set of polynomials applicable to
! successive intervals of the given data points.  This method
! uses third-degree polynomials as the default, but the user has
! an option to use higher-degree polynomial to reduce undulations
! in resulting curves.
!
! This method has the accuracy of a third-degree polynomial if
! the degree of the polynomials for the interpolating function is
! set to three.
!
! The input arguments are
!   NP = degree of the polynomials for the interpolating
!        function,
!   ND = number of input data points
!        (must be equal to 2 or greater),
!   XD = array of dimension ND, containing the abscissas of
!        the input data points
!        (must be in a monotonic increasing order),
!   YD = array of dimension ND, containing the ordinates of
!        the input data points,
!   NI = number of points for which interpolation is desired
!        (must be equal to 1 or greater),
!   XI = array of dimension NI, containing the abscissas of
!        the desired points.
!
! The output argument is
!   YI = array of dimension NI, where the ordinates of the
!        desired points are to be stored.
!
! If an integer value smaller than 3 is given to the NP argument,
! this subroutine assumes NP = 3.
!
! The XI array elements need not be monotonic, but this
! subroutine interpolates faster if the XI array elements are
! given in a monotonic order.
!
! If the XI array element is less than XD(1) or greater than
! XD(ND), this subroutine linearly interpolates the YI value.
!
!
! Specification statement
      REAL, INTENT(IN) :: XD(:), YD(:), XI(:)
      INTEGER, INTENT(IN) :: NP, ND
      REAL, INTENT(OUT) :: YI(:)
! Error check
      IF (ND.LE.1)   GO TO 90
      IF (NI.LE.0)   GO TO 91
      icode = 0
      DO ID=2,ND
        IF (XD(ID).LE.XD(ID-1))     GO TO 92
      ENDDO
! Branches off special cases.
      IF (ND.LE.4)   GO TO 50
! General case  --  Five data points of more
! Calculates some local variables.   20 NP0=MAX(3,NP)
      NP0=MAX(3,NP)
      NPM1=NP0-1
      RENPM1=NPM1
      RENNM2=NP0*(NP0-2)
! Main calculation for the general case
! First (outermost) DO-loop with respect to the desired points
!   30 DO 39  II=1,NI
      DO 39 II=1,NI
        IF (II.EQ.1)      IINTPV=-1
        XII=XI(II)
! Locates the interval that includes the desired point by binary
! search.
        IF (XII.LE.XD(1))  THEN
          IINT=0
        ELSE IF (XII.LT.XD(ND))  THEN
          IDMN=1
          IDMX=ND
          IDMD=(IDMN+IDMX)/2
   31     IF (XII.GE.XD(IDMD))  THEN
            IDMN=IDMD
          ELSE
            IDMX=IDMD
          END IF
          IDMD=(IDMN+IDMX)/2
          IF (IDMD.GT.IDMN)    GO TO 31
          IINT=IDMD
        ELSE
          IINT=ND
        END IF
! End of locating the interval of interest
! Interpolation or extrapolation in one of the three subcases
        IF (IINT.LE.0)  THEN
! Subcase 1  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the first data
!                point or less.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XD(1)
            X1=XD(2)-X0
            X2=XD(3)-X0
            X3=XD(4)-X0
            Y0=YD(1)
            Y1=YD(2)-Y0
            Y2=YD(3)-Y0
            Y3=YD(4)-Y0
            DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
            A1=(((X2*X3)**2)*(X3-X2)*Y1 &
               +((X3*X1)**2)*(X1-X3)*Y2 &
               +((X1*X2)**2)*(X2-X1)*Y3)/DLT
          END IF
! Evaluates the YI value.
          YI(II)=Y0+A1*(XII-X0)
! End of Subcase 1
        ELSE IF (IINT.GE.ND)  THEN
! Subcase 2  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the last data
!                point or greater.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XD(ND)
            X1=XD(ND-1)-X0
            X2=XD(ND-2)-X0
            X3=XD(ND-3)-X0
            Y0=YD(ND)
            Y1=YD(ND-1)-Y0
            Y2=YD(ND-2)-Y0
            Y3=YD(ND-3)-Y0
            DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
            A1=(((X2*X3)**2)*(X3-X2)*Y1 & 
               +((X3*X1)**2)*(X1-X3)*Y2 &
               +((X1*X2)**2)*(X2-X1)*Y3)/DLT
          END IF
! Evaluates the YI value.
          YI(II)=Y0+A1*(XII-X0)
! End of Subcase 2
        ELSE
! Subcase 3  --  Interpolation when the abscissa of the desired
!                point is  between those of the first and last
!                data points.
! Calculates the coefficients of the third-degree polynomial (for
! NP.LE.3) or the factors for the higher-degree polynomials (for
! NP.GT.3), when the interval is not the same as the one for the
! previous desired point.
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
! The second DO-loop with respect to the two endpoints of the
! interval
            DO 37  IEPT=1,2
! Calculates the estimate of the first derivative at an endpoint.
! Initial setting for calculation
              ID0=IINT+IEPT-1
              X0=XD(ID0)
              Y0=YD(ID0)
              SMPEF=0.0
              SMWTF=0.0
              SMPEI=0.0
              SMWTI=0.0
! The third (innermost) DO-loop with respect to the four primary
! estimate of the first derivative
              DO 36  IPE=1,4
! Selects point numbers of four consecutive data points for
! calculating the primary estimate of the first derivative.
                IF (IPE.EQ.1)  THEN
                  ID1=ID0-3
                  ID2=ID0-2
                  ID3=ID0-1
                ELSE IF (IPE.EQ.2)  THEN
                  ID1=ID0+1
                ELSE IF (IPE.EQ.3)  THEN
                  ID2=ID0+2
                ELSE
                  ID3=ID0+3
                END IF
! Checks if any point number falls outside the legitimate range
! (between 1 and ND).  Skips calculation of the primary estimate
! if any does.
                IF (ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR. &
                    ID1.GT.ND.OR.ID2.GT.ND.OR.ID3.GT.ND) &
                      GO TO 36
! Calculates the primary estimate of the first derivative  --
! cf. Equation (8)
                X1=XD(ID1)-X0
                X2=XD(ID2)-X0
                X3=XD(ID3)-X0
                Y1=YD(ID1)-Y0
                Y2=YD(ID2)-Y0
                Y3=YD(ID3)-Y0
                DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
                PE=(((X2*X3)**2)*(X3-X2)*Y1 &
                   +((X3*X1)**2)*(X1-X3)*Y2 &
                   +((X1*X2)**2)*(X2-X1)*Y3)/DLT
! Calculates the volatility factor, VOL, and distance factor,
! SXX, for the primary estimate.  --  cf. Equations (9) and (11)
                SX=X1+X2+X3
                SY=Y1+Y2+Y3
                SXX=X1*X1+X2*X2+X3*X3
                SXY=X1*Y1+X2*Y2+X3*Y3
                DNM=4.0*SXX-SX*SX
                B0=(SXX*SY-SX*SXY)/DNM
                B1=(4.0*SXY-SX*SY)/DNM
                DY0=-B0
                DY1=Y1-(B0+B1*X1)
                DY2=Y2-(B0+B1*X2)
                DY3=Y3-(B0+B1*X3)
                VOL=DY0*DY0+DY1*DY1+DY2*DY2+DY3*DY3
! Calculates the EPSLN value, which is used to decide whether or
! not the volatility factor, VOL, is essentially zero.
                EPSLN=(YD(ID0)**2+YD(ID1)**2 &
                      +YD(ID2)**2+YD(ID3)**2)*1.0E-12
! Accumulates the weighted primary estimates.  --
! cf. Equations (13) and (14)
                IF (VOL.GT.EPSLN)  THEN
! - For finite weight.
                  WT=1.0/(VOL*SXX)
                  SMPEF=SMPEF+PE*WT
                  SMWTF=SMWTF+WT
                ELSE
! - For infinite weight.
                  SMPEI=SMPEI+PE
                  SMWTI=SMWTI+1.0
                END IF
   36         CONTINUE
! End of the third DO-loop
! Calculates the final estimate of the first derivative.  --
! cf. Equation (14)
              IF (SMWTI.LT.0.5)  THEN
! - When no infinite weights exist.
                YP=SMPEF/SMWTF
              ELSE
! - When infinite weights exist.
                YP=SMPEI/SMWTI
              END IF
              IF (IEPT.EQ.1)  THEN
                YP0=YP
              ELSE
                YP1=YP
              END IF
! End of the calculation of the estimate of the first derivative
! at an endpoint
   37       CONTINUE
! End of the second DO-loop
            IF (NP0.LE.3)  THEN
! Calculates the coefficients of the third-degree polynomial
! (when NP.LE.3).  --  cf. Equation (4)
              DX=XD(IINT+1)-XD(IINT)
              DY=YD(IINT+1)-YD(IINT)
              A0=YD(IINT)
              A1=YP0
              YP1=YP1-YP0
              YP0=YP0-DY/DX
              A2=-(3.0*YP0+YP1)/DX
              A3= (2.0*YP0+YP1)/(DX*DX)
            ELSE
! Calculates the factors for the higher-degree polynomials
! (when NP.GT.3).  --  cf. Equation (20)
              DX=XD(IINT+1)-XD(IINT)
              DY=YD(IINT+1)-YD(IINT)
              T0=YP0*DX-DY
              T1=YP1*DX-DY
              AA0= (T0+RENPM1*T1)/RENNM2
              AA1=-(RENPM1*T0+T1)/RENNM2
            END IF
          END IF
! End of the calculation of the coefficients of the third-degree
! polynomial (when NP.LE.3) or the factors for the higher-degree
! polynomials (when NP.GT.3), when the interval is not the same
! as the one for the previous desired point.
! Evaluates the YI value.
          IF (NP0.LE.3)  THEN
! - With a third-degree polynomial.  --  cf. Equation (3)
            XX=XII-XD(IINT)
            YI(II)=A0+XX*(A1+XX*(A2+XX*A3))
          ELSE
! - With a higher-degree polynomial.  --  cf. Equation (19)
            U=(XII-XD(IINT))/DX
            UC=1.0-U
            V=AA0*((U**NP0)-U)+AA1*((UC**NP0)-UC)
            YI(II)=YD(IINT)+DY*U+V
          END IF
! End of Subcase 3
        END IF
   39 CONTINUE
! End of the first DO-loop
! End of general case
      RETURN
! Special cases  --  Four data points or less
! Preliminary processing for the special cases
   50 X0=XD(1)
      Y0=YD(1)
      X1=XD(2)-X0
      Y1=YD(2)-Y0
      IF (ND.EQ.2)   GO TO 60
      X2=XD(3)-X0
      Y2=YD(3)-Y0
      IF (ND.EQ.3)   GO TO 70
      X3=XD(4)-X0
      Y3=YD(4)-Y0
      GO TO 80
! Special Case 1  --  Two data points
! (Linear interpolation and extrapolation)
   60 A1=Y1/X1
      DO 61  II=1,NI
        YI(II)=Y0+A1*(XI(II)-X0)
   61 CONTINUE
! End of Special Case 1
      RETURN
! Special Case 2  --  Three data points
! (Quadratic interpolation and linear extrapolation)
   70 DLT=X1*X2*(X2-X1)
      A1=(X2*X2*Y1-X1*X1*Y2)/DLT
      A2=(X1*Y2-X2*Y1)/DLT
      A12=2.0*A2*X2+A1
      DO 71  II=1,NI
        XX=XI(II)-X0
        IF (XX.LE.0.0)  THEN
          YI(II)=Y0+A1*XX
        ELSE IF (XX.LT.X2) THEN
          YI(II)=Y0+XX*(A1+XX*A2)
        ELSE
          YI(II)=Y0+Y2+A12*(XX-X2)
        END IF
   71 CONTINUE
! End of Special Case 2
      RETURN
! Special Case 3  --  Four data points
! (Cubic interpolation and linear extrapolation)
   80 DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
      A1=(((X2*X3)**2)*(X3-X2)*Y1 &
         +((X3*X1)**2)*(X1-X3)*Y2 &
         +((X1*X2)**2)*(X2-X1)*Y3)/DLT
      A2=(X2*X3*(X2*X2-X3*X3)*Y1 &
         +X3*X1*(X3*X3-X1*X1)*Y2 &
         +X1*X2*(X1*X1-X2*X2)*Y3)/DLT
      A3=(X2*X3*(X3-X2)*Y1 &
         +X3*X1*(X1-X3)*Y2 &
         +X1*X2*(X2-X1)*Y3)/DLT
      A13=(3.0*A3*X3+2.0*A2)*X3+A1
      DO 81  II=1,NI
        XX=XI(II)-X0
        IF (XX.LE.0.0)  THEN
          YI(II)=Y0+A1*XX
        ELSE IF (XX.LT.X3) THEN
          YI(II)=Y0+XX*(A1+XX*(A2+XX*A3))
        ELSE
          YI(II)=Y0+Y3+A13*(XX-X3)
        END IF
   81 CONTINUE
! End of Special Case 3
      RETURN
! Error exit
   90 WRITE (*,99090) ND
      GO TO 99
   91 WRITE (*,99091) NI
      GO TO 99
   92 WRITE (*,99092) ID,XD(ID-1),XD(ID)
      icode = ID
   99 WRITE (*,99099)
      RETURN
! Format statements for error messages
99090 FORMAT (1X/,' ***   Insufficient data points.', &
        7X,'ND =',I3)
99091 FORMAT (1X/,' ***   No desired points.', &
        7X,'NI =',I3)
99092 FORMAT (1X/,' ***   Two data points identical or out of ', &
        'sequence.'/ &
        7X,'ID, XD(ID-1), XD(ID) =',I5,2G20.8)
99099 FORMAT (' Error detected in the UVIP3P subroutine'/)
      END SUBROUTINE UVIP3P

END MODULE INTERPOL
