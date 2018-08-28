! Copyright or © or Copr. Bernard Legras (1984,2006)

! legras@lmd.ens.f

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


module sphereharmSPE

  IMPLICIT NONE 
!===============================================================================
! PARAMETRISATION DE LA TRONCATURE.
!==============================================================================
      INTEGER N_TRUNC, NH_LAT, NH_LONG, N_LAT, N_LONG, N2_LONG
      INTEGER NH1_LONG, N1, N2, NI, NP, NDIM, NCPLX, NB_DEGLIB, N_SIZE, N_C1S
!     PARAMETER ( N_TRUNC = 106, NH_LAT = 80)
!     ERA-40 truncature and regular grid including redundant equator
!     PARAMETER ( N_TRUNC = 159, NH_LAT = 121)  ! T159 regular grid
!     PARAMETER ( N_TRUNC = 159, NH_LAT = 80)   ! TL159
!     PARAMETER ( N_TRUNC = 106, NH_LAT = 256)
!     PARAMETER ( N_TRUNC = 170, NH_LAT = 128)
!     PARAMETER ( N_TRUNC = 213, NH_LAT = 160)
      PARAMETER ( N_TRUNC = 255, NH_LAT = 161)  ! T255 regular grid
!     PARAMETER ( N_TRUNC = 255, NH_LAT = 128)  ! TL255
!     PARAMETER ( N_TRUNC = 319, NH_LAT = 240)
!===============================================================================
! PARAMETRISATION DE LA TRONCATURE (AUTOMATISE).
!==============================================================================
!     PARAMETER ( NH_LONG = 240)
      PARAMETER ( NH_LONG = 256)
      PARAMETER ( N_LAT = 2*NH_LAT, N_LONG = 2*NH_LONG)  
      PARAMETER ( NH1_LONG = NH_LONG + 1 )
      PARAMETER ( N1 = N_TRUNC/2, N2 = (N_TRUNC+1)/2)
      PARAMETER ( NP = (N1+1) * (N2+1), NI = (N1+1) * N2)
      PARAMETER ( NDIM = NP + NI)
      PARAMETER ( NCPLX   =    2, NB_DEGLIB = (N_TRUNC + 1) * (N_TRUNC + 2) / 2)
      PARAMETER ( N2_LONG  = N_LONG + 2, N_SIZE = N2_LONG * N_LAT)
      PARAMETER ( N_C1S = NCPLX)

!
!   --PARAMETRES PRINCIPAUX
!       N_TRUNC : nombre d'onde maximum de la troncature triangulaire
!   --PARAMETRES DERIVES
!       N_LAT   : nombre de latitudes de la grille de collocation
!       NH_LAT  : nombre de latitudes de Gauss dans chaque hemisphere
!                 1/2 N_LAT
!       N_LONG  : nombre de longitudes de la grille de collocation
!       NH_LONG : ordre de la FFT complexe
!       NH1 LONG: longueur des tableaux complexes de la FFT
!       N2_LONG : longueur des tableaux reels de la FFT
!       NB_DEGLIB: nombre de degres de liberte en troncature triangulaire
!                  standard
!       NBE DEGLIB: nombre de degres de libertes en troncature etendue des
!                  fonctions de Legendre associees
!       N_C1S   : increment pour MXMA dans les transformees de Legendre
!
!=========================================================================
! PARAMETRISATION PHYSIQUE ET MATHEMATIQUES (mksa).
! (NB: Taken from Gill ( Atmosphere-Ocean Dynamics) )
!=========================================================================
! MATHEMATICAL CONSTANTS (PI called PIN to avoid interferences)
      REAL PIN, ANGLE, PI2
      PARAMETER ( PIN = 3.141592653589793238 , ANGLE = 2.0 * PIN / N_LONG, &
                  PI2=2.0*PIN) 
! PHYSICAL CONSTANTS
     
! EARTH RADIUS (m)
      REAL RAYTER
      PARAMETER(RAYTER=6.371E6)
! ACCELERATION DUE TO GRAVITATY (m/s2)
      REAL G
      PARAMETER(G=9.81)
! ROTATION RATE OF EARTH (1/s)
      REAL OMEGA, OMEGA2
      PARAMETER(OMEGA = 7.292E-5, OMEGA2 = 2.0*OMEGA)
! GAMMA AND RATIO OF SPECIFIC HEATS FOR DRY AIR
      REAL GAMAIR, RCPAIR
      PARAMETER(GAMAIR = 1.4, RCPAIR = (GAMAIR - 1.0) / GAMAIR)
! EARTH SOLID ROTATION SPEED AT THE EQUATOR (m/s)
      REAL EARTHSP
      PARAMETER(EARTHSP=RAYTER*OMEGA)
! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR (J/kg/K)
      REAL CPAIR, RAIR, CVAIR,T0C,PRESEA, KappaAir
      PARAMETER(CPAIR=1005.)
! THE GAZ CONSTANT FOR DRY AIR (J/kg/K)
      PARAMETER(RAIR=RCPAIR*CPAIR)
! SPECIFIC HEAT AT CONSTANT VOLUME FOR DRY AIR (J/kg/K)
      PARAMETER(CVAIR=CPAIR-RAIR)
! RATIO BETWEEN R AND CP
      PARAMETER(KappaAir=RAIR/CPAIR)
! ZERO DEGREE CELSIUS IN KELVIN (K)
      PARAMETER(T0C=273.15)
! PRESSURE OF 1 BAR (N/m2)
      PARAMETER(PRESEA=1.E5)
! MASS OF EARTH (1.e18 kg)
      REAL MASSEAR,MASSATM,MASSOCE,STEFC
      PARAMETER(MASSEAR=5.977E6)
! MASS OF THE ATMOSPHERE ( 1e18 kg)
      PARAMETER(MASSATM=5.3)
! MASS OF OCEAN (1e18 kg)
      PARAMETER(MASSOCE=1.4E3)
! STEFAN'S CONSTANT (W/m2/K4)
      PARAMETER(STEFC=5.67E-8)

!========================================================================
! COMMON BLOCK POUR LES OPERATIONS DIRECTES ET INVERSES SPECTRALES.
!========================================================================
      REAL RAYDEF_IV
      PARAMETER(RAYDEF_IV = 0.)

! -- CONTROL :
!       M_CTRL  : impression de divers resultats intermediaires
!       M_HEMIS : choix du mode hemispherique (1) ou spherique (2)
!                 initialise par le MASTER
!       I_ERR   : NUMERO D'ERREUR DANS LES ENTREES SORTIES
!       M_MAX   : determine le passage du mode MXMA au mode SDOT dans
!                 les transformees de Legendre.
!                 initialise par le MASTER
!       M_TLG   : choix de la transformee de Legendre dans les tests
!       M_CORIOL:
!
      integer M_CTRL,M_PRINT,M_ESSAI,M_HEMIS,M_INIT,M_TLG,IERR,M_MAX,M_CORIOL

!
! -- TRCVERT
!    Description de la troncature verticale
!
!       MVS_DEB  : index de l'element diagonal d'ordre JMOD (dans un tableau
!                  complet de dimension NB_DEGLIB)
!       MVA_DEB  : index du premier element antisym. d'ordre JMOD  : diag + 1
!       MVS_LG   : nombre d'elements symetriques dans la colonne d'ordre JMOD
!       MVA_LG   : nombre d'elements antisym. dans la colonne d'ordre JMOD
!

      INTEGER MVA_DEB(0:N_TRUNC), MVA_LG(0:N_TRUNC)
      INTEGER MVS_DEB(0:N_TRUNC), MVS_LG(0:N_TRUNC)

      REAL            WORK (4*N_LONG*N_LAT)
      REAL            TRIGS(NH_LONG * 3)
      INTEGER         IFAX(10)

!       PLM : FONCTIONS DE LEGENDRE
      REAL            PLM(NB_DEGLIB,NH_LAT)              

! -- INITIALISE PAR * I N I T   C T R L *
!
!         operateurs de Laplace et de Laplace inverse
!         en composantes spectrales
!        operateur spectral de derivation en longitude

      REAL            OPR_LAP(NB_DEGLIB,2), OPR_LAP_IV
      REAL            OPR_M (NB_DEGLIB), OPRTRUNC(NB_DEGLIB)
  
! temporay fields
! to eliminate after replacement by allocations where they are used

      REAL           Z_FS (NCPLX, 0:NH_LONG, NH_LAT),&
                     Z_FA (NCPLX, 0:NH_LONG, NH_LAT)

  INTERFACE LAPLACE
     MODULE PROCEDURE LAPLACE_R, LAPLACE_C
  END INTERFACE  
  INTERFACE spectogrid
     MODULE PROCEDURE spectogrid_r, spectogrid_c
  END INTERFACE 
 
  CONTAINS

!***********************************************************************
!                                                           INIT CTRL
!***********************************************************************
!
      SUBROUTINE INIT_CTRL 

      INTEGER JMOD

!   --OPTION DE MISE AU POINT ( PAR DEFAUT )
!
!-----------------------------------------------------------------------
!                                                      INIT CTRL.1
!--------------------------------------------INITIALISATION DES TRCCTRLS
!
      MVA_DEB(0) = 1
      MVS_DEB(0) = NI + 1
      MVA_LG(0)  = 1 + (N_TRUNC - 1)/2
      MVS_LG(0)  = 1 + N_TRUNC/2
!
      DO JMOD = 1, N_TRUNC
        MVA_DEB(JMOD) = MVA_DEB(JMOD-1) + MVA_LG(JMOD-1)
        MVS_DEB(JMOD) = MVS_DEB(JMOD-1) + MVS_LG(JMOD-1)
        MVA_LG (JMOD) = 1 + (N_TRUNC - JMOD - 1)/2
        MVS_LG (JMOD) = 1 + (N_TRUNC - JMOD)/2
      ENDDO
      MVA_LG(N_TRUNC) = 0
!
!
!-----------------------------------------------------------------------
!                                                       INIT CTRL.3
!----------------------------------------------INITIALISATION DE CONTROL
      M_CTRL = 0
      M_PRINT = 0
      M_CORIOL = 0
      M_INIT   = 0
      M_HEMIS  = 2
      M_MAX    = 11

      RETURN
      END subroutine init_ctrl

!***********************************************************************
!                                                      INIT TRSF SHORT
!***********************************************************************
!
      SUBROUTINE INIT_TRSF_SHORT() 
    
      INTEGER JMOD, IS, IN, JD, LN, JN
      REAL RAYM1, RAYM2, RAYM3

!*    3.  Initialisation de COMFFT
!     ----------------------------

      call fftfax(n_long,ifax,trigs)

!*    4.  Initialisation de OPERA
!     ---------------------------

      RAYM1 = 1. / RAYTER
      RAYM2 = 1. / ( RAYTER * RAYTER )
      RAYM3 = RAYDEF_IV ** 2
      DO JMOD = 0, N_TRUNC
        IS = MVS_DEB(JMOD)
        IN = MVS_LG (JMOD)
        JD = 0
        IF (JMOD.EQ.0) JD = 1
        DO JN = JD, IN - 1
          LN = JMOD + 2 * JN
          OPR_LAP   (IS + JN, 2) = - FLOAT (LN * (LN+1)) * RAYM2
          OPR_LAP   (IS + JN, 1) = OPR_LAP(IS + JN, 2) - RAYM3
          OPR_M     (IS + JN) = FLOAT(JMOD) * RAYM1
        ENDDO
        IS = MVA_DEB(JMOD)
        IN = MVA_LG (JMOD)
        IF(IN==0) CYCLE
        JD = 0
        DO JN = JD, IN - 1
          LN = JMOD + 2 * JN + 1
          OPR_LAP   (IS + JN, 2) = - FLOAT (LN * (LN+1)) * RAYM2
          OPR_LAP   (IS + JN, 1) = OPR_LAP(IS + JN, 2) - RAYM3
          OPR_M     (IS + JN) = FLOAT(JMOD) * RAYM1
        ENDDO
      ENDDO
      OPR_LAP (NI+1,1) = 1.
      OPR_LAP (NI+1,2) = 1.
      OPR_M (NI+1) = 0.

      RETURN
      END subroutine init_trsf_short

!***********************************************************************
!                                                      SPECTOGRID
!***********************************************************************

      subroutine spectogrid_r(grf,spf)

!  Passage en points de grille spherique

      real, intent(out) :: grf(n2_long,n_lat)
      real, intent(in) ::  spf(ncplx, nb_deglib)

      call tlg_inv_sn (grf,spf,plm,1)
      call rfftmlt(grf,work,trigs,ifax,1,n2_long,n_long,n_lat,1)
       
      return
      end subroutine spectogrid_r
      
      subroutine spectogrid_c(grf,spf)
      real, intent(out) :: grf(n2_long,n_lat)
      complex, intent(in) :: spf(nb_deglib)
      real, allocatable :: spfr(:,:)
      allocate(spfr(ncplx,nb_deglib))
      spfr(1,:)=real(spf(:))
      spfr(2,:)=aimag(spf(:))
      call tlg_inv_sn (grf,spfr,plm,1)
      call rfftmlt(grf,work,trigs,ifax,1,n2_long,n_long,n_lat,1)
      return
      end subroutine spectogrid_c

!******************************************************************************
!                                                           LAPLACE
!******************************************************************************
!
      SUBROUTINE LAPLACE_R (SP_OUT,SP_IN,IDIR) 
!
!**** LAPLACE  laplacien spectral
!
!     B. Legras      LMD        25-05-84
!
!     Objet: calcule le laplacien ou resout l'equation de Poisson en spectral
!
!     Arguments:
!                SP_IN  : champ spectral en entree
!                SP_OUT : champ spectral en sortie, peut etre identique a
!                         SP_IN
!                IDIR   : choix de l'operation
!                         IDIR = 2 : laplacien direct
!                         IDIR =-2 : laplacien inverse(resolution de Poisson)
!                         IDIR = 1 : laplacien modifie direct
!                         IDIR =-1 : laplacien modifie inverse
!
!     Note: Les valeurs propres du laplacien sont rangees dans un tableau
!           de dimension NB_DEGLIB. En mode hemispherique, on n'utilise
!           que la premiere partie des tableaux (de longueur NI).
!
!******************************************************************************

      REAL, INTENT(IN)  ::  SP_IN (NCPLX,NB_DEGLIB)
      REAL, INTENT(OUT) ::  SP_OUT(NCPLX,NB_DEGLIB)
      INTEGER IDIR
      INTEGER JC
!
!*    1.  Laplacien direct
!     --------------------
!
      IF(IDIR.GT.0) THEN
       DO JC = 1, NI + (M_HEMIS-1) * NP
         SP_OUT(1,JC) = SP_IN(1,JC) * OPR_LAP(JC, IDIR)
         SP_OUT(2,JC) = SP_IN(2,JC) * OPR_LAP(JC, IDIR)
       ENDDO
       RETURN
!
!*    2.  Laplacien inverse
!     ---------------------
!
      ELSE
       DO JC = 1, NI + (M_HEMIS-1) * NP
         SP_OUT(1,JC) = SP_IN(1,JC) / OPR_LAP(JC, -IDIR)
         SP_OUT(2,JC) = SP_IN(2,JC) / OPR_LAP(JC, -IDIR)
       ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE LAPLACE_R

      SUBROUTINE LAPLACE_C (SP_OUT,SP_IN,IDIR)
      COMPLEX, INTENT(IN)  :: SP_IN (NB_DEGLIB)
      COMPLEX, INTENT(OUT) :: SP_OUT (NB_DEGLIB)
      INTEGER, INTENT(IN) :: IDIR
      IF(IDIR>0) THEN
        SP_OUT(1:NI+(M_HEMIS-1)*NP)=SP_IN(1:NI+(M_HEMIS-1)*NP) &
                                * OPR_LAP(1:NI+(M_HEMIS-1)*NP,IDIR)
      ELSE
        SP_OUT(1:NI+(M_HEMIS-1)*NP)=SP_IN(1:NI+(M_HEMIS-1)*NP) &
                                / OPR_LAP(1:NI+(M_HEMIS-1)*NP,-IDIR)
      ENDIF
      RETURN
      END SUBROUTINE LAPLACE_C
                                        
!***************************************************************************
!                                                        CALPSI
!***************************************************************************
      SUBROUTINE CALPSI(PSI,KHI,U,V,MTRUNC)

!     --------------------------------------------------
!     DONNE UCOSPHI ET VCOSPHI ET RESSORT
!     PSI ET KHI.
!     ---------------------------------------------------
!     LTF ,LTM NB D'ONDES ZONAUX+1
!     LTN NB D'ONDE TOTAL+1

      COMPLEX, INTENT(IN)  :: U(NI+NP),V(NI+NP)
      INTEGER, INTENT(IN)  :: MTRUNC
      COMPLEX, INTENT(OUT) :: PSI(NI+NP),KHI(NI+NP)
      COMPLEX UM(0:N_TRUNC),VM(0:N_TRUNC),KM(0:N_TRUNC+1),PM(0:N_TRUNC+1)        
      COMPLEX A(2,2, 0:N_TRUNC+1),B(2,0:N_TRUNC+1)                              
      REAL D(0:N_TRUNC+1)
      INTEGER I, M, N, NMAX                                  
!C     BOUCLE SUR M:ON PRESENTE LES COMPO DE U ET V A M DONNE
!      AVEC N CROISSANT
      DO I=1,NI+NP
        KHI(I)=0.
        PSI(I)=0.
      ENDDO
      DO M=0,MTRUNC
        DO N=M,MTRUNC
          IF(MOD(N-M,2).EQ.0)THEN
            UM(N-M)=U(MVS_DEB(M)+(N-M)/2)*rayter
            VM(N-M)=V(MVS_DEB(M)+(N-M)/2)*rayter
          ELSE
            UM(N-M)=U(MVA_DEB(M)+(N-M-1)/2)*rayter
            VM(N-M)=V(MVA_DEB(M)+(N-M-1)/2)*rayter
          ENDIF
        ENDDO
        NMAX=MTRUNC+1-M
        CALL RESOL(UM,VM,PM,KM,D,A,B,M+1,NMAX)                              
!
!     REMISE EN PLACE
        DO N=M,MTRUNC
          IF(MOD(N-M,2).EQ.0)THEN
            PSI(MVS_DEB(M)+(N-M)/2)=PM(N-M)
            KHI(MVS_DEB(M)+(N-M)/2)=KM(N-M)
          ELSE
            PSI(MVA_DEB(M)+(N-M-1)/2)=PM(N-M)
            KHI(MVA_DEB(M)+(N-M-1)/2)=KM(N-M)
          ENDIF
        ENDDO
      ENDDO                                                        
      RETURN
      END SUBROUTINE CALPSI

!*****************************************************************************
!                                                                 CALVENT
!*****************************************************************************
                                                                 
      SUBROUTINE CALVENT(U,V,PSI,KHI,MTRUNC)
      COMPLEX, INTENT(OUT) :: U(NI+NP),V(NI+NP)
      COMPLEX, INTENT(IN)  :: PSI(NI+NP),KHI(NI+NP)
      INTEGER MTRUNC
      INTEGER I, M, N
      REAL D(0:N_TRUNC+1)
      COMPLEX UM(0:N_TRUNC),VM(0:N_TRUNC),PM(0:N_TRUNC+1),KM(0:N_TRUNC+1)

!    ----------------------------------------------
!     CALCULE UCOSTETA ET VCOSTETA PAR
!     UMN=(N-1)DMNPSIMN-1 - (N-2)DMN+1PSIMN+1 + IMKHIMN
!     VMN=-(N-1)DMNKHIMN-1 + (N+2)DMN+1KHIMN+1 + IMPSIMN
!     LE CALCUL DE U A LA TRONCATURE MTRUNC NECESSITE PSI EN MTRUNC+1
!     DONC PAS D'APPEL AVEC MTRUNC.GE.N_TRUNC-1
!     --------------------------------------------------
!     IF(MTRUNC.GE.N_TRUNC)THEN
!       WRITE(*,*)'APPEL CALVENT AVEC MTRUNC.GE.N_TRUNC'
!     STOP
!     ENDIF
      DO I=1,NI+NP
        U(I)=0.
        V(I)=0.
      ENDDO
      DO I=0,N_TRUNC+1
        PM(I)=0.
        KM(I)=0.
        D(I)=0.
      ENDDO
!
!     BOUCLE SUR M
!
      DO M=0,MTRUNC
!       RECOPIE DANS PM ET KM
!
      DO  N=M,MTRUNC
      IF(MOD(N-M,2).EQ.0)THEN
         PM(N)=PSI(MVS_DEB(M)+(N-M)/2)/RAYTER
         KM(N)=KHI(MVS_DEB(M)+(N-M)/2)/RAYTER
      ELSE
         PM(N)=PSI(MVA_DEB(M)+(N-M-1)/2)/RAYTER
         KM(N)=KHI(MVA_DEB(M)+(N-M-1)/2)/RAYTER
      ENDIF
      
      ENDDO
     
      PM(MTRUNC+1)=0.
      KM(MTRUNC+1)=0.
!       CALCUL DES DM
      D(M)=0.
      DO  N=M+1,MTRUNC+1
         D(N)=SQRT((N*N-M*M)/(4.*N*N-1.))
      ENDDO

!     CALCUL DES UM ET VM
      UM(M)=-(M+2)*D(M+1)*PM(M+1)+(0.,1.)*M*KM(M)
      VM(M)= (M+2)*D(M+1)*KM(M+1)+(0.,1.)*M*PM(M)
      IF(MTRUNC < N_TRUNC) THEN
      DO N=M+1,MTRUNC
        UM(N)=(N-1)*D(N)*PM(N-1)-(N+2)*D(N+1)*PM(N+1)+(0.,1.)*M*KM(N)
        VM(N)=-(N-1)*D(N)*KM(N-1)+(N+2)*D(N+1)*KM(N+1)+(0.,1.)*M*PM(N)
      ENDDO
      ELSE
      DO N=M+1,N_TRUNC-1
        UM(N)=(N-1)*D(N)*PM(N-1)-(N+2)*D(N+1)*PM(N+1)+(0.,1.)*M*KM(N)
        VM(N)=-(N-1)*D(N)*KM(N-1)+(N+2)*D(N+1)*KM(N+1)+(0.,1.)*M*PM(N)
      ENDDO
      UM(N_TRUNC)=(N_TRUNC-1)*D(N_TRUNC)*PM(N_TRUNC-1)+(0.,1.)*M*KM(N_TRUNC)
      VM(N_TRUNC)=-(N_TRUNC-1)*D(N_TRUNC)*KM(N_TRUNC-1)+(0.,1.)*M*PM(N_TRUNC)
      ENDIF

!     RECOPIE DES UM ET VM DANS U ET V
      DO N=M,MTRUNC
      IF(MOD(N-M,2).EQ.0)THEN
         U(MVS_DEB(M)+(N-M)/2)=UM(N)
         V(MVS_DEB(M)+(N-M)/2)=VM(N)
      ELSE
         U(MVA_DEB(M)+(N-M-1)/2)=UM(N)
         V(MVA_DEB(M)+(N-M-1)/2)=VM(N)
      ENDIF
    
      ENDDO

      ENDDO
      RETURN
      END SUBROUTINE CALVENT

!****************************************************************************
!                                                          RESOL
!****************************************************************************
      SUBROUTINE RESOL(U,V,P,K,D,A,B,M,N) 

!     --------------------------------------------      
!     Obtention de la fonction de courant
!     et du potentiel de vitesse a partir des
!     composantes de la vitesse en spectral    
!     --------------------------------------------
!     SI HM=M EN NOTATION NORMALE
!     M=HM+1
!     U(1:N)=HU(HM:HMAX)
!     D(L)=HD(HM,HN) OU HN=M+L-2
!     SOIT D(L)=HD(HM,HM+L-1)
!     --------------------------------------------
!     SI ON SUPPOSE U(0:HM),V(0:HM)  CONNUS
!     LES RELATIONS A NB D'ONDE ZONAL FIXE HM
!          UMN= (N-1)DMN PMN-1 + IM KMN -(N-2)DMN+1 PMN+1
!          VMN=-(N-1)DMN KMN-1 + IM PMN + (N-2)DMN+1 KN+1
!     FOURNISSENT LE BON NB D'EQ POUR PSI ET KHI
!     SI ON SUPPOSE PMN+1 ET KMN+1 NULS POUR N=N_TRUNC
!     LA RESOLUTION S'EFFECTUE ALORS PAR DOUBLE BALAYAGE
!     LA MATRICE DU PB EST TRIDIAGONALE PAR BLOCS DE DIM 2
!     1.ON LA REND DIAG SUP
!     2.ELIMINATION
!     --------------------------------------------------
      INTEGER, INTENT(IN) :: M,N 
      COMPLEX, INTENT(IN) :: U(N),V(N)
      COMPLEX, INTENT(OUT) ::  P(N+1),K(N+1),A(2,2,0:N),B(2,0:N)
      REAL, INTENT(OUT) :: D(N+1)
      COMPLEX MAT(2,2),DET,XX,YY,HI                                     
      INTEGER L,I,J
      REAL X,Y,Z
      D(1)=0.                                                           
      DO L=1,N+1                                                      
         D(L)=SQRT(FLOAT((M+L-2)**2-(M-1)**2)/&
              FLOAT(4*(M+L-2)**2-1)) 
      ENDDO      
      IF(M.LE.1) THEN                                                
        K(1)=0.                                                           
        P(1)=0.                                                           
        K(2)=V(1)/D(2)/2.                                                 
        P(2)=-U(1)/D(2)/2.                                                
        DO L=3,N                                                        
          K(L)=( V(L-1)+FLOAT(L-3)*D(L-1)*K(L-2))/D(L)/FLOAT(L)             
          P(L)=(-U(L-1)+FLOAT(L-3)*D(L-1)*P(L-2))/D(L)/FLOAT(L)
        ENDDO             
        RETURN         
      ENDIF                                                   
      HI=(0.,1.)                                                        
      P(N+1)=(0.,0.)                                                    
      K(N+1)=(0.,0.)                                                    
      Z=FLOAT(M-1)                                                      
      DO I=1,2                                                        
        B(I,0)=(0.,0.)                                                    
        DO J=1,2                                                        
          A(I,J,0)=(0.,0.)
        ENDDO
      ENDDO                                                 
      DO L=1,N                                                        
        X=FLOAT(L+M)*D(L+1)                                               
        Y=FLOAT(L+M-3)*D(L)                                               
        XX=U(L)-Y*B(1,L-1)                                                
        YY=V(L)+Y*B(2,L-1)                                                
        MAT(1,1)= A(1,1,L-1)*Y                                            
        MAT(2,1)=-A(2,1,L-1)*Y+HI*Z                                       
        MAT(1,2)= A(1,2,L-1)*Y+HI*Z                                       
        MAT(2,2)=-A(2,2,L-1)*Y                                            
        DET=1./(MAT(1,1)*MAT(2,2)-MAT(1,2)*MAT(2,1))                      
        A(1,1,L)= X*DET* MAT(2,2)                                         
        A(2,1,L)=-X*DET* MAT(2,1)                                         
        A(1,2,L)= X*DET*MAT(1,2)                                          
        A(2,2,L)=-X*DET* MAT(1,1)                                         
        B(1,L)=(XX*MAT(2,2)-MAT(1,2)*YY)*DET                              
        B(2,L)=(MAT(1,1)*YY-XX*MAT(2,1))*DET                               
      ENDDO                                                        
      DO L=N,1,-1                                                     
        P(L)=A(1,1,L)*P(L+1)+A(1,2,L)*K(L+1)+B(1,L)                       
        K(L)=A(2,1,L)*P(L+1)+A(2,2,L)*K(L+1)+B(2,L)                       
      ENDDO                                                         
      RETURN                                                            
      END SUBROUTINE RESOL 

! *****************************************************************************
!                                                           TLG_DIR_SN
! *****************************************************************************
!
      SUBROUTINE TLG_DIR_SN (SPC,F_CHAMP,LEGPOL,PARITE)
!
!**** TLG_DIR_SN  Transformee de Legendre directe
!
!     B. Legras        LMD      creation: 14-03-84
!                               certification: 16-05-84
!                               performance: 0.569 ms pour N_TRUNC = 21
!                                            et M_MAX = 20
!
!     Objet:  Transformee de Legendre directe de la representation
!             semi-Fourier a la representation spectrale.
!             Troncature spherique complete.
!             Description verticale de la troncature.
!             Troncature parametrisable par le COMMON TRCVERT.
!             En standard: troncature triangulaire homogene.
!
!     Arguments:
!             LEGPOL: tableau des fonctions de Legendre associees ou de
!                     leurs pseudo-derivees
!             PARITE: parite de la fonction de transformation contenue
!                     dans LEGPOL
!             F_CHAMP:tableau des coefficients de semi-Fourier en entree.
!                     rangement: l'hemisphere Nord d'abord et dans chaque
!                                hemisphere, de l'Equateur vers le Pole
!             SPC   : tableau des coefficients d'harmoniques spheriques
!                     en sortie
!
!     Methode:Les sommations en latitude sont regroupees et calculees
!             par MXMA pour les premieres colonnes de la troncature.
!             Les dernieres colonnes sont calculees par SDOT.
!             On separe en entree les contributions symetriques et
!             antisymetriques.
!
!     Externes: MATMUL
!               DOT_PRODUCT
!
! ****************************************************************************

      REAL, INTENT(IN) ::    LEGPOL  (NB_DEGLIB,NH_LAT)
      REAL, INTENT(OUT) ::   SPC     (NCPLX,NB_DEGLIB)
      REAL, INTENT(IN) ::    F_CHAMP (NCPLX,NH1_LONG,N_LAT)
      INTEGER :: PARITE
      INTEGER :: I_PLUS, JMOD, JN, ISA, ISS, INA, INS

!*    2. Separation des contributions symetriques et antisymetriques
!      --------------------------------------------------------------
!
      I_PLUS = N2_LONG * NH_LAT
      IF(PARITE==1) THEN
        Z_FS(:,:,1:NH_LAT) = F_CHAMP(:,:,1:NH_LAT) + F_CHAMP(:,:,NH_LAT+1:N_LAT)
        Z_FA(:,:,1:NH_LAT) = F_CHAMP(:,:,1:NH_LAT) - F_CHAMP(:,:,NH_LAT+1:N_LAT)
      ELSE
        Z_FS(:,:,1:NH_LAT) = F_CHAMP(:,:,1:NH_LAT) - F_CHAMP(:,:,NH_LAT+1:N_LAT)
        Z_FA(:,:,1:NH_LAT) = F_CHAMP(:,:,1:NH_LAT) + F_CHAMP(:,:,NH_LAT+1:N_LAT)
      ENDIF
!
!*    3.  Calcul matriciel pour les ordres jusqu'a M_MAX-1
!     ----------------------------------------------------

      ISA = MVA_DEB(0)
      ISS = MVS_DEB(0)
      INA = ISA+MVA_LG (0)-1
      INS = ISS+MVS_LG (0)-1
      SPC(1, ISS:INS) = MATMUL( &
           Z_FS(1, 0, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISS:INS, 1:NH_LAT)))
      SPC(1, ISA:INA) = MATMUL( &
           Z_FA(1, 0, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISA:INA, 1:NH_LAT)))
      SPC(2, ISS:INS) = 0.
      SPC(2, ISA:INA) = 0.

      DO JMOD = 1, M_MAX - 1

!*      3.1 Definition des parametres de boucle
        ISA = MVA_DEB(JMOD)
        ISS = MVS_DEB(JMOD)
        INA = ISA+MVA_LG (JMOD)-1
        INS = ISS+MVS_LG (JMOD)-1

!*      3.2 Transformee symetrique
        SPC(1:2, ISS:INS) = MATMUL( &
           Z_FS(1:2, JMOD, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISS:INS, 1:NH_LAT)))

!*      3.3 Transformee antisymetrique
        SPC(1:2, ISA:INA) = MATMUL( &
           Z_FA(1:2, JMOD, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISA:INA, 1:NH_LAT)))

      ENDDO

!*    4.  Calcul scalaire de l'ordre M_MAX a N_TRUNC
!     ----------------------------------------------
!
      DO JMOD = M_MAX, N_TRUNC-1
        ISA = MVA_DEB(JMOD)
        ISS = MVS_DEB(JMOD)
        INA = ISA+MVA_LG(JMOD)-1
        INS = ISS+MVS_LG(JMOD)-1

!*      4.1 Partie symetrique
        DO JN = ISS, INS
          SPC( 1, JN) = DOT_PRODUCT(LEGPOL( JN , 1:NH_LAT), &
                                    Z_FS  ( 1, JMOD, 1:NH_LAT))
          SPC( 2, JN) = DOT_PRODUCT(LEGPOL( JN , 1:NH_LAT), &
                                    Z_FS  ( 2, JMOD, 1:NH_LAT))
        ENDDO

!*      4.2 Partie antisymetrique
        DO JN = ISA, INA
            SPC( 1, JN ) = DOT_PRODUCT(LEGPOL( JN , 1:NH_LAT), &
                                       Z_FA  ( 1, JMOD, 1:NH_LAT))
            SPC( 2, JN ) = DOT_PRODUCT(LEGPOL( JN , 1:NH_LAT), &
                                       Z_FA  ( 2, JMOD, 1:NH_LAT))
        ENDDO
      ENDDO
      ISS = MVS_DEB(N_TRUNC)
      SPC(1, ISS) = DOT_PRODUCT(LEGPOL( ISS , 1:NH_LAT), &
                                       Z_FS  ( 1, N_TRUNC, 1:NH_LAT))
      SPC(2, ISS) = DOT_PRODUCT(LEGPOL( ISS , 1:NH_LAT), &
                                       Z_FS  ( 2, N_TRUNC, 1:NH_LAT))
      
      RETURN 
      END SUBROUTINE TLG_DIR_SN

! *****************************************************************************
!                                                           TLG_DIR_SN_M
! *****************************************************************************
!
      SUBROUTINE TLG_DIR_SN_M (SPC,F_CHAMP,LEGPOL,PARITE)
!
!**** TLG_DIR_SN  Transformee de Legendre directe
!
!     B. Legras        LMD      creation: 14-03-84
!                               version SN: April 2006
!                               version SN_M: 28-04-2006
!                              
!     Objet:  Transformee de Legendre directe de la representation
!             semi-Fourier a la representation spectrale.
!             Troncature spherique complete.
!             Description verticale de la troncature parametrisee.
!             En standard: troncature triangulaire homogene.
!
!     Arguments:
!             LEGPOL: tableau des fonctions de Legendre associees ou de
!                     leurs pseudo-derivees
!             PARITE: parite de la fonction de transformation contenue
!                     dans LEGPOL
!             F_CHAMP:tableau des coefficients de semi-Fourier en entree.
!                     rangement: l'hemisphere Nord d'abord et dans chaque
!                                hemisphere, de l'Equateur vers le Pole
!             SPC   : tableau des coefficients d'harmoniques spheriques
!                     en sortie
!
!     Multilevel fields
!
!     Methode:Les sommations en latitude sont regroupees et calculees
!             par MATMUL.
!             On separe en entree les contributions symetriques et
!             antisymetriques.
!
!     External: MATMUL
!               
!     Explicit assumed shape interface required for SPC and F_CHAMP
!
! ****************************************************************************

      REAL, INTENT(IN) ::    LEGPOL  (NB_DEGLIB,NH_LAT)
      REAL, INTENT(OUT) ::   SPC     (:,:,:)     ! (nb levels,NB-DEGLIB,N_LAT
      REAL, INTENT(IN) ::    F_CHAMP (:,:,0:,:)   ! (nb levels,NCPLX,0:NH_LONG,N_LAT)
      INTEGER PARITE
      INTEGER JMOD, ISA, ISS, INA, INS, MLEV
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Z_FA, Z_FS

!*    1. Allocate memory
!     ------------------

      MLEV=size(F_CHAMP,1)
      allocate(Z_FA(MLEV*NCPLX,0:NH_LONG,NH_LAT))
      allocate(Z_FS(MLEV*NCPLX,0:NH_LONG,NH_LAT))

!*    2. Separation des contributions symetriques et antisymetriques
!      --------------------------------------------------------------

      if(PARITE==1) then
        Z_FS(:,:,1:NH_LAT) = RESHAPE( &
          F_CHAMP(:,:,:,1:NH_LAT) + F_CHAMP(:,:,:,NH_LAT+1:N_LAT), &
          (/MLEV*NCPLX,NH1_LONG,NH_LAT/))
        Z_FA(:,:,1:NH_LAT) = RESHAPE( &
          F_CHAMP(:,:,:,1:NH_LAT) - F_CHAMP(:,:,:,NH_LAT+1:N_LAT), &
          (/MLEV*NCPLX,NH1_LONG,NH_LAT/))
      ELSE
        Z_FS(:,:,1:NH_LAT) = RESHAPE( &
          F_CHAMP(:,:,:,1:NH_LAT) - F_CHAMP(:,:,:,NH_LAT+1:N_LAT), &
          (/MLEV*NCPLX,NH1_LONG,NH_LAT/))
        Z_FA(:,:,1:NH_LAT) = RESHAPE( &
          F_CHAMP(:,:,:,1:NH_LAT) + F_CHAMP(:,:,:,NH_LAT+1:N_LAT), &
          (/MLEV*NCPLX,NH1_LONG,NH_LAT/))
      ENDIF

!*    3.  Calcul matriciel
!     --------------------

      ISA = MVA_DEB(0)
      ISS = MVS_DEB(0)
      INA = ISA+MVA_LG (0)-1
      INS = ISS+MVS_LG (0)-1
      SPC(1:MLEV,1:1, ISS:INS) = RESHAPE(MATMUL( &
           Z_FS(1:MLEV, 0, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISS:INS, 1:NH_LAT))), &
           (/MLEV,1,MVS_LG(0)/))
      SPC(1:MLEV,1:1, ISA:INA) = RESHAPE(MATMUL( &
           Z_FA(1:MLEV, 0, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISA:INA, 1:NH_LAT))), &
           (/MLEV,1,MVA_LG(0)/))
      SPC(1:MLEV,2, ISS:INS) = 0.
      SPC(1:MLEV,2, ISA:INA) = 0.

      DO JMOD = 1, N_TRUNC

!*      3.1 Definition des parametres de boucle
        ISA = MVA_DEB(JMOD)
        ISS = MVS_DEB(JMOD)
        INA = ISA+MVA_LG (JMOD)-1
        INS = ISS+MVS_LG (JMOD)-1

!*      3.2 Transformee symetrique
        SPC(1:MLEV,1:2, ISS:INS) = RESHAPE(MATMUL( &
           Z_FS(:, JMOD, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISS:INS, 1:NH_LAT))), &
           (/MLEV,NCPLX,MVS_LG(JMOD)/))

!*      3.3 Transformee antisymetrique
        SPC(1:MLEV,1:2, ISA:INA) = RESHAPE(MATMUL( &
           Z_FA(:, JMOD, 1:NH_LAT), &
           TRANSPOSE(LEGPOL(ISA:INA, 1:NH_LAT))), &
           (/MLEV,NCPLX,MVA_LG(JMOD)/))

      ENDDO
 
      ISS = MVS_DEB(N_TRUNC)
      SPC(1:MLEV,1:2, ISS:ISS) = RESHAPE(MATMUL( &
           Z_FS(:, N_TRUNC, 1:NH_LAT), &
           RESHAPE(LEGPOL(ISS, 1:NH_LAT), &
           (/NH_LAT,1/))), &
           (/MLEV,NCPLX,1/))

      deallocate(Z_FA,Z_FS)
      
      RETURN 
      END SUBROUTINE TLG_DIR_SN_M
      
!*****************************************************************************
!                                                        TLG_INV_SN
!*****************************************************************************
!
      SUBROUTINE TLG_INV_SN (F_CHAMP,SPC,LEGPOL,PARITE)
!
!**** TLG_INV_SN  Transformee de Legendre inverse 
!
!     B. Legras        LMD      Creation: 16-04-84
!                               Certification: 16-05-84
!                               Performance: 0.5296 ms pour N_TRUNC = 21
!
!     Objet:  Transformee de Legendre inverse de la representation
!             spectrale a la representation en semi-Fourier.
!             Troncature spherique complete.
!             Description verticale de la troncature parametrisee
!             par le common TRCVERT.
!             En standard: troncature triangulaire homogene.
!
!     Arguments:
!             LEGPOL:  Tableau des fonctions de Legendre associees ou de
!                      leur pseudo-derivees.
!             PARITE:  Parite de la fonction de transformation contenue
!                      dans LEGPOL.
!             SPC   :  Tableau des coefficients d'harmoniques spheriques.
!                      En entree.
!             F_CHAMP: Tableau des coefficients de semi-Fourier.
!                      En sortie.
!                      Rangement: l'hemisphere Nord d'abord et dans chaque
!                      hemisphere, de l'Equateur vers le Pole.
!
!     Methode: Les sommations en l pour les colonnes de la troncature
!             sont regroupees et calculees ensemble pour les differentes
!             latitudes.
!             Les contributions symetriques et antisymetriques sont
!             calculees separement et combinees en sortie.
!
!     Externes: MATMUL
!
! ****************************************************************************

      REAL, INTENT(IN) ::    LEGPOL   (NB_DEGLIB,NH_LAT)
      REAL, INTENT(IN) ::    SPC      (NCPLX,NB_DEGLIB)
      REAL, INTENT(OUT) ::   F_CHAMP  (NCPLX,0:NH_LONG,N_LAT)
      INTEGER PARITE
      INTEGER JC, JMOD, INA, INS, ISA, ISS
 
!*    2. Calcul matriciel de la sommation sur le degre
!     ------------------------------------------------
      ISA = MVA_DEB(0)
      ISS = MVS_DEB(0)
      INA = ISA+MVA_LG (0)-1
      INS = ISS+MVS_LG (0)-1
      
      Z_FS(1, 0, 1:NH_LAT) = MATMUL( &
               SPC(1, ISS:INS), &
               LEGPOL(ISS:INS, 1:NH_LAT))

      Z_FA(1, 0, 1:NH_LAT) = MATMUL( &
               SPC(1, INA:ISA), &
               LEGPOL(INA:ISA, 1:NH_LAT))
      Z_FS(2, 0, 1:NH_LAT) = 0.
      Z_FA(2, 0, 1:NH_LAT) = 0.

      DO JMOD = 0, N_TRUNC-1

!*      2.1 Definition des parametres de boucle
        ISA = MVA_DEB(JMOD)
        ISS = MVS_DEB(JMOD)
        INA = ISA+MVA_LG (JMOD)-1
        INS = ISS+MVS_LG (JMOD)-1

!*      2.2 Transformee symetrique

        Z_FS(1:2, JMOD, 1:NH_LAT) = MATMUL( &
               SPC(1:2, ISS:INS), &
               LEGPOL(ISS:INS, 1:NH_LAT))

!*      2.3 Transformee antisymetrique

        Z_FA(1:2, JMOD, 1:NH_LAT) = MATMUL( &
               SPC(1:2, ISA:INA), &
               LEGPOL(ISA:INA, 1:NH_LAT))

      ENDDO

      ISS = MVS_DEB(N_TRUNC)
      Z_FS(1, N_TRUNC, 1:NH_LAT) = SPC(1, ISS) * LEGPOL(ISS, 1:NH_LAT)
      Z_FS(2, N_TRUNC, 1:NH_LAT) = SPC(2, ISS) * LEGPOL(ISS, 1:NH_LAT)
      Z_FA(1:2, N_TRUNC, 1:NH_LAT) = 0.
      
!*    3.  Combinaison des parties symetriques et antisymetriques
!     ----------------------------------------------------------

      IF (PARITE==1) THEN
      DO JC = 1, NH_LAT
        F_CHAMP(1,0,JC+NH_LAT) = Z_FS(1,0,JC) - Z_FA(1,0,JC)
        F_CHAMP(1,0,JC       ) = Z_FS(1,0,JC) + Z_FA(1,0,JC)
        F_CHAMP(2,0,JC+NH_LAT) = 0.
        F_CHAMP(2,0,JC       ) = 0.
      ENDDO     
      DO JC = 1, NH_LAT
      DO JMOD = 1, N_TRUNC
        F_CHAMP(1,JMOD,JC+NH_LAT) = Z_FS(1,JMOD,JC) - Z_FA(1,JMOD,JC)
        F_CHAMP(1,JMOD,JC       ) = Z_FS(1,JMOD,JC) + Z_FA(1,JMOD,JC)
        F_CHAMP(2,JMOD,JC+NH_LAT) = Z_FS(2,JMOD,JC) - Z_FA(2,JMOD,JC)
        F_CHAMP(2,JMOD,JC       ) = Z_FS(2,JMOD,JC) + Z_FA(2,JMOD,JC)
      ENDDO; ENDDO
!
      ELSE
      DO JC = 1, NH_LAT
        F_CHAMP(1,0,JC+NH_LAT) = Z_FA(1,0,JC) - Z_FS(1,0,JC)
        F_CHAMP(1,0,JC       ) = Z_FA(1,0,JC) + Z_FS(1,0,JC)
        F_CHAMP(2,0,JC+NH_LAT) = 0.
        F_CHAMP(2,0,JC       ) = 0.
      ENDDO
      DO JC = 1, NH_LAT
      DO JMOD = 1, N_TRUNC
        F_CHAMP(1,JMOD,JC+NH_LAT) = Z_FA(1,JMOD,JC) - Z_FS(1,JMOD,JC)
        F_CHAMP(1,JMOD,JC       ) = Z_FA(1,JMOD,JC) + Z_FS(1,JMOD,JC)
        F_CHAMP(2,JMOD,JC+NH_LAT) = Z_FA(2,JMOD,JC) - Z_FS(2,JMOD,JC)
        F_CHAMP(2,JMOD,JC       ) = Z_FA(2,JMOD,JC) + Z_FS(2,JMOD,JC)
      ENDDO; ENDDO
      ENDIF

!*    4.  Mise a zero des modes de Fourier non representes
!     ----------------------------------------------------

      DO JMOD = N_TRUNC+1, NH_LONG
      DO JC = 1, N_LAT
        F_CHAMP(1,JMOD,JC) = 0.
        F_CHAMP(2,JMOD,JC) = 0.
      ENDDO; ENDDO
!
!
      RETURN
      END SUBROUTINE TLG_INV_SN

!*****************************************************************************
!                                                        TLG_INV_SN_M
!*****************************************************************************
!
      SUBROUTINE TLG_INV_SN_M (F_CHAMP,SPC,LEGPOL,PARITE)
!
!**** TLG_INV_SN  Transformee de Legendre inverse 
!
!     B. Legras        LMD      Creation: 16-04-84
!                               version SN: April 2006
!                               version SN_M: 28-04-2006
!
!     Objet:  Transformee de Legendre inverse de la representation
!             spectrale a la representation en semi-Fourier.
!             Troncature spherique complete.
!             Description verticale de la troncature parametrisee
!             En standard: troncature triangulaire homogene.
!
!     Arguments:
!             LEGPOL:  Tableau des fonctions de Legendre associees ou de
!                      leur pseudo-derivees.
!             PARITE:  Parite de la fonction de transformation contenue
!                      dans LEGPOL.
!             SPC   :  Tableau des coefficients d'harmoniques spheriques.
!                      En entree.
!             F_CHAMP: Tableau des coefficients de semi-Fourier.
!                      En sortie.
!                      Rangement: l'hemisphere Nord d'abord et dans chaque
!                      hemisphere, de l'Equateur vers le Pole.
!
!     Methode: Les sommations en l pour les colonnes de la troncature
!             sont regroupees et calculees ensemble pour les differentes
!             latitudes.
!             Les contributions symetriques et antisymetriques sont
!             calculees separement et combinees en sortie.
!
!     Externes: MATMUL
!
!     Explicit assumed shape interface required for SPC and F_CHAMP
!
! ****************************************************************************

      REAL, INTENT(IN) ::    LEGPOL   (NB_DEGLIB,NH_LAT)
      REAL, INTENT(IN) ::    SPC      (:,:,:)   ! (nb levels, NCPLX,NB_DEGLIB)
      REAL, INTENT(OUT) ::   F_CHAMP  (:,:,0:,:) ! (nb levels, NCPLX,0:NH_LONG,N_LAT)
      INTEGER PARITE
      INTEGER JC, JMOD, INA, INS, ISA, ISS, MLEV

      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: Z_FA, Z_FS

!*    1. Allocate memory
!     ------------------
      
      MLEV=size(SPC,1)
      ALLOCATE(Z_FA(MLEV,NCPLX,0:NH_LONG,NH_LAT),Z_FS(MLEV,NCPLX,0:NH_LONG,NH_LAT))
 
!*    2. Calcul matriciel de la sommation sur le degre
!     ------------------------------------------------
      ISA = MVA_DEB(0)
      ISS = MVS_DEB(0)
      INA = ISA+MVA_LG (0)-1
      INS = ISS+MVS_LG (0)-1
      
      Z_FS(1:MLEV, 1, 0, 1:NH_LAT) = MATMUL( &
               SPC(1:MLEV,1, ISS:INS), &
               LEGPOL(ISS:INS, 1:NH_LAT))

      Z_FA(1:MLEV, 1, 0, 1:NH_LAT) = MATMUL( &
               SPC(1:MLEV,1, INA:ISA), &
               LEGPOL(INA:ISA, 1:NH_LAT))
      Z_FS(1:MLEV, 2, 0, 1:NH_LAT) = 0.
      Z_FA(1:MLEV, 2, 0, 1:NH_LAT) = 0.

      DO JMOD = 0, N_TRUNC-1   

!*      2.1 Definition des parametres de boucle
        ISA = MVA_DEB(JMOD)
        ISS = MVS_DEB(JMOD)
        INA = ISA+MVA_LG (JMOD)-1
        INS = ISS+MVS_LG (JMOD)-1

!*      2.2 Transformee symetrique

        Z_FS(1:MLEV, 1:2, JMOD, 1:NH_LAT) = RESHAPE(MATMUL( &
               RESHAPE(SPC(1:MLEV,1:2, ISS:INS), &
               (/2*MLEV,MVS_LG(JMOD)/)), &
               LEGPOL(ISS:INS, 1:NH_LAT)), &
               (/MLEV,NCPLX,NH_LAT/))

!*      2.3 Transformee antisymetrique

        Z_FA(1:MLEV, 1:2, JMOD, 1:NH_LAT) = RESHAPE(MATMUL( &
               RESHAPE(SPC(1:MLEV,1:2, ISA:INA), &
               (/2*MLEV,MVA_LG(JMOD)/)), &
               LEGPOL(ISA:INA, 1:NH_LAT)), &
               (/MLEV,NCPLX,NH_LAT/))

      ENDDO

      ISS = MVS_DEB(N_TRUNC)
      
      do JC=1,MLEV
        Z_FS(JC, 1, N_TRUNC, 1:NH_LAT) = SPC(JC, 1, ISS) * LEGPOL(ISS, 1:NH_LAT)
        Z_FS(JC, 2, N_TRUNC, 1:NH_LAT) = SPC(JC, 2, ISS) * LEGPOL(ISS, 1:NH_LAT)
      enddo
      Z_FA(1:MLEV, 1:2, N_TRUNC, 1:NH_LAT) = 0.
      
!*    3.  Combinaison des parties symetriques et antisymetriques
!     ----------------------------------------------------------

      IF (PARITE==1) THEN
        F_CHAMP(1:MLEV,1,0,1:NH_LAT) = &
                             Z_FS(:,1,0,:) +  Z_FA(:,1,0,:)
        F_CHAMP(1:MLEV,1,0,NH_LAT+1:N_LAT) = &
                             Z_FS(:,1,0,:) -  Z_FA(:,1,0,:)
        F_CHAMP(1:MLEV,2,0,1:N_LAT) = 0.

        F_CHAMP(1:MLEV,1:2,1:N_TRUNC,1:NH_LAT) = &
                             Z_FS(:,:,1:N_TRUNC,:) +  Z_FA(:,:,1:N_TRUNC,:)
        F_CHAMP(1:MLEV,1:2,1:N_TRUNC,NH_LAT+1:N_LAT) = &
                             Z_FS(:,:,1:N_TRUNC,:) -  Z_FA(:,:,1:N_TRUNC,:)
      ELSE
        F_CHAMP(1:MLEV,1,0,1:NH_LAT) = &
                             Z_FA(:,1,0,:) +  Z_FS(:,1,0,:)
        F_CHAMP(1:MLEV,1,0,NH_LAT+1:N_LAT) = &
                             Z_FA(:,1,0,:) -  Z_FS(:,1,0,:)
        F_CHAMP(1:MLEV,2,0,1:N_LAT) = 0.

        F_CHAMP(1:MLEV,1:2,1:N_TRUNC,1:NH_LAT) = &
                             Z_FA(:,:,1:N_TRUNC,:) +  Z_FS(:,:,1:N_TRUNC,:)
        F_CHAMP(1:MLEV,1:2,1:N_TRUNC,NH_LAT+1:N_LAT) = &
                             Z_FA(:,:,1:N_TRUNC,:) -  Z_FS(:,:,1:N_TRUNC,:)
      ENDIF

!*    4.  Mise a zero des modes de Fourier non representes
!     ----------------------------------------------------

      F_CHAMP(1:MLEV,1:2,N_TRUNC+1:NH_LONG,1:N_LAT) = 0.

      DEALLOCATE(Z_FA,Z_FS)

      RETURN
      END SUBROUTINE TLG_INV_SN_M

!*****************************************************************************
!                                                                 BSHARM
!*****************************************************************************

      SUBROUTINE BSHARM(PLM,X,N)
!  Calcul des fonctions de Legendre associees
!  d'apres la routine LGNDRE de l'ECMWF

      REAL, INTENT(IN):: X
      INTEGER, INTENT(IN):: N
      
      REAL (KIND=8) :: SYX, F1M, F2M, RE1, E1, E2
      REAL (KIND=8), ALLOCATABLE :: Y(:)
      integer, allocatable:: mv_long(:),mve_start(:)
      INTEGER :: NI, LM, M1, M, L, JMOD, jl
      integer :: isa, iss, ina, ins, ise
      REAL, INTENT(OUT) :: PLM(:)
      ALLOCATE(Y((N_TRUNC+1) * (N_TRUNC+4) / 2))
      ALLOCATE(mv_long(0:N_TRUNC),mve_start(0:N_TRUNC))
      SYX = DSQRT(1.D0 - X * X)
      NI = N + 1
      LM = 2
      Y(1) = 1.D0
      F1M = DSQRT(3.D0)
      Y(2) = F1M*X

      DO M1 = 1, NI
        M = M1 - 1
        RE1 = DSQRT(M+M+3.D0)
        E1 = 1.D0/RE1
        IF( M > 0 ) THEN
          F2M = F1M * SYX / DSQRT(M+M+0.D0)
          F1M = F2M * RE1
          LM = LM + 1
          Y(LM) = F2M
          LM = LM+1
          Y(LM) = F1M * X
        ENDIF
        IF( M1 < NI ) THEN
          DO L = M+2, NI
            E2 = DSQRT((4.D0*L*L -1.D0)/(L*L - M*M))
            LM = LM + 1
            Y(LM) = E2 * (X * Y(LM-1) - E1 * Y(LM-2))
            E2 = 1.D0 / E2
            E1 = E2
          ENDDO
        ENDIF
      ENDDO
      mve_start(0) = 1
      mv_long  (0) = N_TRUNC + 1

      DO JMOD = 1, N_TRUNC
        mve_start(JMOD) = mve_start(JMOD-1) + mv_long(JMOD-1) + 1
        mv_long  (JMOD) = mv_long  (JMOD-1) - 1
      ENDDO  
 
      DO JMOD = 0, N_TRUNC
          isa = mva_deb(jmod)
          iss = mvs_deb(jmod)
          ina = mva_lg (jmod)
          ins = mvs_lg (jmod)
          ise = mve_start(jmod)
          if(ina.eq.0) goto 542
          do jl = 0, ina - 1
            plm(isa+jl) = y(ise+2*jl+1)
          ENDDO
  542     continue
          do jl = 0, ins - 1
            plm(iss+jl) = y(ise+2*jl)
          enddo
      ENDDO

      DEALLOCATE(Y,mv_long,mve_start)

      RETURN
      END SUBROUTINE BSHARM

END MODULE sphereharmSPE
