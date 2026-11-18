!======================================================================
      SUBROUTINE SIDEWALLBC
      !! Boundary conditions at the "side walls"
        ! JxT : upper wall (y-direction)
        ! JxB : lower wall
        ! KxT : anterior wall (z-direction)
        ! KxB : posterior wall
      !! Options included
        ! y-direction : Dirichlet / Neumann
        ! z-direction : Dirichlet / Neumann (3D)
        !               Periodic (2D,3D)


        ! IPX : periodic in x-direction
        ! ICH : For channel flow
        ! IDPDX : RE_TAU Constant
!======================================================================
       USE MOD_COMMON
       IMPLICIT NONE

        JUT = 1         ! 0: DIRICHLET
        JUB = 1         ! 1: NEUMANN
        JWT = 1
        JWB = 1
        JVT = 0
        !IF IPZ=1, THEN KUT=KUB=KVT=KVB=0
        KUT = 0
        KUB = 0
        KVT = 0
        KVB = 0

        IPX = 0        ! 0 : NOT PREODIC
        IPZ = 1        ! 1 : PREODIC

        ICH = 0        ! 0 : CHANNEL X
                       ! 1 : CHANNEL O (FOR THE CONST. MASS FLOW RATE)

        IDPDX = 0      ! 1 : RETAU CONSTANT

        IF (IPZ.EQ.1) THEN
        IF (KUT.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KUT'
            STOP
        ENDIF
        IF (KUB.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KUB'
            STOP
        ENDIF
        IF (KVT.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KVT'
            STOP
        ENDIF
        IF (KVB.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KVB'
            STOP
        ENDIF
        ENDIF


      WRITE(*,101) JUT,JUB,JWT,JWB
      WRITE(*,102) KUT,KUB,KVT,KVB
      WRITE(*,103) JVT
      WRITE(*,104) IPX,IPZ
      WRITE(*,105) ICH

 101  FORMAT('JUT= ',I4,'  JUB= ',I4,'  JWT= ',I4,'  JWB= ',I4)
 102  FORMAT('KUT= ',I4,'  KUB= ',I4,'  KVT= ',I4,'  KVB= ',I4)
 103  FORMAT('JVT= ',I4)
 104  FORMAT('IPX= ',I4,'  IPZ= ',I4)
 105  FORMAT('ICH= ',I4)

      RETURN
      END

!=======================================================================
        SUBROUTINE GETCDAREA
        ! In 2d-simulation
        ! CDAREA = ZL : d/dz
!=======================================================================
       USE MOD_COMMON
       IMPLICIT NONE

       CDAREA= ZL

      RETURN
      END

!=======================================================================
        SUBROUTINE CONVBC
      !! Boundary conditions at the outlet (at the end of postive x)
      !! Convective boundary condition
      !! IF inlet = 0, you must deactivate this subroutine
!=======================================================================
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8      :: QIN,QOUT,QRATIO
       REAL*8      :: UBAR,UCOEF,VWCOEF


!-----DIRICHLET/NEUMANN CONDITION FOR V AT THE TOP WALL
      IF (JVT.EQ.0) THEN
!$OMP PARALLEL DO
       DO K=1,N3M
       DO I=1,N1M
        VTOP(I,K) = V(I,N2,K)
        DVTOP(I,K)= 0.
       ENDDO
       ENDDO
      ELSE IF (JVT.EQ.1) THEN
!$OMP PARALLEL DO
       DO K=1,N3M
       DO I=1,N1M
        VTOP(I,K) = V(I,N2M,K)
        DVTOP(I,K)= V(I,N2M,K)-V(I,N2,K)
       ENDDO
       ENDDO
      ENDIF

!-----CALCULATE INFLUX
       QIN=0.
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QIN)
       DO K=1,N3M
       DO J=1,N2M
         QIN=U(1,J,K)*SDY(J)*SDZ(K)+QIN
       ENDDO
       ENDDO
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QIN)
       DO K=1,N3M
       DO I=1,N1M
        QIN=(V(I,1,K)-V(I,N2,K))*SDX(I)*SDZ(K)+QIN
       ENDDO
       ENDDO
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QIN)
       DO J=1,N2M
       DO I=1,N1M
        QIN=(W(I,J,1)-W(I,J,N3))*SDX(I)*SDY(J)+QIN
       ENDDO
       ENDDO

!-----CALCULATE CONVECTIVE VELOCITY Uc
       UBAR = 0.
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:UBAR)
       DO K=1,N3M
       DO J=1,N2M
        UBAR=U(N1,J,K)*SDY(J)*SDZ(K)+UBAR
       ENDDO
       ENDDO
       UBAR  = UBAR/(YL*ZL)
       UCOEF = UBAR*DTCONST*SSDX(N1M)
       VWCOEF= UBAR*DTCONST*VVDX(N1)

!-----DVOUT(1,K) IS A DUMMY VARIABLE:
!     NOT NEEDED, NOT AFFECT THE SIMULATION
       QOUT = 0.
!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QOUT)
       DO K=1,N3M
       DO J=1,N2M
        UOUT(J,K)=U(N1,J,K)-UCOEF*(U(N1,J,K)-U(N1M,J,K))
        VOUT(J,K)=V(N1,J,K)-VWCOEF*(V(N1,J,K)-V(N1M,J,K))
        WOUT(J,K)=W(N1,J,K)-VWCOEF*(W(N1,J,K)-W(N1M,J,K))
        QOUT=UOUT(J,K)*SDY(J)*SDZ(K)+QOUT
       ENDDO
       ENDDO
       QRATIO=QIN/QOUT

!-----ADJUST BOUNDARY VELOCITY TO SATISFY GLOBAL MASS CONSERVATION
!$OMP PARALLEL DO
       DO K=1,N3M
       DO J=1,N2M
        UOUT(J,K)=UOUT(J,K)*QRATIO
        VOUT(J,K)=VOUT(J,K)*QRATIO
        WOUT(J,K)=WOUT(J,K)*QRATIO
        DUOUT(J,K)=UOUT(J,K)-U(N1,J,K)
        DVOUT(J,K)=VOUT(J,K)-V(N1,J,K)
        DWOUT(J,K)=WOUT(J,K)-W(N1,J,K)
       ENDDO
       ENDDO

!!!!!   Z PERIODICITY
       IF (IPZ .EQ. 1) THEN
!$OMP PARALLEL DO
       DO J=0,N2
        DUOUT(J,N3)=DUOUT(J,  1)
        DUOUT(J, 0)=DUOUT(J,N3M)
        DVOUT(J,N3)=DVOUT(J,  1)
        DVOUT(J, 0)=DVOUT(J,N3M)
        DWOUT(J,N3)=DWOUT(J,  1)
        DWOUT(J, 0)=DWOUT(J,N3M)
       ENDDO
       ENDIF

      RETURN
      END

!=======================================================================
        SUBROUTINE MAKEFLD
      !! Make an initial field
      !  If you select Dirichlet boundary conditions at the side walls,
      !  (See 'SUBROUTINE SIDEWALLBC')
      !  values applied at this subroutine are imposed
      !  at the side walls until the end of code
!=======================================================================
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W
       IMPLICIT NONE
       INTEGER*8     :: I,J,K

       IHIST= 0
       TIME = 0.

!$OMP PARALLEL DO
        DO 10 K=0,N3
        DO 10 J=0,N2
        DO 10 I=1,N1
 10     U(I,J,K) = 1. ! Inlet boundary condition (Uniform inflow)
                      ! U(I,J,K) = 0. then, deactivate CONVBC
!$OMP PARALLEL DO
        DO 20 K=0,N3
        DO 20 J=1,N2
        DO 20 I=0,N1
 20     V(I,J,K) = 0.

!$OMP PARALLEL DO
        DO 30 K=1,N3
        DO 30 J=0,N2
        DO 30 I=0,N1
 30     W(I,J,K) = 0.

!     Z PERIODICITY
      IF (IPZ .EQ. 1) THEN
!$OMP PARALLEL DO
       DO J=0,N2
       DO I=1,N1
         U(I,J,0) =U(I,J,N3M)
         U(I,J,N3)=U(I,J,1)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO J=1,N2
       DO I=0,N1
          V(I,J,0) =V(I,J,N3M)
          V(I,J,N3)=V(I,J,1)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO J=0,N2
       DO I=0,N1
          W(I,J,0) =W(I,J,N3M)
          W(I,J,N3)=W(I,J,1)
       ENDDO
       ENDDO
      ENDIF

!     X PERIODICITY
      IF (IPX .EQ. 1) THEN
!$OMP PARALLEL DO
       DO K=0,N3
       DO J=0,N2
          U(0 ,J,K)=U(N1M,J,K)
          U(N1,J,K)=U(1  ,J,K)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO K=0,N3
       DO J=1,N2
          V(0 ,J,K)=V(N1M,J,K)
          V(N1,J,K)=V(1  ,J,K)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO K=1,N3
       DO J=0,N2
          W(0 ,J,K)=W(N1M,J,K)
          W(N1,J,K)=W(1  ,J,K)
       ENDDO
       ENDDO
      ENDIF

! Small initial wake perturbation, epsilon streamfunction=0.002 U D.
! The discrete curl cancels exactly in the interior finite-volume divergence.
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1
          U(I,J,K)=U(I,J,K)+0.002D0*( &
            EXP(-((X(I)-1.5D0)**2+Y(J+1)**2)/0.25D0) &
           -EXP(-((X(I)-1.5D0)**2+Y(J  )**2)/0.25D0))/SDY(J)
        ENDDO
        ENDDO
        DO J=1,N2
        DO I=1,N1M
          V(I,J,K)=V(I,J,K)-0.002D0*( &
            EXP(-((X(I+1)-1.5D0)**2+Y(J)**2)/0.25D0) &
           -EXP(-((X(I  )-1.5D0)**2+Y(J)**2)/0.25D0))/SDX(I)
        ENDDO
        ENDDO
        ENDDO

        WRITE(*,*)' '
        WRITE(*,100)
  100   FORMAT('INITIAL FIELD:      MAKING DONE')

      RETURN
      END

!=======================================================================
      SUBROUTINE REAL_TIME(SECONDS)
!=======================================================================
      IMPLICIT NONE
      INTEGER CLOCK_COUNT,CLOCK_MAX,CLOCK_RATE
      REAL SECONDS

      CALL SYSTEM_CLOCK ( CLOCK_COUNT, CLOCK_RATE, CLOCK_MAX )

      SECONDS=REAL(CLOCK_COUNT)/REAL(CLOCK_RATE)

      RETURN
      END

!=======================================================================
      subroutine timestamp ( )
!=======================================================================
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!  Example:
!    31 May 2001   9:45:54.872 AM
!  Modified:
!    06 August 2005
!  Author: John Burkardt
!  Parameters: None
        implicit none
        character ( len = 8 ) ampm
        integer d
        integer h
        integer m
        integer mm
        character ( len = 9 ), parameter, dimension(12) :: month = (/  &
          'January  ', 'February ', 'March    ', 'April    ',          &
          'May      ', 'June     ', 'July     ', 'August   ',          &
          'September', 'October  ', 'November ', 'December ' /)
        integer n
        integer s
        integer values(8)
        integer y

        call date_and_time ( values = values )

        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if ( h < 12 ) then
          ampm = 'AM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
          else
            ampm = 'PM'
          end if
        else
          h = h - 12
          if ( h < 12 ) then
            ampm = 'PM'
          else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
              ampm = 'Midnight'
            else
              ampm = 'AM'
            end if
          end if
        end if

      write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )  &
      d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
      write(*,*)' '

        return
      end

!=======================================================================
       FUNCTION FUNCBODY(X,Y,Z,IS)
       ! IDENTIFICATION OF AN IMMERSED BODY IN A WORKING FLUID
       ! BODY = {(X,Y,Z)|FUNCBODY(X,Y,Z) < 0}
       !
       ! IFUNCorPTS = 1
       ! (PLANFORM PROFILE OF A BODY CAN BE EXPRESSED BY A FUNCTION)
       ! IFUNCorPTS = 2
       ! (LOAD PLANFORM GEOMETRY FILE)
       !
       ! Required condition to use secant method :
       !    1. FUNCBODY(X,Y,Z)<0  inside the body
       !       FUNCBODY(X,Y,Z)=0  at the body
       !       FUNCBODY(X,Y,Z)>0  outside the body
       !    2. FUNCBODY Must be continuous function
       !    *It is sufficient to satisfy only 1. to finding forcing
       !     point but 2. is critical to use secant method
       ! Ex.
       !    FUNCBODY=X**2+Y**2+Z**2-0.5**2       ! for sphere
       !    FUNCBODY=X**2+Y**2-0.5**2            ! for infinite cylinder
       !
       ! Variables
       !   TT         : Time
       !   XFU,YFU,ZFU: Upper boundary of the body
       !   XFL,YFL,ZFL: Lower boundary of the body
       !   XC,YC,ZC   : Moving distance at an instance
       !   WTH        : Half of body thickness (IPZ=0)
       !                2*ZL                   (IPZ=1)
       !
!=======================================================================
       USE MOD_MOVINGIBM
       USE MOD_COMMON, ONLY : N3M,IFUNCorPTS,DT,TIME,SUBDT

       IMPLICIT NONE
       REAL*8     :: X,Y,Z,YC,TT,F1,F2,FF1,FF2
       REAL*8     :: FUNCBODY
       INTEGER*8  :: IS  ! IS : Section variable
       INTEGER*8  :: I

       TT=TIME+SUBDT

       IF    (IFUNCorPTS .EQ. 1) THEN

         FUNCBODY = X**2.+(Y-YC(TT))**2. -0.5D0**2.

       ELSEIF(IFUNCorPTS .EQ. 2) THEN

         IF(N3M .EQ. 1) THEN

           IF ((X .LE. XFU(1)) .OR. (X .GE. XFU(ND1))) THEN
           FUNCBODY=10D0
           GOTO 9999
           ENDIF

         ELSE

           IF ((X .LE. XFU(1)) .OR. (X .GE. XFU(ND1))&
          .OR. (Z .LE. -WTH)   .OR. (Z .GE. WTH)) THEN
           FUNCBODY=10D0
           GOTO 9999
           ENDIF

         ENDIF

         DO I=1,ND1-1
         IF (X.GE.XFU(I) .AND. X.LT.XFU(I+1)) THEN
         F1=(YFU(I+1)-YFU(I))/(XFU(I+1)-XFU(I))*(X-XFU(I))+YFU(I)
         ENDIF
         ENDDO

         DO I=ND2,2,-1
         IF (X.GE.XFL(I) .AND. X.LT.XFL(I-1)) THEN
         F2=(YFL(I-1)-YFL(I))/(XFL(I-1)-XFL(I))*(X-XFL(I))+YFL(I)
         ENDIF
         ENDDO

         IF (F1 .LE. F2) THEN
         FUNCBODY=10D0
         GOTO 9999
         ENDIF

         FF1=F1+YC(TT)
         FF2=F2+YC(TT)

         !--CAUTION--
         IF(N3M .EQ. 1.) THEN
         FUNCBODY=DABS((FF1-Y)*(FF2-Y))
         ELSE
         FUNCBODY=DABS((FF1-Y)*(FF2-Y)*(WTH-Z)*(WTH+Z))
         ENDIF
         !-----------

         IF (Y.LE.FF1 .AND. Y.GE.FF2 .AND. DABS(Z).LE.WTH) THEN
         FUNCBODY=-FUNCBODY
         ENDIF

       ENDIF

9999   RETURN
       END