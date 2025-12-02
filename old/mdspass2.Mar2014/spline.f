!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   SPLINE FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SPLC(X, N, Y, DF, IOPT, C, NC, IER)
!************************************************************************
!*  COMPUTE THE COEFFICIENTS OF THE CUBIC SPLINE.                       *
!*  PARAMETERS                                                          *
!*    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
!*    (2) N: NUMBER OF KNOWN POINTS                                     *
!*    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         *
!*    (4) DF: 1-DIM. ARRAY FOR DIFFERENTIALS AT END POINTS              *
!*    (5) IOPT: 1-DIM. ARRAY SPECIFYING THE CONTENT OF DF               *
!*    (6) C: 2-DIM. WORKING ARRAY                                       *
!*    (7) NC: ROW SIZE OF THE ARRAY (C)                                 *
!*    (8) IER: ERROR CODE                                               *
!*  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  *
!************************************************************************
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION X(N), Y(N), DF(2), IOPT(2), C(NC,3), EC(4), D(2)
!C
      IF (N.LT.2 .OR. NC.LT.N-1 .OR. IOPT(1).LT.1 .OR. IOPT(1).GT.3 
     &    .OR. IOPT(2).LT.1 .OR. IOPT(2).GT.3) THEN
       IER = 2
!       WRITE(*,*) '(SUBR. SPLC) INVALID ARGUMENT.',N,NC,IOPT(1),IOPT(2)
       RETURN
      ENDIF
      DO 5 I=1,N-1
       IF (X(I) .GE. X(I+1)) THEN
        IER = 1
!        WRITE(*,*) '(SUBR. SPLC) X SHOULD SATISFY UPWARD ORDER.'
        RETURN
       ENDIF
    5 CONTINUE
      IER = 0
!C  SET THE END CONDITIONS.
      II = 2
      KS = 1
      KE = MIN0(4,N)
      IDER = 1
      DO 70 I=1,2
       I1 = 2 * I - 1
       I2 = 2 * I
       IB = IOPT(I)
!       GO TO (10, 20, 30), IB
       if (ib==1) then
          go to 10
       else if (ib==2) then
          go to 20
       else
          go to 30
       endif
   10  EC(I1) = 0.0D0
       EC(I2) = 2.0D0 * DF(I)
       GO TO 70
   20  D(I) = DF(I)
   25  IF (I .EQ. 2) II = N
       H = X(II) - X(II-1)
       EC(I1) = 1.0D0
       HY = Y(II) - Y(II-1)
       EC(I2) = 6.0D0 * (HY / H - D(I)) / H
       IF (I .EQ. 2) EC(I2) = - EC(I2)
       GO TO 70
   30  IF (I .NE. 1) THEN
        KS = MAX0(1,N-3)
        KE = N
        IDER = N
       ENDIF
       A2 = 0.0D0
       D(I) = 0.0D0
       DO 60 K=KS,KE
        IF (IDER .NE. K) THEN
         A1 = 1.0D0
         DO 50 J=KS,KE
          IF (J .NE. IDER .AND. J .NE. K) THEN
           X1 = X(IDER) - X(J)
           X2 = X(K) - X(J)
           A1 = A1 * X1 / X2
          ENDIF
   50    CONTINUE
         X3 = X(K) - X(IDER)
         D(I) = D(I) + A1 * Y(K) / X3
         A2 = A2 - 1.0D0 / X3
        ENDIF
   60  CONTINUE
       D(I) = D(I) + Y(IDER) * A2
       GO TO 25
   70 CONTINUE
!C  SET THE ELEMENTS FOR THE SYMMETRIC TRIDIAGONAL EQUATION.
      IF (N .NE. 2) THEN
       H1 = X(2) - X(1)
       Y1 = Y(2) - Y(1)
       DO I=2,N-1
        H2 = X(I+1) - X(I)
        Y2 = Y(I+1) - Y(I)
        HH = H1 + H2
        C(I,1) = H2 / HH
        C(I,2) = 1.0D0 - C(I,1)
        C(I,3) = 6.0D0 * (Y2 / H2 - Y1 / H1) / HH
        H1 = H2
        Y1 = Y2
      enddo
      ENDIF
!C  SOLVE THE EQUATION
      C(1,1) = - EC(1) * 0.5D0
      C(1,2) =   EC(2) * 0.5D0
      IF (N .NE. 2) THEN
       DO K=2,N-1
        PIV = 2.0D0 + C(K,2) * C(K-1,1)
        C(K,1) = - C(K,1) / PIV
        C(K,2) = (C(K,3) - C(K,2) * C(K-1,2)) / PIV
       enddo
      ENDIF
      DY1 = (EC(4) - EC(3) * C(N-1,2)) / (2.0D0 + EC(3) * C(N-1,1))
      DO I=1,N-1
       K = N - I
       DY2 = C(K,1) * DY1 + C(K,2)
       H = X(K+1) - X(K)
       C(K,3) = (DY1 - DY2) / (6.0D0 * H)
       C(K,2) = 0.5D0 * DY2
       C(K,1) = (Y(K+1) - Y(K)) / H - (C(K,2) + C(K,3) * H) * H
       DY1 = DY2
      enddo
!C
      RETURN
      END
  
      SUBROUTINE SPLF(X, N, Y, C, NC, V, M, F, IER)
!************************************************************************
!*  INTERPOLATION BY THE CUBIC SPLINE.                                  *
!*  PARAMETERS                                                          *
!*    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
!*    (2) N: NUMBER OF KNOWN POINTS                                     *
!*    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         *
!*    (4) C: 2-DIM. WORKING ARRAY                                       *
!*    (5) NC: ROW SIZE OF THE ARRAY (C)                                 *
!*    (6) V: 1-DIM. ARRAY FOR POINTS WHICH INTERPOLATION MUST BE MADE   *
!*    (7) M: NUMBER OF POINTS FOR WHICH INTERPOLATION MUST BE MADE      *
!*    (8) F: 1-DIM. WORKING ARRAY                                       *
!*    (9) IER: ERROR CODE                                               *
!*  COPYRIGHT   T. OGUNI   JUNE 30 1989   VERSION 1.0                   *
!************************************************************************
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION X(N), Y(N), C(NC,3), V(M), F(M)
!C
      IF (N .LT. 2 .OR. M .LT. 1 .OR. NC .LT. N-1) THEN
       IER = 2
!       WRITE(*,*) '(SUBR. SPLF) INVALID ARGUMENT. ', N, NC, M
       RETURN
      ENDIF
      IER = 0
      I = 1
      DO 90 K=1,M
       V1 = V(K) - X(I)
!       IF (V1) 10, 30, 40
       if (v1<0.0) then
          go to 10
       else if (v1==0.0) then
          go to 30
       else
          go to 40
       endif
   10  IF (I .GT. 1) GO TO 20
       IER = 1
       GO TO 80
   20  I = I - 1
       V1 = V(K) - X(I)
!       IF (V1) 10, 30, 80
       if (v1<0.0) then
          go to 10
       else if (v1==0.0) then
          go to 30
       else
          go to 80
       endif
   30  F(K) = Y(I)
       GO TO 90
   40  IF (I .LT. N) GO TO 50
       IER = 1
       I = N - 1
       GO TO 80
   50  V2 = V(K) - X(I+1)
!       IF (V2) 80, 60, 70
       if (v2<0.0) then
          go to 80
       else if (v2==0.0) then
          go to 60
       else
          go to 70
       endif
   60  I = I + 1
       GO TO 30
   70  I = I + 1
       V1 = V2
       GO TO 40
   80  F(K) = Y(I) + V1 * (C(I,1) + V1 * (C(I,2) + V1 * C(I,3)))
   90 CONTINUE
!C
      RETURN
      END

      SUBROUTINE SPLD(X, N, C, NC, V, M, D1, D2, IER)
!************************************************************************
!*  DIFFERENTIATION BY THE CUBIC SPLINE.                                *
!*  PARAMETERS                                                          *
!*    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
!*    (2) N: NUMBER OF KNOWN POINTS                                     *
!*    (3) C: 2-DIM. WORKING ARRAY                                       *
!*    (4) NC: ROW SIZE OF THE ARRAY (C)                                 *
!*    (5) V: 1-DIM. ARRAY FOR POINTS WHICH INTERPOLATION MUST BE MADE   *
!*    (6) M: NUMBER OF POINTS FOR WHICH INTERPOLATION MUST BE MADE      *
!*    (7) D1: 1-DIM. ARRAY FOR FIRST ORDER DIFFERENTIALS                *
!*    (8) D2: 1-DIM. ARRAY FOR SECOND ORDER DIFFERENTIALS               *
!*    (9) IER: ERROR CODE                                               *
!*  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  *
!************************************************************************
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION X(N), C(NC,3), V(M), D1(M), D2(M)
!C
      IF (N .LT. 2 .OR. NC .LT. N-1 .OR. M .LT. 1) THEN
       IER = 2
!       WRITE(*,*) '(SUBR. SPLD) INVALID ARGUMENT. ', N, NC, M
       RETURN
      ENDIF
      IER = 0
      I = 1
      DO 90 K=1,M
       V1 = V(K) - X(I)
!       IF (V1) 10, 30, 40
       if (v1<0.0) then
          go to 10
       else if (v1==0.0) then
          go to 30
       else
          go to 40
       endif
   10  IF (I .GT. 1) GO TO 20
       IER = 1
       GO TO 80
   20  I = I - 1
       V1 = V(K) - X(I)
!       IF (V1) 10, 30 ,80
       if (v1<0.0) then
          go to 10
       else if (v1==0.0) then
          go to 30
       else
          go to 80
       endif
   30  D1(K) = C(I,1)
       D2(K) = C(I,2) + C(I,2)
       GO TO 90
   40  IF (I .LT. N) GO TO 50
       IER = 1
       I = N - 1
       GO TO 80
   50  V2 = V(K) - X(I+1)
!       IF (V2) 80, 60, 70
       if (v2<0.0) then
          go to 80
       else if (v2==0.0) then
          go to 60
       else
          go to 70
       endif
   60  IF (I .GE. N-1) GO TO 80
       I = I + 1
       GO TO 30
   70  I = I + 1
       V1 = V2
       GO TO 40
   80  T = C(I,2) + 3.0D0 * C(I,3) * V1
       D1(K) = C(I,1) + (T + C(I,2)) * V1
       D2(K) = T + T
   90 CONTINUE
      RETURN
      END

