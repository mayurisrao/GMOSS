C   IMSL ROUTINE NAME   - MMBSKR
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - MODIFIED BESSEL FUNCTION OF THE SECOND KIND
C                           OF NONNEGATIVE REAL FRACTIONAL ORDER
C                           FOR REAL POSITIVE ARGUMENTS SCALED BY
C                           EXP(ARG)
C
C   USAGE               - CALL MMBSKR (ARG,ORDER,N,BK,IER)
C
C   ARGUMENTS    ARG    - INPUT ARGUMENT. ARG MUST BE TYPED APPRO-
C                           PRIATELY IN THE CALLING PROGRAM. (SEE THE
C                           PRECISION/HARDWARE SECTION.) ARG MUST BE
C                           GREATER THAN ZERO.
C                ORDER  - INPUT VALUE SPECIFYING THE DESIRED ORDER OF
C                           THE BESSEL FUNCTION. ORDER MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) ORDER MUST
C                           BE GREATER THAN OR EQUAL TO ZERO.
C                N       - INPUT VALUE SPECIFYING THE NUMBER OF
C                           FUNCTION VALUES TO BE COMPUTED.
C                BK     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           COMPUTED FUNCTION VALUES. BK MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) BK(1) WILL
C                           CONTAIN THE COMPUTED VALUE FOR THE INPUT
C                           ORDER, BK(2) WILL CONTAIN THE COMPUTED
C                           FUNCTION VALUE FOR ORDER + 1, BK(3) FOR
C                           ORDER + 2, ETC.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ONE OR MORE OF
C                             THE INPUT ARGUMENTS, ARG, ORDER, OR N IS
C                             OUT OF RANGE. BK(I), I=1,N, IS SET TO
C                             MACHINE INFINITY.
C                           IER = 130 INDICATES THAT BK(1) IS GREATER
C                             THAN MACHINE INFINITY. BK(I), I=1,N, IS
C                             SET TO MACHINE INFINITY.
C                           IER = 131 INDICATES THAT BK(N)/BK(N-1) IS
C                             GREATER THAN MACHINE INFINITY, BUT THAT
C                             BK(1) IS COMPUTED CORRECTLY.  BK(I),
C                             I=2,N, IS SET TO MACHINE INFINITY.
C                           IER = 131+J INDICATES THAT BK(J) WOULD HAVE
C                             BEEN GREATER THAN MACHINE INFINITY.
C                             BK(L), L=J,N, IS SET TO MACHINE INFINITY.
C                             BK(L), L=1,J-1, IS COMPUTED CORRECTLY.
C
C   PRECISION/HARDWARE  - DOUBLE/H32,H36
C                       - SINGLE/H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE mmbskr(ARG,ORDER,N,BK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      DOUBLE PRECISION   ARG,ORDER,BK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   DMINF,DNTEN,EX,TEMP,TEMP2
      DOUBLE PRECISION   A,ALPHA,BK1,BK2,C,D,DM,DMAX,DMIN,D1,D2,D3,ENU,
     *                   F0,F1,F2,P(8),P0,Q(7),Q0,R(4),RATIO,S(4),T(7),
     *                   TOL,TWONU,TWOX,T1,T2,X2BY4,XINF
      DATA               A/.1159315156584124D0/
      DATA               D/.7978845608028654D0/
      DATA               DNTEN/38.D0/
      DATA               DMIN/2.938735878D-39/
      DATA               XINF/1.700000000D+38/
      DATA               P(1) /.8056298755690433D0/
      DATA               P(2) /.2040455002053652D02/
      DATA               P(3) /.1577056051066762D03/
      DATA               P(4) /.5366711164692075D03/
      DATA               P(5) /.9003827592912888D03/
      DATA               P(6) /.7309238866506604D03/
      DATA               P(7) /.2292993015094251D03/
      DATA               P(8) /.8224670334241132D00/
      DATA               Q(1) /.2946019862478505D02/
      DATA               Q(2) /.2775778685102212D03/
      DATA               Q(3) /.1206703255910274D04/
      DATA               Q(4) /.2762914441597915D04/
      DATA               Q(5) /.3443740505065646D04/
      DATA               Q(6) /.2210631901133786D04/
      DATA               Q(7) /.5722673383598922D03/
      DATA               R(1) /-.1363145832225935D-02/
      DATA               R(2) /.3314384097399805D-01/
      DATA               R(3) /-.2506099471699884D0/
      DATA               R(4) /.8224670334241132D0/
      DATA               S(1) /.2755746186513600D-02/
      DATA               S(2) /-.6530826876446539D-01/
      DATA               S(3) / .5187115283864121D0/
      DATA               S(4) /-.1456159004551443D01/
      DATA               T(1) /.0779373862881377D-12/
      DATA               T(2) /.1605596120426673D-09/
      DATA               T(3) /.2505214205064134D-07/
      DATA               T(4) /.2755731902076135D-05/
      DATA               T(5) /.1984126984193444D-03/
      DATA               T(6) /.8333333333332269D-02/
      DATA               T(7) /.1666666666666667D0/
      DATA               TOL  /2.775557562D-17/
C                                  A = DLOG(2.D0) - EULER S CONSTANT
C                                  D = DSQRT(2.D0/PI)
C                                  P, Q - APPROXIMATION FOR
C                                    DLOG(DGAMMA(1+ENU))/ENU + EULER S
C                                    CONSTANT
C                                  R, S - APPROXIMATION FOR
C                                    (1-ENU*PI/DSIN(ENU*PI))/(2.D0*ENU)
C                                  T - APPROXIMATION FOR DSINH(Y)/Y
C
C                                  FIRST EXECUTABLE STATEMENT
      EX = ARG
      DMAX = 10.D0**DNTEN
      IER = 0
      IF (ORDER.GE.0.D0.AND.N.GE.1 .AND. EX.GT.0.D0) GO TO 10
      IER = 129
      BK(1) = XINF
      IF (N.LT.2) GO TO 9000
      DO 5 K=2,N
         BK(K) = XINF
    5 CONTINUE
      GO TO 9000
   10 K = ORDER+0.5D0
      ENU = ORDER-FLOAT(K)
      IF (DABS(ENU).LT.1.D-30) ENU = 0.D0
      TWONU = ENU+ENU
      IEND = N+K-1
      C = ENU*ENU
      D3 = -C
      IF (EX.GT.1.D0) GO TO 85
C                                  CALCULATION OF P0, Q0 AND F0
      D1 = 0.D0
      D2 = P(1)
      T1 = 1.D0
      T2 = Q(1)
      DO 15 I=2,7,2
         D1 = C*D1+P(I)
         D2 = C*D2+P(I+1)
         T1 = C*T1+Q(I)
         T2 = C*T2+Q(I+1)
   15 CONTINUE
      D1 = ENU*D1
      T1 = ENU*T1
      F1 = DLOG(ARG)
      F0 = A+ENU*(P(8)-ENU*(D1+D2)/(T1+T2))-F1
      P0 = DEXP(ENU*F0)
      Q0 = DEXP(-ENU*(A-ENU*(P(8)+ENU*(D1-D2)/(T1-T2))-F1))
      D1 = 0.D0
      T1 = 0.D0
      DO 20 I=1,4
         D1 = C*D1+R(I)
         T1 = C*T1+S(I)
   20 CONTINUE
      F1 = ENU*F0
      F1 = F1*F1
      IF (F1.GT.1.3D0) GO TO 30
      D2 = 0.D0
      DO 25 I=1,7
         D2 = F1*D2+T(I)
   25 CONTINUE
      D2 = F0+F0*F1*D2
      GO TO 35
   30 D2 = (DEXP(ENU*F0)-DEXP(-ENU*F0))*.5D0
      D2 = D2/ENU
   35 F0 = D2-ENU*D1/((1.D0+C*T1)*P0)
      IF (EX.GT.1.D-10) GO TO 75
C                                  X.LE.1.0E-10
C                                    CALCULATION OF K(ENU,X) AND
C                                    X*K(ENU+1,X)/K(ENU,X)
      BK(1) = F0+ARG*F0
      RATIO = P0/F0
      C = ARG*DMAX
      IF (K.EQ.0) GO TO 55
C                                  CALCULATION OF K(ORDER,X) AND
C                                    X*K(ORDER+1,X)/K(ORDER,X)
C                                    FOR ORDER .GE. 1/2
      DO 40 I=1,K
         TEMP = C/RATIO
         IF (BK(1).GE.TEMP) GO TO 45
         BK(1) = RATIO*BK(1)/ARG
         TWONU = TWONU+2.D0
         RATIO = TWONU
   40 CONTINUE
      GO TO 55
   45 IER = 130
      DO 50 I=1,N
         BK(I) = XINF
   50 CONTINUE
      GO TO 9000
   55 IF (N.EQ.1) GO TO 9005
C                                  CALCULATION OF
C                                    K(ORDER+L,X)/K(ORDER+L-1,X)
C                                    FOR L = 1, 2,. . ., N
      DO 60 I=2,N
         IF (RATIO.GE.C) GO TO 65
         BK(I) = RATIO/ARG
         TWONU = TWONU+2.D0
         RATIO = TWONU
   60 CONTINUE
      JER = 1
      GO TO 140
   65 DO 70 I=2,N
         BK(I) = XINF
   70 CONTINUE
      IER = 131
      GO TO 9000
C                                  1.0E-10 .LT. X .LE. 1.0
   75 C = 1.D0
      X2BY4 = 0.25D0*ARG*ARG
      P0 = 0.5D0*P0
      Q0 = 0.5D0*Q0
      D1 = -1.D0
      D2 = 0.D0
      BK1 = 0.D0
      BK2 = 0.D0
      F1 = F0
      F2 = P0
   80 D1 = D1+2.D0
      D2 = D2+1.D0
      D3 = D1+D3
      C = X2BY4*C/D2
      F0 = (D2*F0+P0+Q0)/D3
      P0 = P0/(D2-ENU)
      Q0 = Q0/(D2+ENU)
      T1 = C*F0
      T2 = C*(P0-D2*F0)
      BK1 = BK1+T1
      BK2 = BK2+T2
      TEMP = DABS(T1/(F1+BK1))
      TEMP2 = DABS(T2/(F2+BK2))
      IF (TEMP.GT.TOL .OR. TEMP2.GT.TOL) GO TO 80
      D1 = DEXP(-ARG)
      BK1 = (F1+BK1)/D1
      BK2 = 2.D0*(F2+BK2)/(ARG*D1)
      DMINF = 24.67D0*EX+3.57D0
      GO TO 120
C                                  X .GT. 1.0
   85 TWOX = ARG+ARG
      ALPHA = 0.D0
      RATIO = 0.D0
      IF (EX.GT.4.D0) GO TO 105
C                                  1.0 .LT. X .LE. 4.0
C                                    CALCULATION OF K(ENU+1,X)/K(ENU,X)
      M = 52.0583D0/EX+5.7607D0
      D2 = M
      D1 = D2+D2
      D2 = D2-0.5D0
      D2 = D2*D2
      DO 90 I=2,M
         D1 = D1-2.D0
         D2 = D2-D1
         RATIO = (D3+D2)/(TWOX+D1-RATIO)
   90 CONTINUE
C                                  CALCULATION OF I(ABS(ENU),X) AND
C                                    I(ABS(ENU)+1,X) BY BACKWARD
C                                    RECURRENCE AND K(ENU,X) FROM THE
C                                    WRONSKIAN
      M = 2.7782D0*EX+14.4303D0
      C = DABS(ENU)
      D3 = C+C
      D1 = D3-1.D0
      D2 = M
      F1 = DMIN
      F0 = (2.D0*(C+D2)/ARG+0.5D0*ARG/(C+D2+1.D0))*DMIN
      DO 95 I=3,M
         D2 = D2-1.D0
         F2 = (D3+D2+D2)*F0
         ALPHA = (1.D0+D1/D2)*(F2+ALPHA)
         F2 = F2/ARG+F1
         F1 = F0
         F0 = F2
   95 CONTINUE
      F1 = (D3+2.D0)*F0/ARG+F1
      D1 = 0.D0
      T1 = 1.D0
      DO 100 I=1,7
         D1 = C*D1+P(I)
         T1 = C*T1+Q(I)
  100 CONTINUE
      P0 = DEXP(C*(A+C*(P(8)-C*D1/T1)-DLOG(ARG)))/ARG
      F2 = (C+0.5D0-RATIO)*F1/ARG
      BK1 = P0+(D3*F0-F2+F0+ALPHA)/(F2+F1+F0)*P0
      DMINF = 4.263D0*EX+23.977D0
      GO TO 115
C                                  X .GT. 4.0
C                                    CALCULATION OF K(ENU,X) AND
C                                    K(ENU+1,X)/K(ENU,X) BY BACKWARD
C                                    RECURRENCE
  105 M = 185.3004D0/EX+9.3715D0
      DM = M
      D2 = DM-0.5D0
      D2 = D2*D2
      D1 = DM+DM
      DO 110 I=2,M
         DM = DM-1.D0
         D1 = D1-2.D0
         D2 = D2-D1
         RATIO = (D3+D2)/(TWOX+D1-RATIO)
         ALPHA = (RATIO+RATIO*ALPHA)/DM
  110 CONTINUE
      BK1 = 1.D0/((D+D*ALPHA)*DSQRT(ARG))
      DMINF = 0.94219D0*(EX-DABS(EX-20.D0))+52.3363D0
C
C                                  CALCULATION OF K(ENU+1,X) FROM
C                                    K(ENU,X) AND K(ENU+1,X)/K(ENU,X)
  115 BK2 = BK1+BK1*(ENU+0.5D0-RATIO)/ARG
C
C                                  CALCULATION OF  IER , K(ORDER+I,X),
C                                    I = 0, 1, . . ., IER-1, ..
C                                    K(ORDER+I,X)/K(ORDER+I-1,X),
C                                    I = IER, IER+1, . . ., N
  120 BK(1) = BK1
      IF (IEND.EQ.0) GO TO 9005
      J = 2-K
      IF (J.GT.0) BK(J) = BK2
      IF (IEND.EQ.1) GO TO 9005
      M = DMINF-ENU
      IF (M.GT.IEND) M = IEND
      DO 125 I=2,M
         T1 = BK1
         BK1 = BK2
         TWONU = TWONU+2.D0
         BK2 = TWONU/ARG*BK1+T1
         J = J+1
         IF (J.GT.0) BK(J) = BK2
  125 CONTINUE
      IF (M.EQ.IEND) GO TO 9005
      RATIO = BK2/BK1
      MPLUS1 = M+1
      IER = 130
      KK = 1
      DO 135 I=MPLUS1,IEND
         TWONU = TWONU+2.D0
         RATIO = TWONU/ARG+1.D0/RATIO
         J = J+1
         IF (J.LE.1) GO TO 130
         BK(J) = RATIO
         GO TO 135
  130    TEMP = DMAX/RATIO
         IF (BK2.GE.TEMP) GO TO 155
         BK2 = RATIO*BK2
  135 CONTINUE
      JER = MPLUS1-K
      IF (JER.LE.0) JER = 1
      IF (JER.EQ.1) BK(1) = BK2
      IF (N.EQ.1) GO TO 9005
  140 J = JER+1
      DO 145 I=J,N
         IER = 131+I
         TEMP = DMAX/BK(I)
         IF (BK(I-1).GE.TEMP) GO TO 150
         BK(I) = BK(I-1)*BK(I)
  145 CONTINUE
      IER = 0
      GO TO 9005
  150 KK = IER-131
  155 DO 160 I=KK,N
         BK(I) = XINF
  160 CONTINUE
 9000 CONTINUE
C      CALL UERTST(IER,6HMMBSKR)
 9005 RETURN
      END


