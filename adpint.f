
      SUBROUTINE adpint(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
c      IMPLICIT DOUBLE PRECISION (A-H,O,P,R-Z)
      REAL*8 RINT,XL,XU,REPS,AEPS,DIF,F,RM,AEPSL,RL,RU,XU1,FINT,DIF0
      LOGICAL Q
      PARAMETER(IPMAX=100,IFMAX=5,MAXPT=100000)
      EXTERNAL F
      DIMENSION XU1(IPMAX)

      IFAIL=0
      RINT=0.0
      DIF=0.0
      IF(XL.EQ.XU) RETURN
      IF(NMAX.LE.0) NMAX=MAXPT
      AEPSL=AEPS
      IER=0
      NPT=0
      RL=XL
      RU=XU
      IU=0

1000  CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RM=0.5*(RL+RU)
      Q=IU.GE.IPMAX.OR.RM.EQ.RL.OR.RM.EQ.RU
      IF(DIF0.LT.AMAX1(ABS(FINT)*REPS,AEPSL).OR.Q) THEN
        RINT=RINT+FINT
        DIF=DIF+DIF0
        IF(Q.AND.DIF0.GT.AMAX1(ABS(RINT)*REPS,AEPSL)) THEN
          IER=11
          IFAIL=IFAIL+1
          IF(IFAIL.GT.IFMAX) THEN
            IER=22
            AEPSL=DIF*0.5
          ENDIF
        ENDIF
        IF(ABS(RU-XU).LT.REPS.OR.IU.LE.0) RETURN
        RL=RU
        RU=XU1(IU)
        IU=IU-1
      ELSE
        IU=IU+1
        XU1(IU)=RU
        RU=RM
      ENDIF

      IF(NPT.LT.NMAX) GO TO 1000

      IER=13
      RU=XU
      CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

C     ----------------------------------------------------------

      SUBROUTINE KRONRD(RI,A,B,DIF,N,F)
C      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z)
	  REAL*8  AT,BT,A,B,FBT,R1,RI,F1,F2,DIF,F,RIA,RIB
      REAL*8  W7(4),A7(4),WK7(4),WK15(4),AK15(4)


      DATA W7  /0.12948496616886969327D0, 0.27970539148927666790D0,
     *          0.38183005050511894495D0, 0.41795918367346938775D0/
      DATA A7  /0.94910791234275852452D0, 0.74153118559939443986D0,
     *          0.40584515137739716690D0, 0.0/
      DATA WK7 /0.06309209262997855329D0, 0.14065325971552591874D0,
     *          0.19035057806478540991D0, 0.20948214108472782801D0/
      DATA WK15 /0.02293532201052922496D0, 0.10479001032225018383D0,
     *          0.16900472663926790282D0, 0.20443294007529889241D0/
      DATA AK15 /0.99145537112081263920D0, 0.86486442335976907278D0,
     *          0.58608723546769113029D0, 0.20778495500789846760D0/

      AT=(B-A)/2.
      BT=(B+A)/2.
c	  print*,1,BT
      FBT=F(0.0D0+BT)
c      print*,11,FBT
      R1=W7(4)*FBT
      RI=WK7(4)*FBT
      DO 2000 K=1,3
c		print*,2,AT*A7(K)+BT
        F1=F(0.0D0+AT*A7(K)+BT)
c	    print*,22,F1
c	    print*,3,BT-AT*A7(K)
        F2=F(0.0D0+BT-AT*A7(K))
c	    print*,33,F2
        R1=R1+W7(K)*(F1+F2)
        RI=RI+WK7(K)*(F1+F2)
2000  CONTINUE

      DO 2500 K=1,4
C	  print*,45,AT*AK15(K)+BT,BT-AT*AK15(K),
C     *   AT*AK15(K)+BT+BT-AT*AK15(K)
      RIA = F(0.0D0+AT*AK15(K)+BT)
      RIB = F(0.0D0+BT-AT*AK15(K))
      RI = RI+WK15(K)*(RIA + RIB)
C      RI=RI+WK15(K)*(F(0.0D0+AT*AK15(K)+BT)+F(0.0D0+BT-AT*AK15(K)))
C	  print*,4545,RI
2500  CONTINUE

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=15
      END