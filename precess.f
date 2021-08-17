	SUBROUTINE precess(RA,DEC,EPOCH1,EPOCH2)

C**********************************************************************
C
C	REFERENCE  ---> tifr /home/soft/IDL/astron/astro/precess.pro 
C NAME:
C    PRECESS
C PURPOSE:
C    Precess coordinates from EPOCH1 to EPOCH2.
C CALLING SEQUENCE:
C    CALL PRECESS,RA,DEC,EPOCH1,EPOCH2
C INPUTS:
C    RA - Input right ascension in DEGREES
C    DEC - Input declination in DEGREES
C OUTPUTS:
C    RA - Right ascension after precession in DEGREES
C    DEC - Declination after precession in DEGREES
C SIDE EFFECTS:
C    The input RA and DEC are modified to give the values after precession.
C RESTRICTIONS:
C    Accuracy of precession decreases for declination values near 90 
C    degrees.  PRECESS should not be used more than 2.5 centures from
C    1900.    
C PROCEDURE:
C    Algorithm from Computational Spherical Astronomy by Taff (1983), 
C    p. 24. 
C
C	ALL CALCULATIONS ARE DONE IN DOUBLE PRECISION.
C********************************************************************** 

	real*8		CDR,RA_RAD,DEC_RAD,A,CSR,T,ST
	real*8		R(3,3),X(3),X2(3),B,C
	real*8		SINA,SINB,SINC,COSA,COSB,COSC
	real*4		RA,DEC


	CDR = 0.17453292519943D-1

C	Convert to double precision if not already.

	RA_RAD = RA*CDR		
	DEC_RAD = DEC*CDR
	A = DCOS(DEC_RAD)

C	input direction cosines

	X(1) = A*DCOS(RA_RAD)
	X(2) = A*DSIN(RA_RAD)
	X(3) = DSIN(DEC_RAD)
 
	CSR = CDR/3600.D0
	T = 0.001D0*(EPOCH2-EPOCH1)
	ST = 0.001D0*(EPOCH1-1900.D0)

C       Compute 3 rotation angles

	A = CSR*T*(23042.53D0 + ST*(139.75D0 +0.06D0*ST)+
     1		T*(30.23D0 - 0.27D0*ST+18.0D0*T))
	B = CSR*T*T*(79.27D0 + 0.66D0*ST + 0.32D0*T) + A
	C = CSR*T*(20046.85D0 - ST*(85.33D0 + 0.37D0*ST)+
     1		T*(-42.67D0 - 0.37D0*ST -41.8D0*T))
 
	SINA = DSIN(A)
	SINB = DSIN(B)
	SINC = DSIN(C)
	COSA = DCOS(A)
	COSB = DCOS(B)
	COSC = DCOS(C)

	R(1,1) = COSA*COSB*COSC-SINA*SINB
	R(2,1) = SINA*COSB+COSA*SINB*COSC
	R(3,1) = COSA*SINC
	R(1,2) = -COSA*SINB-SINA*COSB*COSC
	R(2,2) = COSA*COSB-SINA*SINB*COSC
	R(3,2) = -SINA*SINC
	R(1,3) = -COSB*SINC
	R(2,3) = -SINB*SINC
	R(3,3) = COSC

C	X2 = R#X	rotate to get output direction cosines

	do i = 1,3
	  X2(i) = 0.D0
	  do j = 1,3
	     X2(i) = X2(i) + R(i,j)*X(j)
	  end do
	end do

	RA_RAD = DATAN2(X2(2),X2(1))

	DEC_RAD = DASIN(X2(3))

	RA = RA_RAD/CDR
	if(RA.LT.0.)RA = RA + 360.
	DEC = DEC_RAD/CDR

	RETURN
	END

