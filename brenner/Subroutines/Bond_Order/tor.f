      SUBROUTINE TOR(XNT1,XNT2,CONJUG,ATOR,DRDL,DRDM,DRDN)
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
c
C TRICUBIC SPLINE FOR TORSIONAL INTERACTION
C
C CONJUG=1: NONCONJUGATED
C       >1: CONJUGATED
C
      ATOR=0.0d0
      DRDL=0.0d0
      DRDM=0.0d0
      DRDN=0.0d0
      IF((XNT1.GE.4.0D0).OR.(XNT2.GE.4.0D0)) RETURN
      L=INT(XNT1)
      M=INT(XNT2)
      N=INT(CONJUG)
      DO 32 J=1,64
           X=TLMN(L,M,N,J)*
     &       (XNT1**IN3(J,1))*(XNT2**IN3(J,2))*(CONJUG**IN3(J,3))
           ATOR=ATOR+X
           DRDL=DRDL+X*IN3(J,1)/XNT1
           DRDM=DRDM+X*IN3(J,2)/XNT2
           DRDN=DRDN+X*IN3(J,3)/CONJUG
32    CONTINUE
C
      return 
      END

