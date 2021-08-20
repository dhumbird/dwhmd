      SUBROUTINE RADIC(KI,KJ,XNT1,XNT2,CONJUG,RAD,DRDL,DRDM,DRDN)
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'
c
C TRICUBIC SPLINE
C
C CONJUG=1: NONCONJUGATED
C       >1: CONJUGATED
C
      L=INT(XNT1)
      M=INT(XNT2)
      N=INT(CONJUG)
      RAD=0.0d0
      DRDL=0.0d0
      DRDM=0.0d0
      DRDN=0.0d0
      KIKJ = KI + KJ - 1
C
      IF(L.GE.4) THEN
           L=4
           XNT1=4.0D0
      ENDIF
C
      IF(M.GE.4) THEN
           M=4
           XNT2=4.0D0
      ENDIF
C
      IF(N.GE.9) THEN
           N=9
           CONJUG=9.0D0
      ENDIF
C
      DO 32 J=1,64
           X=CLMN(KIKJ,L,M,N,J)*
     &       (XNT1**IN3(J,1))*(XNT2**IN3(J,2))*(CONJUG**IN3(J,3))
           RAD=RAD+X
           DRDL=DRDL+X*IN3(J,1)/XNT1
           DRDM=DRDM+X*IN3(J,2)/XNT2
           DRDN=DRDN+X*IN3(J,3)/CONJUG
32    CONTINUE
      RETURN
      END
C

