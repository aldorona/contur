      SUBROUTINE XYZ (XX,YY,YYP,YYPP)
C     COMPUTE Y,Y',Y'' FOR A CURVE DESCRIBED BY CUBIC'S A(5,*)
C     WHERE (1) = X-MAX  (2) = HIGH ORDER COEFFICIENT.
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /COEF/ A(5,200),NA
      DATA ZERO/0.0D+0/
      X=XX
      IF (X .GE. A(1,1)) GO TO 2
1     Y=ZERO
      YP=ZERO
      YPP=ZERO
      GO TO 5
2     DO 3 K=2,200
      IF (X .LE. A(1,K)) GO TO 4
3     CONTINUE
      GO TO 1
4     A3=A(2,K)
      A2=A(3,K)
      A1=A(4,K)
      AZ=A(5,K)
      T=A2+A2
      S=A3*3.0D+0
      R=S+S
      Y=AZ+X*(A1+X*(A2+X*A3))
      YP=A1+X*(T+X*S)
      YPP=T+R*X
5     YY=Y
      YYP=YP
      YYPP=YPP
      RETURN
      END