      SUBROUTINE OFELD (A,B,C,NOCON)
C     TO OBTAIN POINTS IN CHARACTERISTIC NETWORK
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CONTR/ ITLE(3),IE
      DATA ZRO/0.0D+0/,ONE/1.D+0/,TWO/2.D+0/,HALF/5.D-1/
      DIMENSION A(5), B(5), C(5)
      A1=DASIN(ONE/A(3))
      A2=DASIN(ONE/B(3))
      T1=A(5)
      T2=B(5)
      IF (IE .EQ. 0) GO TO 8
      IF (A(2) .EQ. ZRO) GO TO 5
      FSY1=DSIN(A(5))/A(2)/A(3)
      GO TO 6
5     T1=ZRO
      FSY1=A(5)
6     IF (B(2) .EQ. ZRO) GO TO 7
      FSY2=DSIN(B(5))/B(2)/B(3)
      GO TO 8
7     T2=ZRO
      FSY2=B(5)
8     TNI=DTAN(T1-A1)
      IF (B(3) .NE. ONE) TN2=DTAN(T2+A2)
      I=-1
      HDPSI=HALF*(A(4)-B(4))
      HT3=HALF*(T1+T2)+HDPSI
      T3=HT3-HALF*IE*HDPSI
      HPSI3=HALF*(A(4)+B(4)+T1-T2)
      PSI3=HPSI3+HALF*IE*(T1-T2)
      C(3)=FMV(PSI3)
      TOLD=T3
1     I=I+1
      FM3=C(3)
      A3=DASIN(ONE/C(3))
      TNA=HALF*(TNI+DTAN(T3-A3))
      IF (B(3) .NE. ONE) TNB=HALF*(DTAN(T3+A3)+TN2)
      IF (B(3) .EQ. ONE) TNB=TWO*DTAN(T3+A3)
      DTN=TNB-TNA
      X3=(B(1)*TNB-A(1)*TNA+A(2)-B(2))/DTN
      Y3=(A(2)*TNB-B(2)*TNA+(B(1)-A(1))*TNA*TNB)/DTN
      IF ((IE .EQ. 0) .OR. (DABS(Y3) .LT. 1.D-9)) GO TO 4
      FSY3=DSIN(T3)/Y3/FM3
      P1=HALF*(FSY1+FSY3)*(X3-A(1))*DSQRT(ONE+TNA**2)
      P2=HALF*(FSY2+FSY3)*(X3-B(1))*DSQRT(ONE+TNB**2)
      T3=HT3+HALF*(P1-P2)
      PSI3=HPSI3+HALF*(P1+P2)
      C(3)=FMV(PSI3)
      IF (DABS(T3-TOLD) .GT. 1.D-9) GO TO 2
      IF (DABS(C(3)-FM3) .LT. 1.D-9) GO TO 4
2     IF (I .EQ. 40) GO TO 3
      TEMP=T3
      T3=(T3+TOLD)*HALF
      TOLD=TEMP
      GO TO 1
3     NOCON=1
4     C(1)=X3
      C(2)=Y3
      C(4)=PSI3
      C(5)=T3
      RETURN
      END