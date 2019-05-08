      SUBROUTINE SCOND (A,B,C,KING)
C     TO OBTAIN PARABOLIC DERIVATIVE OF CURVE (UNEQUALLY SPACED POINTS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(150), B(150), C(150)
      N=KING-1
      DO 1 K=2,N
      S=A(K)-A(K-1)
      T=A(K+1)-A(K)
1     C(K)=((B(K+1)-B(K))*S*S+(B(K)-B(K-1))*T*T)/(S*S*T+S*T*T)
      SO=A(2)-A(1)
      T0=A(3)-A(2)
      QO=SO+T0
      C(1)=(-T0*(QO+SO)*B(1)+QO*QO*B(2)-SO*SO*B(3))/QO/SO/T0
      SF=A(KING-1)-A(KING-2)
      TF=A(KING)-A(KING-1)
      QF=SF+TF
      QST=QF*SF*TF
      C(KING)=(SF*(QF+TF)*B(KING)-QF*QF*B(KING-1)+TF*TF*B(KING-2))/QST
      RETURN
      END