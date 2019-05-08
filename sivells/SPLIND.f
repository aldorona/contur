      SUBROUTINE SPLIND (X,Y,TN2,TNL,L)
C     COMPUTE CUBIC COEFFICIENTS FOR A CURVE X-Y
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /COEF/ E(5,200),NE
      COMMON /WORK/ A(300),B(300),C(300),D(300),G(300),SB(300),XM(300),D
     1X(300),DY(300)
      DIMENSION X(1), Y(1)
      DATA ZERO/0.0D+0/,ONE/1.D+0/,THR/3.D+0/,SIX/6.D+0/
      CALL OREZ (E,5*200)
      CALL OREZ (A,9*300)
      DX(1)=ZERO
      DY(1)=ZERO
      N=L-1
      DO 1 K=2,L
      DX(K)=X(K)-X(K-1)
1     DY(K)=Y(K)-Y(K-1)
C
      B(1)=DX(2)/THR
      C(1)=DX(2)/SIX
      D(1)=DY(2)/DX(2)-TN2
      A(L)=DX(L)/SIX
      B(L)=DX(L)/THR
      D(L)=TNL-DY(L)/DX(L)
      A(1)=ZERO
      DO 2 K=2,N
      A(K)=DX(K)/SIX
      B(K)=(DX(K)+DX(K+1))/THR
      D(K)=DY(K+1)/DX(K+1)-DY(K)/DX(K)
2     C(K)=DX(K+1)/SIX
      SW=ONE/B(1)
      SB(1)=SW*C(1)
      G(1)=SW*D(1)
      DO 3 K=2,L
      SW=ONE/(B(K)-A(K)*SB(K-1))
      SB(K)=SW*C(K)
3     G(K)=SW*(D(K)-A(K)*G(K-1))
      XM(L)=G(L)
      DO 4 K=1,N
      J=L-K
4     XM(J)=G(J)-SB(J)*XM(J+1)
      DO 5 K=2,L
      DXR=ONE/DX(K)
      Q=DXR/SIX
      P=-XM(K-1)*Q
      Q=Q*XM(K)
      R=DX(K)*XM(K-1)/SIX-DXR*Y(K-1)
      S=Y(K)*DXR-DX(K)*XM(K)/SIX
      XK=X(K)
      PX=XK*P
      PXX=PX*XK
      PXXX=PXX*XK
      XJ=X(K-1)
      QX=XJ*Q
      QXX=QX*XJ
      QXXX=QXX*XJ
      E(2,K)=P+Q
      E(3,K)=-THR*(PX+QX)
      E(4,K)=THR*(PXX+QXX)+R+S
      E(5,K)=-PXXX-QXXX-R*XK-S*XJ
5     CONTINUE
      DO 6 K=2,L
      E(1,K)=X(K)
6     CONTINUE
      E(1,1)=X(1)
      NE=L
      RETURN
      END