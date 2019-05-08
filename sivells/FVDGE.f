      SUBROUTINE FVDGE (X,Y,DS,DY)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(5), Y(5)
      DATA H/0.5D+0/,TWO/2.0D+0/
C
      X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      X5=X(5)
C
      Y1=Y(1)
      Y2=Y(2)
      Y3=Y(3)
      Y4=Y(4)
      Y5=Y(5)
C
C     FIND DELTA-Y
      F1=(X3-X1)*(X3-X2)
      F1=TWO/F1
C
      F2=(X4-X3)*(X3-X2)
      F2=-TWO/F2
C
      F3=(X5-X3)*(X4-X3)
      F3=TWO/F3
C
      Z13=X1+X2+X2-X4-X4-X5
      A1=(X2+X3-X4-X5)/Z13
      A3=(X1+X2-X3-X4)/Z13
C
      YP21=(Y2-Y1)/(X2-X1)
      YP32=(Y3-Y2)/(X3-X2)
      YP43=(Y4-Y3)/(X4-X3)
      YP54=(Y5-Y4)/(X5-X4)
C
      X21=H*(X2+X1)
      X32=H*(X3+X2)
      X43=H*(X4+X3)
      X54=H*(X5+X4)
C
      YPP1=(YP32-YP21)/(X32-X21)
      YPP2=(YP43-YP32)/(X43-X32)
      YPP3=(YP54-YP43)/(X54-X43)
      DS=A1*YPP1+A3*YPP3-YPP2
      FX=F2-A1*F1-A3*F3
      DY=DS/FX
C
      RETURN
      END