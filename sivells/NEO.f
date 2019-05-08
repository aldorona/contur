      SUBROUTINE NEO
C
C     SMOOTH BY LINEAR SECOND DERIVATIVE
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /WORK/ E(400),Z(400),X(400),Y(400),YST(400),WTN(250),WALL(5
     1,200),WAX(200),WAY(200),WAN(200)
      COMMON /CONTR/ ITLE(3),IE,LR,IT,JB,JQ,JX,KAT,KBL,KING,KO,LV,NOCON,
     1IN,MC,MCP,IP,IQ,ISE,JC,M,MP,MQ,N,NP,NR,NUT,NF
      DATA ZERO/0.0D+0/,ONE/1.D+0/,TWO/2.D+0/
      DATA J0/4H  UP/,J1/4HDOWN/
      INTEGER*4 NOUP,NPCT,NODO
C
      CONV=90.D+0/DASIN(ONE)
      TNI=DTAN(WALL(5,1))
C
      IF (JQ.EQ.0.OR.IQ.LT.0) READ (1,14,END=13) NOUP,NPCT,NODO
      NOUP=50
      NPCT=85
      NODO=50
C     IF ((JQ.EQ.0) .OR. (IQ.LT.0)) READ (5,14,END=13)NOUP,NPCT,NODO
      IF (JQ .GT. 0) GO TO 2
      JN=J0
      LIM=NUT
      NOTM=NOUP
      DO 1 J=1,LIM
      X(J+1)=WAX(J)
      Y(J+1)=WAY(J)
1     YST(J+1)=Y(J+1)
      X(1)=TWO*X(2)-X(3)
      Y(1)=Y(3)
      X(LIM+2)=TWO*X(LIM+1)-X(LIM)
      Y(LIM+2)=Y(LIM+1)+TNI*(X(LIM+2)-X(LIM+1))
      GO TO 4
2     LIM=N+NP-1
      NOTM=NODO
      JN=J1
      DO 3 J=1,LIM
      X(J+1)=WALL(1,J)
      Y(J+1)=WALL(2,J)
3     YST(J+1)=Y(J+1)
      X(1)=TWO*X(2)-X(3)
      Y(1)=Y(2)-TNI*(X(2)-X(1))
      X(LIM+2)=TWO*X(LIM+1)-X(LIM)
      Y(LIM+2)=Y(LIM+1)
4     LUS=1+(LIM-3)/6
      IF (NOTM .EQ. 0) RETURN
      YST(1)=Y(1)
      YST(LIM+2)=Y(LIM+2)
      SMP=1.D-2*NPCT
      WRITE (6,16) ITLE,JN,NOTM,SMP
C
      DO 8 M=1,NOTM
      CALL OREZ(E,800)
C
      DO 5 K=3,LIM
      CALL FVDGE (X(K-2),Y(K-2),E(K),Z(K))
5     CONTINUE
      E(1)=ZERO
      E(2)=ZERO
      E(LIM+1)=ZERO
      E(LIM+2)=ZERO
C     SEARCH ARRAY AND FIND MAX ERR
      DO 7 LU=1,LUS
      EMAX=ZERO
      DO 6 K=3,LIM
      TEST=DABS(E(K))
      IF (EMAX .GT. TEST) GO TO 6
      J=K
      EMAX=TEST
6     CONTINUE
C     APPLY CORRECTION
      E(J)=ZERO
      E(J+1)=ZERO
      E(J+2)=ZERO
      E(J-1)=ZERO
      E(J-2)=ZERO
      Y(J)=Y(J)+SMP*Z(J)
7     CONTINUE
8     CONTINUE
C
      ERR=ZERO
      DO 9 J=1,LIM
      K=J+1
      E(K)=Y(K)-YST(K)
      IF (ERR .LT. DABS(E(K))) MAX=J 
      IF (ERR .LT. DABS(E(K))) ERR=DABS(E(K))
      WRITE (6,15) J,X(K),Y(K),YST(K),E(K),J
9     IF (MOD(J,10) .EQ. 0) WRITE (6,17)
      WRITE (6,19) ERR,MAX
C
      LM=LIM-1
      CALL SCOND (X,Y,WTN,LIM+2)
      IF (JQ .EQ. 1) GO TO 11
      DO 10 J=2,LM
      WAY(J)=Y(J+1)
10    WAN(J)=CONV*DATAN(WTN(J+1))
      RETURN
C
11    DO 12 J=2,LM
      WALL(2,J)=Y(J+1)
12    WALL(5,J)=DATAN(WTN(J+1))
      RETURN
C
13    WRITE (6,18)
      STOP
C
14    FORMAT (3I5)
15    FORMAT (1H,20X,I5,2X,0P4F13.7,I8)
16    FORMAT (1H1,3A4,2X,A4,24HSTREAM CONTOUR, SMOOTHED,I5,19H TIMES WIT
     1H FACTOR=,F4.2
     2//34X,1HX,11X,6HY-CALC,7X,4HY-IN,10X,4HDIFF /)
17    FORMAT (1H )
18    FORMAT (1H0,10X,34HCARD NOT AVAILABLE FOR NEGATIVE NF)
19    FORMAT (1H0,26X,21HMAX. ABSOLUTE ERROR=,1PG15.6,10H AT POINT,I5)
      END