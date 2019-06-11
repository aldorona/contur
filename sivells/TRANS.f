      SUBROUTINE TRANS (RTO,TK,WO,AMN,AMP,AMPP,W,AWP,AWPP,CWOPPP,AXN)
C     TO DETERMINE THROAT CHARACTERISTIC
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /GG/ GAM,GM,G1,G2,G3,G4,G5,G6,G7,G8,G9,GA,RGA,QT
      COMMON /CONTR/ ITLE(3),IE,LR
      COMMON /TROAT/ FC(6,51)
      DATA ZRO/0.0D+0/,ONE/1.D+0/,TWO/2.D+0/,SIX/6.D+0/,HALF/5.D-1/
      DATA TRHV/1.5D+0/,THR/3.D+0/,FOUR/4.D+0/,EIT/8.D+0/,TLV/1.2D+1/
      NN=IABS(LR)
      JJ=240/(NN-1)
      IF (MOD(JJ,2).NE.0) JJ=JJ+1
      IF (JJ.LT.10) JJ=10
      KK=JJ*NN-JJ
      GB=IE/EIT
      GK=(GAM*(GAM+2.25D+0*IE-16.5D+0)+2.25D+0*(2+IE))/TLV 
      GU=ONE-GAM/TRHV
      GV=(HALF*(5-3*IE)*GAM+IE)/(9-IE)
      GZ=DSQRT(QT*(GAM+ONE))
      U22=GB+GAM/THR/(3-IE)
      U42=(GAM+(4-IE)*TRHV)/SIX/(3-IE)
      IF (IE.EQ.0) GO TO 1
      GT=(GAM*(GAM*92.D+0+180.D+0)-9.D+0)/1152.D+0
      U23=(GAM*(304.D+0*GAM+255.D+0)-54.D+0)/1728.D+0
      U43=(GAM*(388.D+0*GAM+777.D+0)+153.D+0)/2304.D+0
      U63=(GAM*(556.D+0*GAM+1737.D+0)+3069.D+0)/10368.D+0
      UP0=(GAM*(52.D+0*GAM+75.D+0)-9.D+0)/192.D+0 
      UP2=(GAM*(52.D+0*GAM+51.D+0)+327.D+0)/384.D+0 
      V02=(28.D+0*GAM-15.D+0)/288.D+0 
      V22=(20.D+0*GAM+27.D+0)/96.D+0 
      V42=(GAM/THR+ONE)/THR 
      V03=(GAM*(7100.D+0*GAM+2151.D+0)+2169.D+0)/82944.D+0 
      V23=(GAM*(3424.D+0*GAM+4071.D+0)-972.D+0)/13824.D+0 
      V43=(GAM*(3380.D+0*GAM+7551.D+0)+3771.D+0)/13824.D+0 
      V63=(GAM*(6836.D+0*GAM+23031.D+0)+30627.D+0)/82944.D+0 
      GO TO 2
1     GT=(GAM*(GAM*134.D+0+429.D+0)+123.D+0)/4320.D+0
      U23=(GAM*(854.D+0*GAM+807.D+0)+279.D+0)/12960.D+0
      U43=(GAM*(194.D+0*GAM+549.D+0)-63.D+0)/2592.D+0
      U63=(GAM*(362.D+0*GAM+1449.D+0)+3177.D+0)/12960.D+0
      UP0=(GAM*(26.D+0*GAM+51.D+0)-27.D+0)/144.D+0 
      UP2=(GAM*(26.D+0*GAM+27.D+0)+237.D+0)/288.D+0 
      V02=(34.D+0*GAM-75.D+0)/1080.D+0 
      V22=(10.D+0*GAM+15.D+0)/108.D+0 
      V42=(22.D+0*GAM+75.D+0)/360.D+0 
      V03=(GAM*(7570.D+0*GAM+3087.D+0)+23157.D+0)/544320.D+0 
      V23=(GAM*(5026.D+0*GAM+7551.D+0)-4923.D+0)/77760.D+0 
      V43=(GAM*(2254.D+0*GAM+6153.D+0)+2979.D+0)/25920.D+0 
      V63=(GAM*(6574.D+0*GAM+26481.D+0)+40059.D+0)/181440.D+0
2     WWO=WO+(HALF+(U42-U22+(U63-U43+U23)/RTO)/RTO)/RTO
      WOP=(ONE-(GB-GT/RTO)/RTO)/DSQRT(RTO)
      WOPP=(GU-GV/RTO)/RTO
      HOPPP=GK/RTO/DSQRT(RTO)
      HVPPP=(3*IE-(10-3*IE)*GAM)/FOUR/RTO/DSQRT(RTO)
      AMN=WWO/DSQRT(G7-G8*WWO**2)
      BET=DSQRT(AMN**2-ONE)
      PSI1=G2*DATAN(BET/G2)-DATAN(BET)
      P1=ZRO
      T1=ZRO
      X1=ZRO
      Y1=ONE
      FSY1=ZRO
      TN2=-ONE/BET
      FC(1,NN)=X1
      FC(2,NN)=Y1
      FC(3,NN)=AMN
      FC(4,NN)=PSI1
      FC(5,NN)=ZRO
      FC(6,NN)=ZRO
      BX=ONE
      SUM=ZRO
      FSA=(IE+1)*AMN/(G6+G5*AMN**2)**GA
      DO 8 J=1,KK
      Y=DFLOAT(KK-J)/KK
      IF (IE .EQ. 1) BX=Y+Y
      YY=Y*Y
      TN1=TN2
      VO=(((YY*(YY*(YY*V63-V43)+V23)-V03)/RTO+YY*(YY*V42-V22)+V02)/RTO+H
     1ALF*(YY-ONE)/(3-IE))/RTO
      VP=(ONE+((YY*(TWO*GAM+3*(4-IE))-TWO*GAM-TRHV*IE)/(3-IE)/THR+(YY*(S
     1IX*U63*YY-FOUR*U43)+TWO*U23)/RTO)/RTO)/DSQRT(RTO)
      VPP=TWO*(ONE+(TWO*UP2*YY-UP0)/RTO)/RTO
C     ITERATE FOR X AND MACH NUMBER FROM CHARACTERISTIC EQUATIONS
      DO 4 I=1,10
      TNA=HALF*(TN1+TN2)
      X=X1+(Y-Y1)/TNA
      DXI=DSQRT((Y-Y1)**2+(X-X1)**2)
      XOT=X/GZ
      VY=GZ*(VO+XOT*(VP+XOT*(HALF*VPP+XOT*HVPPP/THR)))/DSQRT(RTO)
      W=AMN/DSQRT(G6+G5*AMN**2)
      T=DASIN(VY*Y/W)
      FSY=IE*VY/W/AMN
      P1=HALF*(FSY1+FSY)*DXI
3     PSI=P1+PSI1+T1-T
      FMA=FMV(PSI)
      IF (DABS(AMN-FMA) .LT. 1.D-10) GO TO 5
      FMU=DASIN(ONE/FMA)
      TN2=DTAN(T-FMU)
      AMN=FMA
4     CONTINUE
C     ITERATION COMPLETE
5     IF (MOD(J,2) .EQ. 0) GO TO 6
      AS=Y1-Y
      FSB=BX/DSIN(FMU-T)/(G6+G5*FMA**2)**GA 
      GO TO 7
6     BS=Y1-Y
      CS=AS+BS
      S1=(TWO-BS/AS)*CS/SIX
      S3=(TWO-AS/BS)*CS/SIX 
      S2=CS-S1-S3
      FSC=BX/DSIN(FMU-T)/(G6+G5*FMA**2)**GA
      ADD=S1*FSA+S2*FSB+S3*FSC
      SUM=ADD+SUM
      FSA=FSC
7     X1=X
      Y1=Y
      T1=T
      FSY1=FSY
      PSI1=PSI
      IF (MOD(J,JJ) .NE. 0) GO TO 8
      K=NN-J/JJ
      FC(1,K)=X
      FC(2,K)=Y
      FC(3,K)=FMA
      FC(4,K)=PSI
      FC(5,K)=T
      FC(6,K)=SUM
8     CONTINUE
      DO 9 J=1,NN
      FC(1,J)=FC(1,J)/TK
      FC(2,J)=FC(2,J)/TK
9     FC(6,J)=ONE-FC(6,J)/SUM
      AXN=FC(1,1)
      AWOP=WOP*TK/GZ
      AWOPP=WOPP*(TK/GZ)**2
      AWOPPP=TWO*HOPPP*(TK/GZ)**3
      CWOPPP=SIX*(W-WO-AXN*(AWOP+AXN*AWOPP/TWO))/AXN**3
      IF (CWOPPP .LT. AWOPPP) CWOPPP=AWOPPP
      AWP=AWOP+AXN*(AWOPP+AXN*CWOPPP/TWO)
      AWPP=AWOPP+AXN*CWOPPP
      AMP=AWP*G7*(AMN/W)**3
      AMPP=AMP*(AWPP/AWP+THR*G5*AMP*W*W/AMN)
      IF (LR .GT. 0) RETURN
      LR=NN
      RC=RTO-ONE
      WRITE (2,12) ITLE,RC,AWOP,AWOPP,AWOPPP
      DO 10 J=1,NN
      Y=DFLOAT(J-1)/(NN-1)
      YY=Y*Y
      Y4=YY**2
      Y6=YY**3
      DUY=(HALF*YY+(U42*Y4-U22*YY+(U63*Y6-U43*Y4+U23*YY)/RTO)/RTO)/RTO
      UY=WO+DUY
      VO=(((YY*(YY*(YY*V63-V43)+V23)-V03)/RTO+YY*(YY*V42-V22)+V02)/RTO+H
     1ALF*(YY-ONE)/(3-IE))/RTO
      VY=GZ*VO*Y/DSQRT(RTO)
      WY=DSQRT(UY**2+VY**2)
      YM=WY/DSQRT(G7-G8*WY**2)
      WRITE (2,13) Y,UY,VY,WY,YM
10    IF (MOD(J,10) .EQ. 0) WRITE (2,14)
      XX1=CUBIC(CWOPPP/SIX,AWOPP/TWO,AWOP,WO-ONE)
      XXI=CUBIC(AWOPPP/SIX,AWOPP/TWO,AWOP,WO-W)
      WRITE (2,15) XX1,XXI,W,CWOPPP,TK
      WRITE (2,16)
      PX=AXN+1.D-1
      DO 11 J=1,11
      X=.1D+0*(J-1)
      XW=WO+X*(AWOP+X*(AWOPP/TWO+X*CWOPPP/SIX))
      XWP=AWOP+X*(AWOPP+X*CWOPPP/TWO)
      XWPP=AWOPP+X*CWOPPP
      XM=XW/DSQRT(G7-G8*XW**2)
      XMP=XWP*G7*(XM/XW)**3
      XMPP=XMP*(XWPP/XWP+THR*G5*XMP*XW*XW/XM)
      IF (X.LT.AXN .OR. X.GT.PX) GO TO 11
      WRITE (2,18) AXN,W,AWP,AWPP,AMN,AMP,AMPP
11    WRITE (2,17) X,XW,XWP,XWPP,XM,XMP,XMPP
      RETURN
C
12    FORMAT (1H1,8X,3A4,39H THROAT VELOCITY DISTRIBUTION, X=O, RC=,F10.
     16//10X,44HDERIVATIVES TAKEN WITH RESPECT TO X/Y*, WOP=,F11.8//10X,
     25HWOPP=,1PE15.7,5X,6HWOPPP=,E15.7//10X,4HY/YO,7X,4HU/A*,10X,4HV/A*
     3,11X,1HW,11X,8HMACH NO. /)
13    FORMAT (1H ,F14.4,4F14.8 )
14    FORMAT (1H )
15    FORMAT (1H0,9X,18HFROM CUBIC, X/Y* =,F11.8,11H FOR W= 1.0 //22X,6H
     1X/Y* =,F11.8,7H FOR W=,F11.8 //10X,16HCORRECTED WOPPP=,1PE15.7 // 
     210X,15HRMASS = Y*/YO =,0PF13.10 //)
16    FORMAT (1H0,9X,32HAXIAL VELOCITY DISTRIBUTION, Y=0 //10X,4HX/Y*,9X
     1,1HW,17X,2HWP,16X,3HWPP,15X,1HM,17X,2HMP,16X,3HMPP /)
17    FORMAT (1H ,F13.3,1P6E18.7 )
18    FORMAT (1H ,F16.8,1PE15.7,5E18.7 )
      END
