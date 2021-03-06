      SUBROUTINE BOUND
C
C     TO OBTAIN THE CORRECTION DUE TO THE TURBULENT BOUNDARY LAYER
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /GG/ GAM,GM,G1,G2,G3,G4,G5,G6,G7,G8,G9,GA,RGA,QT
      COMMON /CORR/ DLA(200),RCO(200),DAX(200),DRX(200),SL(200),DR2
      COMMON /COORD/ S(200),FS(200),WALTAN(200),SD(200),WMN(200),TTR(200
     1),DMDX(200),SPR(200),BTA(200),SREF(200),XBIN,XCIN,GMA,GMB,GMC,GMD
      COMMON /PROP/ AR,ZO,RO,VISC,VISM,SFOA,SBL,CONV
      COMMON /PARAM/ ETAD,RC,AMACH,BMACH,CMACH,EMACH,GMACH,FRC,SF,WWO,WW
     1OP,QM,WE,CBET,XE,ETA,EPSI,BPSI,XO,YO,RRC,SDO,XB,XC,AH,PP,SE,TYE,XA
      COMMON /HTTR/ HAIR,TAW,TWQ,TW,TWAT,QFUN,QFUNW,IPQ,IJ,IV,IW
      COMMON /CONTR/ ITLE(3),IE,LR,IT,JB,JQ,JX,KAT,KBL,KING,KO,LV,NOCON,
     1IN,MC,MCP
      DIMENSION Z(16), D(16), SCV(200), SK(200), CDS(200), RW(200)
      DATA ZRO/0.0D+0/,ONE/1.D+0/,TWO/2.D+0/,SIX/6.D+0/,HALF/5.D-1/
      DATA THR/3.D+0/,FOUR/4.D+0/,TEN/1.D+1/,TLV/1.2D+1/
      DATA CF1/3.865D-2/,CF2/4.561D+0/,CF3/5.46D-1/,FSI/3.17897971D+0/
      DATA LY/4H   Y/,LS/4H   S/,DD/8HD2Y/DX2 /,DK/8H CURV.  /
      DATA Z(1)/.052995325D-1/,Z(4)/.1222977958D+0/,Z(7)/.3591982246D+0/
      DATA Z(2)/.277124885D-1/,Z(5)/.1910618778D+0/,Z(8)/.4524937451D+0/
      DATA Z(3)/.671843988D-1/,Z(6)/.2709916112D+0/
      DATA D(1)/.135762297D-1/,D(2)/.31126762D-1/,D(3)/.475792558D-1/
      DATA D(4)/.623144856D-1/,D(5)/.747979944D-1/,D(6)/.845782597D-1/
      DATA D(7)/.913017075D-1/,D(8)/.947253052D-1/
      DO 1 J=9,16
      D(J)=D(17-J)
1     Z(J)=ONE-Z(17-J)
      DO 2 J=1,KAT
2     SREF(J)=S(J)
      SBIN=XBIN
      SCIN=XCIN
      TRPI=CONV/90.D+0
      FCC=2.05D+0+DLOG(.41D+0)
      CHAIR=GAM*G1*AR/RO/RO/777.64885D+0
      IF (IT .EQ. 0) XBL=SBL
C
3     READ (1,66,END=65) PPQ,T0,TWT,TWAT,QFUN,ALPH,IHT,IR,LD,LV
C
      PPS=PPQ
      RHO=144.D+0*PPS/ZO/AR/T0
      ID=IABS(LD)
      KOR=KO
      IF (IABS(IN) .EQ. 10) KOR=1
      IF (MCP .LT. 0) KOR=KING
      ROY=ONE
      IF (IE .EQ. 0) HW=AH
      IF ((ID .EQ. 0) .OR. (IE .EQ. 1)) HW=ZRO
      IF (HW .EQ. ZRO) YOH=ZRO
      IF (HW .EQ. ZRO) YOHA=ZRO 
      ALF=DABS(ALPH)
      ARC=FRC
      IF (IHT .LT. 0) ARC=FRC**(IE+1)
      IPQ=0
      IW=1
      IF (LV .NE. 0) IW=IABS(LV)
      DO 4 J=1,KAT
      S(J)=SREF(J)
      SL(J)=S(J)
      RW(J)=FS(J)
      RCO(J)=FS(J)
      SCW=DSQRT(ONE+WALTAN(J)**2)
      SK(J)=SD(J)/SCW**3
      IF (KAT .EQ. KING) GO TO 4
      IF (S(J) .LT. SBL) KBL=J+2
4     DRX(J)=WALTAN(J)
      IF (KBL .GT. KAT) KBL=KAT+4
5     DO 58 IV=1,IW
      IF ((IV .GT. 1) .AND. (IV .LT. IW)) GO TO 15
      IF (LD .GE. 0) WRITE (2,80) ITLE,PPS,T0
      IF (ALPH .GT. ZRO) GO TO 6
      ALPHA=ZRO
      IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,71)
      GO TO 7
6     ALPHA=ALPH
      IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,70)
7     IF (IR .EQ. 2) GO TO 13
      IF (ALF .EQ. ONE) GO TO 8
      IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,75)
      GO TO 9
8     IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,72)
9     IF (IR) 10,11,12
10    IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,74)
      GO TO 14
11    IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,73)
      GO TO 14
12    IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,76)
      GO TO 14
13    IF ((LD .GE. 0) .OR. (PPQ .EQ. ZRO)) WRITE (2,77)
14    IF (PPQ .EQ. ZRO) GO TO 60
15    CAPI=.55D+0
      IPP=0
      IJ=1
      DO 56 J=1,KAT
      BET=TTR(J)-ONE
      STR=ONE/TTR(J)
      TE=T0*STR
      RAJ=WMN(J)*(G7*STR)**GA
      IF (IHT .GE. 0) RAJ=RAJ**QT
      SCW=DSQRT(ONE+DRX(J)**2)
      EMU=VISC*TE*DSQRT(TE)/(TE+VISM)
      IF (TE .LT. VISM) EMU=HALF*VISC*TE/DSQRT(VISM)
      IF (VISM .LE. ONE) EMU=VISC*TE**VISM
      TAW=TE*(ONE+RO*BET)
      RHOE=RHO*STR**G1
      VE=WMN(J)*DSQRT(GAM*AR*TE)
      REO=RHOE*VE/EMU/TLV
      IF (HW .GT. ZRO) YOH=FS(J)/HW
      IF (IE .EQ. 0 .AND. HW .GT. ZRO) ROY=(HW/FS(J)+ONE)*TRPI
      K=J
      IF (J .EQ. 1) GO TO 19
      IF (J .GT. KOR) K=J-KOR+1
      IF (K-3) 16,17,18
16    DS=S(J)-S(J-1)
      SMD=HALF*DS
      GO TO 19
17    DT=S(J)-S(J-1)
      DST=DS+DT
      SMA=DST*(TWO-DT/DS)/SIX
      SMC=DST*(TWO-DS/DT)/SIX
      SMB=DST-SMA-SMC
      HB=H
      IF (IV .GT. 1) GO TO 19
      BMA=TWO/DS/DST
      BMB=-TWO/DS/DT
      BMC=TWO/DT/DST
      GO TO 19
18    DU=S(J)-S(J-1)
      DT=S(J-1)-S(J-2) 
      DS=S(J-2)-S(J-3)
      DST=DS+DT
      DSTU=DST+DU
      DTU=DT+DU
      DUT=DU-DT
      DTS=DS-DT
      DTUS=DT+TWO*(DU-DS)
      DTSU=DT+TWO*(DS-DU)
      DSTTU=TWO*(DST+DTU)
      HA=HB
      HB=H
      QMA=HALF*DS*(ONE-DS*(THR+(DTU+DU)/DST)/DSTU/SIX)
      QMB=HALF*DS*(ONE+DS*(TWO+(DST+DT)/DTU)/DT/SIX)
      QMC=-DS**3*(ONE+(DTU+DU)/DST)/DT/DU/TLV
      QMD=DS**3*(DST+DT)/DU/DTU/DSTU/TLV
      SMA=HALF*DS+(DUT*DTU**3/DS-DS*DS*(DS+DSTTU))/DST/DSTU/TLV
      SMB=HALF*DST+(DS*DS*(DSTTU-DS)/DT+DT*DT*DTUS/DS-DU**3*(DSTU+DST)/D
     1S/DT)/DTU/TLV
      SMC=HALF*DTU+(DT*DT*DTSU/DU+DU*DU*(DSTTU-DU)/DT-DS**3*(DSTU+DTU)/D
     1T/DU)/DST/TLV
      SMD=HALF*DU+(DTS*DST**3/DU-DU*DU*(DU+DSTTU))/DTU/DSTU/TLV
19    IF (TWT .NE. ZRO) GO TO 20
      TW=TAW
      GO TO 21
20    TWD=(ARC*RAJ-ONE)*(TWT-TWAT)/(ARC-ONE)
      IF (TWD .LT. ZRO) TWD=ZRO
      TW=TWD+TWAT
21    WMU=VISC*TW*DSQRT(TW)/(TW+VISM)
      IF (VISM .LE. ONE) WMU=VISC*TW**VISM
      DL=TW/TE
      DM=ALPHA*(TAW-TW)/TE
      DN=ONE-DL-DM
      DA=ALF*(TAW-TW)
      DB=DA+TW-TE
      IF (DB) 22,23,24
22    DG=DSQRT(-DB*TE)
      DH=DSQRT(-DB*TW)
      DI=(TWO*(DG+TE-TW)-DA)/(TWO*DH+DA)
      DJ=DLOG(DI)
      TP=-DB/DJ/DJ
      GO TO 25
23    TP=(DSQRT(TE)+DSQRT(TW))**2/FOUR
      GO TO 25
24    DC=DSQRT(DA*DA+FOUR*TW*DB)
      DF=DASIN((DB+TW-TE)/DC)
      DE=DASIN(DA/DC)
      TP=DB/(DF+DE)/(DF+DE)
25    IF (IR) 26,27,28 
26    FRD=TW*EMU/WMU/TP
      GO TO 29
27    FRD=EMU/WMU
      GO TO 29
28    FRD=TE*EMU/WMU/DSQRT(TP*TW)
29    IF (IPP .GT. 0) GO TO 31
      RTHI=1.D-2*REO*FS(1)
      RTII=RTHI
      RDLI=TEN*RTHI
      IF (IR .EQ. 1) GO TO 32
30    RTHG=DLOG10(RTHI)
      CFI=CF1/(RTHG+CF2)/(RTHG-CF3)
31    IF (IR .NE. 2) GO TO 33
      SCFI=DSQRT(CFI)
      TC=TW+17.2D+0*SCFI*DA-305.D+0*CFI*DB
      CMU=VISC*TC*DSQRT(TC)/(TC+VISM)
      IF (VISM .LE. ONE) CMU=VISC*TC**VISM
      TP=TW*CMU/WMU
      FRD=EMU/CMU
      GO TO 33
32    RDLG=DLOG10(RDLI)
      CFI=0.0444D+0/(RDLG+4.6221D+0)/(RDLG-1.4402D+0)
33    CF=CFI*TE/TP
      CFS=CF*SCW
      RTIG=DLOG10(RTII)
      XCF=.41D+0*DSQRT((RTIG+CF2)*(RTIG-CF3)/CF1)
34    C3=TWO+CAPI*(FSI+1.5D+0*CAPI)
      C2=ONE+CAPI
      C1=C2-C3/XCF
      FXCF=XCF+DLOG(C1/RTII)-FCC-TWO*CAPI
      FPCP=(XCF-FSI-THR*CAPI)/XCF/C1-TWO
      CAPI=CAPI-FXCF/FPCP
      IF (DABS(FXCF) .GT. 1.D-8) GO TO 34
      DOTI=XCF/C1
      XN=HALF*(DOTI+DSQRT(DOTI*(DOTI-SIX)+ONE)-THR)
      HI=ONE+TWO/XN
      SUMA=ZRO
      SUMB=ZRO
      SUMC=ZRO
      SUMD=ZRO
      DO 35 L=1,16
      UN=Z(L)**XN
      TR=DL+Z(L)*(DM+Z(L)*DN)
      ADD=D(L)*XN*UN/TR
      BDD=ADD*Z(L)
      CDD=ADD*UN
      DDD=BDD*UN
      SUMA=SUMA+ADD
      SUMB=SUMB+BDD
      SUMC=SUMC+CDD
35    SUMD=SUMD+DDD
      DOT=ONE/(SUMA-SUMB)
      DSOD=ONE-SUMA
      DSM=HALF-SUMC
      THM=SUMC-SUMD
      HU=DSOD*DOT
      IF (IPP .GT. 0) GO TO 36
      H=HU
      DOTR=DOT
36    FMY=(H+TWO-G9*BET)*DMDX(J)*STR/WMN(J)+ID*DRX(J)/(RW(J)+HW)
      IF (J.EQ.1) TH=CFS/FMY
      IF (K.EQ.2) TH=(THA+SMD*(DTHA+CFS))/(ONE+SMD*FMY)
      IF (K.EQ.3) TH=(THA+SMA*DTHA+SMB*DTHB+SMC*CFS)/(ONE+SMC*FMY)
      IF (K.GT.3) TH=(THA+SMA*DTHA+SMB*DTHB+SMC*DTHC+SMD*CFS)/(ONE+SMD*F
     1MY)
      DELST=H*TH
      ASEC=DELST+DSQRT(ID*DELST**2+(FS(J)*SCW*ROY)**2)
      DOR=ID*DOTR*TH/ASEC
      DSROD=DSOD-DOR*DSM
      IPP=1
      DOTR=ONE/(ONE/DOT-THM*DOR)
      HR=DSROD*DOTR
      IF (DABS(H-HR) .LT. 5.D-7) GO TO 37
      H=HR
      GO TO 36
37    DELTA=DOTR*TH
      THU=DELTA/DOT
      DSU=DELTA*DSOD
      RDEL=REO*DELTA
      RTII=RDEL/DOTI
      RDLX=FRD*RDEL
      RTHX=RDLX/DOT
      IF (RTHX .LT. 100.D+0) GO TO 38
      IF (IR .EQ. 1) GO TO 39
      IF (DABS(ONE-RTHX/RTHI) .LT. 1.D-6) GO TO 41
      RTHI=RTHX
      GO TO 30
38    WRITE (2,88) RTHX,REO,FRD,TH,DELTA,DOT
      RETURN
39    IF (DABS(ONE-RDLX/RDLI) .LT. 1.D-6) GO TO 40
      RDLI=RDLX
      GO TO 32
40    RTHG=HALF*(DSQRT((CF2+CF3)**2+FOUR*CF1/CF1)-CF2+CF3)
      RTHX=TEN**RTHG
41    IF (J .GT. 1) GO TO 42
      DTH=ZRO
      HAIR=RHOE*VE*CF*CHAIR
      TAIR=HAIR
      IF (TWAT.EQ.TWT .OR. QFUN.EQ.ZRO) GO TO 46
      TWQ=(HAIR*TAW+QFUN*(TWAT-15.D+0))/(HAIR+QFUN)
      CALL HEAT
      IF (IPQ .GT. 100) GO TO 65
      IF (DABS(TW-TWQ).LT.1.D-2.AND.DABS(QFUN-QFUNW).LT.1.D-5) GO TO 46
      TWT=TWAT+(TWQ-TWAT)*(ARC-ONE)/(ARC*RAJ-ONE)
      QFUN=QFUNW
      GO TO 20
42    DTH=CFS-TH*FMY
      IF (DTH .LT. ZRO) DTH=ZRO
      IF (J .EQ. KOR) GO TO 46
      IF (K-3) 43,45,44
43    DTHB=DTH
      GO TO 47
44    THA=THA+QMA*DTHA+QMB*DTHB+QMC*DTHC+QMD*DTH
      DTHA=DTHB
      DTHB=DTHC
      IF (K .GT. 5) GO TO 45
      SCU=DSQRT(ONE+DRX(J-2)**2)
      DELA=HA*THA
      IF ((IE .EQ. 1).OR.(ID .EQ. 0)) YSEC=FS(J-2)*SCU
      IF ((IE .EQ. 0).AND.(HW .GT. ZRO)) YSEC=SCU*(FS(J-2)+HW)*TRPI
      IF (HW .GT. ZRO) YOHA=FS(J-2)/HW
      ASCA=DELA+DSQRT(ID*DELA**2+YSEC**2)
      RW(J-2)=ASCA/SCU
      DLA(J-2)=SCU*(ASCA-YSEC)*(ONE+YOHA)
      RCO(J-2)=FS(J-2)+DLA(J-2)
45    DTHC=DTH
      GO TO 47
46    THA=TH
      DTHA=DTH
      IF ((IV .GT. 1) .AND. (IV .LT. IW)) GO TO 47
      IF (J .EQ. 1 .AND. LD .GE. 0) WRITE (2,82)
47    CDS(J)=ASEC-SCW*FS(J)*ROY
      DLA(J)=SCW*CDS(J)*(ONE+YOH)
      RCO(J)=FS(J)+DLA(J)
      RW(J)=ASEC/SCW
      IF (IV .LT. IW) GO TO 48
      BTA(J)=-DMDX(J)*DSU/WMN(J)/TTR(J)/SCW/CFI
      IF (J.EQ.1 .OR. J.GT.KO .OR. IHT.EQ.0) GO TO 48
      IF (MOD(J,IHT) .NE. 1) GO TO 48
      IJ=J
      HAIR=RHOE*VE*CF*CHAIR
      CALL HEAT
48    IF (LD.LT.0) GO TO 56
      IF ((IV.GT.1).AND.(IV.LT.IW)) GO TO 56
      CFIK=2000.D+0*CFI
      CFK=2000.D+0*CF
      CFSK=2000.D+0*CFS
      DTHK=1000.D+0*DTH
      CTH=TWO*TH/(ONE+DSQRT(ONE-TWO*TH*ID/ASEC))
      CH=CDS(J)/CTH
      IEO=REO+HALF
      ITHX=RTHX+HALF
      WRITE (2,83) J,TW,TE,TAW,TP,IEO,ITHX,FRD,CFIK,CFK,CFSK,H,HI,FMY,DT
     1HK,TH,DELTA,DELST
      IF (J.LT.KBL-3) GO TO 54
      IF (J-KBL+2) 49,50,51
49    CTHA=CTH
      XNA=XN
      DLTA=DELTA
      REOA=REO
      GO TO 55
50    CTHB=CTH
      XNB=XN
      DLTB=DELTA
      REOB=REO
      GO TO 55
51    IF (J-KBL) 52,53,54
52    CTHC=CTH
      XNC=XN
      DLTC=DELTA
      REOC=REO
      GO TO 55
53    IF (IT.GT.0) GO TO 55
      DLST=GMA*CDS(J-3)+GMB*CDS(J-2)+GMC*CDS(J-1)+GMD*CDS(J)
      THBL=GMA*CTHA+GMB*CTHB+GMC*CTHC+GMD*CTH
      HBL=DLST/THBL
      DLTBL=GMA*DLTA+GMB*DLTB+GMC*DLTC+GMD*DELTA
      REOBL=GMA*REOA+GMB*REOB+GMC*REOC+GMD*REO
      REOFT=TLV*REOBL
      RETH=THBL*REOBL
      REDL=DLTBL*REOBL
      RETHG=DLOG10(RETH)
      REDLG=DLOG10(REDL)
      XNBL=GMA*XNA+GMB*XNB+GMC*XNC+GMD*XN
      GO TO 55
54    IF ((J.GT.3) .AND. (MOD(J,10).NE.0)) GO TO 56
55    WRITE (2,86) S(J),DSU,THU,CTH,HU,H,CH,XN
56    CONTINUE
      RW(1)=RCO(1)
      CALL SCOND(S,DLA,DAX,KAT)
      DO 57 J=1,KAT
57    DRX(J)=WALTAN(J)+DAX(J)
      IF ((IT.GT.0) .OR. (LD.LT.0)) GO TO 58
      IF ((IV.GT.1) .AND. (IV.LT.IW)) GO TO 58
      IF (KBL.LE.KAT) WRITE (2,85) XBL,DLST,THBL,HBL,XNBL,DLTBL,REOFT,RE
     1TH,RETHG,REDL,REDLG
      IF (KBL.LE.KAT) GO TO 58 
      HBL=CDS(KAT)/CTH
      REOFT=TLV*REO
      RETH=CTH*REO
      REDL=DELTA*REO
      RETHG=DLOG10(RETH)
      REDLG=DLOG10(REDL)
      IF (KBL.GT.KAT) WRITE (2,85) S(KAT),CDS(KAT),CTH,HBL,XN,DELTA,REOF
     1T,RETH,RETHG,REDL,REDLG
58    CONTINUE
      DD2=BMA*DLA(1)+BMB*DLA(2)+BMC*DLA(3)
      DR2=SD(1)+DD2
      DXS=DAX(1)/DR2
      XST=S(1)-DXS
      YST=RCO(1)-HALF*DAX(1)**2/DR2
      SCW=DSQRT(ONE+DAX(1)**2)
      DR2=DR2/SCW**3
      RCV=ONE/DR2/YST
      IF (IT.GT.0) XBIN=SBIN-XST
      IF (IT.GT.0) XCIN=SCIN-XST
      WRITE (2,78) ITLE,XBIN,XCIN,SF
      PPQ=ZRO
      WRITE (2,67) RC,ETAD,AMACH,BMACH,CMACH,EMACH,MC,AH
      IF (TWT.NE.ZRO) GO TO 59
      WRITE (2,81) PPS,T0
      GO TO 5
59    WRITE (2,79) PPS,T0,TWT,TWAT,TAIR
      GO TO 5
60    IF (IT.EQ.0) GO TO 63
      DO 61 K=1,KAT
      S(K)=SREF(K)-XST
61    SCV(K)=DSQRT(ONE+DRX(K)**2)
      SCV(1)=ONE
      SL(1)=ZRO
      IM=(KAT-1)/2
      DO 62 I=1,IM
      J=2*I
      SS=S(J)-S(J-1)
      IF (I.EQ.1) SS=S(2)
      TT=S(J+1)-S(J)
      ST=SS+TT
      S1=(TWO-TT/SS)*ST/SIX
      S3=(TWO-SS/TT)*ST/SIX
      S2=ST-S1-S3
      SA=(TWO+TT/ST)*SS/SIX
      SB=(TWO+ST/TT)*SS/SIX
      SC=SS-SA-SB
      SL(J)=SL(J-1)+SA*SCV(J-1)+SB*SCV(J)+SC*SCV(J+1)
62    SL(J+1)=SL(J-1)+S1*SCV(J-1)+S2*SCV(J)+S3*SCV(J+1)
      XST=ZRO
      WRITE (2,68) LS,DK
      WRITE (2,69) (K,S(K),SL(K),DLA(K),RCO(K),WALTAN(K),SK(K),DAX(K),DR
     1X(K),WMN(K),DMDX(K),SPR(K),BTA(K),K=1,KAT)
      IF (KBL.GT.KAT) GO TO 64
      CALL TWIXT (SL,GMA,GMB,GMC,GMD,SBL,KAT,KBL)
      XBL=GMA*S(KBL-3)+GMB*S(KBL-2)+GMC*S(KBL-1)+GMD*S(KBL)
      DLAB=GMA*DLA(KBL-3)+GMB*DLA(KBL-2)+GMC*DLA(KBL-1)+GMD*DLA(KBL)
      RCOB=GMA*RCO(KBL-3)+GMB*RCO(KBL-2)+GMC*RCO(KBL-1)+GMD*RCO(KBL)
      WRITE (2,89) XBL,SBL,DLAB,RCOB,GMA,GMB,GMC,GMD
      GO TO 64
63    WRITE (2,68) LY,DD
      WRITE (2,69) (K,S(K),FS(K),DLA(K),RCO(K),WALTAN(K),SD(K),DAX(K),DR
     1X(K),WMN(K),DMDX(K),SPR(K),BTA(K),K=1,KAT)
      IF (KBL.GT.KAT) GO TO 64
      CALL TWIXT (S,GMA,GMB,GMC,GMD,XBL,KAT,KBL)
      DLAB=GMA*DLA(KBL-3)+GMB*DLA(KBL-2)+GMC*DLA(KBL-1)+GMD*DLA(KBL)
      RCOB=GMA*RCO(KBL-3)+GMB*RCO(KBL-2)+GMC*RCO(KBL-1)+GMD*RCO(KBL)
      YBL=RCOB-DLAB
      WRITE (2,84) XBL,YBL,DLAB,RCOB,GMA,GMB,GMC,GMD
64    WRITE (2,87) XST,YST,DD2,DR2,RCV
      S(1)=XST
      RCO(1)=YST
      DRX(1)=ZRO
      IF (SBL .EQ. 1.D+3) RETURN
      IF (LV .GT. 0) GO TO 3
65    CONTINUE
      IF (J .EQ. 1) WRITE (2,90) IPQ,QFUNW,TWT
      RETURN
C
66    FORMAT (6E10.3,4I5) 
67    FORMAT (1H ,4H RC=,F11.6,3X,5HETAD=,F8.4,4H DEG,3X,6HAMACH=,F10.7,
     &3X,6HBMACH=,F10.7,3X,6HCMACH=,F10.7,3X,6HEMACH=,F10.7,3X,A4,2HH=,F
     &11.7/)
68    FORMAT (1H ,7X,9HSTA(IN)  ,A4,40H(IN)    DELR(IN)    R(IN)    DY/D
     &X     ,A8,50H     DA/DX     DR/DX    MACH NO.    DM/DX    PE/PO,7X
     &,4HBETA/)
69    FORMAT (10(I4,0P2F11.6,2F11.7,4F10.7,F11.7,F10.7,1P2E12.4/))
70    FORMAT (1H+,5X,34HQUADRATIC TEMPERATURE DISTRIBUTION)
71    FORMAT (1H+,5X,34HPARABOLIC TEMPERATURE DISTRIBUTION)
72    FORMAT (1H+,44X,34HSPALDING-CHI REFERENCE TEMPERATURE)
73    FORMAT (1H+,83X,36HVAN DRIEST REFERENCE REYNOLDS NUMBER /)
74    FORMAT (1H+,83X,35HCOLES LAW REFERENCE REYNOLDS NUMBER /)
75    FORMAT (1H+,44X,34HMODIF. SPALDING-CHI REFERENCE TEMP)
76    FORMAT (1H+,83X,40HREFERENCE REYNOLDS NUMBER BASED ON DELTA /)
77    FORMAT (1H+,44X,29HMODIFIED COLES TRANSFORMATION /)
78    FORMAT (1H1,3A4,39HNOZZLE CONTOUR, RADIAL FLOW ENDS AT STA,F12.7,2
     &5H, TEST CONE BEGINS AT STA,F12.7,16H, SCALE FACTOR =,F13.8/)
79    FORMAT (1H ,1X,15HSTAG. PRESSURE=,F5.0,24H PSI, STAG. TEMPERATURE=
     &,F5.0,21H DEG R, THROAT TEMP.=,F5.0,19H DEG R, WALL TEMP.=,F4.0,24
     &H DEG R, THROAT HT COEF.=,F8.5//)
80    FORMAT (1H1,3A4,49HBOUNDARY LAYER CALCULATIONS, STAGNATION PRESSUR
     &E=,F5.0,28HPSI, STAGNATION TEMPERATURE=,F5.0,27H DEG R, N BASED ON
     & RE,DELTA //)
81    FORMAT (1H ,5X,15HSTAG. PRESSURE=,F5.0,24H PSI STAG. TEMPERATURE=,
     &F5.0,34H DEG R ADIABATIC WALL TEMPERATURE//)
82    FORMAT (1H ,5X,38HTW    TE    TAW    TP    RE/IN    RTHI,4X,3HFRD,
     &5X,4HKCF1,4X,3HKCF,5X,4HRCFS,5X,1HH,6X,2HHI,5X,38HFMY     KTHP THE
     &TA-1  DELTA  DELTA*-1 /)
83    FORMAT (1H ,I3,2F6.1,F7.1,F6.1,I9,I7,4F8.5,F8.4,F7.4,2F8.5,F9.6,F7
     &.4,F9.6)
84    FORMAT (1H ,3HSTA,2F11.6,2F11.7,7X,27HINTERPOLATION COEFFICIENTS,,
     & F12.8,1H,,F11.8,1H,,F11.8,1H,,F12.8 /)
85    FORMAT (1H0,5H   X=,F7.3,11H,   DELTA*=,F10.7,10H,   THETA=,F9.7,6
     &H,   H=,F10.6,6H,   N=,F10.7,10H,   DELTA=,F11.7,10H,   RE/FT=,F11
     &.0//35X,9HRE,THETA=,F9.0,8H,   LOG=,F8.5,1H,,16X,9HRE,DELTA=,F11.0
     &,8H,   LOG=,F8.5)
86    FORMAT (1H ,3X,2HX=,F7.3,8H,   DSU=,F8.5,8H,   THU=,F9.7,8H,   CTH
     &=,F9.7,7H,   HU=,F10.6,6H,   H=,F10.6,7H,   CH=,F10.6,6H,   N=,F8.
     &5)
87    FORMAT (1H ,3HSTA,F11.6,9H      Y*=,F11.7,14H,     D2A/DX2=,F12.9,
     &14H,     D2R/DX2=,F12.9,16H,     VISCID RC=,F14.8)
88    FORMAT (1H ,5HRTHX=,1PE12.5,6H, REO=,E12.5,6H, FRO=,0PF8.5,5H, TH=
     &,F8.5,8H, DELTA=,F8.5,6H, DOT=,F9.5)
89    FORMAT (1H ,3HSTA,2F11.6,2F11.7,7X,27HINTERPOLATION COEFFICIENTS,,
     &F12.8,1H,,F11.8,1H,,F11.8,1H,,F12.8 /)  
90    FORMAT (1H0,10H ITERATION,I4,11H,    QFUN =,F8.5,18H,   THROAT TEM
     &P  =,F6.1 /)
      END
