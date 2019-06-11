      subroutine bound
!
!     to obtain the correction due to the turbulent boundary layer
!
      use kinddefine
      use gg, only:gam,g1,g7,g9,ga,qt
      use coord, only:s,fs,waltan,sd,wmn,ttr,dmdx,spr,bta,sref,xbin,xcin
     &,gma,gmb,gmc,gmd
      use corr, only:dla,rco,dax,drx,sl,dr2
      use prop, only:ar,zo,ro,visc,vism,xbl,conv
      use param, only:etad,rc,amach,bmach,cmach,emach,frc,sf,ah
      use httr, only:hair,taw,twq,tw,twat,qfun,qfunw,ipq,ij,iv,iw
      use contr, only:itle,ie,it,kat,kbl,king,ko,lv,in,mc,mcp  
!
      implicit none
!
      interface
       subroutine scond(a,b,c,king)
        use kinddefine
        implicit none
        integer(kind=K4),intent(in) :: king
        real(kind=K8),dimension(150),intent(in) :: a,b
        real(kind=K8),dimension(150),intent(out) :: c
       end subroutine scond
!
       subroutine twixt(s,gma,gmb,gmc,gmd,xbl,kat,kbl)
        use kinddefine
        implicit none
        integer(kind=K4),intent(in) :: kat
        integer(kind=K4),intent(out) :: kbl
        real(kind=K8),intent(out) :: gma,gmb,gmc,gmd
        real(kind=K8),intent(in) :: xbl
        real(kind=K8),dimension(200),intent(in) :: s
       end subroutine twixt
      end interface
!
      integer(kind=K4) :: i,id,ieo,iht,im,ipp,ir,ithx,j,k,kor,l,ld
      real(kind=K8) :: add,alf,alph,alpha,arc,asca,asec
      real(kind=K8) :: bdd,bet,bma,bmb,bmc
      real(kind=K8) :: c1,c2,c3,capi,cdd,ch,chair,cf
      real(kind=K8) :: cf1,cf2,cf3,cfi,cfik,cfk
      real(kind=K8) :: cfs,cfsk,cmu,cth,ctha,cthb,cthc
      real(kind=K8) :: da,db,dc,dd2,ddd,de,dela,delta,delst,df,dg
      real(kind=K8) :: dh,di,dj,dl,dlab,dlst,dlta,dltb,dltbl,dltc
      real(kind=K8) :: dm,dn,dor,dot,doti,dotr
      real(kind=K8) :: ds,dsm,dsod,dsrod,dst,dstu,dsttu,dsu
      real(kind=K8) :: dt,dth,dtha,dthb,dthc,dthk,dts,dtsu,dtu,dtus,du
      real(kind=K8) :: dut,dxs,emu,fcc,fmy,four,fpcp,frd,fsi,fxcf
      real(kind=K8) :: h,ha,half,hb,hbl,hi,hr,hu,hw
      real(kind=K8) :: one,ppq,pps,qma,qmb,qmc,qmd,raj,rcob,rcv,rdel
      real(kind=K8) :: rdlg,rdli,rdlx,redl,redlg,reo,reoa,reob,reobl
      real(kind=K8) :: reoc,reoft,reth,rethg
      real(kind=K8) :: rho,rhoe,roy,rthg,rthx,rtig
      real(kind=K8) :: rthi,rtii,s1,s2,s3,sa,sb,sbin,sc,scfi,scin,scu
      real(kind=K8) :: scw,six,sma,smb,smc,smd,ss
      real(kind=K8) :: suma,sumb,sumc,sumd,st,str
      real(kind=K8) :: t0,tair,tc,te,ten
      real(kind=K8) :: th,tha,thbl,thm,thr,thu,tlv
      real(kind=K8) :: tp,tr,trpi,tt,twd,two,twt,un
      real(kind=K8) :: ve,wmu,xbl1,xcf,xn,xna,xnb,xnbl,xnc,xst
      real(kind=K8) :: ybl,yoh,yoha,ysec,yst,zro
      real(kind=K8),dimension(16) :: d,z
      real(kind=K8),dimension(200) :: cds,rw,scv,sk
      character(len=4,kind=K3) :: ly,ls
      character(len=8,kind=K3) :: dd,dk
!      COMMON /GG/ GAM,GM,G1,G2,G3,G4,G5,G6,G7,G8,G9,GA,RGA,QT
!      COMMON /CORR/ DLA(200),RCO(200),DAX(200),DRX(200),SL(200),DR2
!      COMMON /COORD/ S(200),FS(200),WALTAN(200),SD(200),WMN(200),TTR(200
!     1),DMDX(200),SPR(200),BTA(200),SREF(200),XBIN,XCIN,GMA,GMB,GMC,GMD
!      COMMON /PROP/ AR,ZO,RO,VISC,VISM,SFOA,XBL,CONV
!      COMMON /PARAM/ ETAD,RC,AMACH,BMACH,CMACH,EMACH,GMACH,FRC,SF,WWO,WW
!     1OP,QM,WE,CBET,XE,ETA,EPSI,BPSI,XO,YO,RRC,SDO,XB,XC,AH,PP,SE,TYE,XA
!      COMMON /HTTR/ HAIR,TAW,TWQ,TW,TWAT,QFUN,QFUNW,IPQ,IJ,IV,IW
!      COMMON /CONTR/ ITLE(3),IE,LR,IT,JB,JQ,JX,KAT,KBL,KING,KO,LV,NOCON,
!     1IN,MC,MCP
!      DIMENSION Z(16), D(16), SCV(200), SK(200), CDS(200), RW(200)
      data zro/0.0d+0/,one/1.d+0/,two/2.d+0/,six/6.d+0/,half/5.d-1/
      data thr/3.d+0/,four/4.d+0/,ten/1.d+1/,tlv/1.2d+1/
      data cf1/3.865d-2/,cf2/4.561d+0/,cf3/5.46d-1/,fsi/3.17897971d+0/
      data ly/'   Y'/,ls/'   S'/,dd/'D2Y/DX2 '/,dk/' CURV.  '/
      data z(1)/.052995325d-1/,z(4)/.1222977958d+0/,z(7)/.3591982246d+0/
      data z(2)/.277124885d-1/,z(5)/.1910618778d+0/,z(8)/.4524937451d+0/
      data z(3)/.671843988d-1/,z(6)/.2709916112d+0/
      data d(1)/.135762297d-1/,d(2)/.31126762d-1/,d(3)/.475792558d-1/
      data d(4)/.623144856d-1/,d(5)/.747979944d-1/,d(6)/.845782597d-1/
      data d(7)/.913017075d-1/,d(8)/.947253052d-1/
      do j=9,16
       d(j)=d(17-j)
       z(j)=one-z(17-j)
      enddo
      do j=1,kat
       sref(j)=s(j)
      enddo
      sbin=xbin
      scin=xcin
      trpi=conv/90.d+0
      fcc=2.05d+0+dlog(.41d+0)
      chair=gam*g1*ar/ro/ro/777.64885d+0
      if (it .eq. 0) xbl1=xbl
!
3     read (1,66,end=65) ppq,t0,twt,twat,qfun,alph,iht,ir,ld,lv
!
      pps=ppq
      rho=144.d+0*pps/zo/ar/t0
      id=iabs(ld)
      kor=ko
      if (iabs(in) .eq. 10) kor=1
      if (mcp .lt. 0) kor=king
      roy=one
      if (ie .eq. 0) hw=ah
      if ((id .eq. 0) .or. (ie .eq. 1)) hw=zro
      if (hw .eq. zro) yoh=zro
      if (hw .eq. zro) yoha=zro 
      alf=dabs(alph)
      arc=frc
      if (iht .lt. 0) arc=frc**(ie+1)
      ipq=0
      iw=1
      if (lv .ne. 0) iw=iabs(lv)
      do j=1,kat
       s(j)=sref(j)
       sl(j)=s(j)
       rw(j)=fs(j)
       rco(j)=fs(j)
       scw=dsqrt(one+waltan(j)**2)
       sk(j)=sd(j)/scw**3
       if (kat .eq. king) goto 4
       if (s(j) .lt. xbl) kbl=j+2
4      drx(j)=waltan(j)
      enddo
      if (kbl .gt. kat) kbl=kat+4
5     do iv=1,iw
       if ((iv .gt. 1) .and. (iv .lt. iw)) goto 15
       if (ld .ge. 0) write (2,80) itle,pps,t0
       if (alph .gt. zro) goto 6
       alpha=zro
       if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,71,advance="no")
       goto 7
6      alpha=alph
       if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,70,advance="no")
7      if (ir .eq. 2) goto 13
       if (alf .eq. one) goto 8
       if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,75,advance="no")
       goto 9
8      if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,72,advance="no")
9      if (ir) 10,11,12
10     if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,74)
       goto 14
11     if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,73)
       goto 14
12     if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,76)
       goto 14
13     if ((ld .ge. 0) .or. (ppq .eq. zro)) write (2,77)
14     if (ppq .eq. zro) goto 60
15     capi=.55d+0
       ipp=0
       ij=1
       do j=1,kat
        bet=ttr(j)-one
        str=one/ttr(j)
        te=t0*str
        raj=wmn(j)*(g7*str)**ga
        if (iht .ge. 0) raj=raj**qt
        scw=dsqrt(one+drx(j)**2)
        emu=visc*te*dsqrt(te)/(te+vism)
        if (te .lt. vism) emu=half*visc*te/dsqrt(vism)
        if (vism .le. one) emu=visc*te**vism
        taw=te*(one+ro*bet)
        rhoe=rho*str**g1
        ve=wmn(j)*dsqrt(gam*ar*te)
        reo=rhoe*ve/emu/tlv
        if (hw .gt. zro) yoh=fs(j)/hw
        if (ie .eq. 0 .and. hw .gt. zro) roy=(hw/fs(j)+one)*trpi
        k=j
        if (j .eq. 1) goto 19
        if (j .gt. kor) k=j-kor+1
        if (k-3) 16,17,18
16      ds=s(j)-s(j-1)
        smd=half*ds
        goto 19
17      dt=s(j)-s(j-1)
        dst=ds+dt
        sma=dst*(two-dt/ds)/six
        smc=dst*(two-ds/dt)/six
        smb=dst-sma-smc
        hb=h
        if (iv .gt. 1) goto 19
        bma=two/ds/dst
        bmb=-two/ds/dt
        bmc=two/dt/dst
        goto 19
18      du=s(j)-s(j-1)
        dt=s(j-1)-s(j-2) 
        ds=s(j-2)-s(j-3)
        dst=ds+dt
        dstu=dst+du
        dtu=dt+du
        dut=du-dt
        dts=ds-dt
        dtus=dt+two*(du-ds)
        dtsu=dt+two*(ds-du)
        dsttu=two*(dst+dtu)
        ha=hb
        hb=h
        qma=half*ds*(one-ds*(thr+(dtu+du)/dst)/dstu/six)
        qmb=half*ds*(one+ds*(two+(dst+dt)/dtu)/dt/six)
        qmc=-ds**3*(one+(dtu+du)/dst)/dt/du/tlv
        qmd=ds**3*(dst+dt)/du/dtu/dstu/tlv
        sma=half*ds+(dut*dtu**3/ds-ds*ds*(ds+dsttu))/dst/dstu/tlv
        smb=half*dst+(ds*ds*(dsttu-ds)/dt+dt*dt*dtus/ds-du**3*(dstu+dst)
     &/ds/dt)/dtu/tlv
        smc=half*dtu+(dt*dt*dtsu/du+du*du*(dsttu-du)/dt-ds**3*(dstu+dtu)
     &/dt/du)/dst/tlv
        smd=half*du+(dts*dst**3/du-du*du*(du+dsttu))/dtu/dstu/tlv
19      if (twt .ne. zro) goto 20
        tw=taw
        goto 21
20      twd=(arc*raj-one)*(twt-twat)/(arc-one)
        if (twd .lt. zro) twd=zro
        tw=twd+twat
21      wmu=visc*tw*dsqrt(tw)/(tw+vism)
        if (vism .le. one) wmu=visc*tw**vism
        dl=tw/te
        dm=alpha*(taw-tw)/te
        dn=one-dl-dm
        da=alf*(taw-tw)
        db=da+tw-te
        if (db) 22,23,24
22      dg=dsqrt(-db*te)
        dh=dsqrt(-db*tw)
        di=(two*(dg+te-tw)-da)/(two*dh+da)
        dj=dlog(di)
        tp=-db/dj/dj
        goto 25
23      tp=(dsqrt(te)+dsqrt(tw))**2/four
        goto 25
24      dc=dsqrt(da*da+four*tw*db)
        df=dasin((db+tw-te)/dc)
        de=dasin(da/dc)
        tp=db/(df+de)/(df+de)
25      if (ir) 26,27,28 
26      frd=tw*emu/wmu/tp
        goto 29
27      frd=emu/wmu
        goto 29
28      frd=te*emu/wmu/dsqrt(tp*tw)
29      if (ipp .gt. 0) goto 31
        rthi=1.d-2*reo*fs(1)
        rtii=rthi
        rdli=ten*rthi
        if (ir .eq. 1) goto 32
30      rthg=dlog10(rthi)
        cfi=cf1/(rthg+cf2)/(rthg-cf3)
31      if (ir .ne. 2) goto 33
        scfi=dsqrt(cfi)
        tc=tw+17.2d+0*scfi*da-305.d+0*cfi*db
        cmu=visc*tc*dsqrt(tc)/(tc+vism)
        if (vism .le. one) cmu=visc*tc**vism
        tp=tw*cmu/wmu
        frd=emu/cmu
        goto 33
32      rdlg=dlog10(rdli)
        cfi=0.0444d+0/(rdlg+4.6221d+0)/(rdlg-1.4402d+0)
33      cf=cfi*te/tp
        cfs=cf*scw
        rtig=dlog10(rtii)
        xcf=.41d+0*dsqrt((rtig+cf2)*(rtig-cf3)/cf1)
34      c3=two+capi*(fsi+1.5d+0*capi)
        c2=one+capi
        c1=c2-c3/xcf
        fxcf=xcf+dlog(c1/rtii)-fcc-two*capi
        fpcp=(xcf-fsi-thr*capi)/xcf/c1-two
        capi=capi-fxcf/fpcp
        if (dabs(fxcf) .gt. 1.d-8) goto 34
        doti=xcf/c1
        xn=half*(doti+dsqrt(doti*(doti-six)+one)-thr)
        hi=one+two/xn
        suma=zro
        sumb=zro
        sumc=zro
        sumd=zro
        do l=1,16
         un=z(l)**xn
         tr=dl+z(l)*(dm+z(l)*dn)
         add=d(l)*xn*un/tr
         bdd=add*z(l)
         cdd=add*un
         ddd=bdd*un
         suma=suma+add
         sumb=sumb+bdd
         sumc=sumc+cdd
         sumd=sumd+ddd
        enddo
        dot=one/(suma-sumb)
        dsod=one-suma
        dsm=half-sumc
        thm=sumc-sumd
        hu=dsod*dot
        if (ipp .gt. 0) goto 36
        h=hu
        dotr=dot
36      fmy=(h+two-g9*bet)*dmdx(j)*str/wmn(j)+id*drx(j)/(rw(j)+hw)
        if (j.eq.1) th=cfs/fmy
        if (k.eq.2) th=(tha+smd*(dtha+cfs))/(one+smd*fmy)
        if (k.eq.3) th=(tha+sma*dtha+smb*dthb+smc*cfs)/(one+smc*fmy)
        if (k.gt.3) th=(tha+sma*dtha+smb*dthb+smc*dthc+smd*cfs)/(one+smd
     &*fmy)
        delst=h*th
        asec=delst+dsqrt(id*delst**2+(fs(j)*scw*roy)**2)
        dor=id*dotr*th/asec
        dsrod=dsod-dor*dsm
        ipp=1
        dotr=one/(one/dot-thm*dor)
        hr=dsrod*dotr
        if (dabs(h-hr) .lt. 5.d-7) goto 37
        h=hr
        goto 36
37      delta=dotr*th
        thu=delta/dot
        dsu=delta*dsod
        rdel=reo*delta
        rtii=rdel/doti
        rdlx=frd*rdel
        rthx=rdlx/dot
        if (rthx .lt. 100.d+0) goto 38
        if (ir .eq. 1) goto 39
        if (dabs(one-rthx/rthi) .lt. 1.d-6) goto 41
        rthi=rthx
        goto 30
38      write (2,88) rthx,reo,frd,th,delta,dot
        return
39      if (dabs(one-rdlx/rdli) .lt. 1.d-6) goto 40
        rdli=rdlx
        goto 32
40      rthg=half*(dsqrt((cf2+cf3)**2+four*cf1/cf1)-cf2+cf3)
        rthx=ten**rthg
41      if (j .gt. 1) goto 42
        dth=zro
        hair=rhoe*ve*cf*chair
        tair=hair
        if (twat.eq.twt .or. qfun.eq.zro) goto 46
        twq=(hair*taw+qfun*(twat-15.d+0))/(hair+qfun)
        call heat
        if (ipq .gt. 100) goto 65
        if (dabs(tw-twq).lt.1.d-2.and.dabs(qfun-qfunw).lt.1.d-5) goto 46
        twt=twat+(twq-twat)*(arc-one)/(arc*raj-one)
        qfun=qfunw
        goto 20
42      dth=cfs-th*fmy
        if (dth .lt. zro) dth=zro
        if (j .eq. kor) goto 46
        if (k-3) 43,45,44
43      dthb=dth
        goto 47
44      tha=tha+qma*dtha+qmb*dthb+qmc*dthc+qmd*dth
        dtha=dthb
        dthb=dthc
        if (k .gt. 5) goto 45
        scu=dsqrt(one+drx(j-2)**2)
        dela=ha*tha
        if ((ie .eq. 1).or.(id .eq. 0)) ysec=fs(j-2)*scu
        if ((ie .eq. 0).and.(hw .gt. zro)) ysec=scu*(fs(j-2)+hw)*trpi
        if (hw .gt. zro) yoha=fs(j-2)/hw
        asca=dela+dsqrt(id*dela**2+ysec**2)
        rw(j-2)=asca/scu
        dla(j-2)=scu*(asca-ysec)*(one+yoha)
        rco(j-2)=fs(j-2)+dla(j-2)
45      dthc=dth
        goto 47
46      tha=th
        dtha=dth
        if ((iv .gt. 1) .and. (iv .lt. iw)) goto 47
        if (j .eq. 1 .and. ld .ge. 0) write (2,82)
47      cds(j)=asec-scw*fs(j)*roy
        dla(j)=scw*cds(j)*(one+yoh)
        rco(j)=fs(j)+dla(j)
        rw(j)=asec/scw
        if (iv .lt. iw) goto 48
        bta(j)=-dmdx(j)*dsu/wmn(j)/ttr(j)/scw/cfi
        if (j.eq.1 .or. j.gt.ko .or. iht.eq.0) goto 48
        if (mod(j,iht) .ne. 1) goto 48
        ij=j
        hair=rhoe*ve*cf*chair
        call heat
48      if (ld.lt.0) goto 56
        if ((iv.gt.1).and.(iv.lt.iw)) goto 56
        cfik=2000.d+0*cfi
        cfk=2000.d+0*cf
        cfsk=2000.d+0*cfs
        dthk=1000.d+0*dth
        cth=two*th/(one+dsqrt(one-two*th*id/asec))
        ch=cds(j)/cth
        ieo=int(reo+half)
        ithx=int(rthx+half)
        write (2,83) j,tw,te,taw,tp,ieo,ithx,frd,cfik,cfk,cfsk,h,hi,fmy,
     &dthk,th,delta,delst
        if (j.lt.kbl-3) goto 54
        if (j-kbl+2) 49,50,51
49      ctha=cth
        xna=xn
        dlta=delta
        reoa=reo
        goto 55
50      cthb=cth
        xnb=xn
        dltb=delta
        reob=reo
        goto 55
51      if (j-kbl) 52,53,54
52      cthc=cth
        xnc=xn
        dltc=delta
        reoc=reo
        goto 55
53      if (it.gt.0) goto 55
        dlst=gma*cds(j-3)+gmb*cds(j-2)+gmc*cds(j-1)+gmd*cds(j)
        thbl=gma*ctha+gmb*cthb+gmc*cthc+gmd*cth
        hbl=dlst/thbl
        dltbl=gma*dlta+gmb*dltb+gmc*dltc+gmd*delta
        reobl=gma*reoa+gmb*reob+gmc*reoc+gmd*reo
        reoft=tlv*reobl
        reth=thbl*reobl
        redl=dltbl*reobl
        rethg=dlog10(reth)
        redlg=dlog10(redl)
        xnbl=gma*xna+gmb*xnb+gmc*xnc+gmd*xn
        goto 55
54      if ((j.gt.3) .and. (mod(j,10).ne.0)) goto 56
55      write (2,86) s(j),dsu,thu,cth,hu,h,ch,xn
56      continue
       enddo
       rw(1)=rco(1)
       call scond(s,dla,dax,kat)
       do j=1,kat
        drx(j)=waltan(j)+dax(j)
       enddo
       if ((it.gt.0) .or. (ld.lt.0)) goto 58
       if ((iv.gt.1) .and. (iv.lt.iw)) goto 58
       if (kbl .le. kat) write(2,85) xbl1,dlst,thbl,hbl,xnbl,dltbl,reoft
     &,reth,rethg,redl,redlg
       if (kbl.le.kat) goto 58 
       hbl=cds(kat)/cth
       reoft=tlv*reo
       reth=cth*reo
       redl=delta*reo
       rethg=dlog10(reth)
       redlg=dlog10(redl)
       if(kbl.gt.kat) write(2,85) s(kat),cds(kat),cth,hbl,xn,delta,reoft
     &,reth,rethg,redl,redlg
58     continue
      enddo
      dd2=bma*dla(1)+bmb*dla(2)+bmc*dla(3)
      dr2=sd(1)+dd2
      dxs=dax(1)/dr2
      xst=s(1)-dxs
      yst=rco(1)-half*dax(1)**2/dr2
      scw=dsqrt(one+dax(1)**2)
      dr2=dr2/scw**3
      rcv=one/dr2/yst
      if (it.gt.0) xbin=sbin-xst
      if (it.gt.0) xcin=scin-xst
      write (2,78) itle,xbin,xcin,sf
      ppq=zro
      write (2,67) rc,etad,amach,bmach,cmach,emach,mc,ah
      if (twt.ne.zro) goto 59
      write (2,81) pps,t0
      goto 5
59    write (2,79) pps,t0,twt,twat,tair
      goto 5
60    if (it.eq.0) goto 63
      do k=1,kat
       s(k)=sref(k)-xst
       scv(k)=dsqrt(one+drx(k)**2)
      enddo
      scv(1)=one
      sl(1)=zro
      im=(kat-1)/2
      do i=1,im
       j=2*i
       ss=s(j)-s(j-1)
       if (i.eq.1) ss=s(2)
       tt=s(j+1)-s(j)
       st=ss+tt
       s1=(two-tt/ss)*st/six
       s3=(two-ss/tt)*st/six
       s2=st-s1-s3
       sa=(two+tt/st)*ss/six
       sb=(two+st/tt)*ss/six
       sc=ss-sa-sb
       sl(j)=sl(j-1)+sa*scv(j-1)+sb*scv(j)+sc*scv(j+1)
       sl(j+1)=sl(j-1)+s1*scv(j-1)+s2*scv(j)+s3*scv(j+1)
      enddo
      xst=zro
      write (2,68) ls,dk
      write (2,69) (k,s(k),sl(k),dla(k),rco(k),waltan(k),sk(k),dax(k),dr
     &x(k),wmn(k),dmdx(k),spr(k),bta(k),k=1,kat)
      if (kbl.gt.kat) goto 64
      call twixt (sl,gma,gmb,gmc,gmd,xbl,kat,kbl)
      xbl1=gma*s(kbl-3)+gmb*s(kbl-2)+gmc*s(kbl-1)+gmd*s(kbl)
      dlab=gma*dla(kbl-3)+gmb*dla(kbl-2)+gmc*dla(kbl-1)+gmd*dla(kbl)
      rcob=gma*rco(kbl-3)+gmb*rco(kbl-2)+gmc*rco(kbl-1)+gmd*rco(kbl)
      write (2,89) xbl1,xbl,dlab,rcob,gma,gmb,gmc,gmd
      goto 64
63    write (2,68) ly,dd
      write (2,69) (k,s(k),fs(k),dla(k),rco(k),waltan(k),sd(k),dax(k),dr
     &x(k),wmn(k),dmdx(k),spr(k),bta(k),k=1,kat)
      if (kbl.gt.kat) goto 64
      call twixt (s,gma,gmb,gmc,gmd,xbl1,kat,kbl)
      dlab=gma*dla(kbl-3)+gmb*dla(kbl-2)+gmc*dla(kbl-1)+gmd*dla(kbl)
      rcob=gma*rco(kbl-3)+gmb*rco(kbl-2)+gmc*rco(kbl-1)+gmd*rco(kbl)
      ybl=rcob-dlab
      write (2,84) xbl1,ybl,dlab,rcob,gma,gmb,gmc,gmd
64    write (2,87) xst,yst,dd2,dr2,rcv
      s(1)=xst
      rco(1)=yst
      drx(1)=zro
      if (xbl .eq. 1.d+3) return
      if (lv .gt. 0) goto 3
65    continue
      if (j .eq. 1) write (2,90) ipq,qfunw,twt
      return
!
66    format (6e10.3,4i5) 
67    format (1x,' RC=',f11.6,3x,'ETAD=',f8.4,' DEG',3x,'AMACH=',f10.7,3
     &x,'BMACH=',f10.7,3x,'CMACH=',f10.7,3x,'EMACH=',f10.7,3x,a4,'H=',f1
     &1.7/)
68    format (1x,7x,'STA(IN)  ',a4,'(IN)    DELR(IN)    R(IN)    DY/DX  
     &   ',a8,'     DA/DX     DR/DX    MACH NO.    DM/DX    PE/PO',7x,'B
     &ETA'/)
69    format (10(i4,0p2f11.6,2f11.7,4f10.7,f11.7,f10.7,1p2e12.4/))
70    format (1x,5x,'QUADRATIC TEMPERATURE DISTRIBUTION')
71    format (1x,5x,'PARABOLIC TEMPERATURE DISTRIBUTION')
72    format (1x,5x,'SPALDING-CHI REFERENCE TEMPERATURE')
73    format (1x,5x,'VAN DRIEST REFERENCE REYNOLDS NUMBER'/)
74    format (1x,5x,'COLES LAW REFERENCE REYNOLDS NUMBER'/)
75    format (1x,5x,'MODIF. SPALDING-CHI REFERENCE TEMP')
76    format ('+',83x,'REFERENCE REYNOLDS NUMBER BASED ON DELTA'/)
77    format ('+',44x,'MODIFIED COLES TRANSFORMATION'/)
78    format (1x,3a4,'NOZZLE CONTOUR, RADIAL FLOW ENDS AT STA',f12.7,', 
     &TEST CONE BEGINS AT STA',f12.7,', SCALE FACTOR =',f13.8/)
79    format (1x,1x,'STAG. PRESSURE=',f5.0,' PSI, STAG. TEMPERATURE=',f5
     &.0,' DEG R, THROAT TEMP.=',f5.0,' DEG R, WALL TEMP.=',f4.0,' DEG R
     &, THROAT HT COEF.=',f8.5//)
80    format (1x,3a4,'BOUNDARY LAYER CALCULATIONS, STAGNATION PRESSURE='
     &,f5.0,'PSI, STAGNATION TEMPERATURE=',f5.0,' DEG R, N BASED ON RE,D
     &ELTA'//)
81    format (1x,5x,'STAG. PRESSURE=',f5.0,' PSI STAG. TEMPERATURE=',f5.
     &0,' DEG R ADIABATIC WALL TEMPERATURE'//)
82    format (1x,5x,'TW    TE    TAW    TP    RE/IN    RTHI',4x,'FRD',5x
     &,'KCF1',4x,'KCF',5x,'RCFS',5x,'H',6x,'HI',5x,'FMY     KTHP THETA-1
     &  DELTA  DELTA*-1'/)
83    format (1x,i3,2f6.1,f7.1,f6.1,i9,i7,4f8.5,f8.4,f7.4,2f8.5,f9.6,f7.
     &4,f9.6)
84    format (1x,'STA',2f11.6,2f11.7,7x,'INTERPOLATION COEFFICIENTS,',f1
     &2.8,',',f11.8,',',f11.8,',',f12.8/)
85    format (1x,'   X=',f7.3,',   DELTA*=',f10.7,',   THETA=',f9.7,',  
     & H=',f10.6,',   N=',f10.7,',   DELTA=',f11.7,',   RE/FT=',f11.0//3
     &5x,'RE,THETA=',f9.0,',   LOG=',f8.5,',',16x,'RE,DELTA=',f11.0,',  
     & LOG=',f8.5/)
86    format (1x,3x,'X=',f7.3,',   DSU=',f8.5,',   THU=',f9.7,',   CTH='
     &,f9.7,',   HU=',f10.6,',   H=',f10.6,',   CH=',f10.6,',   N=',f8.5
     &)
87    format (1x,'STA',f11.6,'      Y*=',f11.7,',     D2A/DX2=',f12.9,',
     &     D2R/DX2=',f12.9,',     VISCID RC=',f14.8/)
88    format (1x,'RTHX=',1pe12.5,', REO=',e12.5,', FRO=',0pf8.5,', TH=',
     &f8.5,', DELTA=',f8.5,', DOT=',f9.5)
89    format (1x,'STA',2f11.6,2f11.7,7x,'INTERPOLATION COEFFICIENTS,',f1
     &2.8,',',f11.8,',',f11.8,',',f12.8/)  
90    format ('0',' ITERATION',i4,',    QFUN =',f8.5,',   THROAT TEMP  =
     &',f6.1/)
      end subroutine bound
