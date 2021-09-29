      subroutine axial
!
!     to obtain the axial distributiun of velocity and/or mach number
!
      use kinddefine
      use fg, only:gc,gd,ge,gf,gh,gi,hb,hc,he
      use gg, only:gam,gm,g2,g4,g5,g6,g7,g8,g9,ga,rga,qt
      use cline, only:wip,x1,frip,zonk,seo,cse,axis,taxi
      use prop, only:sfoa,conv
      use param, only: etad,rc,amach,bmach,cmach,emach,gmach,frc,sf,wwo,
     &wwop,qm,we,cbet,xe,eta,epsi,bpsi,xo,yo,rrc,sdo,xb,xc,ah,pp,se,tye,
     &xa
      use contr,only:itle,ie,lr,it,jb,jq,jx,lv,nocon,in,mc,mcp,ip,iq,ise
     &,jc,m,mp,mq,n,np,nf,nr,lc,md,mf,mt,nd,nt
      implicit none
!
      interface
       function cubic(ea,eb,ec,ed)
        use kinddefine
        implicit none
        real(kind=K8) :: cubic
        real(kind=K8),intent(in) :: ea,eb,ec,ed
       end function cubic
!
       function fmv(psi)
        use kinddefine
        implicit none
        real(kind=K8) :: fmv
        real(kind=K8), intent(in) :: psi
       end function fmv
!
       function toric(wip,se)
        use kinddefine
        implicit none
        real(kind=K8) :: toric
        real(kind=K8),intent(in) :: se,wip
       end function toric
!
       subroutine conic(xm,b) 
        use kinddefine
        implicit none
        real(kind=K8),dimension(4),intent(out) :: b
        real(kind=K8),intent(in) :: xm
       end subroutine conic
!
       subroutine scond(a,b,c,king)
        use kinddefine
        implicit none
        integer(kind=K4),intent(in) :: king
        real(kind=K8),dimension(150),intent(in) :: a,b
        real(kind=K8),dimension(150),intent(out) :: c
       end subroutine scond
!
       subroutine sorce(w,b)
        use kinddefine
        implicit none
        real(kind=K8),intent(in) :: w
        real(kind=K8),dimension(4),intent(out) :: b
       end subroutine sorce
!
       subroutine trans(rto,tk,wo,amn,amp,ampp,w,awp,awpp,cwoppp,axn)
        use kinddefine
        implicit none
        real(kind=K8) :: rto,tk,wo,amn,amp,ampp,w,awp,awpp,cwoppp,axn
       end subroutine trans
      end interface
!
      integer(kind=K4) :: ix,j,k,l,n0,n3,n4,n5,nx
      real(kind=K8) :: aa,ab,abcm,aem,am,amp,ampp,amsq
      real(kind=K8) :: apsi,awp,awpp,awppp,bbet,bmp,bmpp,bppp
      real(kind=K8) :: c2,c3,c4,cm,cbm,cmc,cmp,cpp,cppp,dw,dx
      real(kind=K8) :: ea,eb,ebet,ec,ed,eit,eoe,ew
      real(kind=K8) :: fbet,fftn,fiv,fmach,fn,four,fpsi,fq,fxw
      real(kind=K8) :: gj,gk,gmm,gmp,gpsi,gq,gr,gs,gw,gwp,gww
      real(kind=K8) :: h,half,hh,om,one,q
      real(kind=K8) :: ra,rg,rmach,rt,sev,six,sm,smpp,smppp,sxty
      real(kind=K8) :: ten,thr,tk,tlv,trty,two,tyin
      real(kind=K8) :: w,wap,wapp,wb,xbc,wbp,wbpp,wc,wcb
      real(kind=K8) :: wcp,wep,wepp,weppp,whp,whpp
      real(kind=K8) :: wi,wipp,wippp,wm,wo,woppp,wp,wpp,wppp
      real(kind=K8) :: wrppp,wspp,wsppp
      real(kind=K8) :: x1in,xain,xbcm,xbcmn,xbcmx,xbcn,xbet,xbin
      real(kind=K8) :: xcin,xd,xdin,xi,xie
      real(kind=K8) :: xiin,xinch,xj,xm,xmp,xmpp,xmppp,xmw
      real(kind=K8) :: xo1,xoi,xoin,yain,zro
      real(kind=K8),dimension(4) :: d
      real(kind=K8),dimension(6) :: c
      real(kind=K8),dimension(150) :: ax,axm,axmp
      character(len=4,kind=K3) :: iaxis,m1,m2
      data zro/0.0d+0/,one/1.d+0/,two/2.d+0/,six/6.d+0/,half/5.d-1/
      data thr/3.d+0/,four/4.d+0/,fiv/5.d+0/,ten/1.d+1/,tlv/1.2d+1/
      data sev/7.d+0/,eit/8.d+0/,fftn/1.5d+1/,trty/3.d+1/,sxty/6.d+1/
      data m1/'GMAC'/,m2/'2-D '/,iaxis/'AXIS'/
      data n3/4h 3RD/,n4/4h 4TH/,n5/4h 5TH/,n0/4h-DEG/
!      npi=9.d+1/conv
!      if (jq.eq.0.and.jx.eq.0) call orez(axis,2*750)
      if (jq.eq.0.and.jx.eq.0) then
       axis(:,:)=0.0d0
       taxi(:,:)=0.0d0
      endif
      if (jq .gt. 0) goto 50
      if (jx .eq. 0) goto 2
!
!     card used to obtain internal streamlines (jx > 0)
!
      read (1,93,end=91) etad,qm,xj
!
      jx=int(xj,K4)
      if (etad .eq. sxty) goto 1
      eta=etad/conv
      if (ie .eq. 0) se=eta
      if (ie .eq. 1) se=two*dsin(half*eta)
      cse=dcos(eta)
      apsi=bpsi-eta/qt
      amach=fmv(apsi)
      ra=((g6+g5*amach**2)**ga/amach)**qt
      gpsi=epsi+eta/qt
      gmach=fmv(gpsi)
      rg=((g6+g5*gmach**2)**ga/gmach)**qt
      mp=one+thr*(ra-rg)
      goto 14
1     se=qm*seo
      goto 14
!
!     constants used in transonic solution
2     gc=(two*gam/qt-thr)/six/(3+ie)
      ge=(thr*(8+ie)-four*gam/qt)/thr/(7+ie)
      gh=(fftn+(2-6*ie)*gam)/tlv/(5+ie)
      gj=(gam*(gam+9.25d+0*ie-26.5d+0)+.75d+0*(6-ie))/tlv/(3-ie)
      gk=(gam*(gam+2.25d+0*ie-16.5d+0)+2.25d+0*(2+ie))/six
      gr=(fftn-(1+9*ie)*gam)/(15+ie)/18.d+0
      hb=(14.d+0*gam-75.d+0+18*ie)/(270.d+0+18*ie)
      if (ie .eq. 0) goto 3
      gd=(gm*(652.d+0*gm+1319.d+0)+1000.d+0)/6912.d+0
      gf=(3612.d+0+gm*(751.d+0+gm*754.d+0))/2880.d+0
      gi=(909.d+0+gam*(270.d+0+gam*412.d+0))/10368.d+0
      gs=(gam*(gam*2708.d+0+2079.d+0)+2115.d+0)/82944.d+0
      hc=(gam*(2364.d+0*gam-3915.d+0)+14337.d+0)/82944.d+0
      he=(gam*(64.d+0*gam+117.d+0)-1026.d+0)/1152.d+0
      goto 4
!
!     axisym flow, ie=1, qt=0.5, gam=1.4, gc=0.10833333, gd=0.236099537,
!     ge=0.65833333, gf=1.40036111, gh=0.13055556, gi=0.2020177469,
!     gj=0.76833333, gk=-1.87333333, gr=0.003472222, gs=0.1245814043,
!     hb=0.12986111, hc=0.1626331019, he=-0.6395486111
!
3     gd=(gm*(32.d+0*gm-14.d+0)+221.d+0)/1080.d+0
      gf=(4230.d+0+gm*(211.d+0+gm*334.d+0))/3780.d+0
      gi=(738.d+0+gam*(273.d+0-gam*82.d+0))/7560.d+0
      gs=(gam*(gam*782.d+0+3507.d+0)+7767.d+0)/272160.d+0
      hc=(gam*(274.d+0*gam-861.d+0)+4464.d+0)/17010.d+0
      he=(gam*(32.d+0*gam+87.d+0)-561.d+0)/540.d+0
!
!     planar flow, ie=0, qt=1.0,gam=1.4, gc=-0.011111, gd=0.2041851852,
!     ge=0.8761904762, gf=1.155513228, gh=0.29666667, gi=0.1269153439,
!     gj=-0.85111111, gk=-2.7733333, gr=0.05037037037, gs=0.05221017049,
!     hb=-0.2051851852, hc=0.2231416814, he=-0.6971851852
!
!     card used to establish inviscid parameters
!
4     read (1,93,end=91) etad,rc,fmach,bmach,cmc,sf,pp,xc
!
!     card used to control calculations
!
      read (1,92) mt,nt,ix,in,iq,md,nd,nf,mp,mq,jb,jx,jc,it,lr,nx
!
      lc=int(xc,K4)
      if (xc .gt. one) lc=int(xc+one,K4)
      nr=six*rc
      mf=fmach
      if (ie .eq. 1) mc=m1
      if (ie .eq. 0) mc=m2
      nocon=0
      eta=etad/conv
      if (ie .eq. 0) se=eta
      if (ie .eq. 1) se=two*dsin(half*eta)
      if (etad .eq. sxty) se=one
      seo=se
      ise=int(se,K4)
      cse=dcos(eta)
      rt=rc+one
      am=one
      wi=one
      wipp=zro
      mcp=cmc
      cmach=dabs(cmc)
      cbet=dsqrt(cmach*cmach-one)
      frc=((g6+g5*cmach**2)**ga/cmach)**qt
      tye=frc*se
      if (sf .lt. zro) sf=-sf/tye
      if (ise .eq. 0) goto 5
!
!     non-radial flow at inflection point
      iq=1
      amach=cmach
      bmach=cmach
      emach=cmach
      fmach=cmach
      gmach=cmach
      if (ie .eq. 1) am=gmach
      we=g2*emach/dsqrt(emach**2+g9)
      dw=we-wi
      xo=zro
      eoe=zro
      goto 15
!
!     radial flow at inflection point
5     if (in .eq. 0) goto 6
      if ((lc .lt. 0) .and. (in .lt. 0)) in=-1
      if ((lc .eq. 0) .or. (mcp .lt. 0)) in=isign(10,in)
6     bbet=dsqrt(bmach*bmach-one)
      bpsi=g2*datan(g4*bbet)-datan(bbet)
      if (fmach) 9,8,7
7     fbet=dsqrt(fmach*fmach-one)
      fpsi=g2*datan(g4*fbet)-datan(fbet)
      goto 10
8     fmach=-bpsi/eta
      if (bpsi/eta .gt. 7.5d+0) fmach=-7.5d+0
9     fpsi=-fmach*eta
      fmach=fmv(fpsi)
10    epsi=fpsi-two*eta/qt
      emach=fmv(epsi)
      we=g2*emach/dsqrt(emach*emach+g9)
      dw=we-wi
      call sorce(we,d)
      xe=d(1)
      wep=d(2)
      wepp=d(3)
      wrppp=d(4)
      if (nr .ne. 0) goto 15
      if ((lr .ne. 0) .or. (iq .lt. 0)) goto 11
      if (ix .eq. 0) write (2,106) itle,n3
      if (ix .ne. 0) write (2,106) itle,n4
!
!     iteration to determine rc if not specified (nr = 0)
11    ea=wrppp
      eb=-fiv*wepp-wipp
      ec=tlv*wep
      ed=-tlv*dw
      xie=cubic(ea,eb,ec,ed)
      if (xie .le. zro) goto 89
12    wip=two*(we-one)/xie-wep+(wepp-wipp)*xie/six
13    nocon=nocon+1
      if (nocon .gt. 100) goto 90
14    rt=toric(wip,se)
      rc=rt-one
15    tk=(one-g7*(one+(ge+gf/rt)/rt)/rt**2/(15+ie)/thr)**qt
      yo=se/tk
      aa=dsqrt(qt*(gam+one)*rt)
      if (qm .ne. one) goto 19
      whpp=(one-gam/1.5d+0+gj/rt)/(aa*yo)**2
      if ((nr .ne. 0) .or. (ise.eq. 1)) goto 18
      if (dabs(whpp-wipp) .lt. 1.d-10) goto 18
      wipp=whpp
      if (ix) 11,17,16
16    ea=gk/(aa*yo)**3
      eb=thr*(wipp+wepp)
      ec=-tlv*wep
      ed=tlv*dw
      xie=cubic(ea,eb,ec,ed)
      if (xie .le. zro) goto 89
      goto 12
17    h=(eit*wip+sev*wep)/(thr*wipp-two*wepp)
      hh=trty*dw/(thr*wipp-two*wepp)
      xie=hh/(dsqrt(h*h+hh)+h)
      wip=wep-half*xie*(wepp+wipp)
      goto 13
!
!     iteration for rc completed, remainder of transonic values computed
18    wip=(one-(gc-gd/rt)/rt)/yo/aa
      whp=wip
      wipp=whpp
      amp=g7*wip
      ampp=g7*(whpp+thr*g8*wip**2)
19    xoi=yo*dsqrt(g7/two/(9-ie)/rt)*(one+(gh+gi/rt)/rt)
      if (qm .ne. one) goto 21
      if (ise .eq. 1) xi=xoi
      xo1=xoi
      wo=one-(half/(3-ie)+(gr+gs/rt)/rt)/rt
      om=wo/dsqrt(g7-g8*wo**2)
      woppp=gk/(aa*yo)**3
      if (lr .eq. 0) goto 21
!
!     call for throat characteristic values
      call trans (rt,tk,wo,am,amp,ampp,wi,awp,awpp,awppp,xi)
      if ((nx .lt. 0) .and. (nt.lt.0)) goto 87
      if (nx .lt. 0) goto 4
      amp=amp/se
      ampp=ampp/se**2
      wap=awp/se
      wapp=awpp/se**2
      woppp=awppp/se**3
      if (ise .eq. 1) goto 21
      dw=we-wi
      xoi=xi*se
      if (nr .gt. 0) goto 20
      x1=xe-xie
      xo=xe-xie-xo1
      c2=xie*wip
      c3=half*wipp*xie**2
      c4=we-one-c2-c3
      if (ix .ne. 0) c4=four*c4+two*c3+c2-xie*wep
      if (iq .lt. 0) goto 20
      write (2,110) itle,n4,lr
      write (2,96) xie,c2,c3,c4,x1
20    wip=wap
      wipp=wapp
21    wwo=one+(one/(ie+3)-(hb-hc/rt)/rt)/rt
      wwop=(one+(one-ie/eit-he/rt)/rt)/yo/aa
      rrc=one/rc
      sdo=rrc/yo
      zonk=qm+1.0d-03
      np=zonk*(iabs(nf)-1)+1
      if (sf .gt. zro) goto 22
      sf=one/yo
22    if (iq .lt. 0) goto 31
      ip=0
      jq=0
      m=zonk*(mt-1)+1
      n=nt
      if (qm .eq. one) goto 23
      xo=x1-xoi
      return
!23    call orez (c,6)
23    c(:)=0.0d0
      if (ise .eq. 0) goto 31
!
!     length of axial distribution for non-radial flow
      x1=xoi
      aem=emach-am
      c(1)=am
      if (lc) 25,24,27
24    amsq=amp**2+aem*ampp*four/thr
      if (lr .eq. 0) write (2,122) itle,n4,n0
      if (lr .ne. 0) write (2,107) itle,n4,n0,lr
      if (amsq .lt. zro) goto 28
      xie=four*aem/(dsqrt(amsq)+amp)
      xe=xie+xi
      c(5)=thr*aem-amp*xie
      goto 26
25    xie=thr*aem/amp
      xe=xie+xi
      if (lr .eq. 0) write (2,122) itle,n3,n0
      if (lr .ne. 0) write (2,107) itle,n3,n0,lr
26    c(2)=amp*xie
      c(3)=six*aem-thr*c(2)
      c(4)=thr*c(2)-eit*aem
      goto 46
27    if (lc .eq. 1) goto 29
      xe=xc/tk
      xie=fiv*aem/(dsqrt(amp**2+in*aem*ampp/eit)+amp)
      if (xe .gt. xi+xie) xe=xi+xie
      xie=xe-xi
      c(2)=amp*xie
      c(3)=half*in*ampp*xie**2/ten
      c(4)=ten*aem-six*c(2)-thr*c(3)
      c(5)=-fftn*aem+eit*c(2)+thr*c(3)
      c(6)=six*aem-thr*c(2)-c(3)
      if (lr .eq. 0) write (2,122) itle,n5,n0
      if (lr .ne. 0) write (2,107) itle,n5,n0,lr
      goto 46
28    c(2)=two*aem
      c(4)=-c(2)
      c(5)=aem
      xie=two*aem/amp
      xe=xie+xi
      goto 46
!     if xc=1 then read the centerline mach distribution from point b
!     to point c. see explanatory note on xc page 54 in adec-tr-87-63
29    do j=1,nt
       k=nt+1-j
       read(9) ax(k),axm(k),axmp(k)
       if (j .eq. 1) dx=xi-ax(k)
       axis(1,k)=ax(k)+dx
      enddo
      axm(nt)=am
      axmp(nt)=amp
      xe=axis(1,1)
      xie=xe-xi
      if (lr .eq. 0) write (2,122) itle,n5,n0
      if (lr .ne. 0) write (2,107) itle,n5,n0,lr
      goto 46
!
!     length of upstream axial distribution for radial flow
31    if (sfoa .eq. zro) goto 32
      if (lr .eq. 0) write (2,106) itle,n5
      if (lr .ne. 0) write (2,110) itle,n5,lr
      goto 44
32    if (lr .eq. 0) goto 33
      if ((nr .eq. 0) .and. (ix .eq. 0)) goto 41
      if ((nr .eq. 0) .and. (ix .ne. 0)) mf=0
      if (mf .ne. 0) goto 40
      if ((iq .lt. 0) .or. (nr .eq. 0)) goto 35
      if (ix .eq. 0) write (2,110) itle,n3,lr
      if (ix .ne. 0) write (2,110) itle,n4,lr
      goto 35
33    if (mf .eq. 0) goto 34
      if (nr .eq. 0) goto 45
      if (iq .ge. 0) write (2,106) itle,n4
      goto 41
!
!     iteration for emach if not specified (mf = 0)
34    if (iq .lt. 0) goto 35
      if (ix .eq. 0) write (2,106) itle,n3
      if (ix .ne. 0) write (2,106) itle,n4
35    if (nocon .gt. 100) goto 90
      if (ix) 41,36,37
36    xie=six*dw/(dsqrt((wip+wep+wep)**2-six*dw*wepp)+wip+wep+wep)
      fxw=half*xie*(wepp+wipp)/(wep-wip)
      if (fxw .le. zro) ew=we+.1d+0
      if (fxw .le. zro) goto 39
      if (fxw .lt. one) ew=wi+dw*(four+fxw**2)/fiv
      if ((fxw .gt. one) .or. (ie .eq. 0)) ew=wi+dw*(9.d+0+fxw)/ten
      goto 39
37    ea=woppp
      eb=fiv*wipp+wepp
      ec=tlv*wip
      ed=-tlv*dw
      xie=cubic(ea,eb,ec,ed)
      if (xie .gt. zro) goto 38
      ew=we-.1d+0
      if (ew .gt. wi) goto 39
      write (2,113)
      goto 4
38    ew=wi+half*xie*(wip+wep+xie*(wipp-wepp)/six)
39    we=ew
!      if (we .gt. g2) go to 79
      if (we .gt. g2) then
       write (2,119)
       write (2,126)
       stop
      endif
      if (dabs(ew-dw-wi) .lt. 1.d-9) goto 43
      dw=we-wi
      call sorce(we,d)
      xe=d(1)
      wep=d(2)
      wepp=d(3)
      wrppp=d(4)
      nocon=nocon+1
      goto 35
40    if (iq .lt. 0) goto 41
      write (2,110) itle,n4,lr
41    h=thr*(wep+wip)/(wipp-wepp)
      hh=tlv*dw/(wipp-wepp)
      xie=hh/(dsqrt(h*h+hh)+h)
      if (mf) 44,42,45
42    ew=wi+xie*(wip+thr*wep-xie*(wepp-xie*wrppp/six))/four
      goto 39
43    emach=we/dsqrt(g7-g8*we*we)
!
!     iteration for emach completed
      ebet=dsqrt(emach*emach-one)
      epsi=g2*datan(g4*ebet)-datan(ebet)
      fpsi=epsi+two*eta/qt
      fmach=fmv(fpsi)
44    if (bmach .gt. fmach) goto 45
      bmach=fmach
      bpsi=fpsi
      mp=0
45    gpsi=fpsi-eta/qt
      gmach=fmv(gpsi)
      if (ie .eq. 1) ah=gmach
      rg=((g6+g5*gmach**2)**ga/gmach)**qt
      apsi=bpsi-eta/qt
      amach=fmv(apsi)
      ra=((g6+g5*amach**2)**ga/amach)**qt
      xa=ra*cse
      if (sfoa .gt. zro) xie=sfoa/sf+xe-xa-xoi
      if (sfoa .lt. zro) xie=xe-sfoa/sf-rg*cse-xoi
      xi=xe-xie
      xo=xi-xoi
      x1=xo+xo1
      if (iq .lt. 0) goto 48
      xb=((g6+g5*bmach**2)**ga/bmach)**qt
      if (lc .lt. 2) xc=((g6+g5*cmach**2)**ga/cmach)**qt
      c(1)=wi
      c(2)=xie*wip
      c(3)=half*wipp*xie*xie
      c(4)=ten*dw-xie*(four*wep-half*xie*wepp)-six*c(2)-thr*c(3)
      c(5)=xie*(sev*wep+eit*wip-xie*(wepp-thr*wipp/two))-fftn*dw
      c(6)=six*dw-thr*xie*(wep+wip)+half*xie*xie*(wepp-wipp)
      if (mf .eq. 0 .and. ix .eq. 0) c(5)=zro
      if (nr .eq. 0 .and. ix .eq. 0 .and. lr .eq. 0) c(5)=zro
      if (sfoa .eq. zro) c(6)=zro
      eoe=epsi/eta
      wippp=six*c(4)/xie/xie/xie
      weppp=six*(c(4)+four*c(5)+ten*c(6))/xie/xie/xie
46    write (2,99) m,n,eoe,bmach,cmach,gam,etad,rc,sf
      write (2,102) se,tk,wwo,wwop,emach,fmach,mc,ah
      if (lr .ne. 0) write (2,123) wi,wap,wapp,am,amp,ampp
      if (ise.eq.1 .and. lr.eq.0) write (2,123) wi,wip,whpp,am,amp,ampp
      if (ise .eq. 1) goto 47
      write (2,101) wi,wip,wipp,wippp,woppp
      write (2,98) we,wep,wepp,weppp,wrppp
47    write (2,94) c(1),c(2),c(3),c(4),c(5),c(6)
      write (2,95) xoi,xi,xo,yo,xie,xe,nocon
      if (ise .eq. 1) xc=xe
      if (ise .eq. 1) xa=xe+tye*cbet
48    nocon=0
      wip=whp
      if (qm .ne. one) goto 49
      if (pp .lt. zro) frip=zro
      if (pp .eq. zro) frip=-xo*sf
      if (pp .gt. zro) frip=pp-sf*xa
      if (iq .lt. 0) goto 50
      xoin=sf*xo+frip
      x1in=sf*x1+frip
      xiin=sf*xi+frip
      write (2,125) om,xoin,x1in,am,xiin
      if (iq .gt. 0) goto 67
49    if (n) 87,50,68
50    m=zonk*(md-1)+1
      jq=1
      n=nd
      ip=in
      if (qm .ne. one) return
!      call orez(c,6)
      c(:)=0.0d0
      if (iq .lt. 0) goto 51
      if (mq .ge. 0 .and. n .gt. 0) goto 51
      write (2,104)
      goto 52
51    write (2,105)
52    if (ip) 53,67,58
!
!     length of downstream velocity distribution, radial flow
53    wc=g2*cmach/dsqrt(cmach*cmach+g9)
      wb=g2*bmach/dsqrt(bmach*bmach+g9)
      wcb=wc-wb
      call sorce(wb,d)
      xb=d(1)
      wbp=d(2)
      wspp=d(3)
      wsppp=d(4)
      c(1)=wb
      wcp=zro
      if (lc) 54,55,56
54    xbc=thr*wcb/wbp
      wbpp=-two*wbp/xbc
      write (2,109) itle,n3
      goto 57
55    wbpp=wspp
      if (mcp .lt. 0) write (2,109) itle,n3
      if (mcp .lt. 0) xbcn=thr*wcb/wbp
      if (mcp .lt. 0) xbcm=-two*wbp/wbpp
      if (mcp .gt. 0) write (2,109) itle,n4
      if (mcp .gt. 0) xbcn=four*wcb/wbp
      if (mcp .gt. 0) xbcm=-thr*wbp/wbpp
      abcm=one-xbcn/xbcm
      if (abcm .lt. zro) goto 88
      xbc=xbcn/(dsqrt(abcm)+one)
      goto 57
56    wbpp=-wspp*ip/ten
      if (mcp.gt.0) xbcmn=cubic(wsppp/thr,thr*wbpp,tlv*wbp,-two*ten*wcb)
      if (mcp.lt.0) xbcmn=cubic(wsppp/six,wbpp,thr*wbp,-four*wcb)
      xbcmx=fiv*wcb/(dsqrt(wbp**2-ip*wcb*wspp/eit)+wbp)
      if (xc .gt. xb+xbcmx) xc=xb+xbcmx
      if (xc .lt. xb+xbcmn) xc=xb+xbcmn
      xbc=xc-xb
      if (mcp .lt. 0) write (2,109) itle,n4
      if (mcp .gt. 0) write (2,109) itle,n5
57    c(2)=xbc*wbp
      c(3)=half*xbc*xbc*wbpp
      if (mcp .lt. 0) c(4)=four*wcb-thr*c(2)-two*c(3)
      if (mcp .lt. 0) c(5)=-thr*wcb+two*c(2)+c(3)
      if (mcp .gt. 0) c(4)=ten*wcb-six*c(2)-thr*c(3)
      if (mcp .gt. 0) c(5)=-fftn*wcb+eit*c(2)+thr*c(3)
      if (mcp .gt. 0) c(6)=six*wcb-thr*c(2)-c(3)
      if (lc .lt. 0) c(5)=zro
      if (lc .le. 0) c(6)=zro
      xc=xb+xbc
      goto 63
!
!     length of downstream mach no. distribution, radial flow
58    call conic(bmach,d)
      xb=d(1)
      bmp=d(2)
      smpp=d(3)
      smppp=d(4)
      cbm=cmach-bmach
      c(1)=bmach
      bmpp=smpp*ip/ten
      if (lc .ne. 0) goto 59
      if (mcp .lt. 0) write (2,108) itle,n3
      if (mcp .lt. 0) xbcn=thr*cbm/bmp
      if (mcp .lt. 0) xbcm=-two*bmp/bmpp
      if (mcp .gt. 0) write (2,108) itle,n4
      if (mcp .gt. 0) xbcn=four*cbm/bmp
      if (mcp .gt. 0) xbcm=-thr*bmp/bmpp
      abcm=one-xbcn/xbcm
      if (abcm .lt. zro) goto 88
      xbc=xbcn/(dsqrt(abcm)+one)
      xc=xb+xbc
      goto 62
59    if (lc .ne. 1) goto 61
      do k=1,nd
       read (9) ax(k),axm(k),axmp(k)
       if (k .eq. 1) dx=xb-ax(1)
       axis(1,k)=ax(k)+dx
      enddo
      if (axmp(2) .eq. zro) call scond(ax,axm,axmp,nd)
      axm(1)=bmach
      axmp(1)=bmp
      xc=axis(1,nd)
      xbc=xc-xb
      write (2,111) itle
      goto 63
61    if (mcp.gt.0) xbcmn=cubic(smppp/thr,thr*bmpp,tlv*bmp,-two*ten*cbm)
      if (mcp.lt.0) xbcmn=cubic(smppp/six,bmpp,thr*bmp,-four*cbm)
      xbcmx=fiv*cbm/(dsqrt(bmp**2+ip*cbm*smpp/eit)+bmp)
      if (xc .gt. xb+xbcmx) xc=xb+xbcmx
      if (xc .lt. xb+xbcmn) xc=xb+xbcmn
      xbc=xc-xb
      if (mcp .lt. 0) write (2,108) itle,n4
      if (mcp .gt. 0) write (2,108) itle,n5
62    c(2)=xbc*bmp
      c(3)=half*xbc*xbc*bmpp
      if (mcp .lt. 0) c(4)=four*cbm-thr*c(2)-two*c(3)
      if (mcp .lt. 0) c(5)=-thr*cbm+two*c(2)+c(3)
      if (mcp .gt. 0) c(4)=ten*cbm-six*c(2)-thr*c(3)
      if (mcp .gt. 0) c(5)=-fftn*cbm+eit*c(2)+thr*c(3)
      if (mcp .gt. 0) c(6)=six*cbm-thr*c(2)-c(3)
      if (lc .le. 0) c(6)=zro 
63    cpp=zro
      cmp=zro
      if (mcp .lt. 0) cpp=(two*c(3)+six*c(4)+tlv*c(5))/xbc**2
      bppp=six*c(4)/xbc/xbc/xbc
      cppp=six*(c(4)+four*c(5)+ten*c(6))/xbc/xbc/xbc
      xd=xc+tye*cbet
      write (2,100) m,n,np,gam,etad,rc,sf
      if (ip) 64,67,65
64    write (2,116) wb,wbp,wbpp,bppp,wspp,wc,wcp,cpp,cppp,wsppp
      goto 66
65    write (2,117) bmach,bmp,bmpp,bppp,smpp,cmach,cmp,cpp,cppp,smppp
66    write (2,94) c(1),c(2),c(3),c(4),c(5),c(6)
      write (2,118) amach,xa,xb,xbc,xc,xd
      xain=sf*xa+frip
      yain=sf*xa*dtan(eta)
      xbin=sf*xb+frip
      xcin=sf*xc+frip
      xdin=sf*xd+frip
      tyin=sf*tye
      write (2,120) xain,yain,xbin,xcin,xdin,tyin
!     the n=0 goto 4 on next line will fail with a single input card
67    if (n) 87,4,68
68    if (mq .lt. 0) goto 69
!
!     calculate axial distribution
      write (2,103) iaxis
!69    fn=n-1
69    fn=real(n-1,k8)
      l=int((n+40)/41,K4)
      if (ip .ne. 0) xie=xbc
      if (ip .ne. 0) xi=xb
      q=zro
      do k=1,n
       if (ise .eq. 1 .and. lc .eq. 1) goto 72
       if (ip .ne. 0) goto 70
       if (nx .eq. 0) q=((n-k)/fn)**2
       if (nx .ne. 0) q=((n-k)/fn)**(nx*1.d-1)
       goto 71
70     if (lc .eq. 1) goto 72
       q=(k-1)/fn
71     axis(1,k)=xie*q+xi
72     rmach=one
       if (ise .eq. 1) goto 75
       if (axis(1,k) .lt. one+1.d-9) goto 74
       ab=axis(1,k)**(rga/qt)
       if (ab .lt. two) sm=((one+dsqrt(ab*gm-gm))**ga)**2
       if (ab .ge. two) sm=(ab/g5)**g7
!       if (ab .ge. two) sm=(ab/gs)**g7
73     cm=sm**g5
       fq=sm*(g6+g5*sm-cm*ab)/(sm-one)/g5/g6
       sm=sm-fq
       if (dabs(fq) .gt. 1.d-9) goto 73
       rmach=dsqrt(sm)
74     if (ip .lt. 1) goto 78
75     if (lc .eq. 1) goto 76
       xm=c(1)+q*(c(2)+q*(c(3)+q*(c(4)+q*(c(5)+q*c(6)))))
       if (ise .eq. 1 .or. k .eq. 1) goto 77
       if (rmach .lt. xm) write (2,124) k,rmach,xm
       goto 77
76     xm=axm(k)
77     xmp=(c(2)+q*(two*c(3)+q*(thr*c(4)+q*(four*c(5)+q*fiv*c(6)))))/xie
       if (lc .eq. 1) xmp=axmp(k)
       xmpp=two*(c(3)+q*(thr*c(4)+q*(six*c(5)+q*ten*c(6))))/xie/xie
       xmppp=six*(c(4)+q*(four*c(5)+ten*q*c(6)))/xie/xie/xie
       gmm=xm*xm+g9
       gq=dsqrt(gmm)
       w=g2*xm/gq
       wm=g9*g2/gq/gmm
       wp=wm*xmp
       wpp=wm*(xmpp-thr*xm*xmp*xmp/gmm)
       gmp=fiv*xm*xm*xmp*xmp/gmm-thr*xm*xmpp-xmp*xmp
       wppp=wm*(xmppp+thr*xmp*gmp/gmm)
       if (mq .lt. 0) goto 83
       if (mod(k-1,l) .ne. 0) goto 83
       goto 82
78     w=c(1)+q*(c(2)+q*(c(3)+q*(c(4)+q*(c(5)+q*c(6)))))
       wp=(c(2)+q*(two*c(3)+q*(thr*c(4)+q*(four*c(5)+q*fiv*c(6)))))/xie
       wpp=two*(c(3)+q*(thr*c(4)+q*(six*c(5)+q*ten*c(6))))/xie/xie
       wppp=six*(c(4)+q*(four*c(5)+ten*q*c(6)))/xie/xie/xie
       gww=g7-w*w*g8
       if (gww .gt. zro) goto 80
79     write (2,119)
       goto 4
80     gw=dsqrt(gww)
       xm=w/gw
       if (k .eq. 1 .or. k .eq. n) goto 81
       if (ip .eq. 0 .and. rmach .gt. xm) write (2,124) k,rmach,xm
       if (ip .ne. 0 .and. rmach .lt. xm) write (2,124) k,rmach,xm
81     xmw=g7/gw/gww
       xmp=xmw*wp
       xmpp=xmw*(wpp+thr*g8*w*wp*wp/gww)
       gwp=fiv*w*w*wp*wp*g8/gww+thr*w*wpp+wp*wp
       xmppp=xmw*(wppp+thr*wp*g8*gwp/gww)
       if (mq .lt. 0) goto 83
       if (mod(k-1,l) .ne. 0) goto 83
82     xinch=sf*axis(1,k)+frip
       write (2,97) k,axis(1,k),xinch,xm,xmp,xmpp,xmppp,w,wp,wpp,wppp
       if (mod(k+l-1,10*l) .eq. 0) write (2,115)
83     axis(3,k)=xm
       axis(2,k)=zro
       axis(5,k)=ie*half*(xm-one/xm)*wp/w 
       xbet=dsqrt(xm**2-one)
       axis(4,k)=g2*datan(g4*xbet)-datan(xbet)
      enddo
      if (iq .eq. 0 .and. ip .eq. 0 .and. m .le. 0) goto 50
!     the n=0 goto 4 on next line will fail with a single input card
      if (m) 87,4,85
85    if (ip .ne. 0) return
      do k=1,n
       do j=1,5
        taxi(j,k)=axis(j,k)
       enddo
      enddo
      return
87    lv=-1
      return
88    write (2,114)
      goto 4
89    write (2,112)
      goto 4
90    write (2,121) nocon
      goto 4
91    stop
!
92    format (16i5)
93    format (8f10.3)
94    format (1x,9x,'C1=',f11.7,3x,'C2=',f12.8,3x,'C3=',1pe15.7,3x,'C4='
     &,e15.7,3x,'C5=',e15.7,3x,'C6=',e15.7/)
95    format (1x,9x,'XOI=',f12.8,3x,'XI=',f12.8,3x,'XO=',f12.8,3x,'YO=',
     &f12.8,3x,'XIE=',f12.8,3x,'XE=',f12.8,i5,' ITERATIONS'/)
96    format (1x,4x,'CURVE FROM MACH 1,    XIE=',f12.8,'   C2=',f12.8,' 
     &  C3=',1pe15.7,'   C4=',e15.7,'   X1=',0pf12.8/)
97    format (1x,i3,2f10.5,f10.6,1p3e14.6,0pf10.6,1p3e14.6)
98    format (1x,9x,'WE=',f12.8,4x,'WEP=',f12.8,4x,'WEPP=',1pe15.7,4x,'W
     &EPPP=',e15.7,4x,'WRPPP=',e15.7/)
99    format (1x,4x,'NO. OF POINTS ON 1ST CHAR. (M)=',i3,5x,'NO. OF POIN
     &TS ON AXIS (N)=',i3,5x,'EPSI/ETA=',f8.5,4x,'BMACH=',f9.5,4x,'CMACH
     &=',f9.5//5x,'GAMMA=',f7.4,5x,'INFLECTION ANG. (ETA)=',f8.4,2x,'DEG
     &REES',5x,'RAD. OF CURV. (RC)=',f11.6,5x,'SCALE FACTOR (SF)=',f13.8
     &/)
100   format (1x,4x,'NO. OF POINTS ON 1ST CHAR. (M)=',i3,5x,'NO. OF POIN
     &TS ON AXIS (N)=',i3,5x,'NO. OF POINTS ON LAST CHAR. (NP)=',I3//5x,
     &'GAMMA=',f7.4,5x,'INFLECTION ANG. (ETA)=',f8.4,2x,'DEGREES',5x,'RA
     &D. OF CURV. (RC)=',f13.8,5x,'SCALE FACTOR (SF)=',f11.6/)
101   format (1x,9x,'WI=',f12.8,4x,'WIP=',f12.8,4x,'WIPP=',1pe15.7,4x,'W
     &IPPP=',e15.7,4x,'WOPPP=',e15.7/)
102   format (1x,4x,'Y*=',f10.8,4x,'RMASS=',f10.8,4x,'WWO=',f10.7,4x,'WW
     &OP=',f11.8,4x,'EMACH=',f8.5,4x,'FMACH=',f10.7,4x,a4,'H=',f9.5/)
103   format (1x,1x,a4/' POINT',4x,'X',7x,'X(IN)',3x,'MACH NO.',4x,'DM/D
     &X',8x,'D2M/DX2',7x,'D3M/DX3',7x,'W=Q/A*',5x,'DW/DX',8x,'D2W/DX2',7
     &x,'D3W/DX3'/)
104   format ('0',//)
105   format (1x)
106   format (1x,3a4,' THROAT CONTOUR,',a4,'-DEG AXIAL VELOCITY DISTRIBU
     &TION FROM SONIC POINT'/)
107   format ('1',3a4,' INVISCID CONTOUR,',a4,a4,'AXIAL MACH NUMBER DIST
     &RIBUTION FROM THROAT CHARACTERISTIC WHICH HAS',i4,' POINTS'/)
108   format (1x,3a4,' DOWNSTREAM CONTOUR,',a4,'-DEG AXIAL MACH NUMBER D
     &ISTRIBUTION'/)
109   format (1x,3a4,' DOWNSTREAM CONTOUR,',a4,'-DEG AXIAL VELOCITY DIST
     &RIBUTION'/)
110   format (1x,3a4,' THROAT CONTOUR,',a4,'-DEG AXIAL VELOCITY DISTRIBU
     &TION FROM THROAT CHARACTERISTIC WHICH HAS',i4,' POINTS'/)
111   format (1x,3a4,' DOWNSTREAM CONTOUR'/)
112   format (1x,'SOLUTION TO CUBIC EQUATION IS NEGATIVE')
113   format (1x,'RC IS TOO LARGE TO ALLOW A SOLUTION')
114   format (1x,'BMACH IS TOO SMALL TO ALLOW A SOLUTION')
115   format (1x)
116   format ('0',9x,'WB=',f12.8,4x,'WBP=',f12.8,4x,'WBPP=',1pe15.7,4x,'
     &WBPPP=',e15.7,5x,'WSPP=',e15.7//10x,'WC=',0pf12.8,4x,'WCP=',f12.8,
     &4x,'WCPP=',1pe15.7,4x,'WCPPP=',e15.7,4x,'WSPPP=',e15.7)
117   format (1x,9x,'BMACH=',f9.5,4x,'BMP=',f12.8,4x,'BMPP=',1pe15.7,4x,
     &'BMPPP=',e15.7,5x,'SMPP=',e15.7//10x,'CMACH=',0pf9.5,4x,'CMP=',f12
     &.8,4x,'CMPP=',1pe15.7,4x,'CMPPP=',e15.7,4x,'SMPPP=',e15.7/)
118   format (1x,9x,'AMACH=',f11.7,4x,'XA=',f11.7,4x,'XB=',f11.7,4x,'XBC
     &=',f11.7,4x,'XC=',f12.7,4x,'XD=',f12.7/)
119   format ('0','VELOCITY GREATER THAN THEORETICAL MAXIMUM VALUE')
120   format (1x,9x,'XA(IN)=',f11.7,', YA(IN)=',f11.7,', XB(IN)=',f12.7,
     &', XC(IN)=',f12.7,', XD(IN)=',f12.7,', YD(IN)=',f11.7/)
121   format ('1','NO CONVERGENCE IN',i4,'ITERATIONS')
122   format ('1',3a4,' INVISCID CONTOUR,',a4,a4,' AXIAL MACH NUMBER DIS
     &TRIBUTION FROM SONIC POINT'/)
123   format (1x,9x,'WI=',f12.8,4x,'WIP=',f12.8,4x,'WIPP=',1pe15.7,4x,'M
     &I=',0pf12.8,4x,'MIP=',f12.8,4x,'MIPP=',1pe15.7/)
124   format (1x,i3,' RMACH=',2f12.8)
125   format (1x,9x,'MACH',f11.8,' AT',f11.7,' IN.,   MACH 1 AT',f11.7,'
     & IN.,   MACH',f11.8,' AT',f11.7,' IN.'/)
126   format (1x,'CHANGE CARD 3 INPUT MACH NUMBER'/)
      end subroutine axial
