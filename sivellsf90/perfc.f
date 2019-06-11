      subroutine perfc
!
!     to obtain the inviscid contour of the nozzle
!
      use kinddefine
      use gg, only:gam,g1,g2,g4,g5,g6,g7,g8,ga,qt
      use cline, only:axis,taxi,frip,zonk,seo,cse
      use coord, only:s,fs,waltan,sd,wmn,ttr,dmdx,spr,dpx,secd,xbin,xcin
     &,gma,gmb,gmc,gmd
      use work, only:wall,wax,way,wan,a,b,fclast
      use prop, only:xbl,conv
      use param, only:etad,rc,amach,bmach,cmach,emach,gmach,sf,wwo,wwop,
     &qm,cbet,xe,eta,epsi,bpsi,xo,yo,rrc,sdo,xb,xc,ah,se,tye,xa
      use troat
      use contr, only:itle,ie,lr,it,jb,jq,kat,kbl,king,ko,nocon,mc,ip,iq
     &,ise,jc,m,mp,mq,n,np,nf,nut
!
      implicit none
!
      interface
       function fmv(psi)
        use kinddefine
        implicit none
        real(kind=K8) :: fmv
        real(kind=K8), intent(in) :: psi
       end function fmv
!
       subroutine ofeld(a,b,c,nocon)
        use kinddefine
        implicit none
        integer(kind=K4),intent(out) :: nocon
        real(kind=K8),dimension(5),intent(in) :: a,b
        real(kind=K8),dimension(5),intent(out) :: c
       end subroutine ofeld
!
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
      integer(kind=K4) :: i,ib,ichara,icy,inc,iprnt,j,jj,k,kan,kc,kit,kt
      integer(kind=K4) :: kut,l,last,lastp,lin,line,ll,lp,lq,lt,nag,nk
      integer(kind=K4) :: nl,nn
      real(kind=K8),dimension(6,150) :: chara
      real(kind=K8),dimension(150) :: su
      real(kind=K8),dimension(200) :: wdx,wtan,scdf
      real(kind=K8),dimension(100) :: yi
      real(kind=K8) :: add,an,as,bdd,bs,bx,ci,cpsi,cs,cy,cyp
      real(kind=K8) :: del,dhp,dm,dmxx,dn,dpxx,dsx,dydx,em,f,fn,fsx,half
      real(kind=K8) :: hs,la,one,ps,psi,r,rl,s1,s2,s3,sa,sb,sc,sdx
      real(kind=K8) :: six,sne,sprx,ss,st,sum1
      real(kind=K8) :: sumax,t,t1,t2,t3,thr,tne,tt,two,wmnx,x,xbet,xinch
      real(kind=K8) :: xm,xmu,xmur,xs,xs2,xs3,y,ye,yinch,ypx,ys,zro
      character(len=4,kind=K3) :: ifr,iwl,lst,ibl
!      COMMON /GG/ GAM,GM,G1,G2,G3,G4,G5,G6,G7,G8,G9,GA,RGA,QT
!      COMMON /CLINE/ AXIS(5,150),TAXI(5,150),WIP,X1,FRIP,ZONK,SEO,CSE
!      COMMON /COORD/ S(200),FS(200),WALTAN(200),SD(200),WMN(200),TTR(200
!     1),DMDX(200),SPR(200),DPX(200),SECD(200),XBIN,XCIN,GMA,GMB,GMC,GMD
!      COMMON /WORK/ WALL(5,200),WAX(200),WAY(200),WAN(200),A(5,150),B(5,
!     1150),FCLAST(5,150)
!      COMMON /PROP/ AR,ZO,RO,VISC,VISM,SFOA,XBL,CONV
!      COMMON /PARAM/ ETAD,RC,AMACH,BMACH,CMACH,EMACH,GMACH,FRC,SF,WWO,WW
!     1OP,QM,WE,CBET,XE,ETA,EPSI,BPSI,XO,YO,RRC,SDO,XB,XC,AH,PP,SE,TYE,XA
!      COMMON /TROAT/ FC(6,51)
!      COMMON /CONTR/ ITLE(3),IE,LR,IT,JB,JQ,JX,KAT,KBL,KING,KO,LV,NOCON,
!     1IN,MC,MCP,IP,IQ,ISE,JC,M,MP,MQ,N,NP,NF,NUT
!      DIMENSION CHARA(6,150), SU(150), WDX(200), WTAN(200), SCDF(200), YI
!     1(100)
      data zro/0.0d+0/,one/1.d+0/,two/2.d+0/,six/6.d+0/,half/5.d-1/
      data ifr/'FIRS'/,iwl/'WALL'/,lst/'LAST'/,ibl/'    '/,thr/3.d+0/
!      call orez (a,4*750+250)
      a(:,:)=0.0d0
      b(:,:)=0.0d0
      fclast(:,:)=0.0d0
      wall(:,:)=0.0d0
!
      chara(:,:)=0.0d0
      su(:)=0.0d0
      wdx(:)=0.0d0
      wtan(:)=0.0d0
      scdf(:)=0.0d0
      yi(:)=0.0d0
!
      cpsi=g2*datan(g4*cbet)-datan(cbet)
      if (jq.gt.0) goto 6
      if (lr.eq.0) goto 4
!
!     throat characteristic values
      sumax=(se/seo)**(ie+1)
      if (qm.eq.one) sumax=one
      lq=zonk*(lr-1)+1
      nl=n+lq-1
      do j=1,lq
       if (qm.ne.one) goto 1
       fc(1,j)=fc(1,j)*se+xo
       fc(2,j)=fc(2,j)*se
1      fclast(1,j)=fc(1,j)
       fclast(2,j)=fc(2,j)
       fclast(3,j)=fc(3,j)
       fclast(4,j)=fc(4,j)
       fclast(5,j)=fc(5,j)
       if (mq.lt.0) goto 3
       if (j.gt.1) goto 2
       write (2,93) itle
       write (2,99) ibl
2      xmu=conv*dasin(one/fclast(3,j))
       psi=conv*fclast(4,j)
       an=conv*fclast(5,j)
       xinch=sf*fclast(1,j)+frip
       yinch=sf*fclast(2,j)
       write (2,103) j,(fclast(k,j),k=1,3),xmu,psi,an,xinch,yinch
       if (mod(j,10).eq.0) write (2,98)
3      su(j)=fc(6,j)/sumax
      enddo
4     if (ise.eq.0) goto 8
!
!     initial characteristic values if non-radial flow
      do k=1,m
       a(2,k)=(k-1)*tye/(m-1)
       a(1,k)=a(2,k)*cbet+xe
       a(3,k)=cmach
       a(4,k)=cpsi
       a(5,k)=zro
      enddo
      goto 10
!
!     final characteristic values if radial flow
6     nl=n+np-1
      fn=np-1
      do jj=1,np
       if (ie.eq.0) f=(jj-1)/fn
       if (ie.eq.1) f=two*dsin(half*eta*(jj-1)/fn)/se
       fclast(2,jj)=f*tye
       fclast(1,jj)=fclast(2,jj)*cbet+xc
       fclast(3,jj)=cmach
       fclast(4,jj)=cpsi
       fclast(5,jj)=zro
       su(jj)=f**(ie+1)
      enddo
!
!     initial characteristic values if radial flow
8     em=eta/(m-1)
      do k=1,m
       t=(k-1)*em
       if (ip.eq.0) xm=fmv(epsi+t/qt)
       if (ip.ne.0) xm=fmv(bpsi-t/qt)
       r=((g6+g5*xm**2)**ga/xm)**qt
       xbet=dsqrt(xm**2-one)
       a(1,k)=r*dcos(t)
       a(2,k)=r*dsin(t)
       a(3,k)=xm
       a(4,k)=g2*datan(g4*xbet)-datan(xbet)
       a(5,k)=t
      enddo
      if (ie.eq.1 .and. ip.eq.0) a(5,1)=taxi(5,1)
      if (ie.eq.1 .and. ip.ne.0) a(5,1)=axis(5,1)
10    do j=1,5
       wall(j,1)=a(j,m)
      enddo
      line=1
      if (mq.lt.0) goto 14
      if (ise.eq.1) goto 12
      if (jq.eq.0) write (2,91) itle
      if (jq.eq.1) write (2,94) itle
      goto 13
12    write (2,102) itle
13    write (2,106) line
14    su(1)=zro
      if (ie.eq.0) bx=one/se
      nn=1
      do k=1,m
       do j=1,5
        b(j,k)=a(j,k)
       enddo
      enddo
      last=m-1
      goto 20
16    last=m
      line=2
      if (ip.ne.0) goto 38
17    do j=1,5
       b(j,1)=taxi(j,line)
      enddo
      do j=1,last
       k=j
       call ofeld(a(1,k),b(1,k),b(1,k+1),nocon)
       if (nocon.ne.0) goto 83
      enddo
20    lastp=last+1
      if (line.lt.lastp) lp=line
      nk=1+lp/52
      la=conv*dasin(one/b(3,nn))
      iprnt=0
      ichara=0
      if (jc.eq.0) go to 21
      kc=iabs(jc)
      if (jc.gt.0 .and. jq.ne.0) go to 21
      if (jc.lt.0 .and. jq.eq.0) go to 21
      ichara=1
      if (kc.gt.100 .and. kc.lt.101+line) iprnt=1
      if (nn.eq.1 .and. mod(line-1,kc).eq.0) iprnt=1
      if (nn.gt.1 .and. mod(nn-1,kc).eq.0) iprnt=1
21    do j=nn,lastp
       if (ie.eq.1) bx=two*b(2,j)/se**2
       xm=b(3,j)
       xmur=dasin(one/xm)
       xmu=conv*xmur
       psi=b(4,j)*conv
       an=b(5,j)*conv
       if (b(2,j).eq.zro) an=zro
       if (ip.eq.0.or.la.gt.45) goto 22
       s(j)=b(1,nn)-b(1,j)
!     mass integration with  respect  to  x
       dsx=one/dcos(b(5,j)-xmur)
       if (b(2,j).eq.zro) dsx=xm/dsqrt(xm**2-one)
       goto 23
22     s(j)=b(2,j)-b(2,nn)
!     mass integration with respect to y
       if (ip.eq.0) dsx=one/dsin(xmur+b(5,j))
       if (ip.ne.0) dsx=one/dsin(xmur-b(5,j))
       if (b(2,j).eq.zro) dsx=xm
23     if (ichara.eq.0 .or. j.ne.line) goto 24
       chara(1,j)=b(1,j)
       chara(2,j)=b(2,j)
       chara(3,j)=xm
       chara(4,j)=xmu
       chara(5,j)=psi
       chara(6,j)=an
24     fs(j)=dsx*bx/(g6+g5*xm**2)**ga
       if (mq.ge.0 .and. line.eq.1) goto 25
       if (iprnt.eq.0) goto 27
       if (j.gt.nn) goto 25
       if (ip.eq.0) write (2,104) itle
       if (ip.ne.0) write (2,105) itle
       write (2,106) line
25     if ((nk.gt.1) .and. (mod(j,nk).eq.0)) goto 26
       xinch=sf*b(1,j)+frip
       yinch=sf*b(2,j)
       write (2,103) j,b(1,j),b(2,j),xm,xmu,psi,an,xinch,yinch
26     if (mod(j,10*nk).eq.0) write (2,98)
27     continue
      enddo
!
!     integration and interpolation for mass flow
      sa=zro
      sb=zro
      sc=zro
      sum1=su(nn)
      kan=(lastp-nn)/2
      do j=1,kan
       k=nn+2*j
       kt=k
       as=s(k-1)-s(k-2)
       bs=s(k)-s(k-1)
       cs=as+bs
       s1=(two-bs/as)*cs/six
       s3=(two-as/bs)*cs/six
       s2=cs-s1-s3
       add=s1*fs(k-2)+s2*fs(k-1)+s3*fs(k)
       sum1=add+sum1
       if (line.eq.1) goto 28
       del=one-sum1
       if (del) 30,29,28
28     continue
      enddo
      if (line.eq.1) write (2,96) sum1
      if (line.eq.1) goto 16
      bs=s(k+1)-s(k)
      kt=k+1
      dn=two*del/bs
      sc=dn/(fs(k)+dsqrt(fs(k)**2+(fs(kt)-fs(k))*dn))
      sb=one-sc
      goto 34
29    sc=one
      goto 34
30    s2=bs*(two+cs/as)/six
      s3=bs*(two+as/cs)/six
      s1=bs-s2-s3
      bdd=s1*fs(k-2)+s2*fs(k-1)+s3*fs(k)
      if (bdd+del) 31,32,33
31    dn=two*(add+del)/as
      sb=dn/(fs(k-2)+dsqrt(fs(k-2)**2+(fs(k-1)-fs(k-2))*dn))
      sa=one-sb
      go to 34
32    sb=one
      go to 34
33    dn=two*del/bs
      sc=one+dn/(fs(k)+dsqrt(fs(k)**2+(fs(k)-fs(k-1))*dn))
      sb=one-sc
34    do j=1,5
       wall(j,line)=b(j,kt-2)*sa+b(j,kt-1)*sb+b(j,kt)*sc
      enddo
      if (iprnt.eq.1) write (2,107) (wall(j,line),j=1,3)
      last=kt
      if (n-line) 42,41,36
36    line=line+1
      do k=1,5
       do l=1,150
        a(k,l)=b(k,l)
       enddo
      enddo
      if (ip.eq.0) go to 17
38    do j=1,5
       b(j,1)=axis(j,line)
      enddo
      do j=1,last
       k=j
       call ofeld (b(1,k),a(1,k),b(1,k+1),nocon)
       if (nocon.ne.0) goto 83
      enddo
      goto 20
41    if (ip.ne.0) goto 42
      if (lr.eq.0.or.it.ne.0) goto 49
42    if (line.eq.nl-1) goto 48
      nn=nn+1
      line=line+1
      do k=1,5
       do l=1,150
        a(k,l)=b(k,l)
       enddo
      enddo
      do k=1,5
       do l=1,150
        b(k,l)=fclast(k,l)
       enddo
      enddo
      if ((lr.ne.0).and.(jq.eq.0)) goto 46
      do j=nn,last
       k=j
       call ofeld(b(1,k),a(1,k),b(1,k+1),nocon)
       if (nocon.ne.0) goto 83
      enddo
      goto 20
46    do j=nn,last
       k=j
       call ofeld(a(1,k),b(1,k),b(1,k+1),nocon)
       if (nocon.ne.0) go to 83
      enddo
      goto 20
48    if (ip.ne.0) goto 64
!
!     integration of slopes
49    ib=1
      if (iabs(jb).gt.1) ib=2
      lt=0
      if (it.ne.0) lt=ib
      nut=(line-1)/ib+2-lt
      wall(1,line+1)=xo
      wall(5,line+1)=zro
      yi(nut)=wall(2,1)
      y=yi(nut)
      lin=2*((line-lt)/2)
      do j=2,lin,2
       i=nut-j
       ss=wall(1,j)-wall(1,j-1)
       tt=wall(1,j+1)-wall(1,j)
       st=ss+tt
       s1=ss*(two+tt/st)/six
       s2=ss*(two+st/tt)/six
       s3=ss-s1-s2
       t3=tt*(two+ss/st)/six
       t2=tt*(two+st/ss)/six
       t1=tt-t2-t3
       y=y+s1*dtan(wall(5,j-1))+s2*dtan(wall(5,j))+s3*dtan(wall(5,j+1))
       if (ib.eq.1) yi(i+1)=y
       y=y+t1*dtan(wall(5,j-1))+t2*dtan(wall(5,j))+t3*dtan(wall(5,j+1))
       if (ib.eq.1) yi(i)=y
       if (ib.eq.2) yi(i+j/2)=y
      enddo
      if (lr.ne.0.and.line.eq.lin) goto 51
      x=wall(1,line-lt)-xo
      yi(1)=yi(2)-x*(dtan(wall(5,line-lt))+half*x*sdo)/thr
51    do l=2,nut
       jj=1+ib*(nut-l)
       wax(l)=wall(1,jj)
       way(l)=wall(2,jj)
       wmn(l)=wall(3,jj)
       wan(l)=conv*wall(5,jj)
       waltan(l)=dtan(wall(5,jj))
      enddo
      wax(1)=xo
      way(1)=yo
      wan(1)=zro
      wmn(1)=wwo/dsqrt(g7-g8*wwo**2)
      waltan(1)=zro
      if (nf.ge.0) goto 54
!
!     smooth upstream contour if desired
      call neo
      do j=1,nut
       waltan(j)=dtan(wan(j)/conv)
      enddo
54    call scond (wax,waltan,secd,nut)
      secd(1)=sdo
      secd(nut)=zro
      ko=nut+mp
      if (mp.eq.0) goto 56
!
!     radial flow section coordinates
      sne=dsin(eta)
      tne=dtan(eta)
      dm=(amach-gmach)/mp
      do l=1,mp
       ll=nut+l
       wmn(ll)=gmach+l*dm
       rl=((g5*wmn(ll)**2+g6)**ga/wmn(ll))**qt
       wax(ll)=rl*cse
       way(ll)=rl*sne
       wan(ll)=etad
       waltan(ll)=tne
       secd(ll)=zro
      enddo
56    if (mq.lt.0) goto 60
      if (jc.le.0) goto 58
      write (2,105) itle
      write (2,99) lst
      do k=1,lp,nk
       i=(k-1)/nk+1
       xinch=sf*chara(1,k)+frip
       yinch=sf*chara(2,k)
       write (2,103) k,(chara(j,k),j=1,6),xinch,yinch
       if (mod(i,10).eq.0) write (2,98)
      enddo
58    if (ise.eq.0) write (2,91) itle
      if (ise.eq.1) write (2,102) itle
      write (2,84) rc,etad,amach,bmach,cmach,emach,mc,ah
      if (nocon.ne.0) goto 59
      write (2,100) iwl
      write (2,85) (k,wax(k),way(k),wmn(k),wan(k),waltan(k),secd(k),k=1,
     &nut)
      if ((lr.eq.0) .and. (n.lt.42)) goto 59
      if ((lr.ne.0) .and. (n+lr.lt.27)) goto 59
      nocon=1
      goto 58
59    write (2,87)
      nocon=0
!
!     comparison of contour with parabola and hyperbola
60    do j=1,nut
       xs=(wax(j)-xo)/yo
       xs2=xs**2
       xs3=xs**3
       ys=way(j)/yo
       ye=yi(j)/yo
       ps=one+half*xs2*rrc
       dhp=one+xs2*rrc
       hs=dsqrt(dhp)
       if (j.gt.1) goto 61
       if (mq.lt.0) goto 62
       write (2,88) j,xs,ys,ye,ps,hs
       goto 62
61     ypx=waltan(j)/xs
       cy=(ps-ys)/xs3
       ci=(ps-ye)/xs3
       if (j.eq.2) icy=int(1.d+6*(dabs(cy)-dabs(ci)),K4)
       if (mq.lt.0) go to 63
       cyp=(rrc-ypx)/xs/thr
       write (2,88) j,xs,ys,ye,ps,hs,cy,ci,cyp
62     if (mod(j,10).eq.0) write (2,98)
      enddo
63    write (2,97) icy
      if (iq.gt.0) goto 70
      jq=1
      return
64    line=nl
      do j=1,5
       wall(j,nl)=fclast(j,np)
      enddo
!
!     smooth downstream contour if desired
      if (nf.lt.0) call neo
      do j=1,nl
       wdx(j)=wall(1,j)
       wtan(j)=dtan(wall(5,j))
      enddo
      call scond (wdx,wtan,scdf,nl)
      scdf(1)=zro
      scdf(nl)=zro
      if (jc.ge.0) goto 68
      write (2,104) itle
      write (2,99) ifr
      do k=1,lp,nk
       i=(k-1)/nk+1
       xinch=sf*chara(1,k)+frip
       yinch=sf*chara(2,k)
       write (2,103) k,(chara(j,k),j=1,6),xinch,yinch
       if (mod(i,10).eq.0) write (2,98)
      enddo
68    if (iq.lt.0) ko=1
      nag=ko-1
      king=line+nag
      do l=1,line
       wax(nag+l)=wall(1,l)
       way(nag+l)=wall(2,l)
       wmn(nag+l)=wall(3,l)
       wan(nag+l)=conv*wall(5,l)
       waltan(nag+l)=wtan(l)
       secd(nag+l)=scdf(l)
      enddo
      if (mq.lt.0) goto 71
      write (2,94) itle
      write (2,84) rc,etad,amach,bmach,cmach,emach,mc,ah
      write (2,100) iwl
      write (2,85) (k,wax(k),way(k),wmn(k),wan(k),waltan(k),secd(k),k=ko
     &,king)
      goto 71
70    king=ko
!
!     application of scale factor to non-dimensional coordinates
71    do k=1,king
       s(k)=sf*wax(k)+frip
       fs(k)=sf*way(k)
       ttr(k)=one+g8*wmn(k)**2
       spr(k)=one/ttr(k)**(one+g1)
       sd(k)=secd(k)/sf
      enddo
      if (ise.eq.1) xbin=zro
      if (ise.eq.0) xbin=xb*sf+frip
      xcin=xc*sf+frip
      call scond (s,wmn,dmdx,king)
      dmdx(1)=g7*wwop*wmn(1)**3/wwo**3/sf
      if (mp.eq.0 .or. iq.lt.0) go to 74
      do k=nut,ko
       dmdx(k)=wmn(k)*ttr(k)/(wmn(k)**2-one)/qt/sf/wax(k)
      enddo
      goto 75
74    if (ise.eq.0) dmdx(ko)=amach*ttr(ko)/(amach**2-one)/qt/sf/xa
75    if (iq.lt.1 .or. ise.eq.1) dmdx(king)=zro
      do k=1,king
       dpx(k)=-gam*wmn(k)*dmdx(k)*spr(k)/ttr(k)
      enddo
      jq=0
      kat=king
      if (iabs(mq).lt.2) goto 78
!
!     extension of parallel-flow contour
      kit=king+1
      kat=king+iabs(mq)
      kut=s(king)+half
      inc=s(king)-s(king-1)
      if (inc.lt.1) inc=1
      do k=kit,kat
       s(k)=kut+(k-king)*inc
       fs(k)=fs(king)
       wmn(k)=wmn(king)
       ttr(k)=ttr(king)
       spr(k)=spr(king)
       wan(k)=zro
       waltan(k)=zro
       dmdx(k)=zro
       dpx(k)=zro
       sd(k)=zro
      enddo
78    if (xbl.eq.zro) goto 79
      if (s(king-1).lt.xbl) goto 79
!
!     interpolate for values at specified station
      call twixt (s,gma,gmb,gmc,gmd,xbl,king,kbl)
      goto 80
79    kbl=kat+4
80    if (jb.gt.0) return
      if (ise.eq.0) goto 81
      write (2,102) itle
      write (2,92) rc,se,xcin
      goto 82
81    if (iq.gt.0) write (2,91) itle
      if (iq.le.0) write (2,95) itle,xbin,xcin,sf
      write (2,84) rc,etad,amach,bmach,cmach,emach,mc,ah
82    write (2,89)
      write (2,90) (k,s(k),fs(k),waltan(k),sd(k),wmn(k),dmdx(k),spr(k),d
     &px(k),k=1,king)
      if (kbl.gt.kat) return
      j=kbl-1
      fsx=gma*fs(j-2)+gmb*fs(j-1)+gmc*fs(j)+gmd*fs(j+1)
      wmnx=gma*wmn(j-2)+gmb*wmn(j-1)+gmc*wmn(j)+gmd*wmn(j+1)
      dmxx=gma*dmdx(j-2)+gmb*dmdx(j-1)+gmc*dmdx(j)+gmd*dmdx(j+1)
      dydx=gma*waltan(j-2)+gmb*waltan(j-1)+gmc*waltan(j)+gmd*waltan(j+1)
      sdx=gma*sd(j-2)+gmb*sd(j-1)+gmc*sd(j)+gmd*sd(j+1)
      sprx=gma*spr(j-2)+gmb*spr(j-1)+gmc*spr(j)+gmd*spr(j+1)
      dpxx=gma*dpx(j-2)+gmb*dpx(j-1)+gmc*dpx(j)+gmd*dpx(j+1)
      write (2,101) xbl,fsx,dydx,sdx,wmnx,dmxx,sprx,dpxx
      return
83    write (2,86) ip,nn,line,j
      return
!
84    format (1x,' RC=',f11.6,3x,'ETAD=',f8.4,' DEG',3x,'AMACH=',f10.7,3
     &x,'BMACH=',f10.7,3x,'CMACH=',f10.7,3x,'EMACH=',f10.7,3x,A4,'H=',f1
     &1.7/)
85    format (10(8x,i3,2x,1p6e15.7/))
86    format ('0','OFELD,IP=',i3,', NN=',i3,', LINE=',i3,', POINT=',i3)
87    format (1x,9x,'POINT X/YO',8x,'Y/YO',7x,'INT.Y/YO',7x,'PAR/YO',7x,
     &'HYP/YO       C(Y)',11x,'C(YI)',10x,'C(YP)'/)
88    format (1x,9x,i3,5f13.7,1p3e15.6)
89    format (1x,9x,'POINT',7x,'X(IN)',9x,'Y(IN)',9x,'DY/DX',8x,'D2Y/DX2
     &',7x,'MACH NO.',7x,'DM/DX',9x,'PE/PO',11x,'DPR/DX'/)
90    format (10(10x,i3,2x,0p6f14.7,1p2e16.5/))
91    format (1x,3a4,' UPSTREAM CONTOUR'/)
92    format (1x,' RC=',f11.7,',     STREAMLINE RATIO=',f11.8,',   TEST 
     &CONE BEGINS AT',f12.7,' IN.'/)
93    format (1x,3a4,' THROAT CHARACTERISTIC')
94    format (1x,3a4,' DOWNSTREAM CONTOUR'/)
95    format ('1',3a4,' INVISCID NOZZLE CONTOUR, RADIAL FLOW ENDS AT',f1
     &1.6,' IN., TEST CONE BEGINS AT',f11.6,' IN., SCALE FACTOR=',f9.4/)
96    format (1x,8x,'MASS =',f13.10/)
97    format (1x,9x,'ICY =',i13)
98    format (1x)
99    format (1x,8x,a4/8x,'POINT',8x,'X',14x,'Y',10x,'MACH NO.      MACH
     & ANG.(D)     PSI (D)     FLOW ANG.(D)        X(IN)',9x,'Y(IN)'/)
100   format (1x,8x,a4/8x,'POINT',8x,'X',14x,'Y',10x,'MACH NO.      FLOW
     & ANG.(D)     WALTAN',9x,'SECDIF'/)
101   format ('0',14x,6f14.7,1p2e16.5)
102   format ('1',3a4,' INVISCID CONTOUR'/)
103   format (1x,i10,2x,1p6e15.7,0p2f14.7)
104   format (1x,3a4,' INTERMEDIATE LEFT CHARACTERISTIC'/)
105   format (1x,3a4,' INTERMEDIATE RIGHT CHARACTERISTIC'/)
106   format (1x,' CHARACT',i4/8x,'POINT',8x,'X',14x,'Y',10x,'MACH NO.  
     &   MACH ANG.(D)      PSI (D)      FLOW ANG.(D)       X(IN)',9X,'Y(
     &IN)'/)
107   format (1x,'   CONTOUR  ',1p3e15.7/)
      end subroutine perfc
