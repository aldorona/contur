!     main part of
!     program contur(input,output,unit(1)=input,unit(2)=output)
!
!     nozzle contour program vev00028 for axisymmetric or planar flow
!     with radial flow region and/or with center-line velocity or mach
!     number distributions defined by polynomials.
!
!     correction applied for growth of turbulent boundary layer.
!     perfect gas is assumed with constant specific heat ratio, gam,
!     compressibility factor, zo, and recovery factor, ro, as inputs.
!     also input is gas constant, ar, in sq ft per sq second per deg r.
!     if vism is sutherlands temperature, viscosity follows sutherlands
!     law above vism, but is linear with temperature below vism.
!     if (vism .le. 1.d+0) viscosity=visc*temperature**vism mai
!
      program contur
      use kinddefine
      use gg, only:gam,gm,g1,g2,g3,g4,g5,g6,g7,g8,g9,ga,rga,qt
      use coord, only:dmdx,dpx,fs,s,sd,spr,ttr,waltan,wmn
      use corr, only:dr2,rco,drx,sl
      use prop, only:ar,zo,ro,visc,vism,sfoa,xbl,conv
      use param, only:ah,cmach,qm
      use jack, only:aj,sj,xj,yj
      use contr, only:itle,ie,it,jb,jq,jx,kat,kbl,king,ko,lv,nocon
!
      implicit none
!
      interface
       subroutine splind(x,y,tn2,tnl,l)
        use kinddefine
        implicit none
        integer(kind=K4),intent(in) :: l
        real(kind=K8),dimension(1),intent(in) :: x,y
        real(kind=K8),intent(in) :: tn2,tnl
       end subroutine splind
!
       subroutine xyz(xx,yy,yyp,yypp)
        use kinddefine
        implicit none
        real(kind=K8),intent(in) :: xx
        real(kind=K8),intent(out) :: yy,yyp,yypp
       end subroutine xyz
      end interface
!
      integer(kind=K4) :: ipp,jd,jj,k,kap,kup,nc
      real(kind=K8) :: bj,cn,csk,curv
      real(kind=K8) :: dy,one,slong,snk,tmax,two
      real(kind=K8) :: x,xend,xinc,xinc2
      real(kind=K8) :: xlow,xmax,xmid,xst,xx,ymax,yy,yyp,yypp,wang,zro
      character(len=4,kind=K3) :: l1,l2,l3,l4,l5,la,lb
      character(len=8,kind=K3) ::dc1,dc2,dc3,dc4,dc5,dc6,dc7,dca,dcb,dcc
      data zro/0.0d+0/,one/1.d+0/,two/2.d+0/,dc7/'CURVATUR'/
      data dc1/' D2Y/DX2'/,dc2/'        '/,dc3/'   ANGLE'/
      data dc4/'   DY/DX'/,dc5/'   DY/DS'/,dc6/'   DX/DS'/
      data l1/'   X'/,l2/'   Y'/,l3/'   S'/,l4/'    '/,l5/'DIFF'/
      conv=90.d+0/dasin(one)
      it=0
      nc=0
      la=l1
      lb=l4
      dca=dc4
      dcb=dc2
      jj=1000
      dcc=dc1
!
      open(unit=1,file='input.txt',status='old',access='sequential',
     &form='formatted',action='read',blank='null')
      open(unit=2,file='output.txt',status='unknown',access='sequential'
     &,form='formatted',action='write')
1     read (1,30,end=24) itle,jd
      if (itle(1) .eq. l4) goto 24
      ie=1+jd
      qt=one/(1+ie)
!
      read (1,28) gam,ar,zo,ro,visc,vism,sfoa,xbl
!     for gamma=1.4, g9=5, g8=.2, g7=1.2, g6=5/6, g5=1/6, g4=1/sqrt(6),
!     g3=1.5, g2=sqrt(6), g1=2.5
      gm=gam-one
      g1=one/gm
      g9=two*g1
      g8=one/g9
      g7=one+g8
      g6=one/g7
      g5=g8*g6
      rga=two*g5
      ga=one/rga
      g4=dsqrt(g5)
      g3=ga/two
      g2=one/g4
      if (ie .eq. 0) ah=zo
      if (ie .eq. 0) zo=one
      qm=one
      jx=0
2     jq=0
      lv=0
3     call axial
      if (lv .lt. 0) goto 1
      call perfc
      if (nocon .ne. 0) goto 24
      if ((jq .gt. 0) .or. (jx .gt. 0)) goto 3
      if (jb .gt. 0) call bound
      if (xbl .eq. 1.d+3) goto 5
      if (it .lt. 1) goto 4
      la=l3
      dca=dc5
      dcc=dc7
      kup=it
      kap=kup+1
      xend=zro
!
      read (1,28,end=24) (sj(k),k=1,kup),xst
      csk=one/dsqrt(one+drx(kat)**2)
      snk=csk*drx(kat)
      call splind(sl,rco,zro,snk,kat)
      goto 6
4     if (lv .gt. 0) goto 24
      if (jx .lt. 0) goto 1
      goto 2
5     continue
!
      read (1,28,end=24) xst,xlow,xend,xinc,bj,xmid,xinc2,cn
      if (xst .eq. xbl) xst=s(1)
      nc=cn-one
      if (jb .le. 0) call splind(s,fs,waltan(1),waltan(king),king)
      if (jb .gt. 0) call splind(s,rco,drx(1),drx(kat),kat)
      if (xend .gt. zro) goto 6
      lb=l5
      dcb=dc4
6     slong=s(king)-s(1)
      ipp=0
      write (2,25) itle,slong
      write (2,31) la,l2,dca,dc3,dcc,dcb,lb
      if (jb .gt. 0) goto 7
      write (2,26) xst,fs(1),waltan(1),zro,sd(1)
      xmax=slong+xst
      ymax=fs(king)
      tmax=waltan(king)
      goto 8
7     write (2,26) xst,rco(1),drx(1),zro,dr2
      xmax=s(kat)-s(1)+xst
      ymax=rco(kat)
      tmax=drx(kat)
8     if (it .gt. 0) goto 11
      jb=int(bj)
      if (xend .gt. zro) goto 10
      if (jb .lt. 0) goto 9
      kup=king-1
      kap=king-1
      goto 11
9     kup=-jb
      kap=kup+1
!
      read (1,28,end=24) (sj(k),k=1,kup)
      goto 11
10    if (xinc .gt. zro) kup=(xend-xlow)/xinc+1.d-2
      if (xmid .ne. zro) jj=(xmid-xlow)/xinc+1.d-2
      if (xmid .ne. zro) kup=jj+(xend-xmid)/xinc2+1.d-2
      if (jb .gt. 10) kup=jb
      if (jb .gt. 10) xinc=slong/bj
      kap=(xmax-xlow)/xinc+1
      if (xmid .ne. zro) kap=jj+(xmax-xmid)/xinc2+1
11    do k=1,kup
       if (xend .eq. zro) goto 12
       x=xlow+k*xinc
       if (k .gt. jj) x=xmid+(k-jj)*xinc2
       goto 13
12     if ((it .lt. 1) .and. (jb .ge.0)) x=s(k+1)
       if ((it .gt. 0) .or. (jb .lt.0)) x=sj(k)
13     xx=x-xst+s(1)
       if (k .lt. kap) call xyz(xx,yy,yyp,yypp)
       if (k .eq. kap) x=xmax
       if (k .ge. kap) yy=ymax
       if (k .ge. kap) yyp=tmax
       if (k .ge. kap) yypp=zro
       if (it .lt. 1) goto 16
       if (ipp .gt. 0) goto 14
       yj(k)=yy
       aj(k)=dasin(yyp)
       wang=conv*aj(k)
       curv=yypp/dcos(aj(k))
       write (2,26) x,yy,yyp,wang,curv
       goto 18
14     yy=yy-s(1)+xst
       xj(k)=yy
       wang=conv*dacos(yyp)
       write (2,26) x,yy,yyp,wang
       goto 18
16     wang=conv*datan(yyp)
       if ((xend .eq. zro) .and. (jb .ge. 0)) dy=yyp-waltan(k+1)
       if (jb .le. 0) goto 17
       fs(k+1)=yy
       waltan(k+1)=yyp
       sd(k+1)=yypp
17     if (xend.gt.zro.or.jb.lt.0) write (2,26) x,yy,yyp,wang,yypp
       if (xend.eq.zro.and.jb.ge.0) write (2,26) x,yy,yyp,wang,yypp,dy
18     if (mod(k,10) .eq. 0) write (2,29)
       if (mod(k,50) .ne. 0) goto 19
       write (2,25) itle,slong
       write (2,31) la,l2,dca,dc3,dcc,dcb,lb
19     continue
      enddo
      if ((it .gt. 0) .and. (ipp .eq. 1)) call plate
      if (ipp .ge. nc) goto 20
      ipp=ipp+1
      write (2,25) itle,slong
      write (2,31) la,l2,dca,dc3,dcc
      write (2,26) xst,rco(1),drx(1),zro,dr2
      goto 11
20    if ((ipp .gt. 0) .or. (jx .lt. 0)) goto 1
      if (it .eq. 0) goto 21
      ipp=1
      call splind(sl,s,one,csk,kat)
      write (2,29)
      write (2,31) l3,l1,dc6,dc3
      write (2,26) xst,xst,one,zro
      goto 11
21    if (jb) 1,2,22
22    call splind(s,wmn,dmdx(1),dmdx(king),king)
      do k=1,kup
       x=xlow+k*xinc
       if (xend .eq. zro) x=s(k+1)
       xx=x-xst+s(1)
       if (k .lt. kap) call xyz(xx,yy,yyp,yypp)
       if (k .ge. kap) yy=cmach
       if (k .ge. kap) yyp=zro
       s(k+1)=x
       wmn(k+1)=yy
       ttr(k+1)=one+g8*yy**2
       spr(k+1)=one/ttr(k+1)**(one+g1)
       dmdx(k+1)=yyp
       dpx(k+1)=-gam*yy*yyp*spr(k+1)/ttr(k+1)
      enddo
      s(1)=xst
      kat=kup+1
      kbl=kat+4
      ko=1
      call bound
      if (jb .eq. 7) stop
      if (jb .gt. 10) goto 1
      write (2,25) itle,slong
      write (2,31) l1,l2,dc4
      write (2,27) (s(k),rco(k),drx(k),k=1,kat)
      goto 1
24    close(1)
      close(2)
      stop
!
25    format (1x,9x,3a4,'COORDINATES AND DERIVATIVES, LENGTH=',f12.7/)
26    format (1x,8x,2f15.6,1p4e20.8)
27    format (10(9x,0P2f15.6,1pe20.8/))
28    format (8e10.3)
29    format (1x)
30    format (3a4,i3)
31    format (1x,14x,a4,'(IN)',7x,a4,'(IN)',6x,a8,12x,a8,14x,a8,9x,a8,2x
     &,a4 /)
      end program contur
