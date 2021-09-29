      subroutine trans(rto,tk,wo,amn,amp,ampp,w,awp,awpp,cwoppp,axn)
!     to determine throat characteristic
      use kinddefine
      use gg, only:gam,g2,g5,g6,g7,g8,ga,qt
      use contr, only:itle,ie,lr
      use troat
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
       function cubic(ea,eb,ec,ed)
        use kinddefine
        implicit none
        real(kind=K8) :: cubic
        real(kind=K8), intent(in) :: ea,eb,ec,ed
       end function cubic
      end interface
!
      integer(kind=K4) :: i,j,jj,k,kk,nn
      real(kind=K8) :: add,amn,amp,ampp,as,awop,awopp,awoppp,awp,awpp
      real(kind=K8) :: axn,bet,bs,cs,bx,cwoppp,duy,dxi,eit,four,fma,fmu
      real(kind=K8) :: fsa,fsb,fsc,fsy,fsy1,gb,gk,gt,gu,gv,gz,half,hoppp
      real(kind=K8) :: hvppp,one,p1,psi,psi1,px,rc,rto,s1,s2,s3,six,sum
      real(kind=K8) :: t,t1,thr,tk,tlv,tn1,tn2,tna,trhv,two
      real(kind=K8) :: u22,u23,u42,u43,u63,up0,up2,uy
      real(kind=K8) :: v02,v03,v22,v23,v42,v43,v63,vo,vp,vpp,vy
      real(kind=K8) :: w,wo,wop,wopp,wwo,wy,x,x1,xm,xmp,xmpp,xot,xx1,xxi
      real(kind=K8) :: xw,xwp,xwpp,y,y1,y4,y6,ym,yy,zro
      data zro/0.0d+0/,one/1.d+0/,two/2.d+0/,six/6.d+0/,half/5.d-1/
      data trhv/1.5d+0/,thr/3.d+0/,four/4.d+0/,eit/8.d+0/,tlv/1.2d+1/
      nn=iabs(lr)
      jj=240/(nn-1)
      if(mod(jj,2).ne.0) jj=jj+1
      if(jj.lt.10) jj=10
      kk=jj*nn-jj
      gb=ie/eit
      gk=(gam*(gam+2.25d+0*ie-16.5d+0)+2.25d+0*(2+ie))/tlv 
      gu=one-gam/trhv
      gv=(half*(5-3*ie)*gam+ie)/(9-ie)
      gz=dsqrt(qt*(gam+one))
      u22=gb+gam/thr/(3-ie)
      u42=(gam+(4-ie)*trhv)/six/(3-ie)
      if (ie.eq.0) goto 1
      gt=(gam*(gam*92.d+0+180.d+0)-9.d+0)/1152.d+0
      u23=(gam*(304.d+0*gam+255.d+0)-54.d+0)/1728.d+0
      u43=(gam*(388.d+0*gam+777.d+0)+153.d+0)/2304.d+0
      u63=(gam*(556.d+0*gam+1737.d+0)+3069.d+0)/10368.d+0
      up0=(gam*(52.d+0*gam+75.d+0)-9.d+0)/192.d+0 
      up2=(gam*(52.d+0*gam+51.d+0)+327.d+0)/384.d+0 
      v02=(28.d+0*gam-15.d+0)/288.d+0 
      v22=(20.d+0*gam+27.d+0)/96.d+0 
      v42=(gam/thr+one)/thr 
      v03=(gam*(7100.d+0*gam+2151.d+0)+2169.d+0)/82944.d+0 
      v23=(gam*(3424.d+0*gam+4071.d+0)-972.d+0)/13824.d+0 
      v43=(gam*(3380.d+0*gam+7551.d+0)+3771.d+0)/13824.d+0 
      v63=(gam*(6836.d+0*gam+23031.d+0)+30627.d+0)/82944.d+0 
      goto 2
1     gt=(gam*(gam*134.d+0+429.d+0)+123.d+0)/4320.d+0
      u23=(gam*(854.d+0*gam+807.d+0)+279.d+0)/12960.d+0
      u43=(gam*(194.d+0*gam+549.d+0)-63.d+0)/2592.d+0
      u63=(gam*(362.d+0*gam+1449.d+0)+3177.d+0)/12960.d+0
      up0=(gam*(26.d+0*gam+51.d+0)-27.d+0)/144.d+0 
      up2=(gam*(26.d+0*gam+27.d+0)+237.d+0)/288.d+0 
      v02=(34.d+0*gam-75.d+0)/1080.d+0 
      v22=(10.d+0*gam+15.d+0)/108.d+0 
      v42=(22.d+0*gam+75.d+0)/360.d+0 
      v03=(gam*(7570.d+0*gam+3087.d+0)+23157.d+0)/544320.d+0 
      v23=(gam*(5026.d+0*gam+7551.d+0)-4923.d+0)/77760.d+0 
      v43=(gam*(2254.d+0*gam+6153.d+0)+2979.d+0)/25920.d+0 
      v63=(gam*(6574.d+0*gam+26481.d+0)+40059.d+0)/181440.d+0
2     wwo=wo+(half+(u42-u22+(u63-u43+u23)/rto)/rto)/rto
      wop=(one-(gb-gt/rto)/rto)/dsqrt(rto)
      wopp=(gu-gv/rto)/rto
      hoppp=gk/rto/dsqrt(rto)
      hvppp=(3*ie-(10-3*ie)*gam)/four/rto/dsqrt(rto)
      amn=wwo/dsqrt(g7-g8*wwo**2)
      bet=dsqrt(amn**2-one)
      psi1=g2*datan(bet/g2)-datan(bet)
      p1=zro
      t1=zro
      x1=zro
      y1=one
      fsy1=zro
      tn2=-one/bet
      fc(1,nn)=x1
      fc(2,nn)=y1
      fc(3,nn)=amn
      fc(4,nn)=psi1
      fc(5,nn)=zro
      fc(6,nn)=zro
      bx=one
      sum=zro
      fsa=(ie+1)*amn/(g6+g5*amn**2)**ga
      do j=1,kk
       y=dfloat(kk-j)/kk
       if(ie.eq.1) bx=y+y
       yy=y*y
       tn1=tn2
       vo=(((yy*(yy*(yy*v63-v43)+v23)-v03)/rto+yy*(yy*v42-v22)+v02)/rto+
     &half*(yy-one)/(3-ie))/rto
       vp=(one+((yy*(two*gam+3*(4-ie))-two*gam-trhv*ie)/(3-ie)/thr+(yy*(
     &six*u63*yy-four*u43)+two*u23)/rto)/rto)/dsqrt(rto)
       vpp=two*(one+(two*up2*yy-up0)/rto)/rto
!     iterate for x and mach number from characteristic equations
       do i=1,10
        tna=half*(tn1+tn2)
        x=x1+(y-y1)/tna
        dxi=dsqrt((y-y1)**2+(x-x1)**2)
        xot=x/gz
        vy=gz*(vo+xot*(vp+xot*(half*vpp+xot*hvppp/thr)))/dsqrt(rto)
        w=amn/dsqrt(g6+g5*amn**2)
        t=dasin(vy*y/w)
        fsy=ie*vy/w/amn
        p1=half*(fsy1+fsy)*dxi
3       psi=p1+psi1+t1-t
        fma=fmv(psi)
        if(dabs(amn-fma).lt.1.d-10) goto 5
        fmu=dasin(one/fma)
        tn2=dtan(t-fmu)
        amn=fma
       enddo
!      iteration complete
5      if(mod(j,2).eq.0) goto 6
       as=y1-y
       fsb=bx/dsin(fmu-t)/(g6+g5*fma**2)**ga 
       goto 7
6      bs=y1-y
       cs=as+bs
       s1=(two-bs/as)*cs/six
       s3=(two-as/bs)*cs/six 
       s2=cs-s1-s3
       fsc=bx/dsin(fmu-t)/(g6+g5*fma**2)**ga
       add=s1*fsa+s2*fsb+s3*fsc
       sum=add+sum
       fsa=fsc
7      x1=x
       y1=y
       t1=t
       fsy1=fsy
       psi1=psi
       if(mod(j,jj).ne.0) goto 8
       k=nn-j/jj
       fc(1,k)=x
       fc(2,k)=y
       fc(3,k)=fma
       fc(4,k)=psi
       fc(5,k)=t
       fc(6,k)=sum
8     enddo
      do j=1,nn
       fc(1,j)=fc(1,j)/tk
       fc(2,j)=fc(2,j)/tk
       fc(6,j)=one-fc(6,j)/sum
      enddo
      axn=fc(1,1)
      awop=wop*tk/gz
      awopp=wopp*(tk/gz)**2
      awoppp=two*hoppp*(tk/gz)**3
      cwoppp=six*(w-wo-axn*(awop+axn*awopp/two))/axn**3
      if(cwoppp.lt.awoppp) cwoppp=awoppp
      awp=awop+axn*(awopp+axn*cwoppp/two)
      awpp=awopp+axn*cwoppp
      amp=awp*g7*(amn/w)**3
      ampp=amp*(awpp/awp+thr*g5*amp*w*w/amn)
      if(lr.gt.0) return
      lr=nn
      rc=rto-one
      write (2,12) itle,rc,awop,awopp,awoppp
      do j=1,nn
       y=dfloat(j-1)/(nn-1)
       yy=y*y
       y4=yy**2
       y6=yy**3
       duy=(half*yy+(u42*y4-u22*yy+(u63*y6-u43*y4+u23*yy)/rto)/rto)/rto
       uy=wo+duy
       vo=(((yy*(yy*(yy*v63-v43)+v23)-v03)/rto+yy*(yy*v42-v22)+v02)/rto+
     &half*(yy-one)/(3-ie))/rto
       vy=gz*vo*y/dsqrt(rto)
       wy=dsqrt(uy**2+vy**2)
       ym=wy/dsqrt(g7-g8*wy**2)
       write (2,13) y,uy,vy,wy,ym
       if(mod(j,10).eq.0) write (2,14)
      enddo
      xx1=cubic(cwoppp/six,awopp/two,awop,wo-one)
      xxi=cubic(awoppp/six,awopp/two,awop,wo-w)
      write (2,15) xx1,xxi,w,cwoppp,tk
      write (2,16)
      px=axn+1.d-1
      do j=1,11
       x=.1d+0*(j-1)
       xw=wo+x*(awop+x*(awopp/two+x*cwoppp/six))
       xwp=awop+x*(awopp+x*cwoppp/two)
       xwpp=awopp+x*cwoppp
       xm=xw/dsqrt(g7-g8*xw**2)
       xmp=xwp*g7*(xm/xw)**3
       xmpp=xmp*(xwpp/xwp+thr*g5*xmp*xw*xw/xm)
       if(x.lt.axn.or.x.gt.px) goto 11
       write (2,18) axn,w,awp,awpp,amn,amp,ampp
11     write (2,17) x,xw,xwp,xwpp,xm,xmp,xmpp
      enddo
      return
!
12    format (1x,8x,3a4,' THROAT VELOCITY DISTRIBUTION, X=O, RC=',f10.6/
     &/10x,'DERIVATIVES TAKEN WITH RESPECT TO X/Y*, WOP=',f11.8//10x,'WO
     &PP=',1pe15.7,5x,'WOPPP=',e15.7//10x,'Y/YO',7x,'U/A*',10x,'V/A*',11
     &x,'W',11x,'MACH NO.'/)
13    format (1x,f14.4,4f14.8)
14    format (1x)
15    format (1x,9x,'FROM CUBIC, X/Y* =',f11.8,' FOR W= 1.0'//22x,'X/Y* 
     &=',f11.8,' FOR W=',f11.8 //10x,'CORRECTED WOPPP=',1pe15.7//10x,'RM
     &ASS = Y*/YO =',0pf13.10//)
16    format (1x,9x,'AXIAL VELOCITY DISTRIBUTION, Y=0'//10x,'X/Y*',9x,'W
     &',17x,'WP',16x,'WPP',15x,'M',17x,'MP',16x,'MPP'/)
17    format (1x,f13.3,1p6e18.7)
18    format (1x,f16.8,1pe15.7,5e18.7)
      end subroutine trans
