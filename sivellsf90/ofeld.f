      subroutine ofeld(a,b,c,nocon)
!     to obtain points in characteristic network
      use kinddefine
      use contr, only:ie
      implicit none
!
      interface
       function fmv(psi)
        use kinddefine
        implicit none
        real(kind=K8) :: fmv
        real(kind=K8), intent(in) :: psi
       end function fmv
      end interface
!
      integer(kind=K4) :: i
      integer(kind=K4),intent(out) :: nocon
      real(kind=K8),dimension(5),intent(in) :: a,b
      real(kind=K8),dimension(5),intent(out) :: c
      real(kind=K8) :: a1,a2,a3,dtn,fm3,fsy1,fsy2,fsy3,half,hdpsi,hpsi3
      real(kind=K8) :: ht3,one,p1,p2,psi3
      real(kind=K8) :: t1,t2,t3,temp,tn2,tna,tnb,tni,told,two,x3,y3,zro
      data zro/0.0d+0/,one/1.d+0/,two/2.d+0/,half/5.d-1/
      a1=dasin(one/a(3))
      a2=dasin(one/b(3))
      t1=a(5)
      t2=b(5)
      if (ie .eq. 0) goto 8
      if (a(2) .eq. zro) goto 5
      fsy1=dsin(a(5))/a(2)/a(3)
      goto 6
5     t1=zro
      fsy1=a(5)
6     if (b(2) .eq. zro) goto 7
      fsy2=dsin(b(5))/b(2)/b(3)
      goto 8
7     t2=zro
      fsy2=b(5)
8     tni=dtan(t1-a1)
      if (b(3) .ne. one) tn2=dtan(t2+a2)
      i=-1
      hdpsi=half*(a(4)-b(4))
      ht3=half*(t1+t2)+hdpsi
      t3=ht3-half*ie*hdpsi
      hpsi3=half*(a(4)+b(4)+t1-t2)
      psi3=hpsi3+half*ie*(t1-t2)
      c(3)=fmv(psi3)
      told=t3
1     i=i+1
      fm3=c(3)
      a3=dasin(one/c(3))
      tna=half*(tni+dtan(t3-a3))
      if (b(3) .ne. one) tnb=half*(dtan(t3+a3)+tn2)
      if (b(3) .eq. one) tnb=two*dtan(t3+a3)
      dtn=tnb-tna
      x3=(b(1)*tnb-a(1)*tna+a(2)-b(2))/dtn
      y3=(a(2)*tnb-b(2)*tna+(b(1)-a(1))*tna*tnb)/dtn
      if ((ie .eq. 0) .or. (dabs(y3) .lt. 1.d-9)) goto 4
      fsy3=dsin(t3)/y3/fm3
      p1=half*(fsy1+fsy3)*(x3-a(1))*dsqrt(one+tna**2)
      p2=half*(fsy2+fsy3)*(x3-b(1))*dsqrt(one+tnb**2)
      t3=ht3+half*(p1-p2)
      psi3=hpsi3+half*(p1+p2)
      c(3)=fmv(psi3)
      if (dabs(t3-told) .gt. 1.d-9) goto 2
      if (dabs(c(3)-fm3) .lt. 1.d-9) goto 4
2     if (i .eq. 40) goto 3
      temp=t3
      t3=(t3+told)*half
      told=temp
      goto 1
3     nocon=1
4     c(1)=x3
      c(2)=y3
      c(4)=psi3
      c(5)=t3
      end subroutine ofeld
