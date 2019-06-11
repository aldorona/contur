      subroutine xyz(xx,yy,yyp,yypp)
!     compute y,y',y'' for a curve described by cubic's e(5,*)
!     where (1) = x-max  (2) = high order coefficient.
      use kinddefine
      use coef, only:e
      implicit none
      integer(kind=K4) :: k
      real(kind=K8),intent(in) :: xx
      real(kind=K8),intent(out) :: yy,yyp,yypp
      real(kind=K8) :: a1,a2,a3,az,t,s,r,x,y,yp,ypp,zero
      data zero/0.0d+0/
      x=xx
      if(x.ge.e(1,1)) goto 2
1     y=zero
      yp=zero
      ypp=zero
      goto 5
2     do k=2,200
       if(x.le.e(1,k)) goto 4
      enddo
      goto 1
4     a3=e(2,k)
      a2=e(3,k)
      a1=e(4,k)
      az=e(5,k)
      t=a2+a2
      s=a3*3.0d+0
      r=s+s
      y=az+x*(a1+x*(a2+x*a3))
      yp=a1+x*(t+x*s)
      ypp=t+r*x
5     yy=y
      yyp=yp
      yypp=ypp
      end subroutine xyz
