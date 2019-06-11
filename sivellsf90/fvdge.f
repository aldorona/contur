      subroutine fvdge(x,y,ds,dy)
!
      use kinddefine
      implicit none
      real(kind=K8),dimension(5),intent(in) :: x,y
      real(kind=K8),intent(out) :: ds,dy
      real(kind=K8) :: a1,a3,f1,f2,f3,fx,h,two
      real(kind=K8) :: x1,x2,x3,x4,x5,x21,x32,x43,x54
      real(kind=K8) :: y1,y2,y3,y4,y5,yp21,yp32,yp43,yp54,ypp1,ypp2,ypp3
      real(kind=K8) :: z13
      data h/0.5d+0/,two/2.0d+0/
!
      x1=x(1)
      x2=x(2)
      x3=x(3)
      x4=x(4)
      x5=x(5)
!
      y1=y(1)
      y2=y(2)
      y3=y(3)
      y4=y(4)
      y5=y(5)
!
!     find delta-y
      f1=(x3-x1)*(x3-x2)
      f1=two/f1
!
      f2=(x4-x3)*(x3-x2)
      f2=-two/f2
!
      f3=(x5-x3)*(x4-x3)
      f3=two/f3
!
      z13=x1+x2+x2-x4-x4-x5
      a1=(x2+x3-x4-x5)/z13
      a3=(x1+x2-x3-x4)/z13
!
      yp21=(y2-y1)/(x2-x1)
      yp32=(y3-y2)/(x3-x2)
      yp43=(y4-y3)/(x4-x3)
      yp54=(y5-y4)/(x5-x4)
!
      x21=h*(x2+x1)
      x32=h*(x3+x2)
      x43=h*(x4+x3)
      x54=h*(x5+x4)
!
      ypp1=(yp32-yp21)/(x32-x21)
      ypp2=(yp43-yp32)/(x43-x32)
      ypp3=(yp54-yp43)/(x54-x43)
      ds=a1*ypp1+a3*ypp3-ypp2
      fx=f2-a1*f1-a3*f3
      dy=ds/fx
!
      end subroutine fvdge
