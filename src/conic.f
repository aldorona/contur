      subroutine conic(xm,b) 
!     to obtain mach number derivatives in radial flow
      use kinddefine
      use gg, only:g5,g6,g8,ga,qt
      implicit none
      real(kind=K8),dimension(4),intent(out) :: b
      real(kind=K8),intent(in) :: xm
      real(kind=K8) :: area,bmm,c2,c4,cmm,dmm,four,one,thr,two
      real(kind=K8) :: xmm,xmm1,xmm2
      data one/1.d+0/,two/2.d+0/,thr/3.d+0/,four/4.d+0/
      xmm=xm*xm
      xmm1=xmm-one
      xmm2=xmm1**2
      bmm=one+g8*xmm
      area=(g6+g5*xmm)**ga/xm
      b(1)=area**qt
      b(2)=xm*bmm/qt/xmm1/b(1)
      c2=two-(one+thr*g8)/qt
      c4=g8/qt-one
      cmm=xmm*(c2+xmm*c4)-one-one/qt
      b(3)=b(2)*cmm/xmm2/b(1)
      dmm=(four*c4*xmm+two*c2)/cmm-four/xmm1
      b(4)=b(3)*(b(3)/b(2)+xm*b(2)*dmm-one/b(1))
      end subroutine conic
