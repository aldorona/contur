      subroutine sorce(w,b)
!     to obtain velocity derivatives in radial flow
      use kinddefine
      use gg, only:g1,g7,g9,qt
      implicit none
      real(kind=K8) :: al,area,axw,aww,c2,c4,cww,dww,four,one,thr,two,ww
      real(kind=K8) :: ww1
      real(kind=K8),intent(in) :: w
      real(kind=K8),dimension(4),intent(out) :: b
      data one/1.d+0/,two/2.d+0/,thr/3.d+0/,four/4.d+0/
      ww=w*w
      al=g7*g9
      aww=al-ww
      ww1=ww-one
      area=(((al-one)/aww)**g1)/w
      b(1)=area**qt
      axw=al*ww1*b(1)
      b(2)=w*aww/axw/qt
      c2=thr/qt+al*(two-one/qt)
      c4=al+one/qt
      cww=ww*(c2-ww*c4)-al*(one+one/qt)
      b(3)=b(2)*cww/axw/ww1
      dww=(two*c2-four*c4*ww)/cww-four/ww1
      b(4)=b(3)*(b(3)/b(2)+w*b(2)*dww-one/b(1))
      end subroutine sorce
