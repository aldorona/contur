      subroutine scond(a,b,c,king)
!     to obtain parabolic derivative of curve (unequally spaced points)
      use kinddefine
      implicit none
      integer(kind=K4),intent(in) :: king
      integer(kind=K4) :: k,n
      real(kind=K8),dimension(150),intent(in) :: a,b
      real(kind=K8),dimension(150),intent(out) :: c
      real(kind=K8) :: qf,qo,qst,s,sf,so,t,t0,tf
      n=king-1
      do k=2,n
       s=a(k)-a(k-1)
       t=a(k+1)-a(k)
       c(k)=((b(k+1)-b(k))*s*s+(b(k)-b(k-1))*t*t)/(s*s*t+s*t*t)
      enddo
      so=a(2)-a(1)
      t0=a(3)-a(2)
      qo=so+t0
      c(1)=(-t0*(qo+so)*b(1)+qo*qo*b(2)-so*so*b(3))/qo/so/t0
      sf=a(king-1)-a(king-2)
      tf=a(king)-a(king-1)
      qf=sf+tf
      qst=qf*sf*tf
      c(king)=(sf*(qf+tf)*b(king)-qf*qf*b(king-1)+tf*tf*b(king-2))/qst
      end subroutine scond
