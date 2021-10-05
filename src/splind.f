      subroutine splind(x,y,tn2,tnl,l)
!     compute cubic coefficients for a curve x-y
      use kinddefine
      use coef, only:e,ne
      implicit none
      integer(kind=K4),intent(in) :: l
      integer(kind=K4) :: j,k,n
      real(kind=K8),dimension(300) :: a,b,c,d,dx,dy,g,sb,xm
      real(kind=K8),dimension(l),intent(in) :: x,y
      real(kind=K8),intent(in) :: tn2,tnl
      real(kind=K8) :: dxr,one,p,px,pxx,pxxx,q,qx,qxx,qxxx,r,s,six,sw
      real(kind=K8):: thr,xj,xk,zero
      data zero/0.0d+0/,one/1.d+0/,thr/3.d+0/,six/6.d+0/
!      call orez (e,5*200)
!      call orez (a,9*300)
      e(:,:)=0.0d0
      a(:)=0.0d0
      b(:)=0.0d0
      c(:)=0.0d0
      d(:)=0.0d0
      dx(:)=0.0d0
      dy(:)=0.0d0
      g(:)=0.0d0
      sb(:)=0.0d0
      xm(:)=0.0d0
!
      dx(1)=zero
      dy(1)=zero
      n=l-1
      do k=2,l
       dx(k)=x(k)-x(k-1)
       dy(k)=y(k)-y(k-1)
      enddo
!
      b(1)=dx(2)/thr
      c(1)=dx(2)/six
      d(1)=dy(2)/dx(2)-tn2
      a(l)=dx(l)/six
      b(l)=dx(l)/thr
      d(l)=tnl-dy(l)/dx(l)
      a(1)=zero
      do k=2,n
       a(k)=dx(k)/six
       b(k)=(dx(k)+dx(k+1))/thr
       d(k)=dy(k+1)/dx(k+1)-dy(k)/dx(k)
       c(k)=dx(k+1)/six
      enddo
      sw=one/b(1)
      sb(1)=sw*c(1)
      g(1)=sw*d(1)
      do k=2,l
       sw=one/(b(k)-a(k)*sb(k-1))
       sb(k)=sw*c(k)
       g(k)=sw*(d(k)-a(k)*g(k-1))
      enddo
      xm(l)=g(l)
      do k=1,n
       j=l-k
       xm(j)=g(j)-sb(j)*xm(j+1)
      enddo
      do k=2,l
       dxr=one/dx(k)
       q=dxr/six
       p=-xm(k-1)*q
       q=q*xm(k)
       r=dx(k)*xm(k-1)/six-dxr*y(k-1)
       s=y(k)*dxr-dx(k)*xm(k)/six
       xk=x(k)
       px=xk*p
       pxx=px*xk
       pxxx=pxx*xk
       xj=x(k-1)
       qx=xj*q
       qxx=qx*xj
       qxxx=qxx*xj
       e(2,k)=p+q
       e(3,k)=-thr*(px+qx)
       e(4,k)=thr*(pxx+qxx)+r+s
       e(5,k)=-pxxx-qxxx-r*xk-s*xj
      enddo
      do k=2,l
       e(1,k)=x(k)
      enddo
      e(1,1)=x(1)
      ne=l
      end subroutine splind
