      subroutine orez(a,na)
      use kinddefine
      implicit none
      integer(kind=K4),intent(in) :: na
      integer(kind=K4) :: k
      real(kind=K8),dimension(1),intent(out) :: a
      do k=1,na
       a(k)=0.0d+0
      enddo
      end subroutine orez
