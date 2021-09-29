      module work
      use kinddefine
      implicit none
      integer(kind=K8),save :: noup,npct,nodo
      real(kind=K8),dimension(5,200),save :: wall
      real(kind=K8),dimension(200),save :: wax,way,wan
      real(kind=K8),dimension(5,150),save :: a,b,fclast
      real(kind=K8),dimension(400),save :: e,z,x,y,yst
      real(kind=K8),dimension(250),save :: wtn
      end module work

