      function fmv(pma)
!     to obtain mach number from prandtl meyer angle
      use kinddefine
      use gg, only:g2,g7,g9
!
      implicit none
      integer(kind=K4) :: i
      real(kind=K8) :: ang,fmv,one,rem,third,vm,z,zbet
      real(kind=K8),intent(in) :: pma
      one=1.d+0
      third=one/3.d+0
      vm=(dasin(one)*(pma/(g2-one))**2)**third
      z=one+.895d+0*((g7*(g2-one))**2)**third*dtan(vm)
      do i=1,100
       zbet=dsqrt(z*z-one)
       ang=g2*datan(zbet/g2)-datan(zbet)
       rem=(ang-pma)*z*(z*z+g9)/g9/zbet
       if (dabs(rem) .lt. 1.d-10) goto 2
       z=z-rem
      enddo
2     fmv=z-rem
      end function fmv
