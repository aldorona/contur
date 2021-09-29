      subroutine heat
!     dummy to be modified for special calculations of heat transfer
      use httr, only:qfun,qfunw
      implicit none
      qfunw=qfun
      end subroutine heat
