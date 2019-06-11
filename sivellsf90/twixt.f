      subroutine twixt(s,gma,gmb,gmc,gmd,xbl,kat,kbl)
!     to determine interpolation coefficients
      use kinddefine
      implicit none
      integer(kind=K4) :: j,l
      integer(kind=K4),intent(in) :: kat
      integer(kind=K4),intent(out) :: kbl
      real(kind=K8) :: ds,dst,dstu,dt,dtu,du,xbb
      real(kind=K8),intent(out) :: gma,gmb,gmc,gmd
      real(kind=K8),intent(in) :: xbl
      real(kind=K8),dimension(200),intent(in) :: s
      do l=1,kat
       if(s(kat-l).lt.xbl) goto 2
      enddo
2     j=kat-l+1
      xbb=s(j)-xbl
      kbl=j+1
      du=s(j+1)-s(j)
      dt=s(j)-s(j-1)
      ds=s(j-1)-s(j-2) 
      dst=ds+dt 
      dstu=dst+du 
      dtu=dt+du 
      gma=-xbb*(dt-xbb)*(du+xbb)/ds/dst/dstu 
      gmb=xbb*(dst-xbb)*(du+xbb)/ds/dt/dtu 
      gmc=(dst-xbb)*(dt-xbb)*(du+xbb)/dst/dt/du 
      gmd=-xbb*(dst-xbb)*(dt-xbb)/dstu/dtu/du 
      end subroutine twixt
