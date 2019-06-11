      function toric(wip,se)
!     to obtain throat radius of curvature from velocity gradient
      use kinddefine
      use fg, only:gc,gd,ge,gf
      use gg, only:gam,g7,qt
      implicit none
      real(kind=K8) :: ie,ff,fiv,fp,fw,one,thr,tk,toric,tr2,trr
      real(kind=K8),intent(in) :: se,wip
      data one/1.d+0/,thr/3.d+0/,fiv/5.d+0/
      ie=one/qt-one
      fw=wip*se*dsqrt(qt*(gam+one))
      trr=fw*(one+(gc+(thr*gc**2-gd)*fw**2)*fw**2)
1     tr2=trr**2
      tk=(one-g7*(one+(ge+gf*tr2)*tr2)*tr2**2/(45.d+0+3*ie))**qt
      ff=fw/tk-trr*(one-tr2*(gc-gd*tr2))
      fp=one-tr2*(thr*gc-fiv*gd*tr2)
      trr=trr+ff/fp
      if(dabs(ff).gt.1.d-1) go to 1
      toric=one/trr**2
      end function toric
