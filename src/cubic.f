      function cubic(ea,eb,ec,ed)
      use kinddefine
!
      implicit none
      integer(kind=K4) :: ia,ib,ic
      real(kind=K8),intent(in) :: ea,eb,ec,ed
      real(kind=K8) :: a,aa,b,cba,cbb,cbc,cbt1,cbt2,csa,csna,cubic,e3
      real(kind=K8) :: one,q,q1,qr,r1,rq,thr,two,zro
!     to obtain positive real root of cubic equation
      data zro/0.d+0/,one/1.d+0/,two/2.d+0/,thr/3.d+0/
      e3=eb/thr
      q1=ea*ec/thr-e3**2
      r1=ea*(e3*ec-ea*ed)/two-e3**3
      qr=q1**3+r1**2
      rq=dsqrt(dabs(qr))
      q=dsqrt(dabs(q1))
      b=dsign(one,r1)
      cbb=-one
      cbc=-one
      cbt1=zro
      cbt2=zro
      a=zro
      if (qr .gt. zro) goto 1
      if (qr .ne. zro) a=dasin(-rq/q1/q)/thr
      csa=dcos(a)
      csna=dsqrt(thr)*dsin(a)
      cba=(two*b*q*csa-e3)/ea
      cbb=-(b*q*(csa+csna)+e3)/ea
      cbc=-(b*q*(csa-csna)+e3)/ea
      goto 2
1     if (r1+rq .ne. zro) cbt1=dsign(dexp(dlog(dabs(r1+rq))/thr),r1+rq)
      if (r1-rq .ne. zro) cbt2=dsign(dexp(dlog(dabs(r1-rq))/thr),r1-rq)
      cba=(cbt1+cbt2-e3)/ea
2     ia=dsign(one,cba)
      ib=dsign(one,cbb)
      ic=dsign(one,cbc)
      if (ia+ib+ic+1) 11,3,7
3     if (ia .eq. 1) goto 5
      if (ib .eq. 1) goto 6
4     cubic=cbc
      return
5     cubic=cba
      return
6     cubic=cbb
      return
7     if (ia+2*ib+3*ic-2) 8,9,10
8     if (cba .gt. cbb) goto 6
      goto 5
9     if (cba .gt. cbc) goto 4
      goto 5
10    if (cbb .gt. cbc) goto 4
      goto 6
11    aa=a*9.d+1/dasin(one)
      write (2,12) ea,eb,ec,ed,q1,r1,qr,rq,q,aa,cba,cbb,cbc
      cubic=-one
      return
!
12    format ('O','EA=',e14.7,'  EB=',e14.7,'  EC=',e14.7,'  ED=',e14.7,
     &'  Q1=',e14.7,'  R1=',e14.7,'  QR=',e14.7,'  RQ=',e14.7,'   Q=',e1
     &4.7,', AA=',e14.7,',CBA=',e14.7,',CBB=',e14.7,',CBC=',e14.7 /)
      end function cubic
