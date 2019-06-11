      subroutine neo
!
!     smooth by linear second derivative
!
      use kinddefine
      use work, only: wall,wax,way,wan,nodo,noup,npct,e,z,x,y,yst,wtn
      use contr, only: itle,jq,iq,n,np,nut
!
      implicit none
!
      interface
       subroutine fvdge(x,y,ds,dy)
        use kinddefine
        implicit none
        real(kind=K8),dimension(5),intent(in) :: x,y
        real(kind=K8),intent(out) :: ds,dy
       end subroutine fvdge
!
       subroutine scond(a,b,c,king)
        use kinddefine
        implicit none
        integer(kind=K4),intent(in) :: king
        real(kind=K8),dimension(150),intent(in) :: a,b
        real(kind=K8),dimension(150),intent(out) :: c
       end subroutine scond
      end interface
!
      integer(kind=K4) :: lim,lm,lu,lus,j,k,m,maxej
      integer(kind=K4) :: notm
      real(kind=K8) :: conv,emax,error,one,smp,test,tni,two,zero
      character(len=4,kind=K3) :: j0,j1,jn
!      COMMON /WORK/ WALL(5,200),WAX(200),WAY(200),WAN(200),E(400),Z(400)
!     1,X(400),Y(400),YST(400),WTN(250)
!      COMMON /CONTR/ ITLE(3),IE,LR,IT,JB,JQ,JX,KAT,KBL,KING,KO,LV,NOCON,
!     1IN,MC,MCP,IP,IQ,ISE,JC,M,MP,MQ,N,NP,NR,NUT,NF
      data zero/0.0d+0/,one/1.d+0/,two/2.d+0/
      data j0/'  UP'/,j1/'DOWN'/
!
      conv=90.d+0/dasin(one)
      tni=dtan(wall(5,1))
!
      if (jq.eq.0.or.iq.lt.0) read (1,14,end=13) noup,npct,nodo
      if (jq .gt. 0) goto 2
      jn=j0
      lim=nut
      notm=noup
      do j=1,lim
       x(j+1)=wax(j)
       y(j+1)=way(j)
       yst(j+1)=y(j+1)
      enddo
      x(1)=two*x(2)-x(3)
      y(1)=y(3)
      x(lim+2)=two*x(lim+1)-x(lim)
      y(lim+2)=y(lim+1)+tni*(x(lim+2)-x(lim+1))
      goto 4
2     lim=n+np-1
      notm=nodo
      jn=j1
      do j=1,lim
       x(j+1)=wall(1,j)
       y(j+1)=wall(2,j)
       yst(j+1)=y(j+1)
      enddo
      x(1)=two*x(2)-x(3)
      y(1)=y(2)-tni*(x(2)-x(1))
      x(lim+2)=two*x(lim+1)-x(lim)
      y(lim+2)=y(lim+1)
4     lus=1+(lim-3)/6
      if (notm .eq. 0) return
      yst(1)=y(1)
      yst(lim+2)=y(lim+2)
      smp=1.d-2*npct
      write (2,16) itle,jn,notm,smp
!     I HAVE NOT UNDERSTOOD WHY m IS KEPT IN THE COMMON BLOCK contr
      do m=1,notm
!       call orez(e,800)
        e(:)=0.0d0
        z(:)=0.0d0
!
       do k=3,lim
        call fvdge (x(k-2:k+2),y(k-2:k+2),e(k),z(k))
       enddo
       e(1)=zero
       e(2)=zero
       e(lim+1)=zero
       e(lim+2)=zero
!     search array and find max error
       do lu=1,lus
        emax=zero
        do k=3,lim
         test=dabs(e(k))
         if (emax .gt. test) goto 6
         j=k
         emax=test
6        continue
        enddo
!     apply correction
        e(j)=zero
        e(j+1)=zero
        e(j+2)=zero
        e(j-1)=zero
        e(j-2)=zero
        y(j)=y(j)+smp*z(j)
       enddo
      enddo
!
      error=zero
      do j=1,lim
       k=j+1
       e(k)=y(k)-yst(k)
       if (error .lt. dabs(e(k))) maxej=j 
       if (error .lt. dabs(e(k))) error=dabs(e(k))
       write (2,15) j,x(k),y(k),yst(k),e(k),j
       if (mod(j,10) .eq. 0) write (2,17)
      enddo
      write (2,19) error,maxej
!
      lm=lim-1
      call scond (x,y,wtn,lim+2)
      if (jq .eq. 1) goto 11
      do j=2,lm
       way(j)=y(j+1)
       wan(j)=conv*datan(wtn(j+1))
      enddo
      return
!
11    do j=2,lm
       wall(2,j)=y(j+1)
       wall(5,j)=datan(wtn(j+1))
      enddo
      return
!
13    write (2,18)
      stop
!
14    format (3i5)
15    format (1x,20x,i5,2x,0p4f13.7,i8)
16    format (1x,3a4,2x,a4,'STREAM CONTOUR, SMOOTHED',i5,' TIMES WITH FA
     &CTOR=',f4.2//35x,'X',11x,'Y-CALC',7x,'Y-IN',10x,'DIFF' /)
17    format (1x)
18    format ('0',10x,'CARD NOT AVAILABLE FOR NEGATIVE NF')
19    format (1x,26x,'MAX. ABSOLUTE ERROR =',1pg15.6,' AT POINT ',i5/)
      end subroutine neo
