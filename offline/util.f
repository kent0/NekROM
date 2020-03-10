c-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)

      character*1 a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine blank(a,n)

      character*1 a(1)
      character*1 blnk
      save        blnk
      data        blnk /' '/

      do i=1,n
         a(i)=blnk
		enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cfill(a,b,n)
      dimension  a(1)
c
      do 100 i = 1, n
 100     a(i) = b
      return
      end
c-----------------------------------------------------------------------
      subroutine chsign(a,n)
      real a(1)
c
      do 100 i=1,n
         a(i) = -a(i)
 100  continue
      return
      end
c
c-----------------------------------------------------------------------
      function glsum (x,n)
      dimension x(1)
      dimension tmp(1),work(1)
      tsum = 0.
      do 100 i=1,n
         tsum = tsum+x(i)
 100  continue
      tmp(1)=tsum
      CALL GOP(TMP,WORK,'+  ',1)
      glsum = tmp(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine i8copy(a,b,n)
      integer*8 a(1), b(1)
c
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      integer a(1), b(1)
c
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      function iglmin(a,n)
      integer a(1),tmin
      integer tmp(1),work(1)
      tmin=  999999999
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call igop(tmp,work,'m  ',1)
      iglmin=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      function iglmax(a,n)
      integer a(1),tmax
      integer tmp(1),work(1)
      tmax= -999999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call igop(tmp,work,'M  ',1)
      iglmax=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      function iglsum(a,n)
      integer a(1),tsum
      integer tmp(1),work(1)
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call igop(tmp,work,'+  ',1)
      iglsum=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      integer function indx1(s1,s2,l2)
      character*132 s1,s2
c
      n1=132-l2+1
      indx1=0
      if (n1.lt.1) return
c
      do 100 i=1,n1
         i2=i+l2-1
         if (s1(i:i2).eq.s2(1:l2)) then
            indx1=i
            return
         endif
  100 continue
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mxm(a,n1,b,n2,c,n3)

      include 'TIMES'

c     Compute matrix-matrix product C = A*B
c     for contiguously packed matrices A,B, and C.

      real a(n1,n2),b(n2,n3),c(n1,n3)

      tt=dnekclock()
      call mxmf2(a,n1,b,n2,c,n3)
      time_mxm=time_mxm+(dnekclock()-tt)

      return
      end
c-----------------------------------------------------------------------
      subroutine opzero(ux,uy,uz)

      real ux(1),uy(1),uz(1)
c

      call exitti('called unsupported opzero$',1)
c     n = lx1*ly1*lz1*nelfld(ifield)
      call rzero(ux,n)
      call rzero(uy,n)
c     if (if3d) call rzero(uz,n)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine izero(a,n)
      integer a(1)
      do i=1,n
         a(i)=0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rzero(a,n)
      real a(1)
      do i=1,n
         a(i)=0.0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rescale_x(x,x0,x1)

      real x(1)

      write (6,*) 'called rescale_x, unsupported'
      return
c     n = lx1*ly1*lz1*nelt
      xmin = glmin(x,n)
      xmax = glmax(x,n)

      if (xmax.le.xmin) return

      scale = (x1-x0)/(xmax-xmin)
      do i=1,n
         x(i) = x0 + scale*(x(i)-xmin)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function indx2(s1,l1,s2,l2)
      character*132 s1,s2

      n1=l1-l2+1
      indx2=0
      if (n1.lt.1) return

      do i=1,n1
         i2=i+l2-1
         if (s1(i:i2).eq.s2(1:l2)) then
            indx2=i
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function log2(k)
      rk=(k)
      rlog=log10(rk)
      rlog2=log10(2.0)
      rlog=rlog/rlog2+0.5
      log2=int(rlog)
      return
      end
c-----------------------------------------------------------------------
      subroutine ione(a,n)
      integer  a(1)
      do 100 i = 1, n
 100     a(i ) = 1
      return
      end
c-----------------------------------------------------------------------
      subroutine isort(a,ind,n)
c
c     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
c
      integer a(1),ind(1)
      integer aa
c
      do 10 j=1,n
         ind(j)=j
   10 continue
c
      if (n.le.1) return
      l=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         goto 200
         endif
           a(i) = aa
         ind(i) = ii
      goto 100
      end
c-----------------------------------------------------------------------
      subroutine sort(a,ind,n)
c
c     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
c
      real a(1),aa
      integer ind(1)
c
      do 10 j=1,n
         ind(j)=j
   10 continue
c
      if (n.le.1) return
      l=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         goto 200
         endif
           a(i) = aa
         ind(i) = ii
      goto 100
      end
c-----------------------------------------------------------------------
      function ltruncr(string,l)

      character*1 string(l)
      character*1   blnk
      data blnk/' '/

      do i=1,l
         l1=i-1
         if (string(i).eq.blnk) goto 200
      enddo
      l1=0

  200 continue
      ltruncr=l1

      return
      end
c-----------------------------------------------------------------
      subroutine col2(a,b,n)

      real a(1),b(1)

      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine invers2(a,b,n)

      real a(1),b(1)

      do i=1,n
         a(i)=1./b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine vdot2 (dot,u1,u2,v1,v2,n)
C     Compute a Cartesian vector dot product. 2-d version
c
c
      dimension dot(1)
      dimension u1(1),u2(1)
      dimension v1(1),v2(1)
c
c
      do 100 i=1,n
         dot(i) = u1(i)*v1(i) + u2(i)*v2(i)
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine vdot3 (dot,u1,u2,u3,v1,v2,v3,n)
C     Compute a Cartesian vector dot product. 3-d version
c
c
      dimension dot(1)
      dimension u1(1),u2(1),u3(1)
      dimension v1(1),v2(1),v3(1)
c
c
      do 100 i=1,n
         dot(i) = u1(i)*v1(i) + u2(i)*v2(i) + u3(i)*v3(i)
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine vcross (u1,u2,u3,v1,v2,v3,w1,w2,w3,n)
C     Compute a Cartesian vector cross product.
c
c
      dimension u1(1),u2(1),u3(1)
      dimension v1(1),v2(1),v3(1)
      dimension w1(1),w2(1),w3(1)
c
c
      do 100 i=1,n
         u1(i) = v2(i)*w3(i) - v3(i)*w2(i)
         u2(i) = v3(i)*w1(i) - v1(i)*w3(i)
         u3(i) = v1(i)*w2(i) - v2(i)*w1(i)
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine rzero3 (a,b,c,n)
      dimension a(1),b(1),c(1)
      do 100 i=1,n
         a(i)=0.0
         b(i)=0.0
         c(i)=0.0
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine unitvec (x,y,z,n)
      dimension x(1),y(1),z(1)
      do 100 i=1,n
      xlngth = sqrt( x(i)**2 + y(i)**2 + z(i)**2 )
      if (xlngth.ne.0.0) then
         x(i) = x(i)/xlngth
         y(i) = y(i)/xlngth
         z(i) = z(i)/xlngth
      endif
  100 continue
      return
      end
c-----------------------------------------------------------------------
      function vlsc3(x,y,b,n)

c     local inner product, with weight

      real x(1),y(1),b(1)

      vlsc3=0.

      do i=1,n
         vlsc3=vlsc3+x(i)*y(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rone(a,n)
      dimension a(1)
      do 100 i = 1, n
 100     a(i ) = 1.0
      return
      end
c-----------------------------------------------------------------------
      subroutine addcol3(a,b,c,n)

      real a(1),b(1),c(1)

      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine subcol3(a,b,c,n)
      real a(1),b(1),c(1)
c
      do 100 i=1,n
         a(i)=a(i)-b(i)*c(i)
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine col3(a,b,c,n)

      real a(1),b(1),c(1)

      do i=1,n
         a(i)=b(i)*c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      real function vlsc2(x,y,n)

      real x(1),y(1)

      vlsc2=0.
      do i=1,n
         vlsc2=vlsc2+x(i)*y(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine add2(a,b,n)
      real a(1),b(1)
      do i=1,n
         a(i)=a(i)+b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3(ur,us,ut,u,n,e,d,dt)
c     Output: ur,us,ut         Input:u,n,e,d,dt
      real ur(0:n,0:n,0:n),us(0:n,0:n,0:n),ut(0:n,0:n,0:n)
      real u (0:n,0:n,0:n,1)
      real d (0:n,0:n),dt(0:n,0:n)
      integer e
c
      m1 = n+1
      m2 = m1*m1
c
      call mxm(d ,m1,u(0,0,0,e),m1,ur,m2)
      do k=0,n
         call mxm(u(0,0,k,e),m1,dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u(0,0,0,e),m2,dt,m1,ut,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad2(ur,us,u,n,e,d,dt)
c     Output: ur,us         Input:u,n,e,d,dt
      real ur(0:n,0:n),us(0:n,0:n)
      real u (0:n,0:n,1)
      real d (0:n,0:n),dt(0:n,0:n)
      integer e
c
      m1 = n+1
c
      call mxm(d ,m1,u(0,0,e),m1,ur,m1)
      call mxm(u(0,0,e),m1,dt,m1,us,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine transpose(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
c
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cadd(a,const,n)
      real a(1)
c
      include 'OPCTR'
c
      do 100 i=1,n
         a(i)=a(i)+const
 100  continue
      return
      end
c-----------------------------------------------------------------------
      subroutine cmult(a,const,n)
      real a(1)
c
      include 'OPCTR'
c
      do 100 i=1,n
         a(i)=a(i)*const
 100  continue
      return
      end
c-----------------------------------------------------------------------
      real function vlsum(vec,n)
      real vec(1)
      include 'OPCTR'
c
      sum = 0.
c
      do 100 i=1,n
         sum=sum+vec(i)
 100  continue
      vlsum = sum
      return
      end
c-----------------------------------------------------------------------
      subroutine sub2(a,b,n)
      real a(1),b(1)
c
      include 'OPCTR'
c
      do 100 i=1,n
         a(i)=a(i)-b(i)
 100  continue
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat2(a,m,n,k,name)
      include 'SIZE'
      real a(m,n)
      character*4 name
c
      n2 = min(n,8)
      write(6,2) nid,name,m,n,k
      do i=1,m
         write(6,1) nid,name,(a(i,j),j=1,n2)
      enddo
c   1 format(i3,1x,a4,16f6.2)
    1 format(i3,1x,a4,1p8e14.5)
    2 format(/,'Matrix: ',i3,1x,a4,3i8)
      return
      end
c-----------------------------------------------------------------------
      subroutine add2s2(a,b,c1,n)

      real a(1),b(1)
c
      do 100 i=1,n
        a(i)=a(i)+c1*b(i)
  100 continue

      return
      end
c
c-----------------------------------------------------------------------
      subroutine add3(a,b,c,n)
      real a(1),b(1),c(1)

      do i=1,n
         a(i)=b(i)+c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine addcol4(a,b,c,d,n)
      real a(1),b(1),c(1),d(1)

      do 100 i=1,n
         a(i)=a(i)+b(i)*c(i)*d(i)
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine ascol5 (a,b,c,d,e,n)
      real a(1),b(1),c(1),d(1),e(1)

      do 100 i=1,n
         a(i) = b(i)*c(i)-d(i)*e(i)
 100  continue
      return
      end
c-----------------------------------------------------------------------
      subroutine subcol4(a,b,c,d,n)
      real a(1),b(1),c(1),d(1)
c
      do 100 i=1,n
         a(i)=a(i)-b(i)*c(i)*d(i)
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine genevec(vec,val,gram,nsg,ifld)

      include 'LVAR'

      common /scrgvec/ gc(lsg*lsg),wk(lsg*lsg)

      real gram(nsg,nsg),vec(nsg,nsg),val(nsg)

      if (nio.eq.0) write (6,*) 'inside genevec'

      call copy(gc,gram,nsg*nsg)

      call regularev(gc,val,nsg,wk)

      do l=1,nsg
         call copy(vec(1,l),gc(1+(nsg-l)*nsg),nsg)
      enddo

      do l=1,nsg/2
         tmp=val(l)
         val(l)=val(nsg-l+1)
         val(nsg-l+1)=tmp
      enddo

      do i=1,nsg
         if (nio.eq.0) write (6,'(i5,1p1e16.6,3x,a,i1)')
     $      i,val(i),'eval',ifld
      enddo

      if (nio.eq.0) write (6,*) 'exiting genevec'

      return
      end
c-----------------------------------------------------------------------
      subroutine regularev(a,lam,n,wk)
c
c     Solve the eigenvalue problem  A x = lam x
c
c     A -- symmetric matrix
c
c     "SIZE" is included here only to deduce WDSIZE, the working
c     precision, in bytes, so as to know whether dsygv or ssygv
c     should be called.

      include 'LVAR'
      include 'TIMES'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real a(n,n),lam(n),wk(n,n)
      real aa(100)

      tt=dnekclock()
      call dsyev('V','U',n,a,n,lam,wk,n*n,info)
      time_eig=time_eig+(dnekclock()-tt)

      if (info.ne.0) then
         if (mid.eq.0) then
            call outmat2(a  ,n,n,n,'Aeig')
            call outmat2(lam,1,n,n,'Deig')
         endif

         ninf = n-info
         write(6,*) 'Error in regularev, info=',info,n,ninf
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine settranspose(at,a,m,n)

      real a(m,n),at(n,m)

      do j=1,n
         do i=1,m
            at(j,i)=a(i,j)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine opcolv(a1,a2,a3,c,nel)

      include 'LVAR'

      real a1(1),a2(1),a3(1),c(1)

      ntot1=lx1*ly1*lz1*nel

      if (ldim.eq.3) then
         do 100 i=1,ntot1
            a1(i)=a1(i)*c(i)
            a2(i)=a2(i)*c(i)
            a3(i)=a3(i)*c(i)
  100    continue
      else
         do 200 i=1,ntot1
            a1(i)=a1(i)*c(i)
            a2(i)=a2(i)*c(i)
  200    continue
      endif
      return
      end
c-----------------------------------------------------------------------