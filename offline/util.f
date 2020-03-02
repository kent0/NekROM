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
      DIMENSION  A(1)
C
      DO 100 I = 1, N
 100     A(I) = B
      return
      END
c-----------------------------------------------------------------------
      subroutine chsign(a,n)
      REAL A(1)
C
      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
      return
      END
C
c-----------------------------------------------------------------------
      function glsum (x,n)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      TSUM = 0.
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL GOP(TMP,WORK,'+  ',1)
      GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      subroutine i8copy(a,b,n)
      INTEGER*8 A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
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
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*132 S1,S2
C
      N1=132-L2+1
      INDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine mxm(a,n1,b,n2,c,n3)
c
c     Compute matrix-matrix product C = A*B
c     for contiguously packed matrices A,B, and C.
c
      real a(n1,n2),b(n2,n3),c(n1,n3)
c
      call mxmf2(a,n1,b,n2,c,n3)

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
      RK=(K)
      RLOG=LOG10(RK)
      RLOG2=LOG10(2.0)
      RLOG=RLOG/RLOG2+0.5
      LOG2=INT(RLOG)
      return
      END
c-----------------------------------------------------------------------
      subroutine ione(a,n)
      INTEGER  A(1)
      DO 100 I = 1, N
 100     A(I ) = 1
      return
      END
c-----------------------------------------------------------------------
      subroutine isort(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      integer a(1),ind(1)
      integer aa
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
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
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      subroutine sort(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      real a(1),aa
      integer ind(1)
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
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
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
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
C-----------------------------------------------------------------
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
C
C     Compute a Cartesian vector dot product. 2-d version
C
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1)
      DIMENSION V1(1),V2(1)
C
C
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vdot3 (dot,u1,u2,u3,v1,v2,v3,n)
C
C     Compute a Cartesian vector dot product. 3-d version
C
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
C
C
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) + U3(I)*V3(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vcross (u1,u2,u3,v1,v2,v3,w1,w2,w3,n)
C
C     Compute a Cartesian vector cross product.
C
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
      DIMENSION W1(1),W2(1),W3(1)
C
C
      DO 100 I=1,N
         U1(I) = V2(I)*W3(I) - V3(I)*W2(I)
         U2(I) = V3(I)*W1(I) - V1(I)*W3(I)
         U3(I) = V1(I)*W2(I) - V2(I)*W1(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE RZERO3 (A,B,C,N)
      DIMENSION A(1),B(1),C(1)
      DO 100 I=1,N
         A(I)=0.0
         B(I)=0.0
         C(I)=0.0
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UNITVEC (X,Y,Z,N)
      DIMENSION X(1),Y(1),Z(1)
      DO 100 I=1,N
      XLNGTH = SQRT( X(I)**2 + Y(I)**2 + Z(I)**2 )
      IF (XLNGTH.NE.0.0) THEN
         X(I) = X(I)/XLNGTH
         Y(I) = Y(I)/XLNGTH
         Z(I) = Z(I)/XLNGTH
      ENDIF
  100 CONTINUE
      RETURN
      END
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
      DIMENSION A(1)
      DO 100 I = 1, N
 100     A(I ) = 1.0
      return
      END
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
      REAL A(1),B(1),C(1)
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)
  100 CONTINUE
      return
      END
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
      subroutine local_grad3(ur,us,ut,u,N,e,D,Dt)
c     Output: ur,us,ut         Input:u,N,e,D,Dt
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real u (0:N,0:N,0:N,1)
      real D (0:N,0:N),Dt(0:N,0:N)
      integer e
c
      m1 = N+1
      m2 = m1*m1
c
      call mxm(D ,m1,u(0,0,0,e),m1,ur,m2)
      do k=0,N
         call mxm(u(0,0,k,e),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u(0,0,0,e),m2,Dt,m1,ut,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad2(ur,us,u,N,e,D,Dt)
c     Output: ur,us         Input:u,N,e,D,Dt
      real ur(0:N,0:N),us(0:N,0:N)
      real u (0:N,0:N,1)
      real D (0:N,0:N),Dt(0:N,0:N)
      integer e
c
      m1 = N+1
c
      call mxm(D ,m1,u(0,0,e),m1,ur,m1)
      call mxm(u(0,0,e),m1,Dt,m1,us,m1)
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
      REAL A(1)
C
      include 'OPCTR'
C
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine cmult(a,const,n)
      REAL A(1)
C
      include 'OPCTR'
C
      DO 100 I=1,N
         A(I)=A(I)*CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      real function vlsum(vec,n)
      REAL VEC(1)
      include 'OPCTR'
C
      SUM = 0.
C
      DO 100 I=1,N
         SUM=SUM+VEC(I)
 100  CONTINUE
      VLSUM = SUM
      return
      END
c-----------------------------------------------------------------------
      subroutine sub2(a,b,n)
      REAL A(1),B(1)
C
      include 'OPCTR'
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      return
      END
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
C
      DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
  100 CONTINUE

      return
      END
C
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
      REAL A(1),B(1),C(1),D(1)

      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ascol5 (a,b,c,d,e,n)
      REAL A(1),B(1),C(1),D(1),E(1)

      DO 100 I=1,N
         A(I) = B(I)*C(I)-D(I)*E(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine subcol4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
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

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real a(n,n),lam(n),wk(n,n)
      real aa(100)

      call dsyev('V','U',n,a,n,lam,wk,n*n,info)

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