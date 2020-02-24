c-----------------------------------------------------------------------
      program offline

      include 'LVAR'
      include 'OFFLINE'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      character*132 fname
      logical ifavg0,iftherm

      call offline_init(icomm)

      nsg=3
      neg=512
      mel=2

      iftherm=.true.
      ifavg0=.true.

      call gengrams(ga,gb,gc,gt,buf,ieg,indxr,
     $  nsg,mp,neg,mel,ldim,lxyz,iftherm)

      call write_ops(ga,gb,gc,nsg,mid,'  g',.true.)

      call setg(gg,gb,gt,nsg,ifavg0)
      call setq(gvec,gvect,gvecc,gval,gg,nsg,nsc,mp,mid,ifavg0,nsg1)

      call qop2(ga,gt,gvec,gvect,nsg,nsg1)
      call qop2(gb,gt,gvec,gvect,nsg,nsg1)
      call write_ops(ga,gb,gc,nsg1,mid,' ',.false.)

      call qop3(gc,gt,gvec,gvect,gvecc,nsg,nsg1,nsc,mid)

      if (iftherm) then
         call setg(gg,gb((nsg+1)**2+1),gt,nsg,ifavg0)
         call setq(gvec,gvect,gvecc,gval,gg,nsg,nsc,mp,mid,ifavg0,nsg1)

         call qop2(ga((nsg+1)**2+1),gt,gvec,gvect,nsg,nsg1)
         call qop2(gb((nsg+1)**2+1),gt,gvec,gvect,nsg,nsg1)
         call write_ops(ga((nsg+1)**2+1),gb((nsg+1)**2+1),
     $      gc(((nsg-1)/mp+2)*(nsg+1)**2+1),nsg1,mid,' ',.false.)

         call qop3(gc,gt,gvec,gvect,gvecc,nsg,nsg1,nsc,mid)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gengrams(ga,gb,gc,gt,buf,ieg,indxr,
     $   nsg,mp,neg,mel,ldim,lxyz,iftherm)

      integer ieg(1),indxr(1)
      real ga((nsg+1)**2,1),gb((nsg+1)**2,1)
      real gc(((nsg-1)/mp+2)*(nsg+1)**2,1)
      real gt(((nsg-1)/mp+2)*(nsg+1)**2,1)
      real buf(1)

      logical iftherm

      ie=1

      call rzero(ga,nsg*nsg)
      call rzero(gb,nsg*nsg)
      call rzero(gc,nsg*nsg*((nsg-1)/mp+1))

      if (iftherm) then
         call rzero(ga(1,2),nsg*nsg)
         call rzero(gb(1,2),nsg*nsg)
         call rzero(gc(1,2),nsg*nsg*((nsg-1)/mp+1))
      endif

      do while (ie.le.neg)
         ieg(1)=ie
         ieg(2)=min(ie+mel*mp-1,neg)
         call loadsnaps(buf,ieg,indxr,nsg,iftherm)
         nel=ieg(4)
         call setgeom(buf,nel)
         call setops(ga,gb,gc,gt,buf(nel*lxyz*ldim+1),
     $      buf(nel*lxyz*ldim+1),nel,nsg,ldim,ldim)
         if (iftherm) then
            call setops(ga(1,2),gb(1,2),gc(1,2),gt,
     $         buf(nel*lxyz*ldim+nel*lxyz*ldim*nsg+1),
     $         buf(nel*lxyz*ldim+1),nel,nsg,1,ldim)
         endif
         ie=ieg(2)+1
      enddo

      call gop(ga,gt,'+  ',nsg*nsg)
      call gop(gb,gt,'+  ',nsg*nsg)

      if (iftherm) then
         call gop(ga(1,2),gt,'+  ',nsg*nsg)
         call gop(gb(1,2),gt,'+  ',nsg*nsg)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine qop1(g1,gt,gvec,nsg,nsg1)

      real g1(1),gt(1),gvec(1)

      call mxm(g1,1,gvec,nsg,gt,nsg1)
      call copy(g1,gt,nsg1)

      return
      end
c-----------------------------------------------------------------------
      subroutine qop2(ga,gt,gvec,gvect,nsg,nsg1)

      real ga(1),gt(1),gvec(1),gvect(1)

      call mxm(gvect,nsg1,ga,nsg,gt,nsg)
      call mxm(gt,nsg1,gvec,nsg,ga,nsg1)

      return
      end
c-----------------------------------------------------------------------
      subroutine qop3(gc,gt,gvec,gvect,gvecc,nsg,nsg1,nsc,mid)

      real gc(1),gt(1),gvec(1),gvect(1),gvecc(1)

      call mxm(gvect,nsg1,gc,nsg,gt,nsg*nsc)

      do i=1,nsc
         call mxm(gt(1+(i-1)*(nsg1)*nsg),nsg1,
     $        gvec,nsg,gc(1+(i-1)*(nsg1)**2),nsg1)
      enddo

      do k=1,nsg
         call mxm(gc,(nsg1)**2,gvecc(1+(k-1)*nsc),nsc,gt,1)
         call gop(gt,gt((nsg1)**2+1),'+  ',(nsg1)**2)
         do j=0,nsg1-1
         do i=0,nsg1-1
            if (mid.eq.0) write (6,*)
     $         i,j,k-1,gt(1+i+j*(nsg1)),'cc'
         enddo
         enddo
      enddo

      if (mid.eq.0) write (6,*) ' '

      return
      end
c-----------------------------------------------------------------------