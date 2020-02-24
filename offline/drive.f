c-----------------------------------------------------------------------
      program offline

      include 'LVAR'
      include 'OFFLINE'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      character*132 fname
      logical ifavg0

      call offline_init(icomm)

      nsg=3
      neg=512
      mel=2

      call gengrams(ga,gb,gc,gt,buf,ieg,indxr,
     $  nsg,mp,neg,mel,ldim,lxyz)

      call write_ops(ga,gb,gc,nsg,mid,'g',.true.)

      ifavg0=.true.

      call setg(gg,gb,gt,nsg,ifavg0)
      call setq(gvec,gvect,gvecc,gval,gg,nsg,nsc,mp,mid,ifavg0)

      call mxm(gvect,nsg+1,gb,nsg,gt,nsg)
      call mxm(gt,nsg+1,gvec,nsg,gb,nsg+1)

      call mxm(gvect,nsg+1,ga,nsg,gt,nsg)
      call mxm(gt,nsg+1,gvec,nsg,ga,nsg+1)

      call mxm(gvect,nsg+1,gc,nsg,gt,nsg*nsc)

      do i=1,nsc
         call mxm(gt(1+(i-1)*(nsg+1)*nsg),nsg+1,
     $        gvec,nsg,gc(1+(i-1)*(nsg+1)**2),nsg+1)
      enddo

      call write_ops(ga,gb,gc,nsg,mid,' ',.false.)

      do k=1,nsg
         call mxm(gc,(nsg+1)**2,gvecc(1+(k-1)*nsc),nsc,gt,1)
         call gop(gt,gt((nsg+1)**2+1),'+  ',(nsg+1)**2)
         do j=0,nsg-1
         do i=0,nsg-1
            if (mid.eq.0) write (6,*)
     $         i,j,k-1,gt(1+i+j*(nsg+1)),'c'
         enddo
         enddo
      enddo

      if (mid.eq.0) write (6,*) ' '

      return
      end
c-----------------------------------------------------------------------
      subroutine gengrams(ga,gb,gc,gt,buf,ieg,indxr,
     $   nsg,mp,neg,mel,ldim,lxyz)

      integer ieg(1),indxr(1)
      real ga(1),gb(1),gc(1),gt(1),buf(1)

      ie=1

      call rzero(ga,nsg*nsg)
      call rzero(gb,nsg*nsg)
      call rzero(gc,nsg*nsg*((nsg-1)/mp))

      do while (ie.le.neg)
         ieg(1)=ie
         ieg(2)=min(ie+mel*mp-1,neg)
         call loadsnaps(buf,ieg,indxr,nsg)
         nel=ieg(4)
         call setgeom(buf,nel)
         call setops(ga,gb,gc,gt,buf(nel*lxyz*ldim+1),nel,nsg,ldim)
         ie=ieg(2)+1
      enddo

      call gop(ga,gt,'+  ',nsg*nsg)
      call gop(gb,gt,'+  ',nsg*nsg)

      return
      end
c-----------------------------------------------------------------------