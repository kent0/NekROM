c-----------------------------------------------------------------------
      program offline

      include 'LVAR'
      include 'OFFLINE'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      character*132 fname
      logical ifavg0

      call offline_init(icomm)

      nsg=lsg
      nsg=3

      ie=1

      neg=512

      mel=512
      ng=1

      do i=1,ng
         ieg(1)=ie
         ieg(2)=ie+mel-1
         nel=ieg(2)-ieg(1)+1

         call loadsnaps(buf,ieg,indxr,nsg)
         call setgeom(buf,nel)
         call setops(ga,gb,gc,gt,buf(nel*lxyz*ldim+1),nel,nsg,ldim)
         ie=ie+mel
      enddo

      call gop(ga,gt,'+  ',nsg*nsg)
      call gop(gb,gt,'+  ',nsg*nsg)

      ifavg0=.true.

      call setg(gg,gb,gt,nsg,ifavg0)
      call setq(gvec,gval,gg,nsg,ifavg0)

      call settranspose(gvect,gvec,nsg,nsg+1)

      call mxm(gvect,nsg+1,gb,nsg,gt,nsg)
      call mxm(gt,nsg+1,gvec,nsg,gb,nsg+1)

      call mxm(gvect,nsg+1,ga,nsg,gt,nsg)
      call mxm(gt,nsg+1,gvec,nsg,ga,nsg+1)

      call mxm(gvect,nsg+1,gc,nsg,gt,nsg*nsg)
      call mxm(gt,(nsg+1)*nsg,gvec,nsg,gc,nsg+1)
      do i=1,nsg
         call mxm(gc(1+(i-1)*(nsg+1)*nsg),nsg+1,
     $        gvec,nsg,gt(1+(i-1)*(nsg+1)**2),nsg+1)
      enddo

      do j=0,nsg-1
      do i=0,nsg-1
         write (6,*) i,j,ga(1+i+j*(nsg+1)),'ga'
         write (6,*) i,j,gb(1+i+j*(nsg+1)),'gb'
      enddo
      enddo

      do k=0,nsg-1
      do j=0,nsg-1
      do i=0,nsg-1
         write (6,*) i,j,k,gt(1+i+j*(nsg+1)+k*(nsg+1)**2),'gc'
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------