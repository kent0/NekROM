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

      mel=2
      mel=neg
      ng=1

      call rzero(ga,nsg*nsg)
      call rzero(gb,nsg*nsg)
      call rzero(gc,nsg*nsg*((nsg-1)/mp))

      do i=1,ng
         ieg(1)=ie
         ieg(2)=ie+mel-1
         call loadsnaps(buf,ieg,indxr,nsg)
         nel=ieg(4)
         call setgeom(buf,nel)
         call setops(ga,gb,gc,gt,buf(nel*lxyz*ldim+1),nel,nsg,ldim)
         ie=ie+mel
      enddo

      call gop(ga,gt,'+  ',nsg*nsg)
      call gop(gb,gt,'+  ',nsg*nsg)

      do j=0,nsg-1
      do i=0,nsg-1
         if (mid.eq.0) write (6,*) i,j,ga(1+i+j*(nsg)),'ga'
      enddo
      enddo

      do j=0,nsg-1
      do i=0,nsg-1
         if (mid.eq.0) write (6,*) i,j,gb(1+i+j*(nsg)),'gb1'
      enddo
      enddo

      do k=0,nsg-1
      do j=0,nsg-1
      do i=0,nsg-1
         if (mid.eq.0)
     $      write (6,*) i,j,k,gc(1+i+j*(nsg)+k*(nsg)**2),'gc'
      enddo
      enddo
      enddo

      ifavg0=.true.

      call setg(gg,gb,gt,nsg,ifavg0)
      call setq(gvec,gval,gg,nsg,ifavg0)

      call settranspose(gvect,gvec,nsg,nsg+1)

      do i=1,nsg
         if (mod(i-1,mp).eq.mid) nsc=(i-1)/mp+1
      enddo
      write (6,*) 'nsc=',nsc,nsg,mid

      do j=1,nsg+1
      do i=1,nsc
         gvecc(i+(j-1)*nsc)=gvec(1+(i-1)*mp+mid+(j-1)*nsg)
      enddo
      enddo

      call mxm(gvect,nsg+1,gb,nsg,gt,nsg)
      call mxm(gt,nsg+1,gvec,nsg,gb,nsg+1)

      call mxm(gvect,nsg+1,ga,nsg,gt,nsg)
      call mxm(gt,nsg+1,gvec,nsg,ga,nsg+1)

      do k=1,nsc
      do j=1,nsg
      do i=1,nsg
         write (6,*) gc(i+(j-1)*(nsg)+(k-1)*(nsg)*(nsg)),'gc1'
      enddo
      enddo
      enddo

      call mxm(gvect,nsg+1,gc,nsg,gt,nsg*nsc)
      write (6,*) ' '

      do k=1,nsc
      do j=1,nsg
      do i=1,nsg+1
         write (6,*) gt(i+(j-1)*(nsg+1)+(k-1)*(nsg)*(nsg+1)),'gc2'
      enddo
      enddo
      enddo

      do i=1,nsc
         call mxm(gt(1+(i-1)*(nsg+1)*nsg),nsg+1,
     $        gvec,nsg,gc(1+(i-1)*(nsg+1)**2),nsg+1)
      enddo
      write (6,*) ' '

      do k=1,nsc
      do j=1,nsg+1
      do i=1,nsg+1
         write (6,*) gc(i+(j-1)*(nsg+1)+(k-1)*(nsg+1)**2),'gc3'
      enddo
      enddo
      enddo

      write (6,*) ' '

      do k=1,nsg
         call mxm(gc,(nsg+1)**2,gvecc(1+(k-1)*nsc),nsc,gt,1)
         call gop(gt,gt((nsg+1)**2+1),'+  ',(nsg+1)**2)
         do j=0,nsg-1
         do i=0,nsg-1
            if (mid.eq.0) write (6,*)
     $         i,j,k-1,gt(1+i+j*(nsg+1)),'gc4'
         enddo
         enddo
      enddo

      do j=0,nsg-1
      do i=0,nsg-1
         if (mid.eq.0) write (6,*) i,j,ga(1+i+j*(nsg+1)),'ga'
      enddo
      enddo

      do j=0,nsg-1
      do i=0,nsg-1
         if (mid.eq.0) write (6,*) i,j,gb(1+i+j*(nsg+1)),'gb2'
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------