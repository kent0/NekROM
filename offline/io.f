c-----------------------------------------------------------------------
      subroutine rxupt(xyz,upt,ieg,ifread,fname)

      real xyz(1),upt(1)
      integer ieg(1)
      logical ifread(1)
      character*132 fname

      write (6,*) 'starting rxupt'
      call rxupt_open(fname)
      call rxupt_read(xyz,upt,ieg,ifread)
      call rxupt_close(fname)
      write (6,*) 'ending rxupt'

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_open(fname)

      character*132 fname
      integer ieg(1)

      write (6,*) 'starting rxupt_open'
      call my_mfi_prepare(fname)
      write (6,*) 'ending rxupt_open'

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_read(xyz,upt,ieg,ifread)

      ! TODO: add IO without adding OFFLINE

      real xyz(1),upt(1)
      integer ieg(1)
      logical ifread(1)

      common /scrns/ wk(1)

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8

      logical if_byte_sw ! to be removed
      logical ifgetxr ! to be removed
      logical ifgetpr ! to be removed
      logical ifgetur ! to be removed

      write (6,*) 'starting rxupt_read'

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      ic=1
      ie=1

      do while (ic.le.neg)
         mel=min(melt,neg)
         call byte_read(er,mel,ierr)! get element mapping
         if (if_byte_sw) call byte_reverse(er,mel,ierr)
         jc=1
         do while (ie.le.mel)
            if (er(ie).ge.ieg0.and.er(ie).le.ieg1) then
               ieg(er(ie)-ieg0+1)=ie
               jc=jc+1
            endif
            ie=ie+1
         enddo
         ic=ic+mel
      enddo

      iofldsr=0
      if (ifgetxr) iofldsr=ldim

      nelr=ieg1-ieg0+1

      icount=1

      call esort(ieg,ieg1-ieg0+1)

      neltmp=nelr

      iloc=1

      do ig=1,ng
         ie=i3(ig)
         nelr=i4(ig)
         nwk=nelr*ldim*nxyzr8
         offs = offs0 + iofldsr*stride + !ldim*strideB +
     $      ldim*(ie-1)*nxyzr8*wdsizr
         call byte_seek(offs/4,ierr)
         call mfi_getw(wk(iloc),nwk,ifcread)

         iloc=iloc+nwk
      enddo

      nelr=neltmp
      nwk=nelr*ldim*nxyzr8

      if (if_byte_sw) then
      if (wdsizr.eq.8) then
         call byte_reverse8(wk,nwk*2,ierr)
      else
         call byte_reverse(wk,nwk,ierr)
      endif
      endif

      lxyz=lx1*ly1*lz1
      do ie=1,nelr
c        call copy(ux(1,i2(ie)),wk(1+(ie-1)*lxyz*ldim),lxyz)
c        call copy(uy(1,i2(ie)),wk(1+(ie-1)*lxyz*ldim+lxyz),lxyz)
c        if (ldim.eq.3) call
c    $      copy(uz(1,i2(ie)),wk(1+(ie-1)*lxyz*ldim+2*lxyz),lxyz)
      enddo

         iofldsr=iofldsr+ldim
         if (ifgetpr) iofldsr=iofldsr+1
         iloc=1
         do ig=1,ng
            ie=i3(ig)
            nelr=i4(ig)
            nwk=nelr*nxyzr8
            offs = offs0 + iofldsr*stride + !ldim*strideB +
     $         (ie-1)*nxyzr8*wdsizr
            call byte_seek(offs/4,ierr)
            call mfi_getw(wk(iloc),nwk,ifcread)
            iloc=iloc+nwk
         enddo

         nelr=neltmp
         nwk=nelr*nxyzr8

         if (if_byte_sw) then
         if (wdsizr.eq.8) then
            call byte_reverse8(wk,nwk*2,ierr)
         else
            call byte_reverse(wk,nwk,ierr)
         endif
         endif

         do ie=1,nelr
c           call copy(ut(1,i2(ie)),wk(1+(ie-1)*lxyz),lxyz)
         enddo

      write (6,*) 'ending rxupt_read'

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_close(ifcread)

      logical ifcread

      write (6,*) 'starting rxupt_close'
      call byte_close(ierr)
      write (6,*) 'ending rxupt_close'

      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_getw(wk,n)

      real wk(1)

      ! TODO: set n outside to be twice it self if wdsizr = 8

      ierr = 0
      call byte_read(wk,n,ierr)
      ierr=ierr*10*(nid+1)

      return
      end
c-----------------------------------------------------------------------
      subroutine rflist(fnames,ns)

      character*1 fnames(132*ns)
      character*132 fname

      open (unit=2,file='file.list')

      do i=1,ns
         call blank(fnames(1+(i-1)*132),132)
         call blank(fname,132)
         fname=':'
         read (2,1) fname
         if (indx1(fname,':',1).ne.0) goto 3
         call chcopy(fnames(1+(i-1)*132),fname,132)
         itmp=i
      enddo

      close (unit=2)

    3 continue
    1 format(a132)

      return
      end
c-----------------------------------------------------------------------
      subroutine my_mfi_prepare(hname)

      include 'OFFLINE'
      include 'IO'

      logical ifmpiio
      character*132 hname

      common /ipparallel/ nps,lenpb,melt

      integer stride
      character*132 hdr, hname_
      logical if_byte_swap_test,ifcread
      real*4 bytetest

      integer*8 offs0,offs

      ierr = 0

      call chcopy(hname_,hname,132)
      call addfid(hname_,0)
      call byte_open(hname_,ierr)
      if (ierr.ne.0) goto 101

      call blank(hdr,iHeaderSize)
      call byte_read(hdr,iHeaderSize/4,ierr)
      if (ierr.ne.0) goto 101

      call byte_read(bytetest,1,ierr)
      if (ierr.ne.0) goto 101

      if_byte_sw = if_byte_swap_test(bytetest,ierr) ! determine endianess
      if (ierr.ne.0) goto 101

  101 continue

      call mfi_parse_hdr(hdr,ierr)

      if (nelr.gt.lelr) then
         write (6,*) 'ERROR: increase lelr in SIZE!',lelr,nelr
         goto 102
      endif

      stride=np/nfiler
      if (stride.lt.1) then
         write (6,*) nfiler,np,'  TOO MANY FILES, mfi_prepare'
         goto 102
      endif

      pid0r=nid
      pid1r=nid+stride
      fid0r=nid/stride

      call blank(hdr,iHeaderSize)

      call addfid(hname,fid0r)
      call byte_open(hname,ierr)
      if (ierr.ne.0) goto 102

      call byte_read(hdr,iHeaderSize/4,ierr)  
      if (ierr.ne.0) goto 102

      call byte_read(bytetest,1,ierr) 
      if (ierr.ne.0) goto 102

      call mfi_parse_hdr(hdr,ierr) ! replace hdr with correct one

      return
  102 continue

      call byte_close(ierr)
      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_parallel(a,n,fname,nid)

      real a(n)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*128 fname
      character*128 fntrunc

      call blank(fntrunc,128)

      len=ltruncr(fname,128)
      call chcopy(fntrunc,fname,len)

      do id=0,mp-1
         if (id.eq.nid) call dump_parallel_helper(a,n,fntrunc,id.eq.0)
         call nekgsync
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_parallel_helper(a,n,fname,if0)

      real a(n)
      character*128 fname
      logical if0

      if (if0) then
         open (unit=12,file=fname)
      else
         open (unit=12,file=fname,access='append',status='old')
      endif

      do i=1,n
         write (12,1) a(i)
      enddo

      close (unit=12)
    1 format(1pe24.16)

      return
      end
c-----------------------------------------------------------------------
      subroutine esort(ieg,nel)

      integer ieg(1)

      call isort(i1,i2,nel)

      ic=1
      ng=1
c     i3(1)=i1(1)

      call ione(i4,nel)

      do ie=2,nel
c        if (i1(ie)-i1(ie-1).gt.1) then ! if jump
c           ng=ng+1
c           i3(ng)=i1(ie)
c        else
c           i4(ng)=i4(ng)+1
c        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setifread(ifread,ndimt)

      logical ifread(7+ndimt)

      ! xyz

      ifread(1)=.false.
      ifread(2)=.false.
      ifread(3)=.false.

      ! uvw

      ifread(4)=.true.
      ifread(5)=.true.
      ifread(6)=.false.

      ! p

      ifread(7)=.false.

      ! t

      do i=1,ndimt
         ifread(7+i)=.false.
      enddo

      return
      end
c-----------------------------------------------------------------------