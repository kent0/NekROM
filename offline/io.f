c-----------------------------------------------------------------------
      subroutine rxupt(xupt,ieg,indxr,fname)

      real xupt(1)
      integer ieg(1),indxr(1)
      character*132 fname

      write (6,*) 'starting rxupt'
      call rxupt_open(fname)
      call rxupt_read(xupt,ieg,indxr)
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
      subroutine rxupt_read(xupt,ieg,indxr)

      include 'LVAR'
      include 'IO'

      real xupt(1)
      integer ieg(1),indxr(1)

      common /scrread/ wk(lxyz*ldim*lel)

      write (6,*) 'starting rxupt_read'

      ieg(2)=min(ieg(2),nelgr)

      ieg0=ieg(1)
      ieg1=ieg(2)

      ic=1
      ie=1

      mel=ieg1-ieg0+1

      do while (ic.le.mel)
         nel=min(lel,mel-ic+1)
         call byte_read(er,nel,ierr)! get element mapping
         if (if_byte_sw) call byte_reverse(er,nel,ierr)
         jc=1
         do while (ie.le.nel)
            if (er(ie).ge.ieg0.and.er(ie).le.ieg1) then
               ieg(er(ie)-ieg0+1)=ie
               jc=jc+1
            endif
            ie=ie+1
         enddo
         ic=ic+nel
      enddo

      call esort(ieg,mel)

      call rxupt_read_helper1(xupt,wk,
     $   ieg(mel+1),ieg(2*mel+1),ieg(3*mel+1),indxr)

      ieg(1)=ieg0
      ieg(2)=ieg1

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

      call byte_read(wk,n,ierr)

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

      include 'LVAR'
      include 'IO'

      character*132 hdr,hname,hname_
      real*4 bytetest

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

      istride=np/nfiler
      if (istride.lt.1) then
         write (6,*) nfiler,np,'  TOO MANY FILES, mfi_prepare'
         goto 102
      endif

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

      integer ieg(nel,4)

      call isort(ieg,ieg(1,2),nel)

      call izero(ieg(1,3),nel)
      call ione(ieg(1,4),nel)

      ng=1
      ieg(1,3)=ieg(1,1)

      do ie=2,nel
         if (ieg(ie,1)-ieg(ie-1,1).gt.1) then ! jump
            ng=ng+1
            ieg(ng,3)=ieg(ie,1)
         else
            ieg(ng,4)=ieg(ng,4)+1
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setifread(ifread,ndimt)

      logical ifread(3+ndimt)

      ! xup

      ifread(1)=.false.
      ifread(2)=.false.
      ifread(3)=.false.

      ! t

      do i=1,ndimt
         ifread(3+i)=.false.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setindxr(indxr,ndimt)

      integer indxr(3+ndimt)

      ! U only

      indxr(1)=2
      indxr(2)=0

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_read_helper1(xupt,wk,i2,i3,i4,indxr)

      include 'LVAR'

      real xupt(1),wk(1)
      integer i2(1),i3(1),i4(1),indxr(1)

      ner=0
      i=1
      do while (i3(i).ne.0)
         ner=ner+i4(i)
         i=i+1
      enddo

      i=1
      jdim=1
      do while (indxr(i).ne.0)
         call rxupt_read_helper2(wk,ner,i3,i4,indxr(i))
         ndim=1
         if (indxr(i).le.2) ndim=ldim
         do idim=1,ndim
            do ie=1,ner
               call copy(xupt(1+(i2(ie)-1)*lxyz+lxyz*ner*(jdim-1)),
     $            wk(1+(idim-1)*lxyz+(ie-1)*lxyz*ndim),lxyz)
            enddo
            jdim=jdim+1
         enddo
         i=i+1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_read_helper2(wk,ner,i3,i4,ifld)

      include 'LVAR'
      include 'IO'

      real wk(1)
      integer i3(1),i4(1)
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8

      iofldsr=0

      if (ifld.ge.2.and.ifgetxr) iofldsr=iofldsr+ldim
      if (ifld.ge.3.and.ifgetur) iofldsr=iofldsr+ldim
      if (ifld.ge.4.and.ifgetpr) iofldsr=iofldsr+1

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      iloc=1
      ig=1
      do while (i4(ig).ne.0)
         ie=i3(ig)
         nep=i4(ig)
         nwk=nep*ldim*nxyzr8
         if (wdsizr.eq.8) nwk=nwk*2
         offs = offs0 + iofldsr*stride + !ldim*strideB +
     $      ldim*(ie-1)*nxyzr8*wdsizr
         call byte_seek(offs/4,ierr)
         call mfi_getw(wk(iloc),nwk)

c        write (6,*) ig,i3(ig),i4(ig),'i3/i4'
         iloc=iloc+nwk
         ig=ig+1
      enddo
c     call exitt0

      nwk=ner*ldim*nxyzr8

      if (if_byte_sw) then
      if (wdsizr.eq.8) then
         call byte_reverse8(wk,nwk*2,ierr)
      else
         call byte_reverse(wk,nwk,ierr)
      endif
      endif

      return
      end
c-----------------------------------------------------------------------