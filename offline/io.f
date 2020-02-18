c-----------------------------------------------------------------------
      subroutine setxsnap(ixmin,fnames,nsl)

      include 'LVAR'
      include 'IO'

      character*132 fnames(1)
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      write (6,*) 'starting setxsnap'

      nslmax=iglmax(nsl,1)
      ixmin=nslmax*mp+1

      do is=1,nsl
        call rxupt_open(fnames(is))
        call rxupt_close
        if (ifgetxr) then
           ixmin=is+mid*nslmax
           exit
        endif
      enddo

      jxmin=iglmin(ixmin,1)

      if (jxmin.eq.(nslmax*mp+1))
     $   call exitti('could not find any geometry, exiting$',jxmin)

      if (jxmin.eq.ixmin) then
         ixmin=jxmin-mid*nslmax
      else
         ixmin=0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine loadsnaps(buf,ieg,indxr,nsmax)

      include 'LVAR'

      integer*8 ibuf8
      common /iread8/ ibuf8(lxyz*leb*lfld)
      common /iread/ ibuf(lxyz*leb*lfld)

      character*132 fnames
      common /lchr/ fnames(lsg)

      common /isca/ nel,nep,neg,nsl,nsg,ixmin

      common /comm_handles/ ih

      real buf(1)
      integer ieg(1),indxr(1)

      write (6,*) 'starting load_snap'

      nsg=min(nsmax,lsg)

      if (ieg(1).eq.1) then
        call loadflist(fnames,nsg,nsl)
        call setxsnap(ixmin,fnames,nsl)
        call setindxr(indxr)
      endif

      iloc=1
      n=0
      do is=1,nsl
        ieg(4)=is
        if (is.eq.ixmin) then
           indxr(1)=7
        else
           indxr(1)=0
        endif
        call rxupt_open(fnames(is))
        call rxupt_read(buf(iloc),ibuf(iloc),ibuf8(iloc),ieg,indxr)
        call rxupt_close
        iloc=iloc+ieg(3)
        n=n+ieg(3)
      enddo

      m=lxyz*leb*lfld
      do i=1,n
c        write (6,*) i,ibuf(i),ibuf8(i),buf(i),'pre'
      enddo
      call fgslib_crystal_tuple_transfer(ih,n,m,ibuf,1,ibuf8,1,buf,1,1)
      do i=1,n
c        write (6,*) i,ibuf(i),ibuf8(i),buf(i),'mid'
      enddo
      call fgslib_crystal_tuple_sort(ih,n,ibuf,1,ibuf8,1,buf,1,2,1)
      do i=1,n
c        write (6,*) i,ibuf(i),ibuf8(i),buf(i),'post'
      enddo

      write (6,*) 'ending load_snap'

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt(xupt,ieg,indxr,fname)

      include 'LVAR'

      real xupt(1)
      integer ieg(1),indxr(1)
      character*132 fname

      integer*8 ibuf8
      common /iread8/ ibuf8(lxyz*leb*lfld)
      common /iread/ ibuf(lxyz*leb*lfld)

      write (6,*) 'starting rxupt'

      ieg(4)=0
      call rxupt_open(fname)
      call rxupt_read(xupt,ibuf,ibuf8,ieg,indxr)
      call rxupt_close(fname)

      write (6,*) 'ending rxupt'

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_open(fname)

      character*132 fname
      integer ieg(1)

      write (6,*) 'reading ',fname
      call my_mfi_prepare(fname)

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_read(xupt,ibuf,ibuf8,ieg,indxr)

      include 'LVAR'
      include 'IO'

      common /scrread/ wk(lxyz*ldim*lel)

      real xupt(1)
      integer ieg(1),indxr(1),ibuf(1)
      integer*8 ibuf8(1)

      write (6,*) 'starting rxupt_read'

      ieg(2)=min(ieg(2),nelgr)

      ieg0=ieg(1)
      ieg1=ieg(2)
      is=ieg(4)
      jeg=1
      iel=1

      nel=ieg1-ieg0+1

      do while (iel.lt.nel)
         mel=min(lel,nelgr-jeg+1)
         call byte_read(er,mel,ierr) ! get element mapping
         if (if_byte_sw) call byte_reverse(er,mel,ierr)
         do ie=1,mel
            if (er(ie).ge.ieg0.and.er(ie).le.ieg1) then
               ieg(er(ie)-ieg0+1)=jeg
               iel=iel+1
            endif
            jeg=jeg+1
         enddo
      enddo

      call esort(ieg,nel)

      do i=1,nel
         write (6,*) ieg(i),ieg(nel+i),ieg(2*nel+i),ieg(3*nel+i),'ieg'
      enddo

      call rxupt_read_helper(xupt,wk,ibuf,ibuf8,
     $   ieg(nel+1),ieg(2*nel+1),ieg(3*nel+1),indxr,is,n)

      ieg(1)=ieg0
      ieg(2)=ieg1
      ieg(3)=n

      write (6,*) 'ending rxupt_read'

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_close

      write (6,*) 'starting rxupt_close'
      call byte_close(ierr)
      write (6,*) 'ending rxupt_close'

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
      subroutine setindxr(indxr)

      integer indxr(1)

      indxr(1)=0  ! no xyz
      indxr(2)=7  ! all uvw
      indxr(3)=0  ! no p
      indxr(4)=-1 ! no thermal end

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_read_helper(xupt,wk,ibuf,ibuf8,i2,i3,i4,
     $   indxr,is,n)

      include 'LVAR'
      include 'IO'

      common /comm_handles/ ih

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real xupt(1),wk(1)
      integer i2(1),i3(1),i4(1),indxr(1),ibuf(1)
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8,ibuf8(1)

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      ner=0
      i=1
      n=0

      if (is.le.0) then
         js=is
         is=1
      endif

      do while (i3(i).ne.0)
         ner=ner+i4(i)
         i=i+1
      enddo

      if (nzr.eq.1) then
         ndim=2
      else
         ndim=3
      endif

      nfldt=7+ldimt
      nfld=1

      ifld=1
      ifldt=1
      jfld=1
      write (6,*) 'wp 6.2'
      do while (indxr(ifld).ne.-1)
         write (6,*) 'wp 6.9',ifld,ifldt
         if (ifld.le.2) then
             mdim=ndim
         else
             mdim=1
         endif
         indx=indxr(ifld)
         if (indx.ne.0) then
            nfld=1

            iofldsr=0
            if (ifld.ge.2.and.ifgetxr) iofldsr=iofldsr+ndim
            if (ifld.ge.3.and.ifgetur) iofldsr=iofldsr+ndim
            if (ifld.ge.4.and.ifgetpr) iofldsr=iofldsr+1+(ifld-4)
            write (6,*) 'wp 6.3',ifld

            if (ifld.le.2) then
               if (ndim.eq.2) then
                  indx=min(indx,3)
               else
                  indx=min(indx,7)
               endif
               if (indx.eq.3.or.indx.eq.5.or.indx.eq.6) then
                  nfld=2
               else
                  nfld=3
               endif
            else
               indx=1
            endif
            write (6,*) 'wp 6.4',nfld,ndim,mdim,indx
            iloc=1
            ig=1
            do while (i3(ig).ne.0) ! for now, read mdim stuff instead of nfld
                ie=i3(ig)
                nep=i4(ig)
                nwk=nep*mdim*nxyzr8
                if (wdsizr.eq.8) nwk=nwk*2
                offs=offs0+iofldsr*stride+mdim*(ie-1)*nxyzr8*wdsizr
                call byte_seek(offs/4,ierr)
                call byte_read(wk(iloc),nwk,ierr)
                iloc=iloc+nwk
                ig=ig+1
            enddo
            write (6,*) 'wp 6.5',ig,ng

            nwk=ner*ldim*nxyzr8

            if (if_byte_sw) then
            if (wdsizr.eq.8) then
                call byte_reverse8(wk,nwk*2,ierr)
            else
                call byte_reverse(wk,nwk,ierr)
            endif
            endif
            ! TODO; correctly copy in case nfld != mdim
            call copy(xupt(n+1),wk,lxyz*nfld*ner)
            do k=1,ner
            do j=1,nfld
            do i=1,lxyz
               ill=n+i+(j-1)*lxyz+(k-1)*lxyz*nfld
               ibuf8(ill)=
     $            i+(k-1)*lxyz+(is-1)*lxyz*ner+(j+ifldt-2)*lxyz*ner*lsg
               ibuf(ill)=mod(k,mp)
            enddo
            enddo
            enddo
            n=n+lxyz*nfld*ner
            indxr(ifld)=indx
            jfld=jfld+nfld
            write (6,*) 'wp 6.8'
         endif
         ifldt=ifldt+mdim
         ifld=ifld+1
      enddo

      if (js.lt.1) call
     $   fgslib_crystal_tuple_sort(ih,n,ibuf,1,ibuf8,1,xupt,1,2,1)

      if (is.le.0) then
         is=js
         is=1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine loadflist(fnames,nsg,nsl)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*132 fnames(ns),fname

      if (mid.eq.0) open (unit=10,file='file.list')
      ms=1

      nsl=0
      do isg=1,nsg
          call blank(fname,132)
          fname='done '
          if (mid.eq.0) read (10,'(a132)',end=100) fname
          write (6,*) fname,'fname'
  100     call bcast(fname,132)
          if (indx1(fname,'done ',5).eq.0) then
             ms=isg
             if (mod(isg,mp).eq.mid) then
                nsl=nsl+1
                call chcopy(fnames(nsl),fname,132)
             endif
          else
             goto 200
          endif
      enddo

  200 nsg=ms
      if (mid.eq.0) close (unit=10)

      return
      end
c-----------------------------------------------------------------------