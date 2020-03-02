c-----------------------------------------------------------------------
      subroutine setxsnap(ixmin,fnames,nsl)

      include 'LVAR'
      include 'IO'

      character*132 fnames(1)
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

c     write (6,*) 'starting setxsnap'

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
      subroutine loadsnaps(buf,ieg,indxr,nsmax,iftherm)

      include 'LVAR'

      integer*8 ibuf8
      common /iread8/ ibuf8(lxyz*leb*lfld)
      common /iread/ ibuf(lxyz*leb*lfld)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*132 fnames
      common /lchr/ fnames(lsg)

      common /isca/ nel,nep,neg,nsl,nsg,ixmin

      common /comm_handles/ ih

      real buf(1)
      integer ieg(1),indxr(1)

c     write (6,*) 'starting load_snap'

      nsg=min(nsmax,lsg)

      if (ieg(1).eq.1) then
        call loadflist(fnames,nsg,nsl)
        call setxsnap(ixmin,fnames,nsl)
        call setindxr(indxr,iftherm)
      endif

      iloc=1
      n=0
      do is=1,nsl
        ieg(4)=is+mid*((nsg-1)/mp+1)
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

      nel=ieg(4)

      m=lxyz*leb*lfld
      mio=mid
      mio=-1

      do i=1,n
         if (mio.eq.0) write (6,*) i,ibuf(i),ibuf8(i),buf(i),'pre'
      enddo
      call fgslib_crystal_tuple_transfer(ih,n,m,ibuf,1,ibuf8,1,buf,1,1)
      do i=1,n
         if (mio.eq.0) write (6,*) i,ibuf(i),ibuf8(i),buf(i),'mid'
      enddo
      call fgslib_crystal_tuple_sort(ih,n,ibuf,1,ibuf8,1,buf,1,2,1)
      do i=1,n
         if (mio.eq.0) write (6,*) i,ibuf(i),ibuf8(i),buf(i),'post'
      enddo

c     do id=0,mp-1
c        if (id.eq.mid) then
c     nel=ieg(4)
c     nel=ieg(2)-ieg(1)+1
c     do idim=1,(nsg+1)*ldim+1
c     write (6,*) 'fld2',id,idim
c     do iel=1,nel
c     do j=1,lxyz
c        write (6,11)
c    $      j,iel,ibuf(j+(iel-1)*lxyz+(idim-1)*lxyz*nel),
c    $            ibuf8(j+(iel-1)*lxyz+(idim-1)*lxyz*nel),
c    $            buf(j+(iel-1)*lxyz+(idim-1)*lxyz*nel)
c     enddo
c     enddo
c     enddo
c     endif
c        call nekgsync
c     enddo
c     call exitt0

 10   format(i5,i5,i5,i5,'  ',1p1e13.5,'  fld1')
 11   format(i5,i5,i5,i5,'  ',1p1e13.5,'  fld2')

c     write (6,*) 'ending load_snap'

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

c     write (6,*) 'starting rxupt'

      ieg(4)=0
      call rxupt_open(fname)
      call rxupt_read(xupt,ibuf,ibuf8,ieg,indxr)
      call rxupt_close(fname)

c     write (6,*) 'ending rxupt'

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_open(fname)

      character*132 fname
      integer ieg(1)

c     write (6,*) 'reading ',fname
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

      if (ieg(1).gt.nelgr) then
         ieg(2)=nelgr-ieg(1)
         return
      endif

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

      call rxupt_read_helper(xupt,wk,ibuf,ibuf8,
     $   ieg(nel+1),ieg(2*nel+1),ieg(3*nel+1),indxr,is,n,nel)

      ieg(1)=ieg0
      ieg(2)=ieg1
      ieg(3)=n
      ieg(4)=nel

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_close

      call byte_close(ierr)

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
      logical if_byte_swap_test

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
      subroutine setindxr(indxr,iftherm)

      integer indxr(1)
      logical iftherm

      indxr(1)=0  ! no xyz
      indxr(2)=7  ! all uvw
      indxr(3)=0  ! no p
      indxr(4)=-1 ! no thermal end

      if (iftherm) then
         indxr(4)=1
         indxr(5)=-1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rxupt_read_helper(xupt,wk,ibuf,ibuf8,i2,i3,i4,
     $   indxr,is,n,nel)

      include 'LVAR'
      include 'IO'

      common /comm_handles/ ih

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real xupt(1),wk(1)
      integer i2(1),i3(1),i4(1),indxr(1),ibuf(1)
      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8,ibuf8(1)

      logical ifdebug

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      ifdebug=.false.

      if (ifdebug) write (6,*) 'wdsizr',wdsizr

      ner=0
      mel=nel
      i=1
      n=0

      if (is.le.0) then
         js=is
         is=1
      endif

      do while (i3(i).ne.0.and.i.le.mel)
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
      do while (indxr(ifld).ne.-1)
         if (ifdebug) write (6,*) 'rrh 1',ifld,ifldt
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
            if (ifdebug) write (6,*) 'rrh 2',ifld

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
            if (ifdebug) write (6,*) 'rrh 3',nfld,ndim,mdim,indx
            ig=1
            iloc=1
            do while (i3(ig).ne.0.and.ig.le.mel)
                ! for now, read mdim stuff instead of nfld
                ie=i3(ig)
                nep=i4(ig)
                nwk=nep*mdim*nxyzr8
                offs=offs0+iofldsr*stride+mdim*(ie-1)*nxyzr8*wdsizr
                call byte_seek(offs/4,ierr)
                call byte_read(wk(iloc),nwk*(wdsizr/4),ierr)
                iloc=iloc+nwk
                ig=ig+1
            enddo
            if (ifdebug) write (6,*) 'rrh 4',ig,iofldsr,stride,iloc

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
            nel=0
            do k=1,ner
               if (ifdebug) write (6,*)
     $            'rrh 5, k=',k,mod(k,mp),mp,mid,nel
               if (mod(k-1,mp).eq.mid) nel=nel+1
            do j=1,nfld
            do i=1,lxyz
               ill=n+i+(j-1)*lxyz+(k-1)*lxyz*nfld
               ibuf8(ill)=
     $            i+(k-1)*lxyz+(is-1)*lxyz*ner+(j+ifldt-2)*lxyz*ner*lsg
               ibuf(ill)=mod(k-1,mp)
            enddo
            enddo
            enddo
            n=n+lxyz*nfld*ner
            indxr(ifld)=indx
            jfld=jfld+nfld
         endif
         ifldt=ifldt+mdim
         ifld=ifld+1
      enddo

c     if (js.lt.1) call
c    $   fgslib_crystal_tuple_sort(ih,n,ibuf,1,ibuf8,1,xupt,1,2,1)

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
          write (6,*) fname,mid,mp,'fname'
  100     call bcast(fname,132)
          if (indx1(fname,'done ',5).eq.0) then
             ms=isg
             if (mod(isg-1,mp).eq.mid) then
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
      subroutine write_ops(ga,gb,gc,nsg,mid,pfx,ifc)

      real ga(nsg,nsg),gb(nsg,nsg),gc(nsg,nsg,nsg)
      logical ifc
      character*3 pfx

      if (mid.eq.0) then
         do j=1,nsg
         do i=1,nsg
            write (6,*) i,j,ga(i,j),pfx,'a'
         enddo
         enddo

         write (6,*) ' '

         do j=1,nsg
         do i=1,nsg
            write (6,*) i,j,gb(i,j),pfx,'b'
         enddo
         enddo

         write (6,*) ' '

         if (ifc) then
            do k=1,nsg
            do j=1,nsg
            do i=1,nsg
               write (6,*) i,j,k,gc(i,j,k),pfx,'c'
            enddo
            enddo
            enddo

            write (6,*) ' '
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine addfid(fname,fid)

      character*1 fname(132)
      integer fid

      character*8  eight,fmt,s8
      save         eight
      data         eight / "????????" /

      do ipass=1,2      ! 2nd pass, in case 1 file/directory
         do k=8,1,-1
            i1 = indx1(fname,eight,k)
            if (i1.ne.0) then ! found k??? string
               write(fmt,1) k,k
               write(s8,fmt) fid
               call chcopy(fname(i1),s8,k)
               goto 10
            endif
         enddo
   10    continue
      enddo

    1 format('(i',i1,'.',i1,')')

      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest,ierr)

      include 'LVAR'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
 
      real*4 bytetest,test2
      real*4 test_pattern
      save   test_pattern
 
      test_pattern = 6.54321
      eps          = 0.00020
      etest        = abs(test_pattern-bytetest)
      if_byte_swap_test = .true.
      if (etest.le.eps) if_byte_swap_test = .false.
 
      test2 = bytetest
      call byte_reverse(test2,1,ierr)
      if (mid.eq.0 .and. loglevel.gt.2) 
     $   write(6,*) 'byte swap:',if_byte_swap_test,bytetest,test2
      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_parse_hdr(hdr,ierr)

      character*132 hdr

      if (indx2(hdr,132,'#std',4).eq.1) then
          call parse_std_hdr(hdr)
      else
         ierr = 1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_std_hdr(hdr)

      include 'LVAR'
      include 'IO'

      character*132 hdr
      character*4 dummy
      logical if_press_mesh

      p0thr = -1
      if_press_mesh = .false.

      read(hdr,*,iostat=ierr) dummy
     $         ,  wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         ,  ifiler,nfiler
     $         ,  rdcode      ! 74+20=94
     $         ,  p0thr, if_press_mesh

      if (ierr.gt.0) then ! try again without pressure format flag
        read(hdr,*,iostat=ierr) dummy
     $         ,  wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         ,  ifiler,nfiler
     $         ,  rdcode      ! 74+20=94
     $         ,  p0thr
      endif

      if (ierr.gt.0) then ! try again without mean pressure
        read(hdr,*,err=99) dummy
     $         ,  wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         ,  ifiler,nfiler
     $         ,  rdcode      ! 74+20=94
      endif

c     set if_full_pres flag
      if_full_pres = .false.
      if (.not.ifsplit) if_full_pres = if_press_mesh

c      ifgtim  = .true.  ! always get time
      ifgetxr = .false.
      ifgetur = .false.
      ifgetpr = .false.
      ifgettr = .false.
      do k=1,ldimt-1
         ifgtpsr(k) = .false.
      enddo

      NPSR = 0
      do i=1,10
         if (rdcode1(i).eq.'X') ifgetxr = .true.
         if (rdcode1(i).eq.'U') ifgetur = .true.
         if (rdcode1(i).eq.'P') ifgetpr = .true.
         if (rdcode1(i).eq.'T') ifgettr = .true.
         if (rdcode1(i).eq.'S') then
            read(rdcode1(i+1),'(I1)') NPS1
            read(rdcode1(i+2),'(I1)') NPS0
            NPSR = 10*NPS1+NPS0
            NPS  = NPSR
            if(NPSR.gt.ldimt-1) NPS=ldimt-1
            do k=1,NPS
               ifgtpsr(k) = .true.
            enddo
            ! nothing will follow
            GOTO 50
         endif
      enddo

  50  if (NPS.lt.NPSR) then
         if (nid.eq.0) then
c          write(*,'(A,/,A)')
c    &      'WARNING: restart file has a NSPCAL > LDIMT',
c    &      'read only part of the fld-data!'
         endif
      endif

      if (NPS.lt.NPSCAL) then
         if (nid.eq.0) then
c          write(*,'(A,/,A)')
c    &      'WARNING: NPSCAL read from restart file differs from ',
c    &      'currently used NPSCAL!'
         endif
      endif

      p0th = 1
      if (p0thr.gt.0) p0th = p0thr

      return

   99 continue   !  If we got here, then the May 2008 variant of std hdr
                 !  failed and we may have an older input file.

      call parse_std_hdr_2006(hdr,rdcode)  ! try the original header format

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_std_hdr_2006(hdr,rlcode)

      include 'LVAR'
		include 'IO'

      character*132 hdr
      character*1 rlcode(20)

c                4  7  10  13   23    33    53    62     68     74
      read(hdr,1) wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         , ifiler,nfiler
     $         , (rlcode(k),k=1,20)                   ! 74+20=94
    1 format(4x,i2,3i3,2i10,e20.13,i9,2i6,20a1)

      if (nid.eq.0) write(6,*) 'WARNING: reading depreacted header!'

      if (nelr.gt.lelr) then
        write(6,*)nid,nelr,lelr,'parse_std_hdr06: inc. lelr in RESTART'
        call exitt
      endif

c     Assign read conditions, according to rdcode
c     NOTE: In the old hdr format: what you see in file is what you get.
c      ifgtim  = .true.  ! always get time
      ifgetxr = .false.
      ifgetur = .false.
      ifgetpr = .false.
      ifgettr = .false.
      do k=1,npscal
         ifgtpsr(k) = .false.
      enddo

      if (rlcode(1).eq.'X') ifgetxr = .true.
      if (rlcode(2).eq.'U') ifgetur = .true.
      if (rlcode(3).eq.'P') ifgetpr = .true.
      if (rlcode(4).eq.'T') ifgettr = .true.
      do k=1,npscal
         if (rlcode(4+k).ne.' ') ifgtpsr(k) = .true.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flush_io
      return
      end
c-----------------------------------------------------------------------
      subroutine dump_serial(a,n,fname,nid)

      real a(n)

      character*128 fname
      character*128 fntrunc

      write (6,*) 'nid',nid
      if (nid.eq.0) then
         call blank(fntrunc,128)

         len=ltruncr(fname,128)
         call chcopy(fntrunc,fname,len)

         call dump_serial_helper(a,n,fntrunc)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_serial_helper(a,n,fname)

      real a(n)

      character*128 fname

      open (unit=12,file=fname)
      write (6,*) 'writing to ',fname

      do i=1,n
         write (12,1) a(i)
      enddo

      close (unit=12)
    1 format(1pe24.16)

      return
      end
c-----------------------------------------------------------------------