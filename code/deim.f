c-----------------------------------------------------------------------
      subroutine deim_fpts(ipts,irks,flds,nflds)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      logical ifmult,iftmp

      parameter (lt=lx1*ly1*lz1*lelt)

      common /romup/ rom_time
      common /scrns/ err(lt),rpts(lt)

      real flds(lx1,ly1,lz1,lelt,nflds)
      integer ipts(nflds),irks(nflds)

      n=lx1*ly1*lz1*nelt

      call rzero(rpts,n)
      call izero(ipts,nflds)
      call izero(irks,nflds)

      iftmp=ifxyo

      do iter=1,nflds
         ifxyo=(iter.eq.1).or.iftmp
         call copy(err,flds(1,1,1,1,iter),n)
         if (iter.ge.2) then
            call rzero(rtmp1,(iter-1)**2)
            call rzero(rtmp2,iter-1)

            do j=1,iter-1
               do i=1,iter-1
                  if (irks(i).eq.nid)
     $               rtmp1(i+(j-1)*(iter-1),1)=flds(ipts(i),1,1,1,j)
               enddo
               if (irks(j).eq.nid) rtmp2(j,1)=flds(ipts(j),1,1,1,i)
            enddo

            call gop(rtmp1,rtmp3,'+  ',(iter-1)**2)
            call gop(rtmp2,rtmp3,'+  ',iter-1)

            call dgetrf(iter-1,iter-1,rtmp1,iter-1,ipiv,info)
            call dgetrs(
     $         'N',iter-1,1,rtmp1,iter-1,ipiv,rtmp2,iter-1,info)
            do i=1,iter-1
               call add2s2(err,flds(1,1,1,1,i),-rtmp2(i,1),n)
            enddo
         endif

         call glamax_ind(emax,iproc,ind,err,n)
         if (nid.eq.iproc) write (6,*) iter,emax,iproc,ind,'glamax_ind'

         irks(iter)=iproc
         ipts(iter)=ind
         if (irks(iter).eq.nid) rpts(ipts(iter))=rpts(ipts(iter))+1.
         call outpost(err,rpts,vz,pr,t,'err')
      enddo

      ifxyo=iftmp

      return
      end
c-----------------------------------------------------------------------
      subroutine glamin_ind(vmin,iproc,ind,vec,n)

      include 'SIZE'

      common /cube1/ node,pid,np,nullpid,node0

      real vec(1)

      vmin_loc=vlmin(vec,n)
      vmin=glamin(vmin_loc,1)

      iproc=np
      ind=0

      if (vmin.eq.vmin_loc) iproc=nid

      iproc=iglmin(iproc,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine glamax_ind(vmax,iproc,ind,vec,n)

      include 'SIZE'

      real vec(1)

      vmax_loc=vlamax(vec,n)
      vmax=glmax(vmax_loc,1)

      iproc=iglmax(nid,1)+1
      ind=0

      if (vmax.eq.vmax_loc) iproc=nid

      iproc=iglmin(iproc,1)

      if (nid.eq.iproc) then
         do i=1,n
            if (vmax.eq.abs(vec(i))) then
               ind=i
               goto 10
            endif
         enddo
      endif

   10 continue

      return
      end
c-----------------------------------------------------------------------
