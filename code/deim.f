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
      subroutine evalc_pts_deim(peval,u,t,dmat,zmat,nb,ncb,ndim)

      real peval(nb),u(nb),t(nb)
      real dmat(ncb,nb,ndim),zmat(ncb,nb,ndim)

      call rzero(peval,nb)

      do idim=1,ndim
      do i=1,ncb
         dti=0.
         zi=0.
         do j=1,nb
            dti=dti+dmat(i,j,idim)*t(j)
            zi=zi+zmat(i,j,idim)*u(j)
         enddo
         peval(i)=peval(i)+zi*dti
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine deim_check(ifld)

      include 'SIZE'
      include 'MOR'
      include 'SOLN'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /deim_check2/ vxyz(lt,ldim),peval(lb,4),uc(lb,4),
     $                     dxyz(lt,ldim,ldim),
     $                     a1(lb),a2(lb),a3(lb),a4(lb)

      n=lx1*ly1*lz1*nelv

      write (6,*) 'deim_check', ifld
      call reconv(vxyz(1,1),vxyz(1,2),vxyz(1,3),u)
      call opsub2(vxyz(1,1),vxyz(1,2),vxyz(1,3),ub,vb,wb)

      write (6,*) 'wp 0'

      if (ifld.eq.1) then
         ! evaluate points and reconstruct deim field
         call evalc_pts_deim(
     $      peval(1,1),u,u,uxyz_eim,uvw_eim(1,1),nb,ncb,ldim)
         call deim_coeff(uc(1,1),peval(1,1),cbu(1,1,1),
     $      irks_deim(1,1),ipts_deim(1,1),nb,ncb)

         call evalc_pts_deim(
     $      peval(1,2),u,u,vxyz_eim,uvw_eim(1,2),nb,ncb,ldim)
         call deim_coeff(uc(1,2),peval(1,2),cbu(1,1,2),
     $      irks_deim(1,2),ipts_deim(1,2),nb,ncb)

         if (ldim.eq.3) then
            call evalc_pts_deim(
     $         peval(1,3),u,u,wxyz_eim,uvw_eim(1,3),nb,ncb,ldim)
            call deim_coeff(uc(1,3),peval(1,3),cbu(1,1,3),
     $         irks_deim(1,3),ipts_deim(1,3),nb,ncb)
         endif

         call opzero(vxlag,vylag,vzlag)

         do i=1,ncb
            call add2s2(vxlag,cbu(1,i,1),uc(i,1),n)
            call add2s2(vylag,cbu(1,i,2),uc(i,2),n)
            if (ldim.eq.3) call add2s2(vzlag,cbu(1,i,3),uc(i,3),n)
         enddo

         ! reconstruct pod field
         call evalcflds(vxlag(1,1,1,1,2),vxyz,vxyz(1,1),1,1)
         call evalcflds(vylag(1,1,1,1,2),vxyz,vxyz(1,2),1,1)
         if (ldim.eq.3)
     $      call evalcflds(vzlag(1,1,1,1,2),vxyz,vxyz(1,3),1,1)

         ! check point values match
         do i=1,ncb
            if (irks_deim(i,1).eq.nid) then
               v1=peval(i,1)
               v2=vxlag(ipts_deim(1,1),1,1,1,1)
               v3=vxlag(ipts_deim(1,1),1,1,1,2)
               write (6,*)
     $            1,i,v1,v2,v3,v1-v2,v1-v3,v2-v3,'peval-check-u'
            endif
         enddo
         do i=1,ncb
            if (irks_deim(i,2).eq.nid) then
               v1=peval(i,2)
               v2=vylag(ipts_deim(1,2),1,1,1,1)
               v3=vylag(ipts_deim(1,2),1,1,1,2)
               write (6,*)
     $            2,i,v1,v2,v3,v1-v2,v1-v3,v2-v3,'peval-check-v'
            endif
         enddo
         do i=1,ncb
            if (irks_deim(i,3).eq.nid) then
               v1=peval(i,3)
               v2=vzlag(ipts_deim(1,3),1,1,1,1)
               v3=vzlag(ipts_deim(1,3),1,1,1,2)
               write (6,*)
     $            3,i,v1,v2,v3,v1-v2,v1-v3,v2-v3,'peval-check-w'
            endif
         enddo
      else
         call recont(t,ut)
         call sub2(t,tb,n)

         call gradm1(dxyz(1,1,1),dxyz(1,2,1),dxyz(1,3,1),t)

         ! evaluate points and reconstruct deim field
         call evalc_pts_deim(
     $      peval,u(1),ut(1),txyz_eim,uvw_eim(1,4),nb,ncb,ldim)

         call deim_coeff(
     $      uc,peval,cbt,irks_deim(1,4),ipts_deim(1,4),nb,ncb)

         call rzero(tlag,n)
         do i=1,ncb
            call add2s2(tlag,cbt(1,i),uc(i,1),n)
         enddo

         ! reconstruct pod field
         call evalcflds(tlag(1,1,1,1,1,2),vxyz,t,1,1)

         do i=1,ncb
            write (6,*) i,uc(i,1),peval(i,1),'uc'
         enddo

         ! check point values match
         do i=1,ncb
            if (irks_deim(i,4).eq.nid) then
               v1=peval(i,1)
               v2=tlag(ipts_deim(i,4),1,1,1,1,1)
               v3=tlag(ipts_deim(i,4),1,1,1,1,2)
               peval(i,1)=v3
               write (6,*) i,v1,v2,v3,v1-v2,v1-v3,v2-v3,'peval-check1'
            endif
         enddo

         call mxm(txyz_eim(1+0*nb*ncb),ncb,ut(1),nb,a1,1)
         call mxm(txyz_eim(1+1*nb*ncb),ncb,ut(1),nb,a2,1)

         call mxm(uvw_eim(1+0*nb*ncb,4),ncb,u(1),nb,a3,1)
         call mxm(uvw_eim(1+1*nb*ncb,4),ncb,u(1),nb,a4,1)

         ! check convection components
         do i=1,ncb ! at each point
            if (irks_deim(i,4).eq.nid) then
               dx=dxyz(ipts_deim(i,4),1,1)
               dy=dxyz(ipts_deim(i,4),2,1)
               ux=vxyz(ipts_deim(i,4),1)
               uy=vxyz(ipts_deim(i,4),2)
               ct=dx*ux+dy*uy

               write (6,*) i,dx,dy,ux,uy,ct,'conv-eval-check 1'

               dx=a1(i)
               dy=a2(i)
               ux=a3(i)
               uy=a4(i)
               ct=dx*ux+dy*uy

               write (6,*) i,dx,dy,ux,uy,ct,'conv-eval-check 2'
            endif
         enddo

         do i=1,nb ! at each basis
            call gradm1(dxyz(1,1,1),dxyz(1,2,1),dxyz(1,3,1),tb(1,i))
            do j=1,ncb ! at each point
               tx_eim=txyz_eim(j+(i-1)*ncb+0*nb*ncb)
               ty_eim=txyz_eim(j+(i-1)*ncb+1*nb*ncb)
               ux_eim=uvw_eim(j+(i-1)*ncb+0*nb*ncb,4)
               uy_eim=uvw_eim(j+(i-1)*ncb+1*nb*ncb,4)

               tx_bas=dxyz(ipts_deim(j,4),1,1)
               ty_bas=dxyz(ipts_deim(j,4),2,1)
               ux_bas=ub(ipts_deim(j,4),i)
               uy_bas=vb(ipts_deim(j,4),i)

         write (6,*) i,j,tx_eim,tx_bas,ty_eim,ty_bas,'intp-check t'
         write (6,*) i,j,ux_eim,ux_bas,uy_eim,uy_bas,'intp-check u'
            enddo
         enddo

         do i=1,ncb ! at each basis
            write (6,*) 'intp-check'
         enddo

         call outpost(tlag,tlag(1,1,1,1,1,2),vz,pr,t,'ect')

         call deim_coeff(
     $      uc,peval,cbt,irks_deim(1,4),ipts_deim(1,4),nb,ncb)

         call rzero(tlag,n)
         do i=1,ncb
            call add2s2(tlag,cbt(1,i),uc(i,1),n)
         enddo

         do i=1,ncb
            if (irks_deim(i,4).eq.nid) then
               v1=peval(i,1)
               v2=tlag(ipts_deim(i,4),1,1,1,1,1)
               v3=tlag(ipts_deim(i,4),1,1,1,1,2)
               peval(i,1)=v3
               write (6,*) i,v1,v2,v3,v1-v2,v1-v3,v2-v3,'peval-check2'
            endif
         enddo
         call outpost(tlag,tlag(1,1,1,1,1,2),vz,pr,t,'ect')

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine deim_coeff(uc,peval,cb,irks,ipts,nb,ncb)

      include 'SIZE'

      common /deim_coeffr/ tmat(100**2),wk(100**2)
      common /deim_coeffi/ i1(100),i2(100)

      parameter (lt=lx1*ly1*lz1*lelt)

      real uc(ncb),peval(ncb),cb(lt,ncb)
      integer irks(ncb),ipts(ncb)

      call rzero(tmat,ncb*ncb)
      call copy(uc,peval,ncb)

      do j=1,ncb
      do i=1,ncb
         if (irks(i).eq.nid) tmat(i+(j-1)*ncb)=cb(ipts(i),j)
      enddo
      enddo

      call gop(tmat,wk,'+  ',ncb*ncb)

      call dgetrf(ncb,ncb,tmat,ncb,i1,info)
c     write (6,*) info,'info'
c     do j=1,ncb
c     do i=1,ncb
c        write (6,*) i,j,tmat(i+(j-1)*ncb),'rrr'
c     enddo
c     enddo
c     do i=1,ncb
c        write (6,*) i,uc(i),irks(i),ipts(i),i1(i),'sss'
c     enddo
      call dgetrs('N',ncb,1,tmat,ncb,i1,uc,ncb,info)

c     do i=1,ncb
c        write (6,*) i,uc(i),'ttt'
c     enddo

      return
      end
c-----------------------------------------------------------------------
