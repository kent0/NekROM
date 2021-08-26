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
c        call outpost(err,rpts,vz,pr,t,'err')
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

      if (nid.eq.0) write (6,*) 'deim_check', ifld
      call reconv(vxyz(1,1),vxyz(1,2),vxyz(1,3),u)
      call opsub2(vxyz(1,1),vxyz(1,2),vxyz(1,3),ub,vb,wb)

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
            if (nid.eq.0) write (6,*) i,uc(i,1),peval(i,1),'uc'
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

               if (nid.eq.0) then
                  write (6,*)
     $               i,j,tx_eim,tx_bas,ty_eim,ty_bas,'intp-check t'
                  write (6,*)
     $               i,j,ux_eim,ux_bas,uy_eim,uy_bas,'intp-check u'
               endif
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
      subroutine deim_conv_check(ifld)

      include 'SIZE'
      include 'MOR'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /deim_check2/ vxyz(lt,ldim),peval(lb,4),uc(lb,4),
     $                     dxyz(lt,ldim,ldim),s(lt),
     $                     a1(lb),a2(lb),a3(lb),a4(lb)

      n=lx1*ly1*lz1*nelv
      ifxyo=.true.

      jfield=ifield
      ifield=ifld
      ifield=1

      if (nid.eq.0) write (6,*) 'deim_conv_check', ifld,ncb

      do i=0,nb
         call outpost(ub(1,i),vb(1,i),wb(1,i),pr,tb(1,i),'bbb')
      enddo

      call opcopy(vxyz(1,1),vxyz(1,2),vxyz(1,3),
     $   us0(1,1,1),us0(1,2,1),us0(1,3,1))

      call reconv(vxyz(1,1),vxyz(1,2),vxyz(1,3),u)
      call opsub2(vxyz(1,1),vxyz(1,2),vxyz(1,3),ub,vb,wb)

      call recont(s,ut)
      call sub2(s,tb,n)

c     call opcopy(vxyz(1,1),vxyz(1,2),vxyz(1,3),
c    $   us0(1,1,1),us0(1,2,1),us0(1,3,1))

c     call copy(s,ts0,n)

      do i=0,nb
         if (nio.eq.0) write (6,*) i,u(i),'u'
      enddo

      do i=0,nb
         if (nio.eq.0) write (6,*) i,uk(i,1),'uk'
      enddo

      if (ifld.eq.1) then
         ! reconstruct pod field
         call evalcflds(vxlag,vxyz,vxyz(1,1),1,1)
         call evalcflds(vylag,vxyz,vxyz(1,2),1,1)
         if (ldim.eq.3) call evalcflds(vzlag,vxyz,vxyz(1,3),1,1)

         call evalc_pts_deim(
     $      peval(1,1),u(1),u(1),uxyz_eim,uvw_eim(1,1),nb,ncb,ldim)

         call evalc_pts_deim(
     $      peval(1,2),u(1),u(1),vxyz_eim,uvw_eim(1,2),nb,ncb,ldim)

         if (ldim.eq.3) call evalc_pts_deim(
     $      peval(1,3),u(1),u(1),wxyz_eim,uvw_eim(1,3),nb,ncb,ldim)

c        do i=1,ncb
c           peval(i,1)=vxlag(ipts_deim(i,1),1,1,1,1)
c           peval(i,2)=vylag(ipts_deim(i,2),1,1,1,1)
c           if (ldim.eq.3) peval(i,3)=vzlag(ipts_deim(i,3),1,1,1,1)
c        enddo

         do icb=1,ncb
            call rzero(vxlag(1,1,1,1,2),n)
            call rzero(vylag(1,1,1,1,2),n)
            if (ldim.eq.3) call rzero(vzlag(1,1,1,1,2),n)

            call rzero(vxlag(1,1,1,1,3),n)
            call rzero(vylag(1,1,1,1,3),n)
            if (ldim.eq.3) call rzero(vzlag(1,1,1,1,3),n)

            do jcb=1,icb
               cj=glsc3(vxlag,cbu(1,jcb,1),bm1,n)
     $           /glsc3(cbu(1,jcb,1),cbu(1,jcb,1),bm1,n)
               call add2s2(vxlag(1,1,1,1,2),cbu(1,jcb,1),cj,n)

               write (6,*) jcb,icb,cj,'cj-1'

               cj=glsc3(vylag,cbu(1,jcb,2),bm1,n)
     $           /glsc3(cbu(1,jcb,2),cbu(1,jcb,2),bm1,n)
               call add2s2(vylag(1,1,1,1,2),cbu(1,jcb,2),cj,n)

               write (6,*) jcb,icb,cj,'cj-2'

               if (ldim.eq.3) then
                  cj=glsc3(vzlag,cbu(1,jcb,3),bm1,n)
     $              /glsc3(cbu(1,jcb,3),cbu(1,jcb,3),bm1,n)
                  call add2s2(vzlag(1,1,1,1,2),cbu(1,jcb,3),cj,n)
               endif

            enddo

            call deim_coeff(uc(1,1),peval(1,1),cbu(1,1,1),
     $         irks_deim(1,1),ipts_deim(1,1),nb,icb)
            call deim_coeff(uc(1,2),peval(1,2),cbu(1,1,2),
     $         irks_deim(1,2),ipts_deim(1,2),nb,icb)
            if (ldim.eq.3) call deim_coeff(uc(1,3),peval(1,3),
     $         cbu(1,1,3),
     $         irks_deim(1,3),ipts_deim(1,3),nb,icb)

            write (6,*) icb,uc(1,1),uc(1,2),uc(1,3),'uc'

            do i=1,icb
               call add2s2(vxlag(1,1,1,1,3),cbu(1,i,1),uc(i,1),n)
               call add2s2(vylag(1,1,1,1,3),cbu(1,i,2),uc(i,2),n)
               if (ldim.eq.3) call add2s2(
     $            vzlag(1,1,1,1,3),cbu(1,i,3),uc(i,3),n)
            enddo

c           call outpost(vxlag(1,1,1,1,2),vxlag(1,1,1,1,3),
c    $         vz,pr,t,'err')
c           call exitt0

            call sub2(vxlag(1,1,1,1,2),vxlag,n)
            call sub2(vxlag(1,1,1,1,3),vxlag,n)

            call sub2(vylag(1,1,1,1,2),vylag,n)
            call sub2(vylag(1,1,1,1,3),vylag,n)

            if (ldim.eq.3) then
               call sub2(vzlag(1,1,1,1,2),vzlag,n)
               call sub2(vzlag(1,1,1,1,3),vzlag,n)
            endif

            ex1=sqrt(glsc3(vxlag(1,1,1,1,2),vxlag(1,1,1,1,2),bm1,n)
     $              /glsc3(vxlag,vxlag,bm1,n))
            ex2=sqrt(glsc3(vxlag(1,1,1,1,3),vxlag(1,1,1,1,3),bm1,n)
     $              /glsc3(vxlag,vxlag,bm1,n))

            ey1=sqrt(glsc3(vylag(1,1,1,1,2),vylag(1,1,1,1,2),bm1,n)
     $              /glsc3(vylag,vylag,bm1,n))
            ey2=sqrt(glsc3(vylag(1,1,1,1,3),vylag(1,1,1,1,3),bm1,n)
     $              /glsc3(vylag,vylag,bm1,n))

            if (ldim.eq.3) then
               ez1=sqrt(glsc3(vzlag(1,1,1,1,2),vzlag(1,1,1,1,2),bm1,n)
     $                 /glsc3(vzlag,vzlag,bm1,n))
               ez2=sqrt(glsc3(vzlag(1,1,1,1,3),vzlag(1,1,1,1,3),bm1,n)
     $                 /glsc3(vzlag,vzlag,bm1,n))
            endif

            if (nio.eq.0) then
               write (6,*) icb,ex1,ex2,'x-err'
               write (6,*) icb,ey1,ey2,'y-err'
               if (ldim.eq.3) write (6,*) icb,ez1,ez2,'z-err'
            endif
         enddo
      else if (ifld.eq.2) then
         ! reconstruct pod field
         call evalcflds(vxlag,vxyz,s,1,1)

         call outpost(vxlag,vy,vz,pr,t,'ttt')

         call evalc_pts_deim(
     $      peval(1,1),u(1),ut(1),txyz_eim,uvw_eim(1,4),nb,ncb,ldim)

c        do i=1,ncb
c           peval(i,1)=vxlag(ipts_deim(i,4),1,1,1,1)
c        enddo

         do icb=1,ncb
            call rzero(vxlag(1,1,1,1,2),n)
            call rzero(vxlag(1,1,1,1,3),n)

            do jcb=1,icb
               cj=glsc3(vxlag,cbt(1,jcb),bm1,n)
     $           /glsc3(cbt(1,jcb),cbt(1,jcb),bm1,n)
               call add2s2(vxlag(1,1,1,1,2),cbt(1,jcb),cj,n)

               write (6,*) jcb,icb,cj,'cj'
            enddo
c           call deim_coeff_j(uc(1,1),peval(1,1),ju_eim(1,4),icb,ncb)
            call deim_coeff_j(uc(1,1),peval(1,1),jt_eim,icb,ncb)

            do i=1,icb
               write (6,*) i,icb,uc(i,1),'uc_j'
            enddo

c           call deim_coeff(uc(1,1),peval(1,1),cbt(1,1),
c    $         irks_deim(1,4),ipts_deim(1,4),nb,icb)

c           do i=1,icb
c              write (6,*) i,icb,uc(i,1),'uc_std'
c           enddo

            do i=1,icb
               call add2s2(vxlag(1,1,1,1,3),cbt(1,i),uc(i,1),n)
            enddo

            call outpost(vxlag(1,1,1,1,2),vxlag(1,1,1,1,3),
     $         vz,pr,t,'err')

            call sub2(vxlag(1,1,1,1,2),vxlag,n)
            call sub2(vxlag(1,1,1,1,3),vxlag,n)

            et1=sqrt(glsc3(vxlag(1,1,1,1,2),vxlag(1,1,1,1,2),bm1,n)
     $              /glsc3(vxlag,vxlag,bm1,n))
            et2=sqrt(glsc3(vxlag(1,1,1,1,3),vxlag(1,1,1,1,3),bm1,n)
     $              /glsc3(vxlag,vxlag,bm1,n))

            if (nio.eq.0) write (6,*) icb,et1,et2,'t-err'
         enddo
      endif

      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine pod_conv_check(basis,snaps,ms)

      include 'SIZE'
      include 'MOR'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /deim_check2/ vxyz(lt,ldim),peval(lb,4),uc(lb,4),
     $                     dxyz(lt,ldim,ldim),
     $                     a1(lb),a2(lb),a3(lb),a4(lb)

      real basis(lt,ms),snaps(lt,ms)

      n=lx1*ly1*lz1*nelv

      if (nid.eq.0) write (6,*) 'pod_conv_check',ms

      call pod(basis,evec,eval,ug,snaps,1,'L2 ',ms,ms,.true.)

c     do i=1,ms
c        call ps2b1(evec,snaps(1,i),basis)
c        call dgemm('N','N',n,1,ms,1.,basis,lt,evec,ms,0.,vxlag,lt)
c        call sub3(vylag,snaps(1,i),vxlag,n)
c        ifxyo=.true.
c        time=i*1.0
c        call outpost(snaps(1,i),vxlag,vz,pr,vylag,'pcd')
c        call outpost(snaps(1,i),basis(1,i),vz,pr,vylag,'pcb')
c        err=sqrt(glsc3(vylag,vylag,bm1,n)
c    $           /glsc3(snaps(1,i),snaps(1,i),bm1,n))
c        if (nio.eq.0) write (6,*) i,ms,err,'pcc_err'
c     enddo

      do i=1,ms
         call ps2b1(evec,snaps(1,1),basis,i)
         call dgemm('N','N',n,1,i,1.,basis,lt,evec,i,0.,vxlag,lt)
         call sub3(vylag,snaps(1,1),vxlag,n)
         ifxyo=.true.
         time=i*1.0
         call outpost(snaps(1,1),vxlag,vz,pr,vylag,'pcd')
         call outpost(snaps(1,1),basis(1,i),vz,pr,vylag,'pcb')
         err=sqrt(glsc3(vylag,vylag,bm1,n)
     $           /glsc3(snaps(1,1),snaps(1,1),bm1,n))
         if (nio.eq.0) write (6,*) i,ms,err,'pcc_err'
      enddo

      call exitt0

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
      call dgetrs('N',ncb,1,tmat,ncb,i1,uc,ncb,info)

      return
      end
c-----------------------------------------------------------------------
      subroutine deim_coeff_j(uc,peval,jmat,ncb,mcb)

      include 'SIZE'

      common /deim_coeffr/ tmat(100**2),wk(100**2)
      common /deim_coeffi/ i1(100),i2(100)

      parameter (lt=lx1*ly1*lz1*lelt)

      real uc(ncb),peval(ncb),cb(lt,ncb),jmat(mcb,mcb)
      integer irks(ncb),ipts(ncb)

      call copy(uc,peval,ncb)

      do j=1,ncb
      do i=1,ncb
         tmat(i+(j-1)*ncb)=jmat(i,j)
      enddo
      enddo

c     if (mcb.eq.ncb) then
c        do j=1,ncb
c        do i=1,ncb
c           write (6,*) i,j,jmat(i,j),'jmat'
c        enddo
c        enddo
c        do j=1,ncb
c        do i=1,ncb
c           write (6,*) i,j,tmat(i+(j-1)*ncb),'tmat'
c        enddo
c        enddo
c        call exitt0
c     endif

      call dgetrf(ncb,ncb,tmat,ncb,i1,info)
      call dgetrs('N',ncb,1,tmat,ncb,i1,uc,ncb,info)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_c_eim012(c_eim0,c_eim1,c_eim2,ub,vb,wb,tb,
     $   wk,nb,ifcf)

      include 'SIZE'
      include 'INPUT'
      include 'LMOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /morconvect/ c1v(lxd*lyd*lzd*lelm,0:lub),
     $                    c2v(lxd*lyd*lzd*lelm,0:lub),
     $                    c3v(lxd*lyd*lzd*lelm,0:lub),
     $                    u1v(lxd*lyd*lzd*lelm,0:lub),
     $                    u2v(lxd*lyd*lzd*lelm,0:lub),
     $                    u3v(lxd*lyd*lzd*lelm,0:lub)

      real c_eim0(1),c_eim1(1),c_eim2(1)
      real ub(lt,0:nb),vb(lt,0:nb),wb(lt,0:nb)
      real wk(lt,ldim),tb(lt,0:nb)

      logical ifcf

      nv=lx1*ly1*lz1*nelv

      if (.not.ifaxis) then
         do i=0,nb
            if (.not.ifcf) then
               call set_convect_new(c1v(1,i),c2v(1,i),c3v(1,i),
     $                              ub(1,i),vb(1,i),wb(1,i))
            endif
            call intp_rstd_all(u1v(1,i),tb(1,i),nelv)
         enddo
      endif

      if (.not.ifcf) then
         call rzero(c_eim0,nb)
         call rzero(c_eim1,nb*nb)
         call rzero(c_eim2,nb*nb)
      endif

      call cct(wk,0,0)
      do i=1,nb
         c_eim0(i)=c_eim0(i)+glsc2(wk,tb(1,i),nv)
      enddo

      do j=1,nb
         call cct(wk,0,j)
         do i=1,nb
            c_eim1(i+(j-1)*nb)=c_eim1(i+(j-1)*nb)
     $         +glsc2(wk,tb(1,i),nv)
         enddo
      enddo

      do j=1,nb
         call cct(wk,j,0)
         do i=1,nb
            c_eim2(i+(j-1)*nb)=c_eim2(i+(j-1)*nb)
     $         +glsc2(wk,tb(1,i),nv)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_j_eim3(j_eim,w1,ipts,irks,ncb,cb)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real w1(ncb,ncb),cb(lt,1),j_eim(ncb,ncb)
      integer irks(1),ipts(1)
       
      call rzero(j_eim,ncb*ncb)

      do j=1,ncb
      do i=1,ncb
         if (irks(i).eq.nid) j_eim(i,j)=cb(ipts(i),j)
      enddo
      enddo
      
      call gop(j_eim,w1,'+  ',ncb*ncb)

      return
      end
c-----------------------------------------------------------------------
c     subroutine set_c_eim3_(c_eim3_,nb,ncb,ub,cb)
      subroutine set_c_eim3_(c_eim3_,nb,ncb,ub,cb,bm1)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ub(lt,1),cb(lt,1),c_eim3_(nb,ncb),bm1(1)
      integer irks(1),ipts(1),ipiv(1)

      do j=1,ncb
      do i=1,nb
c        c_eim3_(i,j)=glsc2(ub(1,i),cb(1,j),lx1*ly1*lz1*nelv)
         c_eim3_(i,j)=glsc3(ub(1,i),cb(1,j),bm1,lx1*ly1*lz1*nelv)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_c_eim3(c_eim3,c_eim3_,j_eim,w1,ipiv,nb,ncb)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real c_eim3(1),w1(ncb,ncb),c_eim3_(1),j_eim(1)
      integer ipiv(1)
      
      do j=1,ncb
      do i=1,ncb
         w1(i,j)=1.-min(abs(i-j),1)
      enddo
      enddo
      
      call copy(c_eim3,j_eim,ncb*ncb)
      
      call dgetrf(ncb,ncb,c_eim3,ncb,ipiv,info)
      call dgetrs('N',ncb,ncb,c_eim3,ncb,ipiv,w1,ncb,info)
      
      call mxm(c_eim3_,nb,w1,ncb,c_eim3,ncb)

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_eim012_lowmem(rhs,u,t,c_eim0,c_eim1,c_eim2,nb)

      real rhs(nb),u(nb),t(nb)
      real c_eim0(nb),c_eim1(nb,nb),c_eim2(nb,nb)

      call rzero(rhs,nb)

      do i=1,nb
         rhs(i)=rhs(i)+c_eim0(i)
         do j=1,nb
            rhs(i)=rhs(i)+c_eim1(i,j)*t(j)+c_eim2(i,j)*u(j)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_eim3_lowmem(rhs,u,t,c_eim3,dmat,zmat,nb,ncb,ndim)

      real rhs(nb),u(nb),t(nb)
      real c_eim3(nb,ncb),dmat(ncb,nb,ndim),zmat(ncb,nb,ndim)

      do i=1,ncb
         zdti=0.
         do idim=1,ndim
            dti=0.
            zi=0.
            do j=1,nb
               dti=dti+dmat(i,j,idim)*t(j)
               zi=zi+zmat(i,j,idim)*u(j)
            enddo
            zdti=zdti+dti*zi
         enddo
         do j=1,nb
            rhs(j)=rhs(j)+c_eim3(j,i)*zdti
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc_eim(rhs,uu,tt,ifld)

      include 'SIZE'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      real rhs(1),uu(1),tt(1)

      if (icalld.eq.0) then
         evalc_eim_time=0.
         icalld=1
      endif

      stime=dnekclock()

      if (ifld.eq.1) then
         call evalc_eim012_lowmem(rhs(1),uu(1),uu(1),
     $      cu_eim0,cu_eim1,cu_eim2,nb)

         call evalc_eim3_lowmem(rhs(1),uu(1),uu(1),cu_eim3(1,1),
     $      uxyz_eim,uvw_eim(1,1),nb,ncb,ldim)
         call evalc_eim3_lowmem(rhs(1),uu(1),uu(1),cu_eim3(1,2),
     $      vxyz_eim,uvw_eim(1,2),nb,ncb,ldim)
         if (ldim.eq.3) call evalc_eim3_lowmem(rhs(1),uu(1),uu(1),
     $      cu_eim3(1,3),wxyz_eim,uvw_eim(1,3),nb,ncb,ldim)
      else if (ifld.eq.2) then
         call evalc_eim012_lowmem(rhs(1),uu(1),tt(1),
     $      ct_eim0,ct_eim1,ct_eim2,nb)

         call evalc_eim3_lowmem(rhs(1),uu(1),tt(1),ct_eim3,
     $      txyz_eim,uvw_eim(1,4),nb,ncb,ldim)
      endif

      evalc_eim_time=evalc_eim_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
      subroutine set_intp_eim(
     $   uvw,txyz,ub,vb,wb,tb,itmp1,itmp2,wk,nb,ncb)

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ dx(lt),dy(lt),dz(lt)

      real uvw(ncb,nb,ldim)
      real txyz(ncb,nb,ldim),wk(ncb,nb)
      real ub(lt,0:nb),vb(lt,0:nb),wb(lt,0:nb),tb(lt,0:nb)
      integer itmp1(ncb),itmp2(ncb)

      call rzero(uvw,nb*ncb*ldim)
      call rzero(txyz,nb*ncb*ldim)

      if (nid.eq.0) write (6,*) 'start set_intp_eim',nb,ncb
      do j=1,nb
         call gradm1(dx,dy,dz,tb(1,j))
         do i=1,ncb
            if (itmp2(i).eq.nid) then
               uvw(i,j,1)=ub(itmp1(i),j)
               uvw(i,j,2)=vb(itmp1(i),j)
               if (ldim.eq.3) uvw(i,j,3)=wb(itmp1(i),j)

               txyz(i,j,1)=dx(itmp1(i))
               txyz(i,j,2)=dy(itmp1(i))
               if (ldim.eq.3) txyz(i,j,3)=dz(itmp1(i))
            endif
         enddo
      enddo

      do i=1,ldim
         call gop(uvw(1,1,i),wk,'+  ',ncb*nb)
         call gop(txyz(1,1,i),wk,'+  ',ncb*nb)
      enddo

      if (nid.eq.0) write (6,*) 'end set_intp_eim'

      return
      end
c-----------------------------------------------------------------------
