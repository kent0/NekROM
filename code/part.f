c-----------------------------------------------------------------------
      subroutine split_domain(domain,nsplit)

      include 'SIZE'
      include 'GEOM'

      integer domain(lx1,ly1,lz1,lelt)

      lxyz=lx1*ly1*lz1
      n=lxyz*nelv

      call izero(domain,n)

      if (ldim.eq.2) then
         xmax=glmax(x,n)
         xmin=glmax(x,n)

         do i=1,nelv
            isplit=1+nint(nsplit*(xm1(2,2,1,i)-xmin)/(xmax-xmin))
            call ifill(domain(1,1,1,i),isplit,lxyz)
         enddo
      else
         zmax=glmax(z,n)
         zmin=glmax(z,n)

         do i=1,nelv
            isplit=1+nint(nsplit*(zm1(2,2,1,i)-zmin)/(zmax-zmin))
            call ifill(domain(1,1,1,i),isplit,lxyz)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine domain_mask(dmask,isplit,nsplit)

      include 'SIZE'
      include 'GEOM'

      real dmask(lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv
      one=1.
      pi=4.*atan(one)

      if (nsplit.eq.1) then
         if (isplit.eq.1) then
            call rone(dmask,n)
         else
            call rzero(dmask,n)
         endif
         if (nio.eq.0) write (6,*) 'warning: nsplit is 1'
         return
      endif

      if (ldim.eq.2) then
         xmin=glmin(xm1,n)
         xmax=glmax(xm1,n)
      else
         xmin=glmin(zm1,n)
         xmax=glmax(zm1,n)
      endif

      dx=(xmax-xmin)/real(nsplit-1)

      dxi=1./dx
      x0=dx*real(isplit-1)

      do i=1,lx1*ly1*lz1*nelv
         if (ldim.eq.2) then
            x=xm1(i,1,1,1)
         else
            x=zm1(i,1,1,1)
         endif

         if (abs(x-x0).lt.dx) then
            dmask(i)=((1.+cos(pi*(x-x0)*dxi))*.5)
         else
            dmask(i)=0.
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine domain_mask_p(dmask,isplit,nsplit)

      include 'SIZE'
      include 'GEOM'

      real dmask(lx1*ly1*lz1*lelt)

      n=lx1*ly1*lz1*nelv
      one=1.
      pi=4.*atan(one)

      if (nsplit.eq.1) then
         if (isplit.eq.1) then
            call rone(dmask,n)
         else
            call rzero(dmask,n)
         endif
         if (nio.eq.0) write (6,*) 'warning: nsplit is 1'
         return
      endif

      if (ldim.eq.2) then
         xmin=glmin(xm1,n)
         xmax=glmax(xm1,n)
      else
         xmin=glmin(zm1,n)
         xmax=glmax(zm1,n)
      endif

      dx=1./real(nsplit)
      dxi=1./dx
      x0=dx*real(isplit-1)

      do i=1,lx1*ly1*lz1*nelv
         if (ldim.eq.2) then
            x=xm1(i,1,1,1)
         else
            x=zm1(i,1,1,1)
         endif
            
         xp=(x-xmin)/(xmax-xmin)
         d=min(abs(xp-x0),abs(xp-x0-1.))
         if (isplit.eq.1) then
            if (xp.lt.dx) then
               dmask(i)=(1.+cos(pi*xp*dxi))*.5
            elseif (xp.gt.(1.-dx)) then
               dmask(i)=(1.+cos(pi*(1.-xp)*dxi))*.5
            else
               dmask(i)=0.
            endif
         else
            if (abs(x0-xp).lt.dx) then
               d=min(x0-xp,x0+1.-xp)
               dmask(i)=(1.+cos(pi*d*dxi))*.5
            else
               dmask(i)=0.
            endif
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
