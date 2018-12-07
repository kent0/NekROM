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
            isplit=1+integer(nsplit*(xm1(2,2,1,i)-xmin)/(xmax-xmin))
            call ifill(domain(1,1,1,i),isplit,lxyz)
         enddo
      else
         zmax=glmax(z,n)
         zmin=glmax(z,n)

         do i=1,nelv
            isplit=1+integer(nsplit*(zm1(2,2,1,i)-zmin)/(zmax-zmin))
            call ifill(domain(1,1,1,i),isplit,lxyz)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine domain_mask(dmask,domain,isplit)

      include 'SIZE'
      include 'GEOM'

      real dmask(lx1*ly1*lz1*lelt)
      integer domain(lx1*ly1*lz1*lelt)

      do i=1,lx1*ly1*lz1*nelv
         dmask(i)=1-min(abs(domain(i)-isplit),1)
      enddo

      return
      end
c-----------------------------------------------------------------------
