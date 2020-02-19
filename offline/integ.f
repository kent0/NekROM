c-----------------------------------------------------------------------
      subroutine setops(a,b,c,u,nel,nb,ndim)

      real a(1),b(1),c(1),u(1)

      call bip(b,u,u,nel,nb,ndim)

      do j=1,nb
      do i=1,nb
c        write (6,*) i,j,b(i+(j-1)*nb),'b'
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine bip(b,u,v,nel,nb,ndim)

      include 'LVAR'
      include 'INTEG'

      common /tinteg/ tmp(lxyz,lel,ldim)

      real b(nb,nb)
      real u(lxyz,nel,nb,ndim),v(lxyz,nel,nb,ndim)

      do l=1,ndim
      do k=1,nb
         call col3(tmp,u(1,1,k,l),bm1,lxyz*nel)
         do j=1,nb
            b(j,k)=b(j,k)+vlsc2(tmp,v(1,1,j,l),lxyz*nel)
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine aip(a,u,v,nb,n)

      real a(nb,nb)
      real u(n,nb),v(n,nb)

      do k=1,nb
      do j=1,nb
      do i=1,n
         a(j,k)=a(j,k)+u(i,k)*v(i,j)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cip(c,u,v,w,nb,n)

      real c(nb,nb,nb)
      real u(n,nb),v(n,nb),w(n,nb)

      do l=1,nb
      do k=1,nb
      do j=1,nb
      do i=1,n
         c(j,k,l)=c(j,k,l)+u(i,l)*v(i,j)*w(i,k)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------