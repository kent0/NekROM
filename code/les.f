c-----------------------------------------------------------------------
      subroutine set_les_imp(fles1,fles2)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /mytmp/ dfld(lt),t1(lt)

      real fles1(0:nb,nb),fles2(0:nb,nb)

      nv=lx1*ly1*lz1*nelv
      nt=lx1*ly1*lz1*nelv

      call push_sol(vx,vy,vz,pr,t)

      do i=1,nb
         if (ifrom(1)) then
            call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
            call opadd2(vx,vy,vz,ub,vb,wb)
         endif

         if (ifrom(2)) then
            call copy(t,tb(1,i),nt)
            call add2(t,tb,nt)
         endif

         call q_filter(1.)
         if (ifrom(1)) call pv2b(fles1(0,i),vx,vy,vz,ub,vb,wb)
         if (ifrom(2)) call ps2b(fles2(0,i),t,vy,vz,tb)
      enddo

      do i=1,nb
         if (ifrom(1)) fles1(i,i)=fles1(i,i)-1.
         if (ifrom(2)) fles2(i,i)=fles2(i,i)-1.
      enddo

      call pop_sol(vx,vy,vz,pr,t)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_les_exp(fles1,fles2)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /mytmp/ dfld(lt),t1(lt)

      real fles1(0:nb,0:nb),fles2(0:nb,0:nb)

      n=lx1*ly1*lz1*nelv

      if (ifrom(1)) then
         call rzero(fles1,(nb+1)**2)
         do i=0,nb
            call evalnut(dfld,ub(1,i),vb(1,i),wb(1,i),.false.)
            do j=0,nb
               call axhelm(t1,dfld,zeros,1,1)
               fles1(i,j)=fles1(i,j)+glsc2(t1,ub(1,j),n)
               call axhelm(t1,dfld,zeros,1,2)
               fles1(i,j)=fles1(i,j)+glsc2(t1,vb(1,j),n)
               if (ldim.eq.3) then
                  call axhelm(t1,dfld,zeros,1,3)
                  fles1(i,j)=fles1(i,j)+glsc2(t1,wb(1,j),n)
               endif
            enddo
         enddo
      endif

      if (ifrom(2)) then
         call rzero(fles2,(nb+1)**2)
         do i=0,nb
            call evalnut(dfld,ub(1,i),vb(1,i),wb(1,i),.true.)
            do j=0,nb
               call axhelm(t1,dfld,zeros,1,1)
               fles2(i,j)=fles2(i,j)+glsc2(t1,ub(1,j),n)
            enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_les_imp(uu,tt,sig,fles1,fles2,tmp)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /mytmp/ dfld(lt),t1(lt)

      real uu(0:nb),tt(0:nb),fles1(0:nb,nb),fles2(0:nb,nb),tmp(0:nb)

      do i=1,nb
         if (ifrom(1)) then
            call mxm(fles1,nb+1,uu(1),nb,tmp,1)
            tmp(0)=0.
            call add2s1(tmp,uu,-sig,nb+1)
         endif
         if (ifrom(2)) then
            call mxm(fles2,nb+1,tt(1),nb,tmp,1)
            tmp(0)=0.
            call add2s1(tmp,tt,-sig,nb+1)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_les_exp(uu,fles1,fles2)

      return
      end
c-----------------------------------------------------------------------
      subroutine evalnut(nut,u1,u2,u3,ifmol)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      real nut(lt),u1(lt),u2(lt),u3(lt)

      common /scrns/ du(lt),t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)
      common /scruz/ t7(lt),t8(lt),t9(lt)

      logical ifmol

      n=lx1*ly1*lz1*nelv

      call gradm1(t1,t2,t3,u1)
      call gradm1(t4,t5,t6,u2)
      if (ldim.eq.3) call gradm1(t7,t8,t9,u3)

      do i=1,n
         if (ldim.eq.2) then
            du(i)=sqrt(t1(i)*t1(i)+t2(i)*t2(i)+t3(i)*t3(i)
     $                +t4(i)*t4(i)+t5(i)*t5(i)+t6(i)*t6(i))
         else
            du(i)=sqrt(t1(i)*t1(i)+t2(i)*t2(i)+t3(i)*t3(i)
     $                +t4(i)*t4(i)+t5(i)*t5(i)+t6(i)*t6(i)
     $                +t7(i)*t7(i)+t8(i)*t8(i)+t9(i)*t9(i))
         endif
      enddo

      if (ifmol) call cnvlphi(du)

      do i=1,n
         cs=1.
         hk=1.
         nut(i)=(cs*hk)**2*du(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cnvl(du)

      include 'SIZE'
      include 'TOTAL'

      real du(1)

      ! TODO: implement convolutions

      return
      end
c-----------------------------------------------------------------------
