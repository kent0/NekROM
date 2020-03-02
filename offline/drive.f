c-----------------------------------------------------------------------
      program offline

      include 'LVAR'
      include 'OFFLINE'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      character*132 fname
      logical ifavg0,iftherm

      call offline_init(icomm)

      nsg=3
      neg=512
      mel=2

      iftherm=.true.
      ifavg0=.true.

      call gengrams(ga,gb,gc,gt,buf,ieg,indxr,
     $  nsg,mp,neg,mel,ldim,lxyz,iftherm)

      call write_ops(ga,gb,gc,nsg,mid,' gu',.true.)

      if (iftherm) call write_ops(ga((nsg+1)**2+1),gb((nsg+1)**2+1),
     $   gc(((nsg-1)/mp+2)*(nsg+1)**2+1),nsg,mid,' gt',.true.)

      call setg(gg,gb,gt,nsg,ifavg0)
      call setq(gvec,gvect,gvecc,gval,gg,nsg,nsc,mp,mid,ifavg0,nsg1)

      call qop2(ga,gt,gvec,gvect,nsg,nsg1)
      call qop2(gb,gt,gvec,gvect,nsg,nsg1)
      call write_ops(ga,gb,gc,nsg1,mid,'  u',.false.)

      call qop3(gc,gt,gvec,gvect,gvecc,nsg,nsg1,nsc,mid,'  u')

      if (iftherm) then
         call setg(gg,gb((nsg+1)**2+1),gt,nsg,ifavg0)
         call setq(gvec,gvect,gvecc,gval,gg,nsg,nsc,mp,mid,ifavg0,nsg1)

         call qop2(ga((nsg+1)**2+1),gt,gvec,gvect,nsg,nsg1)
         call qop2(gb((nsg+1)**2+1),gt,gvec,gvect,nsg,nsg1)
         call write_ops(ga((nsg+1)**2+1),gb((nsg+1)**2+1),
     $      gc(((nsg-1)/mp+2)*(nsg+1)**2+1),nsg1,mid,'  t',.false.)

         call qop3(gc(((nsg-1)/mp+2)*(nsg+1)**2+1),
     $        gt,gvec,gvect,gvecc,nsg,nsg1,nsc,mid,'  t')
      endif

      return
      end
c-----------------------------------------------------------------------