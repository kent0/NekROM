c-----------------------------------------------------------------------
        program offline

        include 'LVAR'
        include 'OFFLINE'

        common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
        character*132 fname

        call offline_init(icomm)

        nsg=lsg
        nsg=5

        ie=1

        nel=2
        mel=1
        ng=neg/mel
        ng=1

        do i=1,ng
           write (6,*) i,'offline'
           ieg(1)=ie
           ieg(2)=ie+mel-1
           nel=ieg(2)-ieg(1)+1

           call loadsnaps(buf,ieg,indxr,nsg)
           call setgeom(buf,nel)
           call setops(ga,gb,gc,gt,buf(nel*lxyz*ldim+1),nel,nsg,ldim)
           ie=ie+mel
        enddo

        call gop(ga,gt,'+  ',nsg*nsg)
        call gop(gb,gt,'+  ',nsg*nsg)

        call dump_serial(ga,nsg*nsg,'ga ',mid)
        call dump_serial(gb,nsg*nsg,'gb ',mid)

        return
        end
c-----------------------------------------------------------------------