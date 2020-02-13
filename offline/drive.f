c-----------------------------------------------------------------------
        program offline

        include 'LVAR'
        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        nsg=lsg
        nel=512
        call loadflist(fnames,nsg)

        call setindxr(indxr,ldimt)

        do is=1,nsg
           write (6,'(a132)') fnames(is)
        enddo

        ns=1

        mel=nel
        ng=1

        mel=nel/16
        ng=16

        call rzero(gb,ns*ns)

        iloc=1

        is0=1
        is1=nsg

        do ig=1,ng
        do is=is0,is1
           ieg(1)=mel*(ig-1)+1
           ieg(2)=mel*ig
           call rxupt(buf(iloc),ieg,indxr,fnames(is))
           j=1
           nloc=0
        enddo
        enddo

        call setb(gb,buf,buf,ns,mel*lxyz*ldim)
        write (6,*) ig,gb(1),buf(1),'wp deb'

        end
c-----------------------------------------------------------------------