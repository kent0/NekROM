c-----------------------------------------------------------------------
        program offline

        include 'LVAR'
        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        nsg=lsg
        ieg(1)=1
        ieg(2)=10
        nel=ieg(2)-ieg(1)+1

        call loadsnaps(buf,ieg,indxr,nsg)
        call setgeom(buf,nel)
        call geom_check(nel)
        call setops(ga,gb,gc,buf(nel*lxyz*ldim+1),nel,nsg,ldim)

        return
        end
c-----------------------------------------------------------------------