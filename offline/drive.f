c-----------------------------------------------------------------------
        program offline

        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        call blank(fname,132)
        fname='/Users/kaneko/Developer/MOR/offline/r0.f00001 '

        nel=512

        ieg(1)=1
        ieg(2)=nel

        call setifread(ifread,ldimt)
        write (6,*) 'upt(1)',upt(1)
        call rxupt(xyz,upt,ieg,ifread,fname)
        write (6,*) 'upt(1)',upt(1)
        call setb(gb,upt,upt,1,512*lxyz)

        write (6,*) 'gb(1)',gb(1)

        end
c-----------------------------------------------------------------------