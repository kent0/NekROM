      ! Auxiliary

      parameter (lsu=(lcs-1)*lmu+1) ! size of velocity snapshots allocation
      parameter (lsp=(lcs-1)*lmp+1) ! size of pressure snapshots allocation
      parameter (lst=(lcs-1)*lmt+1) ! size of temperature snapshots allocation


      parameter (lub=(lb-1)*lmu+1) ! size of velocity basis allocation
      parameter (lpb=(lb-1)*lmp+1) ! size of pressure basis allocation
      parameter (ltb=(lb-1)*lmt+1) ! size of temperature basis allocation

      parameter (lres=1) ! size of residual storage
      parameter (lres_u=((3*lb+lb**2)-1)*lei+1) ! size of residual storage for vel
      parameter (lres_t=((2*lb+lb**2)-1)*lei+1) ! size of residual storage for temp
      parameter (lelm=lelt)   ! number of local elements
      parameter (ltr=1+0*399) ! max number of tensor rank

      parameter (lcs=ls)  ! max number of coefficient set
