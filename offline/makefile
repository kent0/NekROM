### makefile automatically created by makenek 02/10/2020 16:59:58 ###
MD=/Users/kaneko/Developer/MOR
OD=$(MD)/offline

BINNAME=offline
LIBNAME=libo.a
CASENAME:=h3d
CASEDIR:=$(OD)
S=/Users/kaneko/Developer/Nek5000/core
SCMT=/Users/kaneko/Developer/Nek5000/core/cmt
OPT_INCDIR=./
OBJDIR=oobj

FC=mpif77
CC=mpicc
FFLAGS= -w -std=legacy -fdefault-real-8 -fdefault-double-8 -cpp
CFLAGS= -w
LDFLAGS=

PPPO=
PPS= MPI UNDERSCORE TIMER

USR=io.o integ.o geom.o misc.o util.o
USR_LFLAGS=-mcmodel=medium -L/Users/kaneko/Developer/Nek5000/3rd_party/gslib/lib -lgs

MPI=1
CMT=0
DPROCMAP=0
VISIT=0

VISIT_IFLAGS=
VISIT_LFLAGS=
HYPRE_IFLAGS=
GSLIB_IFLAGS:=-DPREFIX=gslib_ -DFPREFIX=fgslib_ -DGLOBAL_LONG_LONG -I/Users/kaneko/Developer/Nek5000/3rd_party/gslib/include

################################################################################

PPS_F = $(patsubst %,$(PPPO)-D%,$(PPS))
PPS_C = $(patsubst %,-D%,$(PPS))

# CORE##########################################################################
CORE = byte.o byte_mpi.o chelpers.o speclib.o
################################################################################
# MXM
MXM=mxm_std.o
################################################################################
# BLAS
BLAS=blas.o dsygv.o

COMM_MPI := comm_mpi.o
ifeq ($(MPI),0)
  COMM_MPI := ${COMM_MPI} mpi_dummy.o
endif

# CMT append *.o files to CORE
CMTCORE :=

DUMMY:= $(shell cp $S/PARALLEL_ $S/PARALLEL)
ifeq ($(DPROCMAP),1)
	CORE := ${CORE} dprocmap.o
	DUMMY:= $(shell cd $S; patch < PARALLEL.patch)
endif

ifneq ($(VISIT),0)
  VISITNEK_INCLUDES:=-I$S/3rd_party
  VISITO := visit.o
  USR_LFLAGS := ${USR_LFLAGS} $(VISIT_LFLAGS)
endif

ifeq ($(MPI),0)
	DUMMY:= $(shell cp $S/mpi_dummy.h $S/mpif.h)
else
	DUMMY:= $(shell rm -rf $S/mpif.h)
endif

TMP1 = $(CORE) $(MXM) $(USR) $(COMM_MPI) $(VISITO) $(BLAS)

ifeq ($(CMT),1)
TMP1 += $(CMTCORE)
OPT_INCDIR += -I $(SCMT)
endif

OPT_INCDIR += -I $S/experimental

NOBJS = $(patsubst %,$(OBJDIR)/%,$(TMP1))

L0 = $(G) -O0
L2 = $(G) -O2
L3 = $(G) -O3

FL0   = $(L0) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR)
FL2   = $(L2) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR)
FL3   = $(L3) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR)

cFL0   = $(L0) $(CFLAGS) $(PPS_C) -I$S -I$(OPT_INCDIR)
cFL2   = $(L2) $(CFLAGS) $(PPS_C) -I$S -I$(OPT_INCDIR)
cFL3   = $(L3) $(CFLAGS) $(PPS_C) -I$S -I$(OPT_INCDIR)

################################################################################
all : offline

objdir:
	@mkdir $(OBJDIR) 2>/dev/null; cat /dev/null

usrfile:
	@env CASENAME=$(CASENAME) PPS="$(PPS)" $S/mkuserfile

offline: lib usrfile
	$(FC) -I$(OD)/ $(FL2) -o ${BINNAME} drive.f $(CASEDIR)/${CASENAME}.f $(VISITNEK_INCLUDES) $(OBJDIR)/${LIBNAME} $(USR_LFLAGS) $(LDFLAGS)
	@if test -f ${BINNAME}; then \
        echo ""; \
        echo "$(NEK_WARN)"; \
        echo ""; \
	echo "#############################################################"; \
	echo "#                  Compilation successful!                  #"; \
	echo "#############################################################"; \
	size ${BINNAME}; \
	echo ""; \
	else \
	echo -e "\033[1;31;38m" "ERROR: Compilation failed!"; \
	echo -e "\033[0m"; \
	fi
	@rm -f *.i
ifeq ($(MPI),0)
	@rm -rf $S/mpif.h
endif

lib:	objdir $(NOBJS)
	@$(AR) cru $(OBJDIR)/${LIBNAME} $(NOBJS)
	@ranlib $(OBJDIR)/${LIBNAME}

clean:
	@echo "cleaning offline ..."
	@rm -rf ${CASENAME}.f  $(OBJDIR) ${BINNAME}
ifeq ($(MPI),0)
	@rm -rf $S/mpif.h
endif

$(NOBJS) : $(OD)/SIZE

# NEK CORE     #################################################################
$(OBJDIR)/dprocmap.o    :$S/dprocmap.f $S/DPROCMAP;     $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/comm_mpi.o	:$(OD)/comm_mpi.f;	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mpi_dummy.o	:$S/mpi_dummy.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/byte_mpi.o	:$S/byte_mpi.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/speclib.o	:$S/speclib.f;			$(FC) -c $(FL2) $< -o $@

# MXM       ####################################################################
$(OBJDIR)/mxm_wrapper.o	  :$S/mxm_wrapper.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mxm_std.o	  :$S/mxm_std.f;		$(FC) -c $(FL3) $< -o $@
$(OBJDIR)/mxm_bgq.o	  :$S/mxm_bgq.f;		$(FC) -c $(FL3) $< -o $@

# C Files ######################################################################
$(OBJDIR)/nek_comm.o         :$S/nek_comm.c;          $(CC) -c $(cFL2) $< -o $@
$(OBJDIR)/byte.o             :$S/byte.c;              $(CC) -c $(cFL2) $< -o $@
$(OBJDIR)/chelpers.o         :$S/chelpers.c;          $(CC) -c $(cFL2) $< -o $@
$(OBJDIR)/fcrs.o             :$S/fcrs.c;              $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $(HYPRE_IFLAGS) $< -o $@
$(OBJDIR)/crs_xxt.o          :$S/crs_xxt.c;           $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $< -o $@
$(OBJDIR)/sparse_cholesky.o  :$S/sparse_cholesky.c;   $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $< -o $@
$(OBJDIR)/crs_amg.o          :$S/crs_amg.c;           $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $< -o $@
$(OBJDIR)/fem_amg_preco.o    :$S/experimental/fem_amg_preco.c; $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $(HYPRE_IFLAGS) $< -o $@
$(OBJDIR)/crs_hypre.o        :$S/experimental/crs_hypre.c; $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $(HYPRE_IFLAGS) $< -o $@

# 3rd party ####################################################################
$(OBJDIR)/dsygv.o      :$S/3rd_party/dsygv.f;        $(FC) -c $(L0) $(PPS_F) $< -o $@
$(OBJDIR)/blas.o       :$S/3rd_party/blas.f;	      $(FC) -w -c $(L0) $(PPS_F) $< -o $@
$(OBJDIR)/visit.o      :$S/3rd_party/visit.f;	      $(FC) -c $(VISIT_IFLAGS) $(FL2) $< -o $@
$(OBJDIR)/nek_in_situ.o:$S/3rd_party/nek_in_situ.f; $(FC) -c $(VISIT_IFLAGS) $(FL2) $< -o $@
$(OBJDIR)/finiparser.o :$S/3rd_party/finiparser.c;  $(CC)  -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/iniparser.o  :$S/3rd_party/iniparser.c;   $(CC)  -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/dictionary.o :$S/3rd_party/dictionary.c;  $(CC)  -c $(cFL2) $(JL) $< -o $@

# My Files #####

$(OBJDIR)/io.o    :$(OD)/io.f    LVAR1; $(FC) -c -I$(OD)/ $(FL2) $< -o $@
$(OBJDIR)/integ.o :$(OD)/integ.f LVAR1; $(FC) -c -I$(OD)/ $(FL2) $< -o $@
$(OBJDIR)/geom.o  :$(OD)/geom.f  LVAR1; $(FC) -c -I$(OD)/ $(FL2) $< -o $@
$(OBJDIR)/util.o  :$(OD)/util.f  LVAR1; $(FC) -c -I$(OD)/ $(FL2) $< -o $@
$(OBJDIR)/misc.o  :$(OD)/misc.f  LVAR1; $(FC) -c -I$(OD)/ $(FL2) $< -o $@
