D1=$(MOR_DIR)/code/MOR LMOR
D2=$(MOR_DIR)/code/MORDICT $(MOR_DIR)/code/MOR LMOR

F1=-I$(MOR_DIR)/code

$(OBJDIR)/ana.o    :$(MOR_DIR)/code/ana.f    $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/aux.o    :$(MOR_DIR)/code/aux.f    $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/conv.o   :$(MOR_DIR)/code/conv.f   $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/const.o  :$(MOR_DIR)/code/const.f  $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/dump.o   :$(MOR_DIR)/code/dump.f   $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/ei.o     :$(MOR_DIR)/code/ei.f     $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/filter.o :$(MOR_DIR)/code/filter.f $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/lapack.o :$(MOR_DIR)/code/lapack.f $D1; $(FC) -c $(FL0) $F1 $< -o $@
$(OBJDIR)/legacy.o :$(MOR_DIR)/code/legacy.f $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/pod.o    :$(MOR_DIR)/code/pod.f    $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/qoi.o    :$(MOR_DIR)/code/qoi.f    $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/read.o   :$(MOR_DIR)/code/read.f   $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/riesz.o  :$(MOR_DIR)/code/riesz.f  $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/rk.o     :$(MOR_DIR)/code/rk.f     $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/rom.o    :$(MOR_DIR)/code/rom.f    $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/time.o   :$(MOR_DIR)/code/time.f   $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/tensor.o :$(MOR_DIR)/code/tensor.f $D1; $(FC) -c $(FL2) $F1 $< -o $@
$(OBJDIR)/mpar.o   :$(MOR_DIR)/code/mpar.f   $D2; $(FC) -c $(FL2) $F1 $< -o $@

$(OBJDIR)/unit.o   :$(MOR_DIR)/tests/unit.f  $D1; $(FC) -c $(FL2) $F1 $< -o $@
