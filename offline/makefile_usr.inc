MOR_DIR=/Users/kaneko/Developer/MOR
$(OBJDIR)/aux.o :$(MOR_DIR)/offline/aux.f OFFLINE; $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/io.o  :$(MOR_DIR)/offline/io.f  IO OFFLINE; $(FC) -c $(FL2) $< -o $@
