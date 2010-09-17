include ../../Makefile.def

all: 
	@$(CD) $(FE)/../PACKAGES/CompositePackages/changManderConcrete01; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/CompositePackages/sakinoSunConcrete04; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/CompositePackages/shenSteelCCFT; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/CompositePackages/shenSteelRCFT; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/CompositePackages/mixedBeamColumn2d; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/CompositePackages/mixedBeamColumn3d; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(CD) $(FE)/../PACKAGES/CompositePackages/changManderConcrete01; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/CompositePackages/sakinoSunConcrete04; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/CompositePackages/shenSteelCCFT; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/CompositePackages/shenSteelRCFT; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/CompositePackages/mixedBeamColumn2d; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/CompositePackages/mixedBeamColumn3d; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
