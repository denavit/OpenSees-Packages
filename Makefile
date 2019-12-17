include ../../Makefile.def

all: 
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/changManderConcrete01; $(MAKE);
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/sakinoSunConcrete04; $(MAKE);
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/shenSteel01; $(MAKE);
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/multiSurfaceKinematicHardening; $(MAKE);
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/ratchet; $(MAKE);
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/mixedBeamColumn2d; $(MAKE);
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/mixedBeamColumn3d; $(MAKE);
    @$(CD) $(FE)/../DEVELOPER/CompositePackages/inerterTruss2d; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/changManderConcrete01; $(MAKE) wipe;
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/sakinoSunConcrete04; $(MAKE) wipe;
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/shenSteel01; $(MAKE) wipe;
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/multiSurfaceKinematicHardening; $(MAKE) wipe;
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/ratchet; $(MAKE) wipe;
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/mixedBeamColumn2d; $(MAKE) wipe;
	@$(CD) $(FE)/../DEVELOPER/CompositePackages/mixedBeamColumn3d; $(MAKE) wipe;
    @$(CD) $(FE)/../DEVELOPER/CompositePackages/inerterTruss2d; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
