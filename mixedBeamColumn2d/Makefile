include ../../../Makefile.def

OS_FLAG ?= -D_LINUX

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY         = $(COMPOSITE_LIB_BIN)/mixedBeamColumn2d.so

all:   $(LIBRARY)

OBJS = mixedBeamColumn2d.o 

$(LIBRARY): $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

mixedBeamColumn2d.o: mixedBeamColumn2d.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall mixedBeamColumn2d.cpp

else

LIBRARY         = mixedBeamColumn2d.dylib


all: 
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace mixedBeamColumn2d.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)

endif

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
