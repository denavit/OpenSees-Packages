include ../../../Makefile.def


#	$(FE)/matrix/Matrix.o \
#	$(FE)/matrix/ID.o \
#	$(FE)/nDarray/nDarray.o \
#	$(FE)/nDarray/BJtensor.o \
#	$(FE)/nDarray/BJmatrix.o \
#	$(FE)/nDarray/basics.o \
#	$(FE)/material/uniaxial/UniaxialMaterial.o \
#	$(FE)/element/Information.o \
#	$(FE)/material/Material.o \
#	$(FE)/actor/actor/MovableObject.o \
#	$(FE)/recorder/response/MaterialResponse.o \
#	$(FE)/recorder/response/Response.o \
#	$(FE)/tagged/TaggedObject.o

OS_FLAG ?= -D_LINUX

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY         = $(COMPOSITE_LIB_BIN)/changManderConcrete01.so

OUTSIDE_OBJS = 


OBJS = changManderConcrete01.o

all: $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

changManderConcrete01.o: changManderConcrete01.cpp
	$(CC++) -fPIC $(INCLUDES)  -c -Wall changManderConcrete01.cpp

else

LIBRARY         = changManderConcrete01.dylib

OUTSIDE_OBJS = $(FE)/../PACKAGES/dynamicLoadFunctions.o


all: 
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace changManderConcrete01.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)

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
