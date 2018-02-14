# Installing OpenSees-Packages

Download the package source from Github and put it in `OpenSees/DEVELOPER/CompositePackages`.

```sh
cd $HOME/OpenSees/DEVELOPER
git clone https://github.com/denavit/OpenSees-Packages.git CompositePackages
```

You'll need to edit your `Makefile.def` (in the base directory; not the one in `DEVELOPER`) to add the variable `COMPOSITE_LIB_BIN`. This sets where you want the compiled shared libraries to be put. It's common to keep the compiled libraries in the OpenSees `bin` or `lib` folders. For example:

```Makefile
COMPOSITE_LIB_BIN = $(HOME)/bin
```

Run `make` in the `CompositePackages` folder. This should produce the following files in the directory specified by `COMPOSITE_LIB_BIN`:

```
changManderConcrete01.so
mixedBeamColumn2d.so
mixedBeamColumn3d.so
multiSurfaceKinematicHardening.so
ratchet.so
sakinoSunConcrete04.so
shenSteel01.so
```

Finally, you'll need to set your `LD_LIBRARY_PATH` variable to include `$COMPOSITE_LIB_BIN`. For example, add this line to your `~/.profile`:

```sh
export LD_LIBRARY_PATH="$HOME/bin:$LD_LIBRARY_PATH"
```

Note that `LD_LIBRARY_PATH` must be used even though Ubuntu has moved away from using it. Using `/etc/ld.so.conf.d` files will NOT work.

## Tcl scripts

Tcl scripts that provide defaults and useful functions (especially for section creation) are implemented in the Tcl package `OpenSeesComposite`. They are available [here](https://github.com/denavit/OpenSees-Tcl-Scripts). The Tcl scripts are not necessary for using the packages, but are extremely useful.
