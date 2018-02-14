# Installing OpenSeesComposite

There are two steps: installing the C++ code and installing the supporting Tcl scripts.

## C++ packages

Download the package source from Github and put them in the correct directory.

```sh
cd $HOME/OpenSees/DEVELOPER
git clone https://github.com/denavit/OpenSees-Packages.git CompositePackages
```

You'll need to edit the `Makefile.def` in `DEVELOPER` (not the same `Makefile.def` used for compiling OpenSees) to add where you want the compiled shared libraries to be put. It's common to keep these in the OpenSees `bin` or `lib` folders. For example:

```Makefile
COMPOSITE_LIB_BIN = $(HOME)/OpenSees/bin
```

Run `make` in the `CompositePackages` folder.

Finally, you'll need to set your `LD_LIBRARY_PATH` variable to include that same path [1]. For example, add this line to your `~/.profile`:

```sh
export LD_LIBRARY_PATH="$HOME/OpenSees/bin:$LD_LIBRARY_PATH"
```

## Tcl scripts

Tcl scripts that provide defaults and useful functions (especially for section creation) are implemented in the Tcl package `OpenSeesComposite`. They are available [here](https://github.com/denavit/OpenSees-Tcl-Scripts). The Tcl scripts are not necessary for using the packages, but are extremely useful.

[1] `LD_LIBRARY_PATH` must be used even though Ubuntu has moved away from using it. Using `/etc/ld.so.conf.d` files will NOT work.