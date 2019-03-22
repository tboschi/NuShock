# NuShock

NuShock is a flexible framework to study the phenomenology of models with heavy neutrinos.
Given some description of the physics processes (cross-sections or decay rates), 
the neutrino flux, and the detector geometry, it is possible to study the sensitivity 
to a specific process of a neutrino experiment.
Tools for background estimations are also available.

## Installation

The installation requires only the root package (version 5 or newer)
and its environmental variables already set.
The standard c++11 is also preferable, even though compatibility with c++0x should be ok.

```bash
make
```

If your g++ or preferred c++ compiler does not support c++11 please change the Makefile accordingly, 
replacing the flag 
```
-std=c++11
```
with 
```
-std=c++0x
```

Recompile with make everytime you change something in the source code or if you create new executables.

## Usage of the framework

The main executables are placed in the app/ folder.
They can use any of the library inside the include/src/ folder.
It is best to use the header files in include/ to make the code clean and readable.
If you create a new class which is not included in the header files under include/,
please change them or properly include them in the pre-compiler instructions.

At the moment only background, detector, and some tools work.
The rest is still in development as I recently changed a few things.

### Configuration files

Most of the classes, like flux or detector, use description files contained in the config/ subfolder.

The detector configuration, for example ND_Full, described the DUNE near detector, composed of two parts:
a LArTPC and a FGT.
Sizes and dimensions are desribed in the file, as well as resolutions for some particle detections.
At the end there is a list of efficiency files, used to estimate sensitivty including background.

More description inside the configuration file itself.

At the moment the detector class expects a detector made of two parts: a LArTPC in front of a TPC,
but this problem can be overcome by defining a null secondary detector.
The geometry are boxes, so this should be maybe changed.
In the Tracker class, correct material properties should be defined, like interaction length etc.
At the moment olny liquid argon, gaseous argon, iron and lead are defined.
Material description is important for energy loss in material: the length of the track defines
the resolution to use, if it is a contained track or not.
In the case of DUNE ND, a track is contained if 80% of the total length are inside.
This can be changed with the Containment entry in the config file.
