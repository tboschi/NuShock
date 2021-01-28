# NuShock

NuShock is a flexible framework to study the phenomenology of models with heavy neutrinos.
Given some description of the physics processes (cross-sections or decay rates), 
the neutrino flux, and the detector geometry, it is possible to study the sensitivity 
to a specific process of a neutrino experiment.
Tools for background estimations are also available.

## Installation

### Requirments

* make
* C++11 compiler
* [cuba](http://www.feynarts.de/cuba/)
* [LHAPDF](https://lhapdf.hepforge.org/index.html)
* [ROOT](https://root.cern.ch/)

Set the variables:
* `$CUBA` to point to the installation path for cuba
* `$LHAPDF` to point to the installation path for LHAPDF
* `$ROOTSYS` to point to the installation path for ROOT (this already part of the installation process)

If ROOT is properly set, `root-config` should be in the user's `$PATH` and so the Makefile can find the includes and libraries on its own.
Remember to update your `$LD_LIBRARY_PATH`.

### Compilation

Simply doing 
```
make
```
is enough.
This will create binaries in the `bin` folder.

# WARNING: this readme is deprecated and will be updated soon!
## Usage of the framework

The main functions are to be placed in the `app/` folder and the exectuables are moved to the `bin/` folder after compilation.
These are, at the moment


### Configuration files

Most of the classes, like flux or detector, use description files contained in the `config/` subfolder.

The detector configuration, for example `config/detectro.card` describes the DUNE near detector, composed of two parts: a LArTPC and a FGT.
Sizes and dimensions are desribed in the file, whereas the tracking resolutions are in `config/resolution.card`.
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

The flux configuration files expects paths to ROOT files containing TH1D with the fluxes.
As an extra, modifier files for tau flux, created from direct simulation of HNL production from Ds decays.

### Flux

The flux package has a matrioska-like structure.
The lowest level class is `Flux` which copies to memory TH1D from flux configuration file.
Above that the `Driver` class takes care of scaling the light-nu fluxes to HNL fluxes.
The better way to use the `Driver` class is with the `Engine` class, which binds flux, detector and neutrino (HNL) type
together to evaluate the flux AT the ND site and therefore evaluation of number of events is easy.

The number of events formula is split in two parts: flux+energy integration and probability+efficiency.
The latter is done in the `Detector` class (see below).
The flux and integration over energy E is handled by `Engine` which calls the `Detector` at each energy to 
estimate the probability+efficiency of an HNL of given energy E to decay inside the ND.

### Physics

This is the most important package, in which all decay widths and distributions are compute.
It is based on a base class, `Amplitude`, which has definition of the Feynman amplitudes, M2,
of all decay and production modes taken in consideration.
The `Amplitude` class is never called alone, but it is accesed by derived classes:
- `DecayRates`	for computing decay rates and branching ratios
- `Production`	for computing production scale factors
- `PhaseSpace`	to be used for simulation of decays and production 

The derived classes are best implemented in the `Neutrino` class, which in turn is derived from
the `Particle` class (`tools.h` package).
A neutrino object can be constructed on the stack as it doesn't require much memory and has methods that returns
decay widths, branching ratios, scale factors, phase space generation, etc.
When constructing a neutrion object, options can be passed to specify helicity and fermionic nature.
The options are
- `Neutrino::Left`		for a left-helicity particle
- `Neutrino::Right`		for a right-helicity particle
- `Neutrino::Unpolarised`	for an unpolarised particle, so information on helicity is integrated
- `Neutrino::Dirac`		for a Dirac HNL
- `Neutrino::Majorana`		for a Majorana HNL
- `Neutrino::Antiparticle`	for an anti HNL; it makes sense only if HNL is Dirac

It the Antiparticle option is not specified, then the HNL is a particle by default.


In order to call a specific prodution or decay mode, a unique string is passed that identifies the channel of interest.

These strings are for decay channels:
- `nnn`	  	for 3-body decay into 3 light neutrinos
- `nGAMMA`	for radiative decay into a neutrino and a gamma (deprecated now)
- `nEE`		for 3-body decay into nu e e
- `nEM`		for 3-body decay into nu e mu
- `nMM`		for 3-body decay into nu mu mu
- `nET`		for 3-body decay into nu e tau
- `nMT`		for 3-body decay into nu mu tau
- `nPI0`	for 2-body decay into nu pion0
- `EPI`		for 2-body decay into e pion
- `MPI`		for 2-body decay into mu pion
- `TPI`		for 2-body decay into tau pion
- `EKA`		for 2-body decay into e kaon
- `MKA`		for 2-body decay into mu kaon
- `nRHO0`	for 2-body decay into nu rho0
- `ERHO`	for 2-body decay into e rho
- `MRHO`	for 2-body decay into mu rho
- `EKAx`	for 2-body decay into e kaon\*
- `MKAx`	for 2-body decay into mu kaon\*
- `nOMEGA`	for 2-body decay into nu omega
- `nETA`	for 2-body decay into nu eta
- `nETAi`	for 2-body decay into nu eta'
- `nPHI`	for 2-body decay into nu phi
- `ECHARM`	for 2-body decay into e D (charm meson)
- `ExpALL`	for decay with good discovery sensitivity (nEE, nEM, nMM, nPI0, EPI, MPI)

These strings are for production channels:
-  `MuonE`	for 3-body production from muon via Ue mixing
-  `MuonM`	for 3-body production from muon via Um mixing
-  `TauEE`	for 3-body production from tau to electron via Ue mixing
-  `TauET`	for 3-body production from tau to electron via Ut mixing
-  `TauMM`	for 3-body production from tau to muon via Um mixing
-  `TauMT`	for 3-body production from tau to muon via Ut mixing
-  `TauPI`	for 2-body production from tau to pion
-  `Tau2PI`	for 3-body production from tau to 2 pions (only phase space!)
-  `PionE`	for 2-body production from pion to electron
-  `PionM`	for 2-body production from pion to muon
-  `KaonE`	for 2-body production from kaon to electron
-  `KaonM`	for 2-body production from kaon to muon
-  `CharmE`	for 2-body production from D (charm) to electron
-  `CharmM`	for 2-body production from D (charm) to muon
-  `CharmT`	for 2-body production from D (charm) to tau
-  `Kaon0E`	for 3-body production from kaon0 to pion electron
-  `Kaon0M`	for 3-body production from kaon0 to pion muon
-  `KaonCE`	for 3-body production from kaon to pion0 electron
-  `KaonCM`	for 3-body production from kaon to pion0 muon

### Detector

The `Detector` class loads geometry from configuration file.
There are also a `Tracker` and `Efficiency` class, derived from `Detector`, to deal with background generation/smearing
and detector efficiency computation.

### Examples

The executables for production scale and decay branch can be run as

```
./bin/DecayBranch -o outputfile1 -c nEM -E 1.0 -M 0.5
./bin/ProductionScale -o outputfile2 -c Kaon0E -E 1.0 -M 0.5
```

The output files will contain two columns: the first ones being mass of HNL and the second ones the branching ratio of
the decay of HNL into nEM and production scale factor for Kaon0E with mixing ratio 2:1:0.

Another example
```
./bin/Sensitivity -d config/ND_Full_2 -f config/FluxConfig -c EPI -E -o output
./bin/Sensitivity -d config/ND_Full_2 -f config/FluxConfig -c MPI -M -o output
```
which compute the experimental sensitivity to the channels EPI and MPI for zero background expectation.
The lines are computed using bisection search method implemented in the `Exclusion` class from the `analysis` package.
