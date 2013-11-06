Agent-Based Plasmid Simulations
===============================

This is the implementation of the evolutionary agent-based model of
hosts and plasmids used in the following research paper :

K. Kentzoglanakis, D. G. Lopez, S. P. Brown, and R. A. Goldstein. [The
Evolution of Collective Restraint: Policing and Obedience among
Non-conjugative
Plasmids.](http://dx.plos.org/10.1371/journal.pcbi.1003036) PLoS
Comput Biol, 9(4):e1003036+, 2013.



## REQUIREMENTS

A simple `Makefile` is included. Run `make` in the dps directory in
order to compile the code. `dps` depends upon a few external libraries
which should be present in your system along with their respective
header (development) files. The prerequisites for compiling `dps` are
as follows:

1. [GNU GCC](http://gcc.gnu.org) >= 4.4.5

2. [GNU Scientific Library](http://www.gnu.org/software/gsl) >= 1.15

3. [HDF5](http://www.hdfgroup.org/HDF5) 

4. [GLib](https://developer.gnome.org/glib/) >= 2.24.2



## USAGE

Use `dps -h` for a list of options that can be specified in the
command line. The most important options include :

* `mutate` : use any combination of b, k, a corresponding to
activating mutations on β, κ and α respectively. E.g. `--mutate bk
--alpha 1` sets α=1 and activates mutations on β and κ. Defaults: no
mutations and no copy number control (i.e. beta=0.05, kappa=0,
alpha=0).

* `mu` : specifies the probability of mutation per plasmid replication
      event (default: 5e-3).

* `mut_rng` : specifies the width / 2 of the uniform distribution
        around a plasmid's current parameter values used for mutations
        (default: 0.05).

* `pconj` : specifies the probability of a successful horizontal
      transmission event per donor host at a given time step. Use
      `--pconj 0` (default) to switch off conjugation.

* `steps` : how many steps to run the simulation for. The `SIGINT`
      signal (C-c C-c) is taken to signify a user-requested premature
      end to the simulation and is handled gracefully by the program.

* `psize` : specifies the host population size (default: 1000)

* `load_from` : specify the final population of another simulation as
     the initial population for the current simulation.

* `mu` : specify the mutation probability per plasmid replication



## OUTPUT


The simulation's output is stored in an HDF file format; the file name
(typically ending in `.h5`) must be specified as the last argument to
the `dps` command.

The HDF output file has four major sections :

1. `dynamics` : contains the evolutionary dynamics of the simulation
split in four groups

* `counters` : records the numbers of various events (such as
      population size, total copy number, division events etc.) over
      time

	
* `intra`, `inter` and `global` : these contains descriptive
      statistics about evolutionary variables (i.e. means in `M`,
      variances in `V` and pairwise covariances in `C`) at three
      different levels: within hosts (intra), between hosts (inter)
      and across all plasmids regardless of hosts (global)

2. `histograms` : contains the copy number and cell age (i.e. number
of simulation steps required for a host to divide) histograms. The
bins are stored in dataset `x`, whereas the counts in dataset `y`
which has dimensionality `(fparts x max_cn)`, where `fparts` and
`max_cn` are arguments to `dps`. `histograms` also contains all the
joint histograms between β, κ, and α in the respective groups `bk`,
`ba` and `ka`. In the case of the joint distributions, the bins are
stored in datasets `x` and `y`, whereas the counts in dataset `z`
which has dimensionality (`fparts` x `nbins`+1 x `nbins`+1), where
`nbins` is an argument to `dps`.

3. `settings` : provides access to the simulation's parameter values. 

4. `population` : contains the state of the plasmid (group `profiles`)
and host (group `hosts`) population at the end of the simulation as
nested-parentheses strings.




## POST-PROCESSING

R source file `dps.r` is provided as a convenient template for
viewing/processing the output of the simulation (HDF file). This
requires the R package
[`hdf5`](http://cran.r-project.org/web/packages/hdf5/index.html) for
loading/saving HDF files.

In R, load `dps.r` and use

    results <- dps.load("results.h5")

in order to load the simulation results located in file `results.h5`
into R object `results`. 

Having loaded the results, use

    dps.analyze(results)

to plot a basic view of the simulation. 

