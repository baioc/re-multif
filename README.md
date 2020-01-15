# Replicating a Multi-Functional Synthetic Gene Network

[![DOI](https://zenodo.org/badge/214520225.svg)](https://zenodo.org/badge/latestdoi/214520225)

Here we replicate the work of Purcell, di Bernardo, Grierson and Savery on [A Multi-Functional Synthetic Gene Network: A Frequency Multiplier, Oscillator and Switch](https://dx.doi.org/10.1371%2Fjournal.pone.0016140) for [ReScience](https://rescience.github.io/).


## Reproducing our work

The model is built using Octave and programs were tested in version 5.1.0 on a couple of different Linux systems.
No extensions should be needed.

To reproduce our experiments, make sure simulation parameters are correctly configured (in the source code) and then simply run the target script:

```bash
$ octave src/freqdiv/SCRIPT_NAME.m
```

Simulation results are printed to PDF which is put in the same directory the `octave` command is run from.
To crop empty margins in the generated file we have used the `pdfcrop` tool:

```bash
$ pdfcrop uncropped_.pdf cropped.pdf
```

Simply swapping the exported filename extension in the scripts to ".png" is an alternative to this process, but figures may have a lower resolution.


## Study Documentation

The [article](article.pdf) was written over a [ReScience submission template](https://github.com/ReScience/template) with a [custom makefile](doc/Makefile).
Compiling this document may require several TeX packages and a Python installation.

```bash
$ cd doc
$ make
```
