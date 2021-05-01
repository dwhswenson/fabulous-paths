# Alanine dipeptide examples

This example illustrates the use of the OPS CLI and the fabulous-paths plugin.
We start with a file, `ad_setup.db`, which is ready to run TPS on alanine
dipeptide in vacuum. It also includes the additional CV definitions we'll use
in FABULOUS. 

Note that the additional CVs are not calculated during sampling -- the functions
have been defined, but these will only be calculated when processing the data
for FABULOUS.

## Requirements

In addition to fabulous-paths and its requirements (such as FABULOUS and
OpenPathSampling), you will need to install OpenMM and OpenMMTools. They can be
installed via conda with:

```bash
conda install -c conda-forge openmm openmmtools
```

## Contents of this directory

* `README.md`: This file. Instructions for the tutorial.
* `ad_setup.db`: A OpenPathSampling file with information about how to run TPS,
  as well as initial conditions.
* `ad_setup.py`: A Python script that generates `ad_setup.db`. This does not
  need to be run as part of the tutorial; it is included for informational
  purposes.
* `ad_fabulous.yml`: FABULOUS parameters for this system
* `ad.pdb`: A PDB file describing the topology of this system.

## Running TPS

To run the TPS for 5000 TPS MC steps, run the following command:

```bash
openpathsampling pathsampling ad_setup.db -o tps.db --nsteps 5000
```

This may take about an hour or so.


## Modeling the results with FABULOUS

First, you'll need to copy the CV definitions from `ad_setup.db` into `tps.db`.
Use the `openpathsampling append` command for that:

```bash
openpathsampling append ad_setup.db -a tps.db --cv omega --cv theta --cv NCaR \
  --cv end_Ca_end --cv dOO
```

You can check the contents of a file with the `openpathsampling contents`
command, e.g., `openpathsampling contents tps.db` to see what is contained in
that file.

To perform the analysis with FABULOUS, run the following command:

```bash
openpathsampling fabulous tps.db --cv phi --cv psi --cv omega --cv theta --cv NCaR \
  --cv end_Ca_end --cv dOO --ref c7ax_input.pdb --keep-atoms backbone \
  -c ad_fabulous.yml --results results -n 50
```

The initial stages of extracting the data will take about 20 minutes or so.
Then the data is passed from `fabulous-paths` over the `FABULOUS` for the real
analysis, which can take up to ??6 hours??

## Analyzing the ML models

## Alternate ways to run the analysis

### `--keep-atoms`

In the example above, we used the MDTraj atom selection language to select only
the backbone atoms. However, you can also explicitly provide a list. That list
can be provided as comma-separated integers on the command line (no spaces!),
as a NumPy array stored in a .npy file, or as a JSON array stored in a file
ending in .json or .jsn.

### `--cv`

Rather than explicitly mentioning each CV, you can use the `--cv-list` argument. This 
