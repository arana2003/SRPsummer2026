########################################################################
#
# Pipeline Automation for MLCP                                         #
## Kaiwan Bilal, Phillip Thomas, Alexis Rana                          ##
#
########################################################################

### The Pipeline

Multi-Layer Cannonical Polyadic, or MLCP, eigensolver (as created by Phillip Thomas 
at Lawrence-Berkeley National Lab) enables computational researchers to employ a fully variational 
method of modeling molecular vibrational spectra without requiring an obscene amount of computer resources.

On its own, MLCP is a fully variational eigensolver that can take anharmonic constants well beyond the fourth order, 
but potential energy surfaces (PES) are commonly available as quartic force fields. 
Generating a viable PES usually calls another pre-existing electronic structure package and generally
only yields harmonic force constants, but MLCP is built to consider fourth-order anharmonicity. 

**That's where this pipeline comes in!**

Our pipeline automates this collective process from start to finish. Calling Push-Button Quartic Force Fields,
MOPAC, NWChem, and MLCP, users provide a "master" input file (as shown in `inputs/sample.inp`), which is a combined
input document merging PBQFF, NWChem, and MLCP into one joint input file for the pipeline.

### Installation

The entire pipeline is contained within a multi-stage Podman/Docker Container for ease of distribution.
Considering such, there are two options for installation:

## Downloading container as image
If your CPU/GPU architecture matches what we used to build (LBL/NERSC Perlmutter):

1) podman-hpc pull docker.io/library/ubuntu:latest

2) podman-hpc migrate docker.io/library/ubuntu:latest

## Building image from Containerfile
1) Download the container and cd into the root directory `interface`

2) Build the container set (replace `podman-hpc` with respective container handler)

```bash
podman-hpc build -f ./containers/Containerfile_env -t env .
podman-hpc build -f ./containers/Containerfile_working -t working .

podman-hpc migrate working # If needed, for HPC centers
```

### Usage

This entire `interface` directory gets mounted into the container, so users can 
simply add/create/modify files within the interface directory for pipeline use.

1) Create a plain-text input file following this model (complete version in `inputs/sample.inp`)

```bash
# System Name
system_name = h2o

# Geometry
 O       0.0  0.0  1.0  
...
--

# PBQFF Input
charge = 0
optimize = true
...
-- 
# Intder.in**
# INTDER ##########################
    3    3    3    0    0    3    0    0    0    1    0    0    0    1    1    0
STRE     1    2
...
--

# NWChem Input
BASIS  
 H library sto-3g  
 O library sto-3g  
...
--

# MLCP Input
$control
system='h2o'
rs='0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
...
--
```

If you've worked with any of these models, you'll notice that each section (with the exception 
of a shared geometry input) is simply a complete input file for running standalone calculations 
with package. Each section is copied directly into individual input files, allowing users to 
fully tune all available parameters granted by each standalone component.

**Note: Intder.in section only necessary depending on PBQFF options. For ease, setting
PBQFF coordinates to "normal" eliminates this dependency.

```bash
# PBQFF Input
# ...
coord_type = "normal"
```

2) Edit system variables in `checkpoint.sh`

```bash
export MLCPPIPE="/global/cfs/cdirs/m5128/kbilal/interface"     # absolute path to directory
export INPUT="${MLCPPIPE}/inputs/sample.inp"                   # relative path to working input file
export PROJ="m5128"                                            # Slurm project name
export EXTEND_PES="false"                                      # false: use NWChem optimization. true: use "Morsified" long-range multi-fidelity PES
export SYS_NAME="h2o"                                          # system name of choice
```

3) Edit `interface/slurm_jobs/control.slurm`\
This Slurm cron script is the root job for the pipeline. Make sure
everything looks correct (namely, the -A field and absolute path to `interface` directory)

4) Run the pipeline

Submit the cron Slurm job from your terminal:

```bash
scrontab interface/slurm_jobs/control.slurm
```

**Note: when re-running the same simulation, be sure to
delete any old `pipeline.checkpoint` files or irrelevant
`system.vars` variables