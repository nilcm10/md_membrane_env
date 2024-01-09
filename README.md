# MD in a membrane environment

![](practical/files/images/box+water.png)

Welcome to the repository for the 11/01/2024 and 16/01/2024 lectures on MD of bilayers and membrane proteins.

## Distribution

Under [theory/](theory/) you have the theoretical presentation to the topic. Under [practical/](practical/) you have a [protocol](practical/README.md) with the different steps to be followed to:
- Build, minimize, equilibrate, simulate and analyze a POPC membrane system (under [just_popc/](practical/just_popc)).
- Build and analyze a POPC+CHL membrane system (under [popc+chl/](practical/popc+chl)).
- Build and analyze a POPC+CHL membrane system with an embedded membrane protein (under [membrane_protein/](practical/membrane_protein)).

## Main software

- [PACKMOL-memgen](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00269) to build the systems.
- [GROMACS](https://manual.gromacs.org/) to prepare, simulate and analyze the simulations.
- [FATSLiM](http://fatslim.github.io/) to analyze the membranes in the simulations.

## Setting up

In this repo we assume the [Anaconda](https://www.anaconda.com/) package manager is installed in your computer (you can check that by typing `conda` on your terminal session). If you do, continue from [Creating an environment](./README.md#creating-an-environment). If not, go on with [Installing miniconda](./README.md#installing-miniconda).

### Installing miniconda

A lighter version of [Anaconda](https://www.anaconda.com/) can be installed, and that is [miniconda](https://docs.conda.io/en/latest/miniconda.html). To download and install it:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

Follow the instructions.

Once you have it installed, you probably have to restart your terminal or do `source ~/.bashrc`, and then `conda init bash`.

### Creating an environment

Now we're going to create an environment and we're going to call it `membranes`:

```
conda create -n membranes python=3.6
conda activate membranes
```

And install all the packages needed through `conda install`:

```
conda install -c conda-forge numpy pandas matplotlib ambertools acpype lipyphilic
```

Say `Y`.

### Installing the rest of the software

We might need to install the molecular dynamics machine as well. First check if it's installed in your system by typing `gmx help` or `which gmx` in your shell. If there's no prompt or the software is not installed, do:

```
sudo apt install gromacs
```

To install FATSLiM and test if it works:

```
pip install fatslim
fatslim self-test
```

To start working with the practical part, `git clone` this repo anywhere you want in your local machine and follow the steps of the [protocol](practical/README.md).