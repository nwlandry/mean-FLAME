# mean-FLAME

This code accompanies the preprint "Stochastic diffusion with approximate master equations with mean-field limits" by Laurent Hébert-Dufresne, Matthew M. Kling, Samuel F. Rosenblatt, Stephanie N.
Miller, P. Alexander Burnham, Nicholas W. Landry, Nicholas J. Gotelli, and Brian J.
McGill.

This repository provides the C++ code necessary to generate the results presented in this preprint. Below we detail the 

## Repository organization

The repository contains several directories:
* `life-cycle`: contains all of the code used for simulating an ecological model using the mean-FLAME model approach.
* `Lokta-Volterra`: contains all of the code used for simulating the Lokta-Volterra model using the mean-FLAME model approach.
* `SIRS`: contains all of the code used for simulating the SIRS model using the mean-FLAME model approach.
* `tutorials`: contains a Jupyter notebook which demonstrates simulating a birth-death process in Python using the mean-FLAME modeling approach.

## Google Colab

In addition to the notebook in the [tutorials](/tutorials/), there is also a [Google Colab notebook](https://bit.ly/meanFLAME).

## Installation

This library uses the Boost and GSL libraries.

Compiling for Mac: Using HomeBrew (after installing, of course), run the following in Terminal
```
brew install cmake
brew install boost
brew install gsl
```

For each of the folders (`life-cycle`, `Lokta-Volterra`, and `SIRS`), build the executable by doing the following:

From the top-level directory, run the following commands:
```
cd <folder>
mkdir build
cd build
cmake ..
cmake --build .
```
where <folder> is one of the three listed sub-directories.

## life-cycle
This model tracks the population of trees over four stages of their life-cycle: seeds, seedlings, saplings, and adult trees. The two dimensions of interest are the number of adult trees and the number of young trees. The seed and seedling dimensions are captured as a fully mean-field description (no explicit master equation state). Across these states, we capture arrival of new saplings and death of saplings , growth of saplings into adult trees, and death of adult trees.

The installation process creates an executable called `life-cycle` that is located in `life-cycle/build`. To run this executable, from the main directory run
```
life-cycle/build/life-cycle <tree2seed> <seed2seedling> <seed_death> <seedling2sapling> <seedling_death> <sapling2tree> <sapling_death> <tree_death> <n_l> <n_me1> <n_me2>
```
If you want to use the parameters in Fig. 8, the command to run is
```
life-cycle/build/life-cycle 1000.0 0.5 20.0 0.2 1.0 0.1 0.5 0.05 50 5 5
```

## Lokta-Volterra
Species interactions, modeled with Lokta-Volterra systems of equations, are one of the most ubiquitous examples of mean-field models in ecology. They capture a simplified version of predator-prey or competition dynamics. This system contains no spatial aspects or coupling between subsystems to demonstrate how to generalize the internal dynamics of mean-FLAME models to systems with more than one degree of freedom (or variable of interest).

The installation process creates an executable called `lokta-volterra` that is located in `life-cycle/build`. To run this executable, from the main directory run
```
Lokta-Volterra/build/lokta-volterra <beta> <mu> <K> <nu> <n_l> <n_me1> <n_me2>
```
If you want to use the parameters in Fig. 5, the command to run is
```
Lokta-Volterra/build/lokta-volterra 0.1 0.005 20 0.1 1 3 3
```

## SIRS
In this example we model an epidemic process where individuals can take one of three states: susceptible ($s$), infected ($i$), or recovered ($r$). Individuals can cycle through all states and back based on a general set of mechanisms referred to as the SIRS process. Susceptible individuals become infected by contact with infected individuals at a infection rate $\beta$. Infected individuals recover and gain immunity at a recovery rate $\alpha$. And recovered individuals lose their immunity and become susceptible again according to a waning rate $\gamma$.

The installation process creates an executable called `SIRS` that is located in `SIRS/build`. To run this executable, from the main directory run
```
SIRS/build/SIRS <beta> <alpha> <gamma> <n_l> <n_me1> <n_me2>
```
If you want to use the parameters in Fig. 6, the command to run is
```
SIRS/build/SIRS 0.025 1.0 0.01 50 10 10
```