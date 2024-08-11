# mean-FLAME

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
What does this do??

The installation process creates an executable called `life-cycle` that is located in `life-cycle/build`. To run this executable, run
```
./life-cycle <beta> <seed_rate> <>

## Lokta-Volterra
What does this do??

The installation process creates an executable called `lokta-volterra` that is located in `life-cycle/build`. To run this executable, run

## SIRS