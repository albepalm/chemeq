# Chemical Equilibrium Calculator for homogeneous gaseous mixtures 

This simple tool evaluates the equilibrium chemical composition of different terrestrial air models (N2, O2, NOx, Ar, and related dissociated and ionized species) at different mixture pressures and temperatures. It solves a system of nonlinear equations, taking into account the variation of the thermodynamic properties of gases as a function of temperature [[1]](#1) [[2]](#2).
The code was written for part of my master's thesis work [[3]](#3).

## Installation

Clone the repo: `git clone https://github.com/albepalm/chemeq.git`
Run the make utility inside the folder: `make`

## Input file

In the `./input` folder can be found two examples of input file. Input files are structured as follows:

```
[m/t] [p_min] [p_max] [p_step] [t_min] [t_max] [t_step] [output_name] [# Nitrogen-only equations] [type of reactions] [# Oxygen-only equations] [type of reactions] [# Nitrogen monoxide-only equations] [type of reactions] [# Argon-only equations] [type of reactions] [ionization] [Nitrogen ratio] [Oxygen ratio] [Argon ratio]
```
- `m/t`: `m` (matrix format, temperature-pressure map), `t` (pressure is fixed);
- pressure varies in log10 scale, temperature in linear scale;
- ```output_name``` less than 256 chars;
- type of reactions (for Nitrogen and Oxygen): 0-bit X2, 1-bit X, 2-bit X+, 3-bit X2+;
- type of reactions (for Nitrogen monoxide): 0-bit NO, 1-bit NO+;
- type of reactions (for Argon): 0-bit Ar, 1-bit Ar+;

The bit of the reaction is `1` if it occurs, `0` otherwise.

## Output file
In the `./output` file are stored the results as tab-separeated values:
- the temperature is in the first column;
- the pressure is in the second (only if `m` option is specified;
- the others contain the molar fractions (or partial pressures, `pp` string is appended at the end of the file name).

## chemeq.config
This file contains the settings of the Newton-Raphson method:
- maximum iterations;
- relative tolerance;
- cli output (0 less verbose, 1, 2 extremely verbose).

## Some post-processed results

![Nitrogen model at 1 bar](/examples/n2.png)

![Oxygen model at 1 bar](/examples/o2.png)

![11 species model at 1 bar](/examples/11s.png)

![Plasma temperature/pressure map](/examples/plasma13s.png)

## References

<a id="1">[1]</a> John D. Anderson. Hypersonic and High-Temperature Gas Dynamics, Third Edition. American Institute of Aeronautics and Astronautics, Inc., Jan. 2019. doi: 10.2514/4.105142.

<a id="2">[2]</a> Bonnie J. McBride, M.J. Zehe, and S. Gordon. NASA Glenn coefficients for calculating thermodynamic properties of individual species NASA. tech. rep. TP-2002,2002.

<a id="3">[3]</a> Palmieri, Alberto. Parametrized numerical analysis of hypersonic flow field. Diss. Politecnico di Torino, 2023.
