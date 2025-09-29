# GrossPitaevskiiSCF

Numerical study of the discretization of the Gross-Pitaevskii eigenvalue problem in Julia

clone as:

```bash
git clone https://github.com/djuliannader/GrossPitaevskiiSCF.git
```

## Data 


This repository provides the data and code associated with the paper Critical Energies and Wigner Functions of the Stationary States of Bose–Einstein Condensates in a Double-Well Trap by D. J. Nader and E. Serrano-Ensástiga, published open access in Advanced Quantum Technologies.

The repository is still under development, so we kindly ask for your patience until its completion. It currently includes a Mathematica notebook that can be used to reproduce the figures from the article.

## Environment Requirements  

To run this repository, please make sure the following environment is available:

- LinearAlgebra (stdlib)  
- PyPlot


## Usage

- For the stationary states, navigate to the /src manipulate the parameters from the file input.dat and execute as:

```bash
julia main.jl
```

Cite as:

D. J. Nader, E. Serrano-Ensástiga, Critical Energies and Wigner Functions of the Stationary States of the Bose Einstein Condensates in a Double-Well Trap. Adv Quantum Technol. 2025, 2400451. https://doi.org/10.1002/qute.202400451

- For the quench dynamics, navigate to the /src manipulate the parameters from the file input_quenchdynamics.dat and execute as:

```bash
julia main_quenchdynamics.jl
```

Cite as:

D. J. Nader, Quench-Induced Dynamical Phase Transitions in BECs: A Gross–Pitaevskii Approach, (In progress )


