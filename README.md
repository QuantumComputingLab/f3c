# Fast free fermion compiler (FFFC or F3C)

<p align="center"><img src="doc/doxygen/F3C.png?raw=true" /></p>

Fast free fermion compiler is an application-specific quantum circuit compiler for
time-evolution circuits of spin Hamiltonian systems that can be mapped to free fermions.
We support the following types of Hamiltonians:
- (classical) Ising model
- Kitaev chains
- XY, XZ, and YZ models
- TFIM 
- TFXY, TFXZ, and TFYZ models


### How to run? ###

The F3C compiler is implemented using MATLAB object oriented functionalities
and is compatible with MATLAB R2018a or newer. It uses the 
[QCLAB toolbox](https://github.com/QuantumComputingLab/qclab).

0. Install [QCLAB](https://github.com/QuantumComputingLab/qclab)

1. Clone repository:

        git clone https://github.com/QuantumComputingLab/f3c.git

2. Add all files from qclab directory to MATLAB path to install F3C:

		addpath(f3croot);
		savepath;

3. Run tests in MATLAB:
		
		cd test/
		runTests.m
 
4. Generate documentation with doxygen. Requires [doxygen](https://www.doxygen.nl/index.html) and [doxymatlab](https://github.com/simgunz/doxymatlab). Adjust tags `FILTER_PATTERNS` and `FILTER_SOURCE_PATTERNS`  in `doxygen/Doxyfile.dox` to local `m2cpp.pl` script.
	
		cd doc/doxygen/
		doxygen Doxyfile.dox


## References
The F3C compiler is based on:
- An Algebraic Quantum Circuit Compression Algorithm for Hamiltonian Simulation, Daan Camps, 
Efekan K&ouml;ck&uuml;, Lindsay Bassman, Wibe A. de Jong, Alexander F. Kemper, and 
Roel Van Beeumen (2021) 
- Algebraic Compression of Quantum Circuits for Hamiltonian Evolution, Efekan K&ouml;ck&uuml;,
Daan Camps, Lindsay Bassman, J. K. Freericks, Wibe A. de Jong, Roel Van Beeumen, and 
Alexander F. Kemper (2021) 

## Developers - Lawrence Berkeley National Laboratory
- [Daan Camps](http://campsd.github.io/) - dcamps@lbl.gov
- [Roel Van Beeumen](http://www.roelvanbeeumen.be/) - rvanbeeumen@lbl.gov


## Funding
The F3C project is supported by the Laboratory Directed Research and
Development Program of Lawrence Berkeley National Laboratory under U.S.
Department of Energy Contract No. DE-AC02-05CH11231.

## About
F3C Copyright (c) 2021, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy). All rights 
reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.
