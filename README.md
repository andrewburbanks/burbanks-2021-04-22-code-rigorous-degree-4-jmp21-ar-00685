# burbanks-2021-04-22-code-rigorous-degree-4-jmp21-ar-00685

## Source code to reproduce the rigorous results and figures of the preprint jmp21-ar-00685

### Preprint

**REFERENCE**

Published Manuscript doi: <a href="https://doi.org/10.1063/5.0054823">10.1063/5.0054823</a>

Submitted Manuscript #JMP21-AR-00685

**TITLE**

Rigorous computer-assisted bounds on the period doubling renormalisation fixed point and eigenfunctions in maps with critical point of degree 4.

**ABSTRACT**

We gain tight rigorous bounds on the renormalisation fixed point for period doubling in families of unimodal maps with degree 4 critical point. We use a contraction mapping argument to bound essential eigenfunctions and eigenvalues for the linearisation of the operator and for the operator controlling the scaling of added noise. Multi-precision arithmetic with rigorous directed rounding is used to bound operations in a space of analytic functions yielding tight bounds on power series coefficients and universal constants to over 320 significant figures. 

**AUTHORS**

* Andrew D Burbanks.
* Judi A Thurlby.
* Andrew H Osbaldestin.

**SUBMITTED**

2021-04-22

**ACCEPTED**

2021-10-18

## Code

### 1 Quick start

#### 1.1 Python version (requires Python 3)

1. Download the repository.
2. Change directory into the /python directory.
3. Launch Jupyter Notebook (Python 3 required) and load file 2_feigenbaum_existence_2019-10-11-1815_180.ipynb.

The file contains the results of a previous run and can also be run to reproduce the results (in particular, the proof at truncation degree 40) from the paper.

#### 1.2 Julia version (requires Julia 1.3.1)

1. Download the repository.
2. Change directory to julia/ directory.
3. Launch Jupyter Notebook (Julia kernel required) and load file 2 feigenbaum_degree4_2020-02-26-2357_40.ipynb.

The file contains the results of a previous run and can also be run to reproduce the results (in particular, the proof at truncation degree 40) from the paper.

### 2 Manifest

* /python --- Python code, using decimal arithmetic with rigorous directed rounding.
* /julia --- Julia code, using binary arithmetic with rigorous directed rounding.

Each folder contains Jupyter notebooks.

Documentation (in the files indicated) is in LaTeX and is embedded in the files. The files contain documentation (where indicated), source code, results from the code, and relevant figures generated directly by the code.

### 3 Python code

The python code imports a library of rigorous functions (also provided).  It is based on the 2-variable formulation due to [Eckmann, Koch, Wittwer, 1982 and 1984]. The python implementation of the rigorous framework relies on the decimal module, which conforms to the relevant industry standards, in particular the directed rounding modes. The pathos multiprocessing suite is used to distribute some tasks between processors. The context in the decimal module is per-process so the relevant precision and rounding modes are safe to use in a multiprocessing environment.

#### 3.1 Files

All files are provided as the original Jupyter notebooks and also for convenience of viewing on systems that don't have Jupyter installed, as html files that can be downloaded and viewed in any browser.

1. nonrigorous_starting_points_2019-12-06-0951 (nonrigorous computations to get approximate fixed points)
2. feigenbaum_existence_2019-10-11-1815_180 (proofs; with thorough LaTeX documentation)
3. Rigorous_Framework (the rigorous framework and tests that are imported into the above)
4. rigorous/ (contains the files of the rigorous framework and the unit tests)
5. dat/ (data files generated by 1 (nonrigorous starting points) required for the proof)
6. fig/ (figures generated)
7. plotting.py (convenience plotting functions)
8. renorm.py (some common functions defining the renormalisation operators)

#### 3.2 Some dependencies

* https://www.python.org/ Python 3.
* https://jupyter.org/ Jupyter Notebook (with Python 3 kernel).
* https://docs.python.org/3/library/decimal.html Decimal library (built-in).
* https://pypi.org/project/pathos/ Pathos may be used for multiprocessing.

### 4 Julia code

The Julia code has the rigorous code embedded into the notebook itself in order to reduce dependencies. Const globals are used in various places where this results in increased efficiency.

The rigorous framework relies on the BigFloat module, which conforms to the relevant industry standards, in particular the directed rounding modes. The Distributed package, in particular the @everywhere macro and careful use of pmap parallel mapping are used to distribute some tasks over processors. The context in the BigFloat module is per-process so the relevant precision and rounding modes are safe to use in a multiprocessing environment. At the highest truncation degrees, memory usage becomes significant so the code serialises objects at various checkpoints and deserialises again, allowing experiments with different numbers of processors for the various parts of the computation.

Separate, datestamped, Julia notebooks contain the results of runs with different truncation degrees.

#### 4.1 Files

All files are provided as the original Jupyter notebooks and also for convenience of viewing on systems that don't have Jupyter installed, as html files that can be downloaded and viewed in any browser.

1. bifurcations_degree4 (nonrigorous location of period-doubling accumulation)
2. feigenbaum_degree4_2020-02-26-2357_40 (proof; thorough LaTeX documentation)

The remaining files provide the code and outputs for truncation degrees 80, 160, 320, 480, and 640. These should be used to view the output from the runs at higher truncation degree.

3. feigenbaum_degree4_2020-02-26-2333_80
4. feigenbaum_degree4_2020-02-27-0041_160
5. feigenbaum_degree4_2020-02-26-2257_320
6. feigenbaum_degree4_2020-02-27-0225_480
7. feigenbaum_degree4_2020-05-19-2248_parallel_41
8. feigenbaum_degree4_2020-05-19-2248_parallel_480
9. feigenbaum_degree4_2020-05-19-2248_parallel_640

#### 4.2 Some dependencies

* https://julialang.org/ Julia 1.3.1.
* https://jupyter.org/ Jupyter Notebook (with Julia kernel).
* https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/ BigFloat (Built-in).
* https://github.com/RalphAS/GenericSchur.jl GenericSchur.
* https://docs.julialang.org/en/v1/stdlib/Distributed/ Distributed (Standard Library).

The authors plan to distribute a more general framework in a more widely accessible form, most likely via GitHub, at a later date.
