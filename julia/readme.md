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
