IBMF — Individual-Based Mean Field (IBMF)
=====================================

Overview
--------
This folder contains a C++ program implementing the Individual-Based Mean Field
(IBMF) stationary solution for the Generalised Lotka–Volterra model on sparse graphs.
The implementation supports both zero-temperature (T=0) and finite-temperature (T>0)
solutions, with optional immigration (lambda), graph generation, and multiple
random initialization sequences.

Files
-----
- `IBMF_LV_sequential.cpp`  — main program (single translation unit; includes the headers below)
- `IBMF_common.h`           — common types, graph generation, utilities
- `IBMF_convergence_T0_seq.h` — convergence routine for T=0
- `IBMF_convergence_finite_T_seq.h` — convergence routine for finite T (uses GSL special functions)

Prerequisites
-------------
- A C++ compiler (g++ recommended)
- GNU Scientific Library (GSL) installed (for random numbers and special functions)

On Debian/Ubuntu-like systems you can install GSL with:

```bash
sudo apt-get update
sudo apt-get install build-essential libgsl-dev
```

Compile
-------
The program is contained in a single `.cpp` file that includes the headers.
Use the following command to compile (example with optimization):

```bash
cd /path/to/GenLotkaVolterra_SparseGraphs/IBMF
g++ -std=c++11 -O3 IBMF_LV_sequential.cpp -lgsl -lgslcblas -lm -o IBMF_LV_sequential
```

If your compiler is in a different location or you need to specify include/library
paths, add `-I` and `-L` flags as appropriate.

Run
---
The program reads an interaction graph from standard input unless `--gr_inside` is
provided (in which case the program generates a graph internally). See command-line
options below for control.

Basic usage (prints help):

```bash
./IBMF_LV_sequential --help
```

Example: run zero-temperature (T=0) with an internally generated random regular graph (RRG)

```bash
./IBMF_LV_sequential --gr_inside --graph_type RRG -N 1024 -c 3 --eps 0.0 --mu 0.2 --sigma 0.0 --seed_graph 1 -T 0
```

Example: run finite-temperature with immigration

```bash
./IBMF_LV_sequential --gr_inside --graph_type ER -N 1024 -c 3 --mu 0.2 --sigma 0.1 --seed_graph 2 -T 0.01 --lambda 1e-6
```

Input graph format (when not using `--gr_inside`)
-------------------------------------------------
If you do not use `--gr_inside`, the program expects the graph on standard input in the following format:

First line:
```
N M
```
where `N` is the number of nodes (species) and `M` the number of edges.

Then `M` lines, each with:
```
i j aij aji
```
where `i` and `j` are node indices (0-based) and `aij` is the interaction from i->j and
`aji` is the interaction from j->i. If your input uses the reversed order, use
`--alpha_inverse` to flip the parsing.

Key command-line arguments
--------------------------
The program supports many options. The most relevant are listed below (use `--help` to see the full list):

- `--avn_0` (double, default 0.08) : initial average abundance
- `--random_init dn id_0 num_init_conds` : random initial conditions in [avn_0-dn, avn_0+dn]
- `-T` or `--temp` (double, default 0.01) : temperature; set `T=0` for zero-temperature routine
- `--lambda` (double, default 1e-6) : immigration rate
- `--tol` (double) : tolerance for convergence (e.g. `1e-6`)
- `--max_iter` (int) : maximum number of iterations
- `--seed_seq` (unsigned long) : seed for generating update sequences
- `--num_seq` (int) : number of sequences to try
- `--tol_fp` (double) : tolerance to decide if two fixed points are the same
- `--damping` (double) : damping factor for updates (1.0 means no damping)
- `--print_avgs` : if present, save individual average abundances per run
- `--print_only_last` : print results only for the last attempted sequence
- `--gr_inside` : generate the interaction graph internally instead of reading stdin
- `--graph_type` : `RRG` or `ER` (used with `--gr_inside`)
- `--eps` : degree of symmetry for generated edges (probability that aij == aji)
- `--mu` `--sigma` : mean and std.dev. of interaction strengths (used with `--gr_inside`)
- `-N` or `--size` : number of species
- `-c` or `--connect` : degree (RRG) or average connectivity (ER)
- `--alpha_inverse` : read input edges in reversed interaction order

Output
------
The program prints a tab-separated summary per run (iterations, converged/diverged status,
average abundance, fluctuations, number of unconverged nodes, number of extinct nodes,
seeds used, whether fixed points matched, runtime in seconds). When `--print_avgs` is set,
individual abundances are written to files named like `<fileout_base>_seedseq_<s>_seedinit_<i>.txt`.

Notes and troubleshooting
-------------------------
- The finite-temperature implementation uses GSL special functions (gamma, confluent hypergeometric).
  For extreme parameter ranges these functions may overflow or become numerically unstable; the
  implementation includes asymptotic fallbacks but watch the program stderr for warnings about
  divergence or gamma-function approximations.
- If you get linker errors about GSL, ensure `libgsl` and `libgslcblas` are installed and the
  `-lgsl -lgslcblas -lm` flags are present.

Authors / Contact
-----------------
This code is part of the `GenLotkaVolterra_SparseGraphs` repository.
If you need help or wish to improve documentation, open an issue or contact the maintainer.

License
-------
Check the repository root for license information; this folder does not contain a separate license file.
