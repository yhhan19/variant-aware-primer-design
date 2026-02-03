# ∆-PRO (Delta-PRO): Risk & Dimer Optimization for Multiplex PCR Design

This repository contains a C++ command-line tool for optimizing primer design regions (PDRs) and solving a primer–primer dimer minimization subproblem.

From the provided `main.cpp`, the program supports two families of workflows:

- **Risk optimization**: select PDRs to minimize a risk/loss objective under amplicon-length constraints.
  - `risk-d1`: DP-based relaxed convex optimizer (strategy s1)
  - `risk-d2`: DP-based relaxed convex optimizer (strategy s2, with a min/max amplicon length range)
  - `risk-h`: greedy/random heuristic (Olivar-style)
- **Dimer minimization**:
  - `dimer-h`: fast solver on a k-partite graph formulation

## Build

You can compile with any C++17-capable compiler.

```bash
g++ -O3 -std=c++17 -o delta-pro *.cpp
