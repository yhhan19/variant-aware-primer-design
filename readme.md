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

## Usage

Compile the program into a single executable (for example, `delta-pro`), then run:

```bash
./delta-pro -i <input_file> -o <output_file> -x <mode> [options]
```

---

### Required arguments

| Flag | Description |
|------|-------------|
| `-i` | Input file |
| `-o` | Output file |
| `-x` | Execution mode |

---

### Modes

#### `risk-d1`
DP-based relaxed convex optimizer (strategy s1) using a fixed amplicon length.

```bash
./delta-pro -i rate.txt -o pdr.txt -x risk-d1
```

---

#### `risk-d2`
DP-based relaxed convex optimizer (strategy s2) using a minimum and maximum amplicon length.

```bash
./delta-pro -i rate.txt -o pdr.txt -x risk-d2
```

---

#### `risk-h`
Greedy / random heuristic optimizer (Olivar-style baseline).

```bash
./delta-pro -i rate.txt -o pdr.txt -x risk-h
```

---

#### `dimer-h`
Fast heuristic solver for the k-partite graph dimer minimization problem.

```bash
./delta-pro -i graph.txt -o solution.txt -x dimer-h
```

---

### Optional parameters

| Flag | Description | Default |
|------|-------------|---------|
| `-S` | Random seed | `42` |
| `-Ln` | Maximum (or nominal) amplicon length | `420` |
| `-Lx` | Minimum amplicon length | `252` |
| `-Lp` | PDR length | `40` |
| `-Ux` | Upper bound for convex search parameter | `10000` |
| `-Un` | Lower bound for convex search parameter | `0.1` |
| `-I` | Iterations for `dimer-h` | `1000` |

---

### Output

- **Risk modes** (`risk-d1`, `risk-d2`, `risk-h`)  
  Write an optimized Primer Design Region (PDR) to the output file.

- **Dimer mode** (`dimer-h`)  
  Writes a single line of space-separated indices representing the selected solution.

All modes print the final loss value and total runtime in seconds.
