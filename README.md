# Elliptic PDE Solver

## 1. Overview

This project implements a numerical solver for the elliptic partial differential equation

    -a * u_xx - b * u_yy = f(x,y)

on the unit square [0,1] x [0,1] with Dirichlet boundary conditions.  
The solution is obtained using the **Seidel method**. Additionally, the largest and smallest eigenvalues of the discrete Laplace operator are computed using the **power method** and the **shifted power method**.

The code is written in C++ with OpenMP for parallelization. A Python script (using numpy and matplotlib) is provided for visualization of the results.

## 2. Building and Running

### Requirements

- C++17
- OpenMP support (/openmp, enabled by default in .sln project for Release)
- Python 3.x with packages: numpy, matplotlib

### Build Instructions

1. Open the solution file `Elliptic_PDE_Solver.sln` in Visual Studio.
2. Select the desired configuration (Release or Debug).
3. Build the solution (Build → Build Solution).

### Running the Program

- Ensure that an `output` folder exists in the root directory (next to `src`). If not, create it manually.
- Run the executable (no command‑line arguments required).  
- Upon completion, the following files will be created in the `output` folder:
  - `solution.bin` – binary file containing the grid and the numerical solution.
  - `comparison.png` – 3D plots of the exact solution, numerical solution, and absolute error.
  - `error_map.png` – 2D map of the absolute error.

### Configuration

All adjustable parameters are located in `src/main.cpp`. Example:
```cpp
    const DataType a = 0.9, b = 1.1; // equation coefficients
    const DataType h = 0.1; // desired grid step (will be adjusted)
    const DataType tol_power = 1e-6; // tolerance for the power method
    const DataType tol_seidel = 1e-6; // tolerance for the Seidel method
    const size_t max_iter_power = 100000; // max power‑method iterations
    const size_t max_iter_seidel = 100000; // max Seidel iterations
    const std::string output_dir = "../output"; // output directory
```

Modify these values as needed and rebuild the project.

## 3. Generating Documentation

The project uses **Doxygen** to generate both HTML and LaTeX (PDF) documentation.  
Configuration is stored in `docs/Doxyfile`.

### Prerequisites

- Install Doxygen (from [doxygen.nl](https://www.doxygen.nl/))
- For PDF generation: install a LaTeX distribution (e.g., MiKTeX, TeX Live)

### Steps

1. Open a terminal in the `docs` folder.
2. Run Doxygen:
   
       doxygen Doxyfile
   
3. After successful execution:
   - HTML documentation is located in `docs/html/index.html`.
   - LaTeX files are placed in `docs/latex`. To generate a PDF, navigate to that folder and run:
     
           cd latex
           make pdf          # or run pdflatex refman.tex several times
     
     The resulting PDF is named `refman.pdf`.

### Important Doxyfile Settings

- `INPUT = ../src ./` – processes source files from `../src` and the current folder (for `mainpage.dox`).
- `USE_MATHJAX = YES` – renders formulas with MathJax (no need for external tools).
- `GENERATE_LATEX = YES` – creates LaTeX output.
- `EXTENSION_MAPPING = py=Python` – ensure proper handling of Python files.

### Documentation Structure

- `mainpage.dox` contains the main page content.
- All source files are documented using Doxygen comments (`///` or `/** ... */` for C++, `"""! ... """` for Python).

After modifying any documentation, re‑run Doxygen to update the generated files.