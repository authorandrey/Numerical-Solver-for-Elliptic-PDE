# -*- coding: utf-8 -*-
## @file plot_full.py

"""!
@file plot_full.py
@brief Visualization of the numerical solution and comparison with the exact one.
@details The script reads a binary file solution.bin containing the grid and numerical solution values,
         creates 3D plots of the exact solution, numerical solution, and absolute error, as well as a 2D error map.
         Uses matplotlib and numpy.
"""

import numpy as np
import matplotlib.pyplot as plt
import struct
import sys
import os
from mpl_toolkits.mplot3d import Axes3D
from typing import Tuple, Union

def u_exact(x: Union[float, np.ndarray], y: Union[float, np.ndarray]) -> np.ndarray:
    """!
    @brief Exact analytical solution of the problem.
    @param x x-coordinate (scalar or numpy array).
    @param y y-coordinate (scalar or numpy array).
    @return u(x,y) = x^2 + cos^2(y*x).
    """
    return x**2 + np.cos(y * x)**2

def read_full_solution(filename: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """!
    @brief Reads a binary solution file created by OutputHandler::save_full_solution.
    @param filename File name (e.g., 'solution.bin').
    @return Tuple (X, Y, u_num), where X, Y are coordinate matrices (meshgrid),
            u_num is a 2D array of the numerical solution.
    @throws ValueError if the file is too short or has an invalid format.
    @details File format:
             - first 8 bytes: number of nodes along x (nx, size_t)
             - next 8 bytes: number of nodes along y (ny, size_t)
             - then nx*ny double values in row-major order (j=0..ny-1, i=0..nx-1)
    """
    with open(filename, 'rb') as f:
        nx_bytes = f.read(8)
        ny_bytes = f.read(8)
        if len(nx_bytes) < 8 or len(ny_bytes) < 8:
            raise ValueError("File too short: cannot read dimensions")
        nx = struct.unpack('Q', nx_bytes)[0]
        ny = struct.unpack('Q', ny_bytes)[0]

        total = nx * ny
        data_bytes = f.read(total * 8)
        if len(data_bytes) < total * 8:
            raise ValueError("File too short: cannot read all data")

        data = struct.unpack('d' * total, data_bytes)

    u_num = np.array(data).reshape((ny, nx))
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    X, Y = np.meshgrid(x, y)
    return X, Y, u_num

def main() -> None:
    """!
    @brief Main function of the script.
    @details Performs reading of the file, computation of exact solution and error,
             creates and saves plots (comparison.png, error_map.png),
             and displays them on the screen.
             Expects the output directory as first command-line argument (default '.').
    """
    if len(sys.argv) > 1:
        output_dir = sys.argv[1]
    else:
        output_dir = '.'

    input_file = os.path.join(output_dir, 'solution.bin')
    comparison_file = os.path.join(output_dir, 'comparison.png')
    error_map_file = os.path.join(output_dir, 'error_map.png')

    try:
        X, Y, u_num = read_full_solution(input_file)
    except Exception as e:
        print(f"Error reading binary file: {e}")
        sys.exit(1)

    u_ex = u_exact(X, Y)
    error = np.abs(u_ex - u_num)

    fig = plt.figure(figsize=(18, 5))

    # Exact solution
    ax1 = fig.add_subplot(131, projection='3d')
    surf1 = ax1.plot_surface(X, Y, u_ex, cmap='viridis', linewidth=0, antialiased=True)
    ax1.set_title('Exact solution')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('u')
    fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)

    # Numerical solution
    ax2 = fig.add_subplot(132, projection='3d')
    surf2 = ax2.plot_surface(X, Y, u_num, cmap='viridis', linewidth=0, antialiased=True)
    ax2.set_title('Numerical solution')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('u')
    fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)

    # Absolute error
    ax3 = fig.add_subplot(133, projection='3d')
    surf3 = ax3.plot_surface(X, Y, error, cmap='hot', linewidth=0, antialiased=True)
    ax3.set_title('Absolute error')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_zlabel('error')
    fig.colorbar(surf3, ax=ax3, shrink=0.5, aspect=10)

    plt.tight_layout()
    plt.savefig(comparison_file, dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved as {comparison_file}")

    # Error map (2D)
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, error, shading='auto', cmap='hot')
    plt.colorbar(label='Absolute error')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Error map |exact - numerical|')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(error_map_file, dpi=300)
    print(f"Error map saved as {error_map_file}")

    plt.show()

if __name__ == "__main__":
    main()