import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------
# USER INPUT
# ----------------------------------
nx = int(input("Enter number of grid points in x-direction: "))
ny = int(input("Enter number of grid points in y-direction: "))

length_x = float(input("Enter plate length in x (default 1): ") or 1)
length_y = float(input("Enter plate length in y (default 1): ") or 1)

# Boundary temperatures (all fixed)
T_L = float(input("Left boundary temperature: "))
T_R = float(input("Right boundary temperature: "))
T_B = float(input("Bottom boundary temperature: "))
T_T = float(input("Top boundary temperature: "))

# ----------------------------------
# GRID SETUP
# ----------------------------------
dx = length_x / (nx - 1)
dy = length_y / (ny - 1)

total_nodes = nx * ny
A = np.zeros((total_nodes, total_nodes))
b = np.zeros(total_nodes)

def node(i, j):
    return i * ny + j

# ----------------------------------
# BUILDING FDM EQUATIONS
# ----------------------------------
for i in range(nx):
    for j in range(ny):

        p = node(i, j)

        # interior nodes
        if (i > 0 and i < nx-1) and (j > 0 and j < ny-1):

            A[p, node(i, j)] = -2 * (1/dx**2 + 1/dy**2)
            A[p, node(i+1, j)] = 1/dx**2
            A[p, node(i-1, j)] = 1/dx**2
            A[p, node(i, j+1)] = 1/dy**2
            A[p, node(i, j-1)] = 1/dy**2

            b[p] = 0

        # left boundary
        elif i == 0:
            A[p, p] = 1
            b[p] = T_L

        # right boundary
        elif i == nx - 1:
            A[p, p] = 1
            b[p] = T_R

        # bottom boundary
        elif j == 0:
            A[p, p] = 1
            b[p] = T_B

        # top boundary
        elif j == ny - 1:
            A[p, p] = 1
            b[p] = T_T

# ----------------------------------
# SOLVE SYSTEM
# ----------------------------------
print("\nSolving...\n")
T_sol = np.linalg.solve(A, b)

T_num = T_sol.reshape((nx, ny))

# ----------------------------------
# ANALYTICAL SOLUTION (valid here)
# ----------------------------------
def exact_temp(x, y, n_terms=50):
    val = 0.0
    for n in range(1, n_terms, 2):
        val += (1/n) * np.sinh(n*np.pi*y) / np.sinh(n*np.pi) * np.sin(n*np.pi*x)
    return (4 * T_L / np.pi) * val

x_vals = np.linspace(0, 1, nx)
y_vals = np.linspace(0, 1, ny)

T_exact = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        T_exact[i, j] = exact_temp(x_vals[i], y_vals[j])

# ----------------------------------
# ERROR
# ----------------------------------
error = np.abs(T_num - T_exact)
print("Maximum error =", np.max(error))

# ----------------------------------
# PLOTTING
# ----------------------------------
plt.figure(figsize=(12, 5))

plt.subplot(1, 3, 1)
plt.title("Numerical")
plt.contourf(T_num, 20)
plt.colorbar()

plt.subplot(1, 3, 2)
plt.title("Analytical")
plt.contourf(T_exact, 20)
plt.colorbar()

plt.subplot(1, 3, 3)
plt.title("Error")
plt.contourf(error, 20)
plt.colorbar()

plt.tight_layout()
plt.show()