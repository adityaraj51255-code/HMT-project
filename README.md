# 2D Steady-State Heat Conduction using Finite Difference Method (FDM)

## 📌 Overview
This project solves the **2D steady-state heat conduction equation** in a rectangular plate using the **Finite Difference Method (FDM)**.

The temperature distribution is obtained numerically and compared with an analytical solution derived using separation of variables.

---

## 🧮 Governing Equation

The steady-state heat conduction equation (Laplace equation):

\[
\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0
\]

### Assumptions:
- Steady-state condition  
- No internal heat generation  
- Constant thermal conductivity  

---

## ⚙️ Methodology

### 1. Domain Discretization
The rectangular plate is divided into a grid:

- Grid size: `nx × ny`
- Grid spacing:
\[
dx = \frac{L_x}{nx-1}, \quad dy = \frac{L_y}{ny-1}
\]

---

### 2. Finite Difference Approximation

The governing equation is approximated as:

\[
\frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{dx^2} +
\frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{dy^2} = 0
\]

This results in a system of linear equations:

\[
A \cdot T = b
\]

---

### 3. Boundary Conditions

All four edges are maintained at fixed temperatures (Dirichlet conditions):

- Left boundary → \(T_L\)  
- Right boundary → \(T_R\)  
- Bottom boundary → \(T_B\)  
- Top boundary → \(T_T\)  

---

### 4. Numerical Solution

The system is solved using:

```python
np.linalg.solve(A, b)
