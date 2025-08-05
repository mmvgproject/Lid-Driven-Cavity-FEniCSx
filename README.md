# **Lid-Driven Cavity Flow Simulation Using FEniCS**

## **1. Project Title**
**Lid-Driven Cavity Flow Simulation Using the Incremental Pressure Correction Scheme (IPCS) in FEniCS**

---

## **2. Project and Program Specifications**
### **Programming Language & Libraries Used**
- **Python** (Primary language)
- **FEniCS** (Finite element solver)
- **NumPy** (Numerical computations)
- **Matplotlib** (Visualization)
- **ImageIO** (GIF generation)

### **Computational Method**
- **Finite Element Method (FEM)**  
- **Incremental Pressure Correction Scheme (IPCS)** (Time-stepping scheme for incompressible Navier-Stokes)

### **Simulation Parameters**
| Parameter       | Value   | Description |
|----------------|--------|-------------|
| `T`            | 4.0    | Total simulation time |
| `num_steps`    | 500    | Number of time steps |
| `dt`           | 0.008  | Time step size (`T / num_steps`) |
| `Re`           | 400    | Reynolds number |
| `nu`           | 0.0025 | Kinematic viscosity (`1.0 / Re`) |

### **Mesh & Function Spaces**
- **Mesh:** Unit square (`40×40` elements)
- **Velocity Space:** Quadratic Lagrange (`P2`)
- **Pressure Space:** Linear Lagrange (`P1`)

---

## **3. Project Goals**
The primary objectives of this project are:
1. **Simulate lid-driven cavity flow** using the incompressible Navier-Stokes equations.
2. **Implement the IPCS scheme** for stable time-stepping.
3. **Visualize velocity, pressure, and streamlines** to analyze flow patterns.
4. **Generate an animated GIF** showing the evolution of the flow.
5. **Validate the results** against known benchmark solutions for cavity flow.

---

## **4. Logic and Methodology**
### **Governing Equations**
The **incompressible Navier-Stokes equations** are solved:
\[
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\nabla p + \nu \nabla^2 \mathbf{u}
\]
\[
\nabla \cdot \mathbf{u} = 0
\]
where:
- \(\mathbf{u}\) = velocity field
- \(p\) = pressure
- \(\nu\) = kinematic viscosity

### **Numerical Scheme (IPCS)**
1. **Step 1:** Solve for intermediate velocity (ignoring pressure):
   \[
   \frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} + (\mathbf{u}^n \cdot \nabla) \mathbf{u}^n = \nu \nabla^2 \mathbf{u}^*
   \]
2. **Step 2:** Solve for pressure correction:
   \[
   \nabla^2 p^{n+1} = \frac{\nabla \cdot \mathbf{u}^*}{\Delta t}
   \]
3. **Step 3:** Correct velocity to enforce incompressibility:
   \[
   \mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t \nabla p^{n+1}
   \]

### **Boundary Conditions**
- **Lid (Top Boundary):** \(\mathbf{u} = (1, 0)\) (moving wall)
- **Other Walls:** \(\mathbf{u} = (0, 0)\) (no-slip condition)

---

## **5. Expected Results**
1. **Velocity Field:** Formation of a primary vortex and secondary vortices at the corners.
2. **Pressure Distribution:** High pressure at the top-right corner due to flow turning.
3. **Streamlines:** Closed loops indicating recirculation zones.
4. **Animated GIF:** Shows the transient evolution of the flow.

### **Benchmark Validation**
- Compare with **Ghia et al. (1982)** benchmark data for Re=400.
- Expected agreement in vortex locations and velocity profiles.

---

## **6. Suggestions for Improvement**
### **Numerical Enhancements**
- **Adaptive Mesh Refinement (AMR):** Improve resolution near vortices.
- **Higher-Order Elements:** Use `P3` for velocity for better accuracy.
- **Parallel Computing:** Use MPI for large-scale simulations.

### **Physics Extensions**
- **Turbulence Modeling:** Add LES or RANS for higher Re.
- **Temperature Effects:** Include Boussinesq approximation for natural convection.

### **Visualization Improvements**
- **3D Visualization:** Use Paraview for 3D post-processing.
- **Quantitative Analysis:** Plot velocity profiles along centerlines.

### **Code Optimization**
- **Preconditioners:** Use PETSc for faster linear solves.
- **Memory Efficiency:** Store only key time steps for large simulations.

---

## **7. Conclusion**
This project successfully implements a **lid-driven cavity flow solver** using FEniCS. The IPCS scheme ensures stability, and the results align with expected fluid dynamics behavior. Future work could extend the solver to more complex geometries or turbulent flows.

---

## **8. How to Run the Code**
1. **Install Dependencies:**
   ```bash
   pip install fenics matplotlib numpy imageio
   ```
2. **Run the Simulation:**
   ```bash
   python lid_driven_cavity.py
   ```
3. **View Results:**
   - Check `results_cavity/` for output images and GIF.

---

## **9. References**
- Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method.* Journal of Computational Physics.
- Logg, A., Mardal, K. A., & Wells, G. N. (2012). *Automated Solution of Differential Equations by the Finite Element Method.* Springer.

---

### **License**
This project is open-source under the **MIT License**. Contributions are welcome!  

🚀 **Happy Coding!** 🚀
