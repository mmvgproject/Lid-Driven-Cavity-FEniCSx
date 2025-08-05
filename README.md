# Lid-Driven Cavity Simulation using FEniCS and IPCS Scheme

## Project Title
**Lid-Driven Cavity Simulation with FEniCS and IPCS Scheme**

This project implements a numerical simulation of the lid-driven cavity flow, a benchmark problem in computational fluid dynamics (CFD), using the **FEniCS** finite element library. The simulation employs the **Incremental Pressure Correction Scheme (IPCS)** to solve the 2D incompressible Navier-Stokes equations for a square cavity with a moving lid. The code generates velocity, pressure, and stream function fields, along with visualizations and an animated GIF to illustrate the transient flow evolution.

---

## Project and Program Specifications

### Overview
The lid-driven cavity problem models the flow of an incompressible fluid in a unit square domain where the top boundary (lid) moves at a constant velocity, while the other walls are stationary. This creates a complex flow with vortices driven by viscous effects and pressure gradients, making it a standard test case for CFD solvers.

### Program Specifications
- **Programming Language**: Python 3.x
- **Key Libraries**:
  - **FEniCS (DOLFIN)**: For solving partial differential equations (PDEs) using the finite element method.
  - **NumPy**: For numerical operations and array handling.
  - **Matplotlib**: For plotting velocity, pressure, and streamline fields.
  - **ImageIO**: For generating animated GIFs of the solution evolution.
- **Algorithm**: Incremental Pressure Correction Scheme (IPCS).
- **Numerical Method**: Finite element method with P2 elements for velocity and P1 elements for pressure.
- **Boundary Conditions**:
  - Top wall (lid): Velocity \( u = (1.0, 0.0) \).
  - Other walls (left, right, bottom): No-slip condition \( u = (0.0, 0.0) \).
  - Stream function: Zero on all boundaries.
- **Physical Parameters**:
  - Domain: Unit square (\( 1.0 \times 1.0 \)).
  - Final time: \( T = 4.0 \, \text{s} \).
  - Reynolds number: \( \text{Re} = 400 \).
  - Kinematic viscosity: \( \nu = \frac{1}{\text{Re}} \).
- **Numerical Parameters**:
  - Mesh: \( 40 \times 40 \) elements.
  - Time steps: 500 (\( \Delta t = \frac{T}{\text{num_steps}} = 0.008 \, \text{s} \)).
- **Output**:
  - Contour plots for \( u \)-velocity, \( v \)-velocity, pressure, stream function, and streamlines.
  - Animated GIF showing the evolution of velocity and pressure fields over time.
  - Results saved in the `results_cavity/` directory.

### Dependencies
To run the code, install the required Python packages and FEniCS:
```bash
pip install fenics numpy matplotlib imageio
```
Note: FEniCS installation may vary by platform. Refer to the [FEniCS documentation](https://fenicsproject.org/) for detailed instructions (e.g., using Docker, Ubuntu packages, or Conda).

---

## Project Goals
The primary goals of this project are:
1. **Accurate Simulation**: Solve the transient 2D incompressible Navier-Stokes equations for the lid-driven cavity problem using the IPCS scheme.
2. **Visualization**: Generate high-quality plots of velocity, pressure, and stream function fields to analyze flow behavior.
3. **Educational Tool**: Provide a clear, well-documented implementation of the IPCS scheme in FEniCS for academic and research purposes.
4. **Transient Analysis**: Capture the time-dependent evolution of the flow through an animated GIF.
5. **Reproducibility**: Ensure the code is modular, portable, and suitable for CFD studies.

---

## Logic Used in the Project

### Overview of the IPCS Scheme
The Incremental Pressure Correction Scheme (IPCS) is a fractional-step method for solving the incompressible Navier-Stokes equations. It decouples the velocity and pressure computations to improve efficiency. The key steps are:
1. **Tentative Velocity Step**: Solve the momentum equation with the previous pressure to compute an intermediate velocity.
2. **Pressure Correction Step**: Solve a Poisson equation for the pressure to enforce continuity (\( \nabla \cdot \mathbf{u} = 0 \)).
3. **Velocity Update**: Correct the velocity using the updated pressure field.

### Code Structure
The code is organized into 11 main sections:
1. **Simulation Parameters**: Defines time, time steps, Reynolds number, and viscosity.
2. **Mesh and Function Spaces**: Creates a \( 40 \times 40 \) mesh and defines P2 (velocity) and P1 (pressure) function spaces.
3. **Boundary Conditions**: Applies the lid velocity and no-slip conditions using `SubDomain` classes.
4. **Trial and Test Functions**: Defines trial and test functions for velocity (\( u, v \)) and pressure (\( p, q \)).
5. **Solution Functions**: Initializes functions for current and previous time steps (\( u_n, u_, p_n, p_ \)).
6. **Variational Problem (IPCS)**: Sets up the weak forms for the tentative velocity, pressure correction, and velocity update steps.
7. **Assemble Matrices**: Pre-assembles the left-hand side matrices for efficiency.
8. **Time-Stepping Loop**: Advances the solution in time, solving the three IPCS steps and saving frames for visualization.
9. **GIF Creation**: Combines saved frames into an animated GIF.
10. **Post-Processing**: Generates final plots for velocity and pressure components.
11. **Stream Function Calculation**: Computes and plots the stream function and streamlines using the vorticity.

### Key Implementation Details
- **Finite Element Method**: FEniCS uses P2 elements for velocity and P1 elements for pressure to satisfy the Ladyzhenskaya-Babuška-Brezzi (LBB) condition for stability.
- **IPCS Scheme**: Splits the Navier-Stokes equations into three linear systems, solved sequentially at each time step.
- **Progress Bar**: Uses FEniCS’s `Progress` class to display simulation progress (corrected to use `update()`).
- **Stream Function**: Computed by solving a Poisson equation with the vorticity (\( \omega = \nabla \times \mathbf{u} \)) as the source term.
- **Visualization**: Saves frames every 5 time steps for the GIF and final plots for detailed analysis.

---

## Theory Used in the Project

### Governing Equations
The lid-driven cavity flow is governed by the 2D incompressible Navier-Stokes equations:
1. **Continuity Equation**:
   \[
   \nabla \cdot \mathbf{u} = 0
   \]
   where \( \mathbf{u} = (u, v) \) is the velocity vector.
2. **Momentum Equations**:
   \[
   \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\nabla p + \nu \nabla^2 \mathbf{u}
   \]
   where \( p \) is the pressure divided by density, and \( \nu \) is the kinematic viscosity.

### IPCS Scheme Theory
The IPCS scheme is a time-stepping method that splits the Navier-Stokes equations into:
1. **Tentative Velocity**:
   \[
   \frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} + (\mathbf{u}^n \cdot \nabla) \mathbf{u}^n = \nu \nabla^2 \left( \frac{\mathbf{u}^* + \mathbf{u}^n}{2} \right) - \nabla p^n
   \]
2. **Pressure Correction**:
   \[
   \nabla^2 p^{n+1} = \nabla^2 p^n - \frac{1}{\Delta t} \nabla \cdot \mathbf{u}^*
   \]
3. **Velocity Update**:
   \[
   \mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t \nabla p^{n+1}
   \]
This approach ensures the velocity field is divergence-free at each time step.

### Stream Function and Vorticity
The stream function (\( \psi \)) is defined such that:
\[
u = \frac{\partial \psi}{\partial y}, \quad v = -\frac{\partial \psi}{\partial x}
\]
The vorticity (\( \omega = \nabla \times \mathbf{u} \)) is computed using the `curl` operator, and the stream function satisfies:
\[
\nabla^2 \psi = -\omega
\]
This Poisson equation is solved with zero boundary conditions to visualize the flow structure.

### Reynolds Number
The Reynolds number (\( \text{Re} = \frac{U_{\text{lid}} L}{\nu} = 400 \)) indicates a laminar flow regime with a primary vortex and secondary vortices in the corners, which the simulation captures over time.

---

## Expected Results
The simulation produces the following outputs:
1. **Velocity Fields**:
   - \( u \)-velocity: Strong horizontal flow near the lid, diminishing toward the bottom.
   - \( v \)-velocity: Vertical components forming a circulatory pattern.
2. **Pressure Field**: Smooth distribution with gradients driving the flow.
3. **Stream Function and Streamlines**: A primary vortex in the cavity center, with secondary vortices in the lower corners for \( \text{Re} = 400 \).
4. **Animated GIF**: Shows the transient evolution of velocity and pressure fields over \( t = 0 \) to \( t = 4.0 \, \text{s} \).
5. **Final Plots**: Detailed visualizations of \( u \)-velocity, \( v \)-velocity, pressure, stream function, and streamlines at \( t = 4.0 \, \text{s} \).

The results should exhibit the characteristic vortex structure of the lid-driven cavity flow, consistent with benchmark solutions (e.g., Ghia et al., 1982).

---

## Suggestions for Improving the Project
1. **Mesh Refinement**:
   - Increase the mesh resolution (e.g., \( 80 \times 80 \)) to capture finer flow details.
   - Use adaptive mesh refinement near the walls to resolve boundary layers.
2. **Higher Reynolds Numbers**:
   - Test higher \( \text{Re} \) (e.g., 1000, 3200) to study transitional or unsteady flows, adjusting the time step for stability.
3. **Alternative Schemes**:
   - Implement other fractional-step methods (e.g., Chorin’s projection method) for comparison.
   - Use higher-order time integration (e.g., Adams-Bashforth for convection).
4. **Parallel Computing**:
   - Leverage FEniCS’s MPI support for parallel simulations on larger meshes.
5. **Additional Visualizations**:
   - Add vorticity plots to highlight rotational structures.
   - Include velocity vector plots alongside streamlines.
6. **Quantitative Validation**:
   - Compare velocity profiles along the cavity centerlines with benchmark data (e.g., Ghia et al., 1982).
7. **Error Handling**:
   - Add checks for numerical stability (e.g., CFL condition) and solver convergence.
8. **Interactive Interface**:
   - Integrate with Jupyter Notebook or a GUI (e.g., ParaView) for interactive parameter tuning and visualization.
9. **File Management**:
   - Add cleanup of temporary frame files after GIF creation.
   - Save solution data in XDMF/HDF5 format for post-processing in ParaView.

---

## Installation and Usage
1. **Clone the Repository**:
   ```bash
   git clone <repository_url>
   cd lid-driven-cavity-fenics
   ```
2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
   Create a `requirements.txt` with:
   ```
   fenics==2019.1.0
   numpy==1.26.4
   matplotlib==3.8.3
   imageio==2.34.0
   ```
   Install FEniCS separately if needed (e.g., via Docker: `docker pull fenicsproject/stable`).
3. **Run the Simulation**:
   ```bash
   python lid_driven_cavity_ipcs.py
   ```
4. **Outputs**:
   - Directory: `results_cavity/`
   - Plots: `u_velocity.png`, `v_velocity.png`, `stream_function.png`, `stream_lines.png`.
   - Animated GIF: `lid_driven_cavity.gif`.

---

## Contributing
Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make changes and commit (`git commit -m "Add feature"`).
4. Push to the branch (`git push origin feature-branch`).
5. Open a pull request.

---

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Acknowledgments
- **FEniCS Project**: For providing a powerful finite element framework.
- **CFD Community**: For benchmark data and literature on the lid-driven cavity problem.
- **Ghia et al. (1982)**: For reference solutions used in validation.

---

## Contact
For questions or suggestions, please open an issue on GitHub .
