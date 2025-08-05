import os
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import imageio

# --- 1. Simulation Parameters ---
T = 4.0             # Final time
num_steps = 500     # Number of time steps
dt = T / num_steps  # Time step size
Re = 400.0          # Reynolds number
nu = 1.0 / Re       # Kinematic viscosity

# Create a directory to save GIF frames
if not os.path.exists('results_cavity'):
    os.makedirs('results_cavity')
    
# --- 2. Mesh and Function Spaces ---
mesh = UnitSquareMesh(40, 40)
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

# --- 3. Boundary Conditions ---
class Lid(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and on_boundary

class Walls(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[1], 0.0) or near(x[0], 0.0) or near(x[0], 1.0)) and on_boundary

lid_velocity = Constant((1.0, 0.0))
noslip_velocity = Constant((0.0, 0.0))

bc_lid = DirichletBC(V, lid_velocity, Lid())
bc_walls = DirichletBC(V, noslip_velocity, Walls())
bcu = [bc_walls, bc_lid]

# --- 4. Define Trial and Test Functions ---
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# --- 5. Define Functions for Solutions ---
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# --- 6. Define Variational Problem (IPCS Scheme) ---
U = 0.5 * (u + u_n)
F1 = (1/dt)*inner(u - u_n, v)*dx \
   + inner(dot(u_n, nabla_grad(u_n)), v)*dx \
   + nu*inner(grad(U), grad(v))*dx \
   - inner(p_n, div(v))*dx
a1 = lhs(F1)
L1 = rhs(F1)

a2 = inner(grad(p), grad(q))*dx
L2 = -(1/dt)*div(u_)*q*dx

a3 = inner(u, v)*dx
L3 = inner(u_, v)*dx - dt*inner(grad(p_), v)*dx

# --- 7. Assemble Matrices ---
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)
[bc.apply(A1) for bc in bcu]

# --- 8. Time-Stepping Loop ---
t = 0
frame_count = 0
filenames = []

# --- PROGRESS BAR CORRECTION ---
# 1. Initialize with the total number of steps
progress = Progress('Time-stepping', num_steps)
set_log_level(LogLevel.PROGRESS)

for n in range(num_steps):
    t += dt
    
    # Solve system
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1)
    
    b2 = assemble(L2)
    solve(A2, p_.vector(), b2)
    
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3)
    
    u_n.assign(u_)
    p_n.assign(p_)

    # Save frame for GIF
    if n % 5 == 0:
        fig = plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plot(u_, title=f'Velocity Field at t={t:.2f}s')
        plt.xlabel('x'); plt.ylabel('y')
        plt.subplot(1, 2, 2)
        c = plot(p_, title=f'Pressure Field at t={t:.2f}s')
        plt.colorbar(c)
        plt.xlabel('x'); plt.ylabel('y')
        plt.tight_layout()
        filename = f'results_cavity/frame_{frame_count:03d}.png'
        filenames.append(filename)
        plt.savefig(filename)
        plt.close(fig)
        frame_count += 1
        
    # --- PROGRESS BAR CORRECTION ---
    # 2. Call update() with no arguments to increment
    #progress.update()

# --- 9. Create GIF ---
print("\nCreating GIF animation...")
with imageio.get_writer('results_cavity/lid_driven_cavity.gif', mode='I', duration=0.1) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
print("GIF saved as 'results_cavity/lid_driven_cavity.gif'")

# --- 10. Post-processing and Final Plots ---
# (The rest of the script is unchanged and should work as intended)
print("Generating final plots...")
u_sol, p_sol = u_.split(deepcopy=True)

plt.figure(figsize=(6, 5))
c = plot(u_sol, title=f'u-velocity at t={T}s (Re={Re})')
plt.colorbar(c, label='u')
plt.xlabel('x'); plt.ylabel('y')
plt.savefig('results_cavity/u_velocity.png')

plt.figure(figsize=(6, 5))
c = plot(p_sol, title=f'v-velocity at t={T}s (Re={Re})')
plt.colorbar(c, label='v')
plt.xlabel('x'); plt.ylabel('y')
plt.savefig('results_cavity/v_velocity.png')

# --- 11. Stream Function Calculation and Plotting ---
W = FunctionSpace(mesh, 'P', 1)
w = project(curl(u_), W)
psi = TrialFunction(W)
phi = TestFunction(W)
bc_psi = DirichletBC(W, Constant(0.0), 'on_boundary')
a_psi = inner(grad(psi), grad(phi))*dx
L_psi = -w*phi*dx
psi_sol = Function(W)
solve(a_psi == L_psi, psi_sol, bc_psi)

plt.figure(figsize=(6, 5))
c = plot(psi_sol, title=f'Stream Function ψ (Re={Re})')
plt.colorbar(c, label='ψ')
plt.xlabel('x'); plt.ylabel('y')
plt.savefig('results_cavity/stream_function.png')

plt.figure(figsize=(6, 5))
plt.title(f'Stream Lines (Re={Re})')
contour = plot(psi_sol, 20)
plt.colorbar(contour, label='ψ (Streamlines)')
plt.xlabel('x'); plt.ylabel('y')
plt.savefig('results_cavity/stream_lines.png')

plt.show()
print("\nSimulation and visualization complete.")
