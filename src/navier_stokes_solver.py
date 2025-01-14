#!/usr/bin/env python3

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation

from argparse import ArgumentParser # for type hint
import argument_parsing


def central_diff_x(f, dx):
    cdx = np.zeros_like(f)
    cdx[1:-1, 1:-1] = (f[1:-1, 2:] - f[1:-1, 0:-2]) / (2 * dx)
    return cdx


def central_diff_y(f, dy):
    cdy = np.zeros_like(f)
    cdy[1:-1, 1:-1] = (f[2:, 1:-1] - f[0:-2, 1:-1]) / (2 * dy)
    return cdy


def laplace(f, dx, dy):
    """five point stencil"""

    _laplace = np.zeros_like(f)
    _laplace[1:-1, 1:-1] = (
          f[1:-1, 0:-2]     # f(x - dx, y)
        + f[1:-1, 2:  ]     # f(x + dx, y)
        + f[0:-2, 1:-1]     # f(x, y - dy)
        + f[2:  , 1:-1]     # f(x, y + dy)
        - 4 * f[1:-1, 1:-1] # f(x, y)
    ) / (dx * dy)
    return _laplace


def open_top_lid_BC_velocity(u, v):
    _u = u.copy()
    _u[ -1, : ] = 1.0
    return _u, v


def left_right_fight_BC_velocity(u, v):
    _u = u.copy()
    _u[ : , 0 ] =  1.0 # left
    _u[ : , -1] = -1.0 # right
    return _u, v


def homogeneous_boundary_condition_velocity(u, v):
    _u, _v = u.copy(), v.copy()
    _u[ 0 , : ] = 0.0
    _u[-1 , : ] = 0.0
    _u[ : ,  0] = 0.0
    _u[ : , -1] = 0.0
    _v[ 0 , : ] = 0.0
    _v[-1 , : ] = 0.0
    _v[ : ,  0] = 0.0
    _v[ : , -1] = 0.0
    return _u, _v


# BOUNDARY_CONDITION = left_right_fight_BC_velocity
BOUNDARY_CONDITION = open_top_lid_BC_velocity

def main():
    # Command-line argument parser
    argparser: ArgumentParser = argument_parsing.create_argparser()

    args: dict[str, int | float] = argument_parsing.parse_args_with_defaults(argparser)

    # shortcuts
    dx, dy = args["dx"], args["dy"]
    T, DT = args["T"], args.get("DT", None)
    Lx, Ly = args["Lx"], args["Ly"]
    Nx, Ny = args["Nx"], args["Ny"]
    N_ITERATIONS = args.get("N_ITERATIONS", None)
    KINEMATIC_VISCOSITY = args["KINEMATIC_VISCOSITY"]
    DENSITY = args["DENSITY"]

    # time discretization
    assert T is not None, "T should not be None"
    if T and DT and N_ITERATIONS:
        raise AttributeError("T, DT and N_ITERATIONS cannot be all set at the same time.")
    elif T and DT:
        assert N_ITERATIONS is None
        N_ITERATIONS: int = int(T / DT)
    elif T and N_ITERATIONS:
        assert DT is None
        DT: float = T / N_ITERATIONS

    # space discretization
    dx, dy = Lx / (Nx - 1), Ly / (Ny - 1) # force adequate values for dx, dy
    if dx != args["dx"]:
        raise ValueError(f"dx value error: found specified dx={args['dx']} but needs dx={dx} (Lx / (Nx - 1))")
    if dy != args["dy"]:
        raise ValueError(f"dy value error: found specified dy={args['dy']} but needs dy={dy} (Ly / (Ny - 1))")
    xs = np.linspace(0.0, Lx, Nx) # [ 0.0 ... Lx ]
    ys = np.linspace(0.0, Ly, Ny) # [ 0.0 ... Ly ]
    X, Y = np.meshgrid(xs, ys)

    # Initial Conditions
    u_prev = np.zeros_like(X)
    v_prev = np.zeros_like(X)
    p_prev = np.zeros_like(X)

    u_initial = u_prev
    v_initial = v_prev
    p_initial = p_prev

    # used to store state variables during simulation
    u_list, v_list, p_list = [], [], []

    # See [Chorin's projection method](https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics))
    for step in range(1, N_ITERATIONS + 1):
        if DT > 0.5 * dx ** 2 / KINEMATIC_VISCOSITY:
            raise ValueError("time step too large (DT > 0.5 * dx^2 * KINEMATIC_VISCOSITY)")

        d_u_prev__d_x = central_diff_x(u_prev, dx)
        d_u_prev__d_y = central_diff_y(u_prev, dy)
        d_v_prev__d_x = central_diff_x(v_prev, dx)
        d_v_prev__d_y = central_diff_y(v_prev, dy)
        laplace_u_prev = laplace(u_prev, dx, dy)
        laplace_v_prev = laplace(v_prev, dx, dy)

        # (1) Prediction step: guess velocity using viscous force
        u_star = u_prev + DT * ( - ( (u_prev * d_u_prev__d_x) - (v_prev * d_v_prev__d_y) )
                                 + KINEMATIC_VISCOSITY * laplace_u_prev
                               )
        v_star = v_prev + DT * ( - ( (u_prev * d_u_prev__d_x) - (v_prev * d_v_prev__d_y) )
                                 + KINEMATIC_VISCOSITY * laplace_v_prev
                               )

        # Homogeneous Dirichlet Boundary conditions
        u_star, v_star = homogeneous_boundary_condition_velocity(u_star, v_star)
        u_star, v_star = BOUNDARY_CONDITION(u_star, v_star)

        # (2) Solving pressure at n + 1 using 5 points star stencil diffusion approximation
        d_u_star__d_x = central_diff_x(u_star, dx)
        d_v_star__d_y = central_diff_y(v_star, dy)
        rhs = KINEMATIC_VISCOSITY / DT * ( d_u_star__d_x + d_v_star__d_y ) # rhs for diffusion

        p_prev_copy = p_prev.copy()
        for _ in range(args["N_PRESSURE_POISSON_ITERATIONS"]):
            p_next = np.zeros_like(p_prev)
            p_next[1:-1, 1:-1] = 1 / 4 * (
                  p_prev[1:-1, 0:-2] # left
                + p_prev[1:-1, 2:  ] # right
                + p_prev[0:-2, 1:-1] # top
                + p_prev[2:  , 1:-1] # bottom
                - dx * dy * rhs[1:-1, 1:-1]
            )

            # Neumann Boundary Conditions for pressure
            p_next[ : , -1] = p_next[ : , -2]
            p_next[ : ,  0] = p_next[ : ,  1]
            p_next[-1 , : ] = p_next[-2 , : ]
            p_next[ 0 , : ] = p_next[ 1 , : ]

            p_prev = p_next
        # print(f"delta poisson ({N_PRESSURE_POISSON_ITERATIONS} iters) = {np.linalg.norm(p_next - p_prev_copy, 2)}")

        # (3) Projection step:
        d_p_next__d_x = central_diff_x(p_next, dx)
        d_p_next__d_y = central_diff_y(p_next, dy)

        u_next = u_star - DT / DENSITY * d_p_next__d_x
        v_next = v_star - DT / DENSITY * d_p_next__d_y

        # Enforce boundary conditions
        u_next, v_next = homogeneous_boundary_condition_velocity(u_next, v_next)
        u_next, v_next = BOUNDARY_CONDITION(u_next, v_next)

        # Adaptative time stepping
        # dist = np.linalg.norm(u_next - u_prev, 2)
        # DT = max(0.5 * dx ** 2 / KINEMATIC_VISCOSITY, 0.1 - dist)
        # print("\t", DT)

        u_prev = u_next
        v_prev = v_next
        p_prev = p_next

        assert N_ITERATIONS > 10, "N_ITERATIONS must be greater than 10."
        # if step % 10 == 0:
        print(f"  │ step {step:>5}/{N_ITERATIONS}")

        u_list.append(u_next.copy())
        v_list.append(v_next.copy())
        p_list.append(p_next.copy())


    if args["plot"]:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=100)

        # Initial state plot
        ax1.set_title("Champ de vitesse et pression initiaux")
        im1 = ax1.imshow(p_initial, extent=[0, Lx, 0, Ly], origin='lower', cmap='viridis', alpha=0.8)
        # quiver1 = ax1.quiver(X[::2], Y[::2], u_initial[::2], v_initial[::2], scale=1, scale_units='xy', color='white')
        ax1.streamplot(X, Y, u_initial, v_initial, color="black")
        fig.colorbar(im1, ax=ax1, label="Pression")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")

        # Final state plot
        ax2.set_title("Champ de vitesse et pression finaux")
        im2 = ax2.imshow(p_next, extent=[0, Lx, 0, Ly], origin='lower', cmap='viridis', alpha=0.8)
        # quiver2 = ax2.quiver(X[::2], Y[::2], u_next[::2], v_next[::2], scale=5, scale_units='xy', color='white')
        ax2.streamplot(X, Y, u_next, v_next, color="black")
        fig.colorbar(im2, ax=ax2, label="Pression")
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")

        # Adjust layout and display
        plt.tight_layout()
        # plt.show()
        filename: str = args["plot_name"](N_ITERATIONS)
        plt.savefig(filename)

        print(f"Genereted initial / final states plot in '{filename}'")


    if args["gif"]:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        # Initial pressure field plot
        pressure_im = ax.imshow(p_list[0], extent=[0, args["Lx"], 0, args["Ly"]], origin='lower', cmap='viridis', alpha=0.8)
        plt.colorbar(pressure_im, ax=ax, label="Pression")

        # Initial velocity field plot (quiver)
        quiver = ax.quiver(X[::2], Y[::2], u_list[0][::2], v_list[0][::2], scale=1, scale_units='xy', color='black')

        # Animation function
        def animate(i):
            # Update the quiver plot
            quiver.set_UVC(u_list[i][::2], v_list[i][::2])

            # Update the pressure heatmap
            # print(f"Pressure at step {i}: {p_list[i]}")
            pressure_im.set_array(p_list[i])  # Update the pressure field array

            # Update the title
            ax.set_title(f"Champ de vitesse et pression - Iter {i}")

            return quiver, pressure_im  # Return both for blitting

        # Create and save the animation
        anim = FuncAnimation(fig, animate, frames=len(u_list), blit=True, interval=100)  # Add interval for speed
        filename: str = args["gif_name"](N_ITERATIONS)
        anim.save(filename, writer="pillow")

        print(f"Generated gif with name '{filename}'")

    print("End")



if __name__ == "__main__":
    main()

