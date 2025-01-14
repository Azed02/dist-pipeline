from argparse import ArgumentParser, Namespace
from datetime import datetime # default gif name


_DEFAULTS = {
    "Nx": 41,
    "Ny": 41,
    "Lx": 1.0,
    "Ly": 1.0,
    "dx": 0.025,  # Lx / (N_POINTS_X - 1)
    "dy": 0.025,  # Ly / (N_POINTS_Y - 1)
    "KINEMATIC_VISCOSITY": 0.1,
    "DENSITY": 1.0,
    # "N_ITERATIONS": 40,
    "N_PRESSURE_POISSON_ITERATIONS": 30,
    "gif": False,
    "gif_name": lambda steps: f"anim_{datetime.today().strftime('%y%m%d_%H%M')}_steps-{steps}.gif",
    "plot": False,
    "plot_name": lambda steps: f"plot_{datetime.today().strftime('%y%m%d_%H%M')}_steps-{steps}.png",
}


def create_argparser() -> ArgumentParser:
    argparser = ArgumentParser(
        prog="navier_stokes_solver",
        description="Implementation of a global NS 2D solver using Chorin's projection method",
    )
    argparser.add_argument("T", type=float, help="Ending time of the simulation")
    argparser.add_argument(
        "DT",
        nargs="?",
        type=float,
        help="Stepping time of the simulation (/!\\ careful with this one)",
    )
    argparser.add_argument("--Lx", type=float, help="Width of the simulation")
    argparser.add_argument("--Ly", type=float, help="Height of the simulation")
    argparser.add_argument(
        "--Nx",
        type=int,
        help="Number of points used in the x-axis space discretization",
    )
    argparser.add_argument(
        "--Ny",
        type=int,
        help="Number of points used in the y-axis space discretization",
    )
    argparser.add_argument(
        "--dx",
        type=float,
        help="Step used in the x-axis space discretization"
    )
    argparser.add_argument(
        "--dy",
        type=float,
        help="Step used in the y-axis space discretization"
    )
    argparser.add_argument(
        "-I",
        "--iter",
        dest="N_ITERATIONS",
        type=int,
        help="Number of simulation's iterations",
    )
    argparser.add_argument(
        "-p",
        "--pressure-iter",
        dest="N_PRESSURE_POISSON_ITERATIONS",
        type=int,
        help="Number of iterations to resolve next step pressure",
    )
    argparser.add_argument(
        "--nu",
        "--kinematic-viscosity",
        "--viscosity",
        dest="KINEMATIC_VISCOSITY",
        type=float,
        help="Fluid's kinematic viscosity",
    )
    argparser.add_argument(
        "--rho", "--density", dest="DENSITY", type=float, help="Fluid's density"
    )
    argparser.add_argument("--gif", action="store_true", help="Enable gif creation")
    argparser.add_argument(
        "--gif_name",
        help="Gif name (default: 'anim_{date}_{hour}_steps-{N_ITERATIONS}.gif')",
    )
    argparser.add_argument("--plot", action="store_true", help="Enable plot creation of initial and final states")
    argparser.add_argument("--plot_name", help="Plot name (default: plot_{date}_{hour}_steps-{N_ITERATIONS}.gif)")

    return argparser


def parse_args_with_defaults(argparser: ArgumentParser) -> dict[str, int | float]:

    cmd_line_args: Namespace = argparser.parse_args()

    args_with_default: dict = {
        key: vars(cmd_line_args)[key]  # if specified on command-line
        or _DEFAULTS[key]  # else, default value
        for key in _DEFAULTS.keys()
    }
    args_with_default["T"] = cmd_line_args.T
    args_with_default["DT"] = cmd_line_args.DT or None

    return args_with_default


if __name__ == "__main__":
    argparser = create_argparser()

    args = parse_args_with_defaults(argparser)

    print(args)

