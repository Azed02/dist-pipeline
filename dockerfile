# Utiliser une image Ubuntu avec Python préinstallé
FROM python:3.10-slim

# Installer les dépendances système nécessaires
RUN apt-get update && apt-get install -y \
    build-essential \
    libopenmpi-dev \
    openmpi-bin \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Installer les bibliothèques Python requises
RUN pip install --no-cache-dir numpy scipy matplotlib mpi4py

# Définir le répertoire de travail
WORKDIR /app

# Copier les fichiers de l'application dans l'image
COPY src/argument_parsing.py src/navier_stokes_solver.py /app/

# Définir la commande par défaut pour exécuter l'application avec MPI
CMD ["mpirun", "--allow-run-as-root", "-np", "1", "python3", "navier_stokes_solver.py", "2.0", "0.01", "--Lx", "1.0", "--Ly", "1.0", "--Nx", "50", "--Ny", "50", "--plot"]
