# ===== Base image Micromamba =====
FROM mambaorg/micromamba:1.5.8

# ===== Installer les paquets système =====
USER root
RUN apt-get update -y --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
        wget unzip git curl build-essential default-jdk less vim && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ===== Définir le répertoire de travail =====
WORKDIR /app

# ===== Copier environment.yml et créer l'environnement =====
COPY environment.yml /app/environment.yml
RUN micromamba create -y -n pipeline_env -f /app/environment.yml && \
    micromamba clean --all --yes

# ===== Copier tout le pipeline =====
COPY . /app

# ===== Créer dossiers nécessaires et donner les droits =====
RUN mkdir -p /app/output /app/logs /app/data && \
    chmod -R 777 /app/output /app/logs /app/data

# ===== Utiliser l'environnement conda pour la suite =====
SHELL ["micromamba", "run", "-n", "pipeline_env", "/bin/bash", "-c"]

# ===== Point d'entrée =====
ENTRYPOINT ["micromamba", "run", "-n", "pipeline_env", "python", "pipeline_python.py"]

