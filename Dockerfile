FROM mambaorg/micromamba:1.5.8

USER root
RUN apt-get update -y --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
        wget unzip git curl build-essential default-jdk less vim && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n pipeline_env -f /tmp/environment.yml && \
    micromamba clean --all --yes

# utilisateur non-root
RUN useradd -m -u 1000 pipelineuser
USER pipelineuser

# Executer le pipeline
CMD ["micromamba", "run", "-n", "pipeline_env", "python", "pipeline_python.py"]
