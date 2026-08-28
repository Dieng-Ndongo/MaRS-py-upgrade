FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update -y --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
        wget unzip git curl build-essential default-jdk less vim && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /pipeline

COPY pf_3D7_Ref /pipeline/pf_3D7_Ref
COPY pf_3D7_snpEff_db /pipeline/pf_3D7_snpEff_db
COPY bin /pipeline/bin
COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n pipeline_env --strict-channel-priority -f /tmp/environment.yml && \ 
    micromamba clean --all --yes


# IMPORTANT: dossier de travail writable
RUN mkdir -p /data /data/output /data/logs /data/data && \
    chmod -R 777 /data

# utilisateur non-root — UID/GID paramétrables pour matcher le compte hôte
# qui possède les répertoires bind-montés (ex: build --build-arg PUID=$(id -u mars) PGID=$(id -g mars))
ARG PUID=1000
ARG PGID=1000
RUN groupadd -g ${PGID} pipelineuser && useradd -m -u ${PUID} -g ${PGID} pipelineuser

USER pipelineuser


CMD ["micromamba", "run", "-n", "pipeline_env", "python", "/app/pipeline_python.py"]
