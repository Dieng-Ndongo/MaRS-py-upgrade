FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update -y --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
        wget unzip git curl build-essential less vim && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /pipeline

COPY environment.yml /tmp/environment.yml

RUN micromamba create -y \
        -n pipeline_env \
        --strict-channel-priority \
        -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Make var2vcf_valid.pl available in the Conda environment PATH
RUN VAR2VCF=$(find /opt/conda/envs/pipeline_env/share \
        -maxdepth 2 \
        -iname 'var2vcf_valid.pl' \
        -print -quit) && \
    test -n "$VAR2VCF" && \
    chmod +x "$VAR2VCF" && \
    ln -sf "$VAR2VCF" \
        /opt/conda/envs/pipeline_env/bin/var2vcf_valid.pl

# Conda environment executables
ENV PATH="/opt/conda/envs/pipeline_env/bin:${PATH}"

# Host/container UID/GID matching
ARG PUID=1000
ARG PGID=1000

RUN groupadd -g ${PGID} pipelineuser && \
    useradd -m -u ${PUID} -g ${PGID} pipelineuser && \
    mkdir -p /pipeline/output /pipeline/logs && \
    chown -R ${PUID}:${PGID} /pipeline

USER pipelineuser

CMD ["micromamba", "run", "-n", "pipeline_env", "python", "/app/pipeline_python.py"]
