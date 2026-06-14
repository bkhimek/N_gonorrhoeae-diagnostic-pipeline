# ─── N. gonorrhoeae AMR pipeline — single container image ────────────────────
# Tools: AMRFinderPlus 4.0.23 | BLAST+ | MAFFT | samtools | Python 3.10
#
# AMRFinder database is NOT bundled (too large; mount from host at runtime).
# See README for build and run instructions.

FROM mambaorg/micromamba:1.5.8

LABEL maintainer="bkhimek"
LABEL description="N. gonorrhoeae AMR profiling pipeline tools"
LABEL version="1.0"

USER root

# Minimal system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy conda environment spec and install all tools in one layer
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"
ENV LC_ALL=C

WORKDIR /data
CMD ["/bin/bash"]
