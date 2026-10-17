# One image per conda environment in envs/, built from the same file the conda
# profile uses, so both ways of running the pipeline get the same packages.
# There are two because only the sgRNA steps need Java (for FlashFry) and
# pyfaidx; the GWAS, LD and overlap steps run in a smaller image without them.
#
#   docker build -f containers/Dockerfile --build-arg ENV=gwas  -t cause-nf/gwas:0.4.0  envs
#   docker build -f containers/Dockerfile --build-arg ENV=sgrna -t cause-nf/sgrna:0.4.0 envs
#
# The build context is envs/, so the large work/ and genome_cache/ folders are
# never sent to Docker.
FROM mambaorg/micromamba:1.5.10

ARG ENV
COPY --chown=$MAMBA_USER:$MAMBA_USER ${ENV}.yml /tmp/env.yml
# procps gives Nextflow the ps command it reads task CPU and memory use with
RUN micromamba install -y -n base -f /tmp/env.yml procps-ng && micromamba clean -a -y

# Nextflow runs tasks with bash and bypasses the image's entrypoint, so the
# environment's tools have to be on PATH without activating it
ENV PATH=/opt/conda/bin:$PATH
