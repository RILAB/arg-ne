FROM condaforge/mambaforge:latest

COPY argprep.yml /opt/argprep.yml
RUN mamba env create -f /opt/argprep.yml && mamba clean -afy

ENV PATH="/opt/conda/envs/argprep/bin:$PATH"

COPY . /opt/argprep
WORKDIR /opt/argprep

RUN pytest -q -p no:cacheprovider

ENTRYPOINT ["snakemake"]
