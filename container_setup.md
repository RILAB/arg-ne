## Container Setup (alternative)

Instead of a local Conda environment, you can use Singularity or Docker. Both
containers bundle all dependencies, the pipeline code, example data, and the
full test suite. Tests are run automatically during the build so the image is
verified before it is produced.

### Singularity

```bash
singularity build --fakeroot argprep.sif singularity.def
```

Run the pipeline on your data (from your project directory):

```bash
singularity exec argprep.sif snakemake -j 8 \
  --snakefile /opt/argprep/Snakefile \
  --configfile options.yaml
```

Run the bundled example:

```bash
singularity exec argprep.sif snakemake -j 4 \
  --snakefile /opt/argprep/Snakefile \
  --configfile /opt/argprep/example_data/options.yaml \
  --config results_dir=example_results
```

Run tests inside the container:

```bash
singularity exec argprep.sif bash -c 'cd /opt/argprep && pytest -q -p no:cacheprovider'
```

Our singularity setup defaults to `/tmp` but for HPC users that want to use some variable made available by their HPC admins like `$SCRATCH` then you could re-write your ENV variable:

```bash
singularity exec --env XDG_CACHE_HOME=$SCRATCH/argprep_cache argprep.sif
```

### Docker

```bash
docker build -t argprep .
```

Run the pipeline (mount your data directory):

```bash
docker run --rm -v $(pwd):/data -w /data argprep \
  -j 8 --snakefile /opt/argprep/Snakefile --configfile options.yaml
```

Run tests:

```bash
docker run --rm --entrypoint pytest argprep -q -p no:cacheprovider
```
