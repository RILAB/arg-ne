#!/bin/bash
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --job-name=argprep-controller
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=logs/slurm/controller-%j.out
# NOTE: --error is intentionally unset so stderr merges into the --output file.
#
# Snakemake CONTROLLER driver job.
# ---------------------------------
# Runs the long-lived snakemake controller as its own SLURM job instead of on
# the head node, then submits the per-rule jobs from there.
#
# Why run the controller on a NON-PREEMPTABLE partition:
#   The controller must outlive every rule job. If it ran on a preemptable
#   partition (e.g. `low`) and got preempted, the whole run would die. The
#   per-rule jobs it submits go to whatever `slurm_partition` you set in your
#   --configfile (set that to `low` to use the cheap preemptable queue); a
#   preempted rule job is auto-requeued by SLURM and the status command
#   (profiles/slurm/status-sacct.sh) keeps the controller waiting for it.
#
# CLUSTER-SPECIFIC DEFAULTS — the ONLY two values here that are not portable are
#   --partition=high   (a non-preemptable partition on the UCD farm cluster)
#   --account=jrigrp
# Override them for another cluster without editing this file, e.g.:
#   sbatch --partition=<your-nonpreemptable> --account=<your-acct> \
#          profiles/slurm/run-controller.sbatch config/options.yaml
# (sbatch CLI flags and the SBATCH_PARTITION / SBATCH_ACCOUNT env vars both
# override the #SBATCH lines above.)
#
# Usage (from the repo root):
#   mkdir -p logs/slurm
#   sbatch profiles/slurm/run-controller.sbatch [CONFIGFILE] [EXTRA snakemake args...]
# Defaults to options.yaml. --rerun-incomplete is always applied.
# Examples:
#   sbatch profiles/slurm/run-controller.sbatch admix/options.yaml
#   sbatch profiles/slurm/run-controller.sbatch admix/options.yaml -n   # dry-run
#
# Watch progress:  tail -f logs/slurm/controller-<jobid>.out

set -euo pipefail

mkdir -p logs/slurm

# First positional arg is the configfile, but only if it isn't itself a flag
# (so `sbatch ... -n` still works without a configfile). Everything left over
# is forwarded to snakemake verbatim.
configfile="options.yaml"
if [[ "${1:-}" && "${1:-}" != -* ]]; then
    configfile="$1"
    shift
fi

# Make snakemake available. On the farm cluster it lives in a conda env; this
# block is best-effort and portable: if snakemake is already on PATH it does
# nothing. Set ARGPREP_CONDA_ENV to force a specific env, or rely on an
# argprep.yml / environment.yml `name:` field if one is present.
if ! command -v snakemake >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load conda 2>/dev/null || true
    env_name="${ARGPREP_CONDA_ENV:-}"
    if [[ -z "$env_name" && -f argprep.yml ]]; then
        env_name="$(awk -F': *' '/^name:/{print $2; exit}' argprep.yml)"
    fi
    if [[ -z "$env_name" && -f environment.yml ]]; then
        env_name="$(awk -F': *' '/^name:/{print $2; exit}' environment.yml)"
    fi
    if [[ -n "$env_name" ]]; then
        conda activate "$env_name"
    fi
fi
if ! command -v snakemake >/dev/null 2>&1; then
    echo "ERROR: snakemake not found on PATH. Activate its conda env first or set ARGPREP_CONDA_ENV." >&2
    exit 1
fi

echo "=== argprep snakemake controller ==="
echo "host:       $(hostname)"
echo "slurm job:  ${SLURM_JOB_ID:-<none>}"
echo "configfile: $configfile"
echo "extra args: $*"
echo "===================================="

exec snakemake \
    --profile profiles/slurm \
    --configfile "$configfile" \
    --rerun-incomplete \
    "$@"
