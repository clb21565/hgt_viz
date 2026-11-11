#!/bin/bash
# One-command launcher for colleagues: reads site settings from slurm/site.env
# (next to this script; not in git), submits slurm/hgt_viz.sbatch and prints how to
# connect once the job runs.
#
#   bash /abs/path/to/hgt_viz/slurm/launch.sh            # start the app
#   PORT=8642 bash /abs/path/to/hgt_viz/slurm/launch.sh  # pick another port if 8501 is taken
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(dirname "$HERE")
SITE=${SITE:-$HERE/site.env}
[[ -f $SITE ]] || { echo "Missing $SITE: copy slurm/site.env.example and fill in the paths." >&2; exit 1; }
PORT_OVERRIDE=${PORT:-}
# shellcheck disable=SC1090
source "$SITE"
PORT=${PORT_OVERRIDE:-${PORT:-8501}}
: "${STORE:?STORE not set in $SITE}"
: "${RESULTS:?RESULTS not set in $SITE}"
: "${ENV:?ENV not set in $SITE}"
[[ -f $STORE/meta.json ]] || { echo "No hgt_viz store at $STORE (no meta.json)." >&2; exit 1; }

LOGDIR=${LOGDIR:-$HOME/hgt_viz_logs}
mkdir -p "$LOGDIR"
JOB=$(sbatch --parsable --chdir="$LOGDIR" ${PARTITION:+--partition="$PARTITION"} \
      --export=ALL,REPO="$REPO",RESULTS="$RESULTS",STORE="$STORE",ENV="$ENV",PORT="$PORT",LOGIN_HOST="${LOGIN_HOST:-}" \
      "$HERE/hgt_viz.sbatch")
echo "Submitted job $JOB. Waiting for it to start…"
for _ in $(seq 1 120); do
    NODE=$(squeue -h -j "$JOB" -o "%N" 2>/dev/null || true)
    STATE=$(squeue -h -j "$JOB" -o "%T" 2>/dev/null || true)
    [[ $STATE == RUNNING && -n $NODE ]] && break
    sleep 5
done
if [[ ${STATE:-} != RUNNING ]]; then
    echo "Job $JOB is not running yet (state: ${STATE:-unknown}). Check later with: squeue -j $JOB"
    echo "Log: $LOGDIR/hgt_viz_$JOB.log"
    exit 0
fi
cat <<MSG

The app is starting on $NODE (allow ~1 min to load the store).

On your own computer, run:
    ssh -N -L $PORT:$NODE:$PORT $USER@${LOGIN_HOST:-<login-host>}
then open http://localhost:$PORT in a browser.

Stop the app when done:   scancel $JOB
Log:                      $LOGDIR/hgt_viz_$JOB.log
MSG
