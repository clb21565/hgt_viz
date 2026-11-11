#!/bin/bash
#SBATCH --job-name=hgt_viz
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=hgt_viz_%j.log
#
# Build an hgt_viz store from a kairos-lr results/ directory (if it does not exist
# yet) and serve the app on this node.
#
#   sbatch --export=ALL,REPO=/path/to/hgt_viz,RESULTS=/path/to/kairos-lr/results \
#          [,STORE=...][,CARD_INDEX=/path/to/aro_index.tsv][,PORT=8501][,REBUILD=1] \
#          [,ENV=name-or-path][,LOGIN_HOST=...] \
#          slurm/hgt_viz.sbatch
#
# Then follow the "Tunnel:" line in hgt_viz_<jobid>.log and open http://localhost:PORT.
set -euo pipefail

: "${REPO:?set REPO=/path/to/hgt_viz}"
: "${RESULTS:?set RESULTS=/path/to/kairos-lr/results}"
STORE=${STORE:-$(dirname "$RESULTS")/hgt_viz_store}
PORT=${PORT:-8501}
ENV=${ENV:-hgt_viz}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"
cd "$REPO"

if [[ ! -f "$STORE/meta.json" || -n "${REBUILD:-}" ]]; then
    echo "[$(date +%T)] building store $STORE from $RESULTS"
    /usr/bin/time -v python -m hgtviz.build --results "$RESULTS" --out "$STORE" \
        ${CARD_INDEX:+--card-index "$CARD_INDEX"} 2> >(tee build_time_${SLURM_JOB_ID}.txt >&2)
    du -sh "$STORE"
fi

echo "Tunnel: ssh -N -L ${PORT}:$(hostname):${PORT} ${USER}@${LOGIN_HOST:-<login-host>}"
HGT_VIZ_STORE="$STORE" exec streamlit run app.py \
    --server.headless true --server.port "$PORT" --browser.gatherUsageStats false
