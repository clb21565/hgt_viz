# HGT Viz

Interactive Streamlit dashboard for visualizing horizontal gene transfer (HGT) analysis results from metagenome assemblies.

## Pages

| Page | Description |
|------|-------------|
| Overview | Key metrics, detection level counts, region type distribution, group scatter |
| Transfer Flows | Sankey diagram and heatmap of donor→recipient flows at any taxonomic level |
| Groups | Treemap and scatter explorer of HGT groups |
| Mismatch Heatmap | Region taxon vs. contig taxon matrix |
| Contig View | Linear gene diagrams per contig or HGT group; HGT2 synteny ribbons |

## Expected data format

The app expects output files from an HGT pipeline with a shared assembly prefix:

```
<data_dir>/
  <assembly>.allHGT1.tsv
  <assembly>.allHGT2.tsv
  <assembly>.allDR.tsv
  <assembly>.hgt_groups.tsv
  <assembly>.hgt_report.tsv
```

The assembly prefix is auto-detected from any `*.allHGT1.tsv` file in the data directory.

## Setup

```bash
git clone https://github.com/clb21565/hgt_viz
cd hgt_viz

conda create -n hgt_viz python=3.11 -y
conda activate hgt_viz
pip install -r requirements.txt
```

## Running on HPC

### Login node (quick interactive use)

```bash
conda activate hgt_viz
HGT_DATA_DIR=/path/to/your/data streamlit run app.py \
    --server.headless true \
    --server.port 8501
```

Then on your local machine, open an SSH tunnel:

```bash
ssh -L 8501:localhost:8501 youruser@hpc.host
```

Open `http://localhost:8501` in your browser.

### Compute node

If running inside a job on a compute node, you need a two-hop tunnel. In your job script, print the node name:

```bash
echo "Compute node: $(hostname)"
HGT_DATA_DIR=/path/to/your/data streamlit run app.py \
    --server.headless true \
    --server.port 8501
```

Then on your local machine:

```bash
ssh -L 8501:<compute-node-name>:8501 youruser@hpc.host
```

### SLURM job script example

```bash
#!/bin/bash
#SBATCH --job-name=hgt_viz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=hgt_viz_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate hgt_viz

echo "Compute node: $(hostname)"
echo "Port: 8501"
echo "Tunnel: ssh -L 8501:$(hostname):8501 ${USER}@hpc.host"

HGT_DATA_DIR=/path/to/your/data \
    streamlit run /path/to/hgt_viz/app.py \
    --server.headless true \
    --server.port 8501
```

## Environment variables

| Variable | Description | Default |
|----------|-------------|---------|
| `HGT_DATA_DIR` | Path to directory containing HGT output files | Current working directory |
| `HGT_ASSEMBLY` | Assembly filename prefix | Auto-detected from `*.allHGT1.tsv` |

If neither variable is set, the sidebar controls can be used to set the data directory and assembly prefix interactively after launch.

## Switching datasets

The **Dataset** panel in the sidebar lets you change the data directory and assembly prefix without restarting the app. Each unique combination is cached independently.
