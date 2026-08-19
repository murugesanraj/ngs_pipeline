# Shiny results portal

The portal is deliberately read-only. It displays only files explicitly listed
in each pipeline-generated `portal/manifest.json`; it does not accept arbitrary
paths or launch analysis jobs.

## Install

From an interactive Hellbender allocation:

```bash
bin/bootstrap_hellbender.sh --with-shiny
```

## Launch through Open OnDemand

1. Sign in to `https://ondemand.rnet.missouri.edu`.
2. Start an RStudio Server interactive session.
3. In its terminal, activate the Shiny environment.
4. Set the approved portal root and run the app.

```bash
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$HOME/data/conda/envs/hellbender-ngs-shiny"
export NGS_PORTAL_ROOT=/cluster/pixstor/GROUP/ngs_portal
Rscript /path/to/hellbender-ngs-pipeline/app/run_app.R
```

Each client runs the app as their own Hellbender account. A project appears only
when its manifest and listed files are readable through that account's existing
POSIX group/ACL permissions.

The default host is `127.0.0.1`. Do not set `NGS_PORTAL_HOST=0.0.0.0` to expose
a compute-node port. A central service needs an approved reverse proxy, SSO,
TLS, logging, and an IT security review.

