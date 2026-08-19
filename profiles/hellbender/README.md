# Hellbender profiles

`config.v9+.yaml` targets the public `general` partition and deliberately uses
`slurm-no-account`. It assumes Snakemake 9 and the official SLURM executor
plugin.

If your lab has priority resources, copy `priority.config.example.yaml` to
`../hellbender-priority/config.v9+.yaml`, replace both placeholders, and keep
that directory untracked. Confirm the live values with your PI/HPC allocation
manager or:

```bash
sacctmgr show assoc where user="$USER" format=Account,Partition,QOS
scontrol show partition PARTITION_NAME
```

Resource overrides can be passed without editing rules, for example:

```bash
snakemake --profile profiles/hellbender \
  --set-resources star_align:mem_mb=96000 star_align:runtime=960
```

