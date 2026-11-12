# Segmentation workflows

`spatialhub` supports segmentation through different workflows.


## SOPA image segmentation

### Setup

The setup task creates a minimal copy of the original zarr stack, and ensure all elements in this object match expected SOPA conventions. Note that this if you modify your zarr stack upstream, Snakemake will not detect it and not re-run this step by default. To force re-execution of this step, delete the sentinel files it creates.

### StarDist

Before running the pipeline, make sure to download the model of choice and noting its path by running the following script in an interactive session (replacing "2D_versatile_fluo" with your StarDist model of choice):

```
from stardist.models import StarDist2D
model = StarDist2D.from_pretrained("2D_versatile_fluo")
print(model.logdir)
```

### Registration

To force image registration to re-run (e.g. after switching maks name overwrote from false to true), delete the corresponding sentinel in `sopa.segment.dir/patches/P6/stardistResolve_stardist_fluo_<HASH>.sentinel`

If you wish to register an old mask which wasn't the latest created, you must force the re-run of the SOPA patch-resolution first. This can be achieved by removing 2 sentinel files:

```
rm sopa.image.dir/patches/P6/stardistResolve_stardist_fluo_<HASH>.sentinel
rm sopa.image.dir/logs/P6_sopaRegister_stardist_<HASH>.sentinel
```