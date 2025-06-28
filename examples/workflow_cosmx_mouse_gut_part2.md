# CosMx mouse gut analysis workflow - Part 2

The following markdown provides a step-by-step workflow to reproduce the **re-segmentation of a CosMx study of the mouse gut in an IBD model** using the `spatialhub` pipeline. (See Part 3 for the next steps of the analysis, and Part 1 for the required preliminary steps.)

> [!IMPORTANT]
> At this stage (June 2025), this part is very much work in progress, as we haven't been able to generate a segmentation mask which is convincingly better than the default `AtoMx` one. This document thus provides an overview of what has been tested.

> [!TIP]
> Regardless of the segmenation tool used, we recommend first running a pilot `spatialhub` run on a subset of samples (this can be achieved by pointing to a subset `samples.tsv` table in the YAML file, instead of the `samples.tsv` table covering the full study dataset). **This piloting is essential since most segmentation tools rely on a lot of parameters, which will need to be optimized for each dataset.**  


* * *

## Study overview

The experiment covers . . .

### Data storage on BMRC

Raw and processsed data alongside key scripts to reproduce the analysis on the BMRC cluster in this directory: `/users/sansom/tme871/work/cosmx_mouse`

Here is an overview of where to find key files:

- Original raw data, as exported from `AtoMx SIP` through their user interface: `/users/sansom/tme871/work/cosmx_mouse/data/raw/atomx`

> [!NOTE]
> At this stage of the analysis (June 2025), we are confidently moving away from `AtoMx` default analyses tools, so these data files could be deleted (knowing they can always be downloaded again from our `AtoMx SIP`). I have kept them as an example of an `AtoMx` export, but the minimal relevant files are copied in the directories described in the bullet point below.

- **Raw data, re-organized for analysis with `spatialhub`**: `/users/sansom/tme871/work/cosmx_mouse/data/raw/${SLIDE_ID}`, with the following slide IDs available: NL4S3a, NL4S3b, NL5H1a, NL5H1b ('a' indicates a KIR run, 'b' indicates a CAMS run - see experimental design section above)
- **Metadata**, including `spatialhub` samples.tsv table (full and filtered for samples passing QC, see details below) and lists of CosMx probes (standard 'LBL-11176-05-Mouse-Universal-Cell-Characterization-Gene-List.txt' and custom add-on panel '1K_Oxford_Sansom_add-on_list.csv')

> [!WARNING]
> **In this study, some of the custom probes use special characters (e.g. *SiglecF (170)* instead of *Siglecf*) and/or non-conventional gene symbols (e.g. *Ly6G* instead of *Ly6g*; *Cd64* instead of *Fcgr1*). In addition, some probes in the standard panel are cross-reactive and match multiple genes. Whenever comparing a CosMx probe panel to reference scRNA-seq datasets/pathway libraries/etc., it is essential to ensure that gene names are properly encoded for the intended purpose.** This point is highlighted again where most relevant in the workflow description.

- **IBD_mouse_scRNA-seq_atlas**: Directory (and [GitHub repo](https://github.com/sansomlab/IBD_mouse_scRNA-seq_atlas)) collecting data and scripts used to generate a scRNA-seq atlas fit for the purpose of annotating this specific CosMx dataset.
- **spatialhub**: Main analysis folder, containing YAML, log out output files/directories from steps run with the `spatialhub` pipeline, as well as additional scripts for steps requiring more customization.
- Other directories contain pilot scripts (irrelevant to reproduce the final analysis), used to test different spatial transcriptomics tools before rolling them out to the dataset as a whole and/or adding them to the `spatialhub` suite of pipelines, if applicable

* * *


## Option 1: Using image-based tools

For re-segmentation using image-based tools, **`spatialhub` currently builds on top of the [`SOPA`](https://gustaveroussy.github.io/sopa/) implementation of [`cellpose`](https://cellpose.readthedocs.io/en/latest/). We thus recommend reading the documentation for both tools for additional details.**


> [!NOTE]
> `SOPA` is well suited to enable fast re-segmentation, since it patches the image and parallelize the segmentation to run on each patch before (properly) stitching the segmentation mask pieces into one (unlike `AtoMx` tools, which just append pieces of the segmentation mask next to each other without resolving conflicts at the boundaries). At this stage, we have only implemented the `cellpose` tool from `SOPA` into `spatialhub`, but the `sopa_segment_img` pipeline could be expanded to cover more tools (please refer to `SOPA` Read The Docs for available image-based segmentation tools).

After running `spatialhub ashlar`, you can readily run the `sopa_segment_img` within the same directory, since the input files will be fetched from the `zarr.dir` directory created by Ashlar, using the following commands:

```
spatialhub sopa_segment_img config
spatialhub sopa_segment_img make full -v5 -p20
```

The YAML file covers major parameters from the [`sopa.segmenation.cellpose` function](https://gustaveroussy.github.io/sopa/api/segmentation/#sopa.segmentation.cellpose). How those parameters impact the resulting segmentation mask are further described in the [`cellpose` Read The Docs](https://cellpose.readthedocs.io/en/latest/settings.html). Beyond running `cellpose`, `sopa.segmenation.cellpose` includes some helpful segmentation mask clean-up functionalities, such as the ability to remove cells under a certain area (again, see `SOPA` for details).

The output segmentation mask will be automatically added to the `zarr.dir` directory, using the `mask_key` specified in the YAML file.

> [!TIP]
> The best results with `cellpose` are likely to be achieved when running it based on a custom, re-trained model (see [`cellpose` tutorial here](https://cellpose.readthedocs.io/en/latest/gui.html#training-your-own-cellpose-model)). If using a pre-trained model, we've had most success with the `nuclei`, `cyto3` and `tissuenet_cp3` models (see [here for a full list of available `cellpose` pre-trained modes](https://cellpose.readthedocs.io/en/v3.1.1.1/models.html#full-built-in-models).

* * *


## Option 2: Using transcripts-based tools

```
spatialhub sopa_segment_tx config
```

```
spatialhub sopa_segment_tx make full -v5 -p20
```
