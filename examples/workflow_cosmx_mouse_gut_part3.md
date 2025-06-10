# CosMx mouse gut analysis workflow - Part 3

The following markdown provides a step-by-step workflow to reproduce the **cell type annotation** using the `spatialhub` pipeline and [`scANVI`](https://docs.scvi-tools.org/en/1.3.0/user_guide/models/scanvi.html). (See Parts 1-2 for the required preliminary steps.)

> [!CAUTION]
> At this stage (June 2025), the `spatialhub` pipeline includes an `annot` pipeline, meant to streamline the steps covered in the [matching `jupyter` notebook](../notebooks/TO_BE_ADDED). **While this pipeline currently can run to completion without throwing an error, it produces incorrect outputs and needs to be fixed before it can be used.**


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

## Step 1: Prepare a reference scRNA-seq atlas

This step is highly study-specific, although it involves some common and essential data pre-processing steps outlined below. For this mouse gut study, please refer to the dedicated [GitHub repo](https://github.com/sansomlab/IBD_mouse_scRNA-seq_atlas) to reproduce the atlas creation steps in full.


### 1.1. Conversion of gene names in the scRNA-seq atlas to CosMx probe names

Regardless of the atlas used, the conversion of gene names (symbols) used in the input scRNA-seq dataset to match CosMx probe names is essential because: 

1. Some CosMx probes (and thus, some features in the counts matrix) are annotated with non-conventional gene symbols (which may happen in particular when using a custom panel, as in the case of this study);
2. Some CosMx probes (including some in the standard panel) cover mutiple genes (e.g. _Ccl21a_, _Ccl21b_ and _Ccl21c_ all matching by the same _Ccl21a/b/c_ CosMx probe). 

If you don't update your atlas gene names, you may therefore miss some essential information when performing label transfer.

> [!NOTE]
> While the current version of the `spatialhub_annot` pipeline needs fixing, we point the user to the [`annot_prepare_reference.R` script](../R/annot_prepare_reference.R) for suggested ways of handling these pre-processing steps, or to `atlas_combine_datasets.py` in this [GitHub repo](https://github.com/sansomlab/IBD_mouse_scRNA-seq_atlas) for an alternative approach (in Python).


### 1.2. Feature selection

`scANVI` relies on training a model to annotate cell types using a selected set of features. **We've found that the way these features are selected can have a large impact on the resulting predicted cell types.** When selecting features, we thus want to find a balance between:

- Retaining enough of the highly variable genes (HVGs) from the original scRNA-seq atlas, to make sure that cell types as defined from this dataset are still properly differentiated when training the model;
- Retaining all (most) of the probes in CosMx panel (for small-scale panels, e.g. 1k);
- Ensuring that the number of HVGs that are not in the CosMx panel is not too large, otherwise predictions will be derived from too little information.

> [!TIP]
> As a rule of thumb, we've found that a total number of HVGs largely exceeding 3x the size of the original CosMx panel is likely to lead to inaccurate predictions (as flagged by the `scANVI` pipeline). In doubt, default to using CosMx probes only as HVGs for integration, and double-check that the trained model accurately re-predicts the original annotations for a subset of the atlas (see [`jupyter` notebook](../notebooks/TO_BE_ADDED)]). In addition, note that this step is likely to be iterative, and you may need to run the full label transfer pipeline to completion before you can decide which HVGs to use.

In this study, we determined the top 800 HVGs (_after_ aggregating features at the CosMx probe level) for each of the 5 study datasets making up the combined atlas, and used the union of these 5 sets + all CosMx probes as our set of features (which resulted in ~3000 HVGs, i.e. ~3x the 1k panel). 


### 1.3. Batch and covariate selection

While it is possible to run `scANVI` without re-integrating the scRNA-seq atlas, [the developers point out that a `scANVI` model is likely to perform best if it is derived from a `scVI` model](https://docs.scvi-tools.org/en/1.3.0/user_guide/models/scanvi.html). Therefore, **you may need to determine a `batch_key` and covariates to integrate your scRNA-seq dataset across (even if it has been integrated in other ways for the purpose of your original analysis)**. 

Again, this choice will be highly specific to your study and is likely to be iterative; however, we've found that it is important to let yourself be guided by the structure of both your scRNA-seq atlas and your CosMx dataset.

In this study, we used `dataset` as a `batch_key`, where `dataset` corresponds to one of 5 studies included in the reference atlas, and corresponds to old vs. new sample blocks in our CosMx study. In addition, we used . . . as covariates. 

* * *


## Step 2: Integrate 

TBC

* * *


## Step 2: Review outcomes

TBC