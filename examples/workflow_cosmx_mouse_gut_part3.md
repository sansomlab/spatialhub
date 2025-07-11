# CosMx mouse gut analysis workflow - Part 3

The following markdown provides a step-by-step workflow to reproduce the **cell type annotation** using the `spatialhub` pipeline and [`scANVI`](https://docs.scvi-tools.org/en/1.3.0/user_guide/models/scanvi.html). (See Parts 1-2 for the required preliminary steps.)

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

This step is highly study-specific, although it involves some common and essential data pre-processing steps, as outlined below. For this specific mouse gut study, please refer to the dedicated [GitHub repo](https://github.com/sansomlab/IBD_mouse_scRNA-seq_atlas) to reproduce the atlas creation steps in full.

> [!IMPORTANT]
> **Even if you have a pre-existing scRNA-seq dataset at hand, it is likely to require some additional post-processing before it is suitable for the purpose of transferring labels (cell type annotations) to a CosMx dataset**. Namely, you should:
> 1. (ESSENTIAL) Ensure that your reference counts matrix uses gene names that match your CosMx probe panel;
> 2. Split, regroup and/or exclude some of the cell types covered in your reference atlas, or look for further data to cover missing cell types, in order to best match the cell types that are expected in your CosMx dataset;
> 3. Reconsider your strategy for HVG selection in light of the genes in your CosMx panel and cell types of interest;
> 4. Reconsider your strategy for integration (batch covariates) in light of the batch effects found in _both_ your reference atlas and your CosMx study.
> 
> **Note that the decisions made in Steps 2 to 4 are not independent and highly iterative.**


### 1.1. ESSENTIAL: Gene names conversion

Regardless of the atlas used, the conversion of gene names (symbols) used in the input scRNA-seq dataset to match CosMx probe names is essential because: 

1. Some CosMx probes (and thus, some features in the counts matrix) are annotated with non-conventional gene symbols (which may happen in particular when using a custom panel, as in the case of this study);
2. Some CosMx probes (including some in the standard panel) cover mutiple genes (e.g. _Ccl21a_, _Ccl21b_ and _Ccl21c_ all matching by the same _Ccl21a/b/c_ CosMx probe). 

**If you don't update your atlas gene names, you may therefore miss some essential information when performing label transfer.**

First, we recommend creating a probe name to gene name conversion table. This can be achieved by first running [`_retrieve_AnnotHub_db.R`](../R/_retrieve_AnnotHub_db.R) to retrieve Ensembl annotations for the Ensembl version(s) used in your reference atlas(es). Next, you can adapat [`_map_CosMx_panels.R`](../R/_map_CosMx_panels.R) to the case of your panel (especially if it is custom), reviewing databases like [GeneCards](https://www.genecards.org/) to resolve gene synonyms at best.

With this table at hand, you can easily convert the gene names in your reference atlas, and **aggregate your counts matrix at a level that is appropriate for your CosMx panel** - see `atlas_combine_datasets.py` in this [GitHub repo](https://github.com/sansomlab/IBD_mouse_scRNA-seq_atlas) for an example script (in Python).

> [!NOTE]
> While the current version of the `spatialhub_annot` pipeline produces unexpected results and needs fixing (hence currently stored in [`archive`](../archive/)), we also point the user to the underlying [`annot_prepare_reference.R` script](../R/annot_prepare_reference.R) for suggested ways of handling this pre-processing step (in R).


### 1.2. Label definition at the 'right' resolution

Finding the correct resolution for cell type annotation can be difficult, but needs to be considered carefully. The better your reference atlas matches cell types (and experimental design) of your CosMx study, the more accurate the label transfer is likely to be.

> [!NOTE]
> Different resolutions of cell type annotation (e.g. lineage-level "T cell" vs. detailed cell type "FOXP3+ Treg cell") are likely to require different sets of HVGs, likely computed from differernt reference (sub-) atlases, possibly affected by different batch effects, etc. In this example mouse gut study, we opted for a two-stage label transfer approach, with a first round to define mid-level lineages, and a second round on each lineage to define detailed cell types (which required going through all of the steps outlined in this tutorial a second time for each lineage).

See also Bruker NanoString recommendations on the matter [here](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/cell-typing-strategies/cell-typing-strategies.html).


### 1.3. Feature selection

`scANVI` relies on training a model to annotate cell types using a selected set of features. **We've found that the way these features are selected can have a large impact on the resulting predicted cell types.** When selecting features, we thus want to find a balance between:

- Retaining enough of the highly variable genes (HVGs) from the original scRNA-seq atlas, to make sure that cell types as defined from this dataset are still properly differentiated when training the model;
- Retaining all (most) of the probes in CosMx panel (for small-scale panels, e.g. 1k) - unless focusing on very specific cell types;
- Ensuring that the number of HVGs that are not in the CosMx panel is not too large, otherwise the mismatch between the information available when training the model on the reference vs. when predicting labels on the CosMx dataset is too large for accurate predictions.

> [!TIP]
> As a rule of thumb, we've found that a total number of HVGs largely exceeding 3x the size of the original CosMx panel is likely to lead to inaccurate predictions (as flagged by the `scANVI` pipeline). In doubt, default to using CosMx probes only (or a subset of these for large CosMx panels) as HVGs for integration, and double-check that the trained model accurately re-predicts the original annotations for a subset of the atlas (see [`jupyter` notebook](../notebooks/perform_label-transfer_scANVI_example.ipynb)). In addition, remember that this step is likely to be iterative, and you may need to run the full label transfer pipeline to completion before you can decide which HVGs to use.

In this study, we determined the top 800 HVGs (_after_ aggregating features at the CosMx probe level) for each of the 5 study datasets making up the combined atlas, and used the union of these 5 sets + all CosMx probes as our set of features (which resulted in ~3000 HVGs, i.e. ~3x the 1k panel) - see `atlas_combine_datasets.py` in this [GitHub repo](https://github.com/sansomlab/IBD_mouse_scRNA-seq_atlas) for details (lines 282-385). 


### 1.4. Batch and covariate selection for integration

While it is possible to run `scANVI` without re-integrating the scRNA-seq atlas, [the developers point out that a `scANVI` model is likely to perform best if it is derived from a `scVI` model](https://docs.scvi-tools.org/en/1.3.0/user_guide/models/scanvi.html). Therefore, **you may need to determine a `batch_key` and covariates to integrate your scRNA-seq dataset across (even if it has been integrated in other ways for the purpose of your original analysis)**. 

Again, this choice will be highly specific to your study and is likely to be iterative; however, we've found that it is important to let yourself be guided by the structure of both your scRNA-seq atlas and your CosMx dataset. For example, in this mouse gut study, we used `dataset` as a `batch_key`, where `dataset` corresponds to one of 5 studies included in the reference atlas, and reflects whether a sample was cut from 'old' vs. 'new' histology blocks (regardless of the microscopy slide it then ended up on) in our CosMx study.

* * *


## Step 2: Train model and transfer labels

Once you've carefully considered the steps above (which is likely to involve running the whole cell type annotation workflow through to the end several times over), adapt the Jupyter notebooks provided in the `notebooks` folder on a GPU to perform label transfer and a preliminary assessment of validity.

* * *


## Step 2: Review cell annotations

Work in progress
