<div class="admonition spinner" style="text-align: center;">
  <p class="admonition-title">
    <span style="display: inline-block; animation: pulse 2s ease-in-out infinite;">🚧</span>
    Work in Progress
  </p>
  <p>SpatialHub is currently in pre-release and under active development. Although we have taken great care 
     to ensure its quality and reliability, bugs and breaking changes may still occur. Please use it with caution 
    and report any issues by opening a  <a href="https://github.com/sansomlab/spatialhub/issues">GitHub issue</a>.</p>
</div>

<style>
@keyframes pulse {
  0%, 100% { opacity: 1; }
  50% { opacity: 0.4; }
}
</style>

# SpatialHub

`spatialhub` is a collection of modular pipelines for the processing of spatial transcriptomics data, compatible with the [scverse](https://scverse.org/) ecosystem.

Currently supports analysis of CosMx and VisiumHD data. Support for Xenium coming soon.

<div class="sh-grid" markdown>

<div class="sh-card" markdown>

:lucide-book-open:{ .sh-icon }

**Overview**
{ .sh-title }

The supported platforms, the data formats used throughout, and how the pipelines connect from raw output to analysis.

[Read the overview :lucide-arrow-right:](overview.md){ .sh-link }

</div>

<div class="sh-card" markdown>

:lucide-download:{ .sh-icon }

**Installation**
{ .sh-title }

Installing SpatialHub on the BMRC cluster, configuring the environment and verifying the installation.

[Install SpatialHub :lucide-arrow-right:](installation.md){ .sh-link }

</div>

<details class="sh-card sh-pipelines" markdown>
<summary markdown="span">:lucide-workflow:{ .sh-icon } <span class="sh-summary-text"><span class="sh-title">Pipelines</span><span class="sh-meta">Seven pipelines covering data preparation, conversion and cell-type deconvolution</span></span> <span class="sh-toggle"><span class="sh-toggle-label"></span>:lucide-chevron-down:</span></summary>

<div class="sh-stages" markdown>

<div class="sh-stage" markdown>

**Build SpatialData objects**
{ .sh-stage-title }

- [`cosmx_makeZarr`](pipelines.md#cosmx_makezarr)
  One Zarr store per CosMx sample, with or without Ashlar stitching
- [`visiumhd_makeZarr`](pipelines.md#visiumhd_makezarr)
  One Zarr store per Visium HD capture area, via Space Ranger count
- [`visiumhd_prepareBin2Cell`](pipelines.md#visiumhd_preparebin2cell)
  Bin2Cell intermediate files up to the StarDist step

</div>

<div class="sh-stage" markdown>

**Create AnnData objects**
{ .sh-stage-title }

- [`extractH5AD`](pipelines.md#extracth5ad-aggregateh5ad)
  Extract an existing table to H5AD, typically for Visium HD
- [`aggregateH5AD`](pipelines.md#extracth5ad-aggregateh5ad)
  Aggregate points and shapes into H5AD, typically for CosMx

</div>

<div class="sh-stage" markdown>

**Cell-type deconvolution**
{ .sh-stage-title }

- [`runRCTD`](pipelines.md#runrctd-runcell2location)
  Deconvolution with RCTD against a single-cell reference
- [`runCell2Location`](pipelines.md#runrctd-runcell2location)
  Deconvolution with Cell2Location against a single-cell reference

</div>

</div>

</details>

</div>