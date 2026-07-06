import os
import numpy as np
import pandas as pd
import scanpy as sc
import cell2location as c2l

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, GREEN, RESET, print_arguments


def main():
    p = AP(description="Deconvolute cell abundances with cell2location.")
    p.add_argument("outdir", help="Path to output directory.")
    p.add_argument("--ref", required=True, help="H5AD for the single-cell reference.")
    p.add_argument("--qry", required=True, help="H5AD for the spatial data in query.")
    p.add_argument("--label-key", required=True, help="Labels key (e.g., cell_type).")
    p.add_argument("--batch-key", default=None, help="Batch key (e.g., sample_id).")
    p.add_argument("--cat-keys", default=None, help="Comma-separated categorical keys.")
    p.add_argument("--ctn-keys", default=None, help="Comma-separated continuous keys.")
    p.add_argument(
        "--max-epochs-reg",
        type=int,
        default=400,
        help="Max epoch number for regression.",
    )
    p.add_argument(
        "--qt-lst",
        default="q05,q50,q95,q0001",
        help="Comma-separated quantiles to export.",
    )
    p.add_argument(
        "--n-cellperloc", type=int, default=50, help="Number of cells per location."
    )
    p.add_argument(
        "--detection-alpha",
        type=float,
        default=20.0,
        help="Detection alpha for cell2location.",
    )
    p.add_argument(
        "--max-epochs-c2l",
        type=int,
        default=30000,
        help="Max number of epochs for cell2location.",
    )
    args = p.parse_args()
    print_arguments(args)

    qt_lst = [q.strip() for q in args.qt_lst.split(",") if q.strip()]

    out_dict = {
        "mdl.reg": os.path.join(args.outdir, "mdl.reg.pt"),
        "inf_aver": os.path.join(args.outdir, "inf_aver.csv"),
        "mdl.c2l": os.path.join(args.outdir, "mdl.c2l.pt"),
    } | {
        f"cellabd_{qt}": os.path.join(args.outdir, f"{qt}_cell_abundance_w_sf.csv")
        for qt in qt_lst
    }
    for key, path in out_dict.items():
        if os.path.exists(path):
            raise FileExistsError(f"Output file {path} already exists for {key}.")
    os.makedirs(args.outdir, exist_ok=True)

    cat_cov_lst = args.cat_keys.split(",") if args.cat_keys else None
    ctn_cov_lst = args.ctn_keys.split(",") if args.ctn_keys else None

    # Train regression model on reference data, export posterior, extract inferred average expression per cell type
    rdata = sc.read_h5ad(args.ref)
    print(rdata, "\n", rdata.X[rdata.X > 5], "\n")
    if "counts" not in rdata.layers:
        raise ValueError("Reference data must contain a 'counts' layer.")

    required = {args.label_key}
    if args.batch_key:
        required |= {args.batch_key}
    if cat_cov_lst:
        required |= set(cat_cov_lst)
    if ctn_cov_lst:
        required |= set(ctn_cov_lst)
    missing = required - set(rdata.obs.columns)
    if missing:
        raise ValueError(f"Missing columns in reference data: {missing}")

    c2l.models.RegressionModel.setup_anndata(
        adata=rdata,
        layer="counts",
        batch_key=args.batch_key,
        labels_key=args.label_key,
        categorical_covariate_keys=cat_cov_lst,
        continuous_covariate_keys=ctn_cov_lst,
    )
    mod = c2l.models.RegressionModel(rdata)
    mod.train(max_epochs=args.max_epochs_reg, batch_size=2500, train_size=1)
    rdata = mod.export_posterior(
        rdata, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
    )
    rdata = mod.export_posterior(
        rdata,
        use_quantiles=True,
        add_to_varm=qt_lst,
        sample_kwargs={"batch_size": 2500},
    )
    mod.save(out_dict["mdl.reg"], overwrite=False)

    inf_aver = rdata.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in rdata.uns["mod"]["factor_names"]]
    ].copy()
    inf_aver.columns = rdata.uns["mod"]["factor_names"]
    inf_aver.to_csv(out_dict["inf_aver"])

    # Run cell2location on query data
    qdata = sc.read_h5ad(args.qry)
    print(qdata, "\n", qdata.X[qdata.X > 5], "\n")
    if "counts" not in qdata.layers:
        raise ValueError("Query data must contain a 'counts' layer.")

    intersect = np.intersect1d(qdata.var_names, inf_aver.index)
    print(f"{len(intersect)} shared genes ({len(intersect)/len(qdata.var_names):.1%})")

    qdata = qdata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    c2l.models.Cell2location.setup_anndata(
        adata=qdata,
        layer="counts",
        batch_key=args.batch_key,
        categorical_covariate_keys=cat_cov_lst,
        continuous_covariate_keys=ctn_cov_lst,
    )
    mod = c2l.models.Cell2location(
        qdata,
        cell_state_df=inf_aver,
        N_cells_per_location=args.n_cellperloc,
        detection_alpha=args.detection_alpha,
    )
    mod.train(max_epochs=args.max_epochs_c2l, batch_size=None, train_size=1)
    qdata = mod.export_posterior(
        qdata,
        sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs},
    )
    mod.save(out_dict["mdl.c2l"], overwrite=False)

    for qt in qt_lst:
        pd.DataFrame(
            qdata.obsm[f"{qt}_cell_abundance_w_sf"],
            index=qdata.obs_names,
            columns=qdata.uns["mod"]["factor_names"],
        ).to_csv(out_dict[f"cellabd_{qt}"])

    print(f"{GREEN}Cell2location completed successfully.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to run Cell2Location.")
