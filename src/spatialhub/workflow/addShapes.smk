import os


configfile: "addShapes.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

# Base resources, in native-Slurm dialect (slurm_partition / slurm_account).
# gres_cmd defaults to "" so the CLI's {resources.gres_cmd} template substitutes
# cleanly on CPU rules (DRMAA requires the key present on every rule).
RESOURCES = {
    "threads": 4,
    "mem_mb": 16000,
    "time": "04:00:00",
    "slurm_partition": "short",
    "gres_cmd": "",
}
RESOURCES.update(config.get("resources", {}))
RESOURCES.setdefault("gres_cmd", "")

# Transitional shim: the shared DRMAA submission string still references the
# legacy {resources.partition}/{resources.account} names. Populate those from a
# single native-dialect source so one config value drives both, with no
# duplicate key in the YAML to drift. After migration to the native executor,
# it reads slurm_partition / slurm_account directly and these two lines go away.
RESOURCES["partition"] = RESOURCES["slurm_partition"]
RESOURCES["account"] = RESOURCES.get("slurm_account", "")

# GPU resources: inherit the base, override partition (both names), add GRES.
_gpu_cfg = config.get("gpu_resources", {})
GPU_RESOURCES = dict(RESOURCES)
_gpu_part = _gpu_cfg.get("slurm_partition", RESOURCES["slurm_partition"])
GPU_RESOURCES["slurm_partition"] = _gpu_part
GPU_RESOURCES["partition"] = _gpu_part
_n_gpus = _gpu_cfg.get("n_gpus", 1)
GPU_RESOURCES["gres_cmd"] = _gpu_cfg.get("gres_cmd", f"--gres=gpu:{_n_gpus}")

FILE_LST = [f for f in os.listdir(config["zarr_dir"]) if f.endswith(".zarr")]

KEY_DICT = {
    "cellpose": "cellpose_segments",
    "FOV": "fov_boxes",
}
SHAPES_TO_ADD = config["shapes_to_add"]


rule full:
    input:
        expand(
            os.path.join(config["zarr_dir"], "{fname}", "shapes", "{shape}"),
            fname=FILE_LST,
            shape=[KEY_DICT[shape] for shape in SHAPES_TO_ADD.keys()],
        ),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        flst=[os.path.join(config["zarr_dir"], fname) for fname in FILE_LST],
    shell:
        """
        if [ {params.lock} = true ]; then
            for fname in {params.flst}; do
                chmod -R a-w $fname
            done
        fi
        echo "=======> All done! <======="
        """


rule add_cellpose:
    input:
        zarr_in=os.path.join(config["zarr_dir"], "{fname}"),
    output:
        directory(
            os.path.join(config["zarr_dir"], "{fname}", "shapes", "cellpose_segments")
        ),
    log:
        os.path.join(config["zarr_dir"], "{fname}.addShapeCellpose.log"),
    resources:
        **GPU_RESOURCES,
    params:
        zarr_dir=config["zarr_dir"],
        ch_nuc=SHAPES_TO_ADD["cellpose"]["ch_nuc"],
        ch_cyto_cmd=(
            f"--ch-cyto {SHAPES_TO_ADD['cellpose']['ch_cyto']}"
            if SHAPES_TO_ADD["cellpose"].get("ch_cyto")
            else ""
        ),
        mdl_path=SHAPES_TO_ADD["cellpose"]["mdl_path"],
        dp=SHAPES_TO_ADD["cellpose"]["dp"],
        gpu_cmd="--gpu" if SHAPES_TO_ADD["cellpose"].get("gpu") else "",
        patch_width=SHAPES_TO_ADD["cellpose"].get("patch_width", None),
        patch_overlap=SHAPES_TO_ADD["cellpose"].get("patch_overlap", None),
        img_key=SHAPES_TO_ADD["cellpose"].get("img_key", "image"),
        flow_thr=SHAPES_TO_ADD["cellpose"].get("flow_threshold", 2),
        pb_thr=SHAPES_TO_ADD["cellpose"].get("cellprob_threshold", -6),
        clip_limit=SHAPES_TO_ADD["cellpose"].get("clip_limit", 0.2),
        gaussian_sigma=SHAPES_TO_ADD["cellpose"].get("gaussian_sigma", 1),
    shell:
        """
        if [ {params.patch_width} = None ]; then
            pwidth_cmd=""
        else
            pwidth_cmd="--pwidth {params.patch_width}"
        fi
        if [ {params.patch_overlap} = None ]; then
            povlp_cmd=""
        else
            povlp_cmd="--povlp {params.patch_overlap}"
        fi

        chmod +w {params.zarr_dir} {input.zarr_in} {input.zarr_in}/shapes
        python -m spatialhub.scripts.addShapeCellpose \
            --zarr-in {input.zarr_in} \
            --ch-nuc {params.ch_nuc} \
            {params.ch_cyto_cmd} \
            --mdl-path {params.mdl_path} \
            --dp {params.dp} \
            {params.gpu_cmd} \
            $pwidth_cmd \
            $povlp_cmd \
            --img-key {params.img_key} \
            --flow-thr {params.flow_thr} \
            --pb-thr {params.pb_thr} \
            --clip-limit {params.clip_limit} \
            --gaussian-sigma {params.gaussian_sigma} \
            >{log} 2>&1
        """
