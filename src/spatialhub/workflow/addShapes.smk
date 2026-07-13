import os


configfile: "addShapes.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

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
        os.path.join(config["zarr_dir"], "{fname}", "shapes", "cellpose_segments"),
    log:
        os.path.join(config["zarr_dir"], "{fname}.addShapeCellpose.log"),
    resources:
        **RESOURCES,
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
        img_key=SHAPES_TO_ADD["cellpose"].get("img_key", "image"),
        flow_thr=SHAPES_TO_ADD["cellpose"].get("flow_thr", 0.4),
        pb_thr=SHAPES_TO_ADD["cellpose"].get("pb_thr", 0.0),
        clip_limit=SHAPES_TO_ADD["cellpose"].get("clip_limit", 0.2),
        gaussian_sigma=SHAPES_TO_ADD["cellpose"].get("gaussian_sigma", 1),
    shell:
        """
        chmod +w {params.zarr_dir} {input.zarr_in} {input.zarr_in}/shapes
        python -m spatialhub.scripts.addShapeCellpose \
            --zarr-in {input.zarr_in} \
            --ch-nuc {params.ch_nuc} \
            {params.ch_cyto_cmd} \
            --mdl-path {params.mdl_path} \
            --dp {params.dp} \
            {params.gpu_cmd} \
            --img-key {params.img_key} \
            --flow-thr {params.flow_thr} \
            --pb-thr {params.pb_thr} \
            --clip-limit {params.clip_limit} \
            --gaussian-sigma {params.gaussian_sigma} \
            >{log} 2>&1
        """
