"""
Structure of config["hists"]:
    [
        {"datasets": [...], "tree_name": "...", "hists_1D": [{"model": (<name>, <xlabel>, <bins>, <low>, <high>), "expression": ""}, {}, ...], "hists_2D": [...]},
        {"datasets": [...], "tree_name": "...", "hists_1D": [...], "hists_2D": [...]},
        ...
    ]
"""
hists = config["hists"]
output_file_path_prefix = config["output_file_path_prefix"]
output_file_path_suffix = config["output_file_path_suffix"]


def get_input(hists):
    return [dataset for plot_info in hists for dataset in plot_info["datasets"]]


def get_output(hists, prefix, suffix):
    hist_paths_1D = [
        prefix + hist_1D["model"][0] + suffix
        for plot_info in hists
        for hist_1D in plot_info.get("hists_1D", [])
    ]
    hist_paths_2D = [
        prefix + hist_2D["model"][0] + suffix
        for plot_info in hists
        for hist_2D in plot_info.get("hists_2D", [])
    ]
    return {"hist_paths_1D": hist_paths_1D, "hist_paths_2D": hist_paths_2D}


rule plot:
    input:
        *get_input(hists),
    output:
        **get_output(hists, output_file_path_prefix, output_file_path_suffix),
    params:
        hists=hists,
        output_file_path_prefix=output_file_path_prefix,
        output_file_path_suffix=output_file_path_suffix,
    threads: 32
    conda:
        "../envs/plot.yml"
    script:
        "../scripts/plot.py"
