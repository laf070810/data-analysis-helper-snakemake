__author__ = "Anfeng Li"
__copyright__ = "Copyright 2024, Anfeng Li"
__email__ = "anfeng.li@cern.ch"
__license__ = "MIT"


import math
import pickle
from itertools import product

import pandas as pd
import ROOT
from pytee2 import Tee

tee = Tee(output_filepath=snakemake.log[0])
tee.start()

ROOT.EnableImplicitMT(snakemake.threads)

rdf_signal = ROOT.RDataFrame(
    snakemake.params["signal_tree_name"], snakemake.input["signal_filepaths"]
)
rdf_background = ROOT.RDataFrame(
    snakemake.params["background_tree_name"], snakemake.input["background_filepaths"]
)

n_total_signal = rdf_signal.Count().GetValue()
n_total_background = rdf_background.Count().GetValue()

fom_result = {"signal_eff": {}, "background_eff": {}, "S": {}, "B": {}, "FoM": {}}
for bin_indices in product(*snakemake.params["cut_points"].values()):
    print(f"calculating {bin_indices}...")

    criteria = []
    for i, variable_name in enumerate(snakemake.params["cut_points"].keys()):
        criteria.append(f"({variable_name} > {bin_indices[i]})")

    fom_result["signal_eff"][bin_indices] = (
        rdf_signal.Filter(" && ".join(criteria)).Count().GetValue() / n_total_signal
    )
    fom_result["background_eff"][bin_indices] = (
        rdf_background.Filter(" && ".join(criteria)).Count().GetValue()
        / n_total_background
    )
    fom_result["S"][bin_indices] = (
        snakemake.params["nsignal"] * fom_result["signal_eff"][bin_indices]
    )
    fom_result["B"][bin_indices] = (
        snakemake.params["nbackground"] * fom_result["background_eff"][bin_indices]
    )
    fom_result["FoM"][bin_indices] = fom_result["S"][bin_indices] / math.sqrt(
        fom_result["S"][bin_indices] + fom_result["B"][bin_indices]
    )
df = pd.DataFrame(fom_result)

with open(snakemake.output["fom_filepath"], "wb") as f:
    pickle.dump(
        {
            "fom_result": df,
            "cut_points_best": {
                variable_name: df["FoM"].idxmax()[i]
                for i, variable_name in enumerate(snakemake.params["cut_points"].keys())
            },
        },
        f,
    )

tee.stop()
