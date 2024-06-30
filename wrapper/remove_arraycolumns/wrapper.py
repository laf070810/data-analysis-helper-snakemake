__author__ = "Anfeng Li"
__copyright__ = "Copyright 2024, Anfeng Li"
__email__ = "anfeng.li@cern.ch"
__license__ = "MIT"


import ROOT
from pytee2 import Tee

tee = Tee(output_filepath=snakemake.log[0])
tee.start()

ROOT.EnableImplicitMT(snakemake.threads)

df = ROOT.RDataFrame(snakemake.params["input_tree_name"], snakemake.input[0])

# -------- convert DecayTreeFitter variables to non-array --------
dtf_column_names = [
    column_name
    for column_name in df.GetColumnNames()
    if (snakemake.params["dtfvar_namepattern"] in str(column_name))
    and ("RVec" in df.GetColumnType(column_name))
]
for dtf_column_name in dtf_column_names:
    df = df.Redefine(dtf_column_name, f"{dtf_column_name}[0]")

# -------- save non-array variables --------
branches_to_save = [
    column_name
    for column_name in df.GetColumnNames()
    if "RVec" not in df.GetColumnType(column_name)
]
df.Snapshot(
    snakemake.params["output_tree_name"],
    snakemake.output[0],
    branches_to_save,
)

tee.stop()
