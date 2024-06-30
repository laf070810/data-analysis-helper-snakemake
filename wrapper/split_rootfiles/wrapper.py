__author__ = "Anfeng Li"
__copyright__ = "Copyright 2024, Anfeng Li"
__email__ = "anfeng.li@cern.ch"
__license__ = "MIT"


import ROOT
from pytee2 import Tee

tee = Tee(output_filepath=snakemake.log[0])
tee.start()

ipartition = snakemake.params["ipartition"]
npartitions = snakemake.params["npartitions"]

rdf = ROOT.RDataFrame(snakemake.params["input_tree_name"], snakemake.input[0])
num_entries = rdf.Count().GetValue()

rdf.Range(
    round(ipartition * num_entries / npartitions),
    round((ipartition + 1) * num_entries / npartitions),
).Snapshot(snakemake.params["output_tree_name"], snakemake.output[0])


tee.stop()
