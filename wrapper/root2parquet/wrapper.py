__author__ = "Anfeng Li"
__copyright__ = "Copyright 2024, Anfeng Li"
__email__ = "anfeng.li@cern.ch"
__license__ = "MIT"


import pandas as pd
import ROOT
from pytee2 import Tee

tee = Tee(output_filepath=snakemake.log[0])
tee.start()

rdf = ROOT.RDataFrame(snakemake.params["input_tree_name"], snakemake.input[0])
pd.DataFrame(
    {str(key): val for key, val in rdf.AsNumpy(list(rdf.GetColumnNames())).items()}
).to_parquet(snakemake.output[0], compression="zstd")

tee.stop()
