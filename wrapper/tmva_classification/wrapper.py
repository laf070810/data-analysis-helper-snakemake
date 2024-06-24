__author__ = "Anfeng Li"
__copyright__ = "Copyright 2024, Anfeng Li"
__email__ = "anfeng.li@cern.ch"
__license__ = "MIT"


import os

import ROOT
from pytee2 import Tee

tee = Tee(output_filepath=snakemake.log[0])
tee.start()

signal_chain = ROOT.TChain(snakemake.params.signal_tree_name)
for signal_filepath in snakemake.input.signal_filepaths:
    signal_chain.Add(signal_filepath)

background_chain = ROOT.TChain(snakemake.params.background_tree_name)
for background_filepath in snakemake.input.background_filepaths:
    background_chain.Add(background_filepath)

dataloader = ROOT.TMVA.DataLoader("TMVAdataset")
for training_variable_name in snakemake.params.training_variables:
    dataloader.AddVariable(training_variable_name, "F")
for spectator_variable_name in snakemake.params.get("spectator_variables", []):
    dataloader.AddSpectator(spectator_variable_name)
dataloader.AddSignalTree(signal_chain, 1.0)
dataloader.AddBackgroundTree(background_chain, 1.0)
dataloader.PrepareTrainingAndTestTree(
    "",
    "",
    snakemake.params.dataset_options,
)

output_file = ROOT.TFile(snakemake.output.output_filepath, "recreate")

ROOT.TMVA.gConfig().GetIONames().fWeightFileDirPrefix = os.path.join(
    os.path.dirname(snakemake.output.BDTG_weights_filepath), "../.."
)

if snakemake.params.get("use_cv", False):
    cv = ROOT.TMVA.CrossValidation(
        "TMVAClassification", dataloader, output_file, snakemake.params.factory_options
    )
    cv.BookMethod(
        ROOT.TMVA.Types.kBDT,
        "BDTG",
        snakemake.params.method_options,
    )
    cv.Evaluate()
else:
    factory = ROOT.TMVA.Factory(
        "TMVAClassification",
        output_file,
        snakemake.params.factory_options,
    )
    factory.BookMethod(
        dataloader,
        ROOT.TMVA.Types.kBDT,
        "BDTG",
        snakemake.params.method_options,
    )
    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()

output_file.Close()

tee.stop()
