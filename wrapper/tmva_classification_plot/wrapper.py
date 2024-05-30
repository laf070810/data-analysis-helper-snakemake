__author__ = "Anfeng Li"
__copyright__ = "Copyright 2024, Anfeng Li"
__email__ = "anfeng.li@cern.ch"
__license__ = "MIT"


import os

import ROOT
from pytee2 import Tee

tee = Tee(output_filepath=snakemake.log[0])
tee.start()

input_absfilepath = os.path.abspath(snakemake.input.input_filepath)
os.chdir(os.path.dirname(os.path.dirname(os.path.dirname(snakemake.output[0]))))

ROOT.TMVA.variables("TMVAdataset", input_absfilepath)
ROOT.TMVA.correlations("TMVAdataset", input_absfilepath)
ROOT.TMVA.mvas("TMVAdataset", input_absfilepath, ROOT.TMVA.kMVAType)
ROOT.TMVA.mvas("TMVAdataset", input_absfilepath, ROOT.TMVA.kCompareType)
ROOT.TMVA.mvaeffs("TMVAdataset", input_absfilepath)

tee.stop()
