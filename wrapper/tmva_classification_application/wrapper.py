__author__ = "Anfeng Li"
__copyright__ = "Copyright 2024, Anfeng Li"
__email__ = "anfeng.li@cern.ch"
__license__ = "MIT"


import os
import tempfile

from pytee2 import Tee

tee = Tee(output_filepath=snakemake.log[0])
tee.start()


root_codes = f"""
using namespace TMVA::Experimental;
ROOT::EnableImplicitMT({snakemake.threads});
RReader model("{os.path.abspath(snakemake.input.mva_weights_filepath)}");
auto variables = model.GetVariableNames();
ROOT::RDataFrame rdf("{snakemake.params.input_tree_name}", "{os.path.abspath(snakemake.input.input_filepath)}");
auto rdf2 = rdf.Define("{snakemake.params.mva_response_name}", Compute<{snakemake.params.mva_variable_num}, float>(model), variables);
rdf2.Snapshot("{snakemake.params.output_tree_name}", "{os.path.abspath(snakemake.output.output_filepath)}");
"""


with tempfile.NamedTemporaryFile(delete=False) as fp:
    fp.write(bytearray(f"void {os.path.basename(fp.name)}() {{", encoding="utf-8"))
    fp.write(bytearray(root_codes, encoding="utf-8"))
    fp.write(bytearray(f"}}", encoding="utf-8"))
    fp.close()
    os.system(f"root -l -q {fp.name}")

tee.stop()
