import os

# Import the apd tools, vcersions customized for Snakemake
from apd.snakemake import get_analysis_data, remote

# Get the APD dataset for my analysis
dataset = get_analysis_data(config["wg"], config["analysis"])


rule fetch_apdata:
    input:
        lambda w: dataset(
            eventtype=w.eventtype,
            datatype=w.datatype,
            polarity=w.polarity,
        ),
    output:
        os.path.join(
            config.get("output_filepath_prefix", "data"),
            "{eventtype}_{datatype}_{polarity}.root",
        ),
    log:
        os.path.join(
            config.get("log_filepath_prefix", "logs"),
            "fetch_apdata_{eventtype}_{datatype}_{polarity}.log",
        ),
    threads: 8
    retries: 3
    wrapper:
        "v3.9.0/phys/root/hadd"


rule gather_apdata:
    input:
        data=expand(
            os.path.join(
                config.get("output_filepath_prefix", "data"),
                "{eventtype}_{datatype}_{polarity}.root",
            ),
            eventtype=config["eventtype"],
            datatype=config["datatype"],
            polarity=config["polarity"],
        ),
