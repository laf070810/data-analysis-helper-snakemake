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
        "data/input/{eventtype}_{datatype}_{polarity}.root",
    log:
        "logs/fetch_apdata_{eventtype}_{datatype}_{polarity}.log",
    threads: 8
    wrapper:
        "v3.9.0/phys/root/hadd"


rule gather_apdata:
    input:
        data=expand(
            "data/input/{eventtype}_{datatype}_{polarity}.root",
            eventtype=config["eventtype"],
            datatype=config["datatype"],
            polarity=config["polarity"],
        ),
