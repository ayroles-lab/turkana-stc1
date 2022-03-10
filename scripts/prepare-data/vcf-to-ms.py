sys.path.append(snakemake.scriptdir + "/../..")
from utils.conversion import vcf_to_ms
vcf_to_ms(snakemake.input["data"], snakemake.output["ms"])
