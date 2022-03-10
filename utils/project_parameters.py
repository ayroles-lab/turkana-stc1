"""Parameters of the current project and analysis."""

from .popgen_summary_statistics import (
    NumberOfSNPs,
    Pi,
    NumberOfHaplotypes,
    HaplotypeHomozygosity,
    H12,
    H2byH1,
    TajimasD,
)

default_summary_statistics = [
    NumberOfSNPs(),
    Pi(fraction=True),
    NumberOfHaplotypes(),
    HaplotypeHomozygosity(),
    H12(),
    H2byH1(),
    TajimasD(),
]

summary_statistic_order = ["pi", "num_snps", "num_haps", "H1", "H12", "H2/H1", "taj_D"]

locus_size = 1_074_835 
data_dimension = 13
smallest_window = 10000

default_simulation_parameters = {
    # Biology
    "locus-size": locus_size,
    "sample-size": 216,
    # SLiM options
    "last-generation": 5000,
    "simplification-interval": 500,
    "restarts-limit": 1000,
    # Feature calculation and normalizing
    "data-dimension": data_dimension,
    "smallest-window": smallest_window,
    # For recapitation demography
    # How many generations ago was there a size change?
    "size-change-generation": 7000,
    # What is N_ancestral/N_current?
    "size-change-factor": 1/2,
}

neural_network_batch_size = 64
