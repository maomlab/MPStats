library(MPStats)

# Load known compound mechanism of action
# and merge it with the well scores


load("intermediate_data/well_scores.Rdata")

compound_moa <- MPStats::read_compound_moa("raw_data/AGM_moa.csv")
save(compound_moa, file="intermediate_data/compound_moa.Rdata")


