
library(MPStats)

cat("Reading in Zprime by plate data\n")

Zprime_by_plate <- MPStats::read_Zprime_by_plate("raw_data/master55.csv")
save(compound_moa, file="intermediate_data/Zprime_by_plate.Rdata")