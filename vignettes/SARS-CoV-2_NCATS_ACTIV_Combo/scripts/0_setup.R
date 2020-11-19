

library(batchtools)


if (!dir.exists(paths="/scratch/maom_root/maom99/maom/SARS-CoV-2_iAEC2_Combo")) {
    cat("Creating work direcotry '/scratch/maom_root/maom99/maom/SARS-CoV-2_iAEC2_Combo'\n")
    dir.create("/scratch/maom_root/maom99/maom/SARS-CoV-2_iAEC2_Combo")
}

batchtools_registry <- batchtools::makeRegistry(
    file.dir = "intermediate_data/batchtools_registry",
    work.dir = "/scratch/maom_root/maom99/maom/SARS-CoV-2_iAEC2_Combo",
    seed = 22336)

# https://mllg.github.io/batchtools/reference/makeClusterFunctionsSlurm
batchtools_registry$cluster.functions = batchtools::makeClusterFunctionsSlurm(
  template = "scripts/batchtools_slurm-greatlakes.tmpl",
  array.jobs = TRUE,
  nodename = "localhost",
  scheduler.latency = 1,
  fs.latency = 65)


batchtools::saveRegistry(
    reg = batchtools_registry)
