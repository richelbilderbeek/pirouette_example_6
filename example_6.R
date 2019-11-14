# Code of example 6
#
# Works under Linux and MacOS only
#
#
#

# Set the RNG seed
rng_seed <- 314
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  arg <- suppressWarnings(as.numeric(args[1]))
  if (is.na(arg)) {
    stop(
      "Please supply a numerical value for the RNG seed. \n",
      "Actual value: ", args[1]
    )
  }
  rng_seed <- arg
  if (rng_seed < 1) {
    stop("Please supply an RNG seed with a positive non-zero value")
  }
}
if (length(args) > 1) {
  stop(
    "Please supply only 1 argument for the RNG seed. \n",
    "Number of arguments given: ", length(args) - 1
  )
}

library(pirouette)
suppressMessages(library(ggplot2))

root_folder <- getwd()
example_no <- 6
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())
phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)

alignment_params <- create_alignment_params(
  root_sequence = create_blocked_dna(length = 1000),
  rng_seed = rng_seed
)


# JC69, strict, Yule
generative_experiment <- create_gen_experiment()
check_experiment(generative_experiment)

# All non-Yule tree priors
candidate_experiments <- create_all_experiments(
  exclude_model = generative_experiment$inference_model
)
check_experiments(candidate_experiments)

experiments <- c(list(generative_experiment), candidate_experiments)

check_experiments(experiments)

# Set the RNG seed
for (i in seq_along(experiments)) {
  experiments[[i]]$beast2_options$rng_seed <- rng_seed
}

check_experiments(experiments)

# Testing
if (beastier::is_on_ci()) {
  experiments <- experiments[1:3]
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc <- create_mcmc(chain_length = 3000, store_every = 1000)
    experiments[[i]]$est_evidence_mcmc <- create_mcmc_nested_sampling(
      chain_length = 3000,
      store_every = 1000,
      epsilon = 100.0
    )
  }
}
check_experiments(experiments)

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments,
  twinning_params = create_twinning_params(
    rng_seed_twin_tree = rng_seed,
    rng_seed_twin_alignment = rng_seed
  )
)

check_pir_params(pir_params)

# Make Peregrine friendly
if (1 == 2) {
  pir_params <- peregrine::to_pff_pir_params(pir_params)
}
check_pir_params(pir_params)

rm_pir_param_files(pir_params)

errors <- pir_run(
  phylogeny,
  pir_params = pir_params
)

utils::write.csv(
  x = errors,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plot(errors) +
  ggsave(file.path(example_folder, "errors.png"))

pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)

pir_to_tables(
  pir_params = pir_params,
  folder = example_folder
)
