# Code of example 3
#
# Works under Linux and MacOS only
library(pirouette)
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))

root_folder <- getwd()
example_no <- 3
example_folder <- file.path(root_folder, paste0("example_", example_no))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(314)
testit::assert(is_beast2_installed())
phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)

alignment_params <- create_alignment_params(
  root_sequence = create_blocked_dna(length = 1000),
  mutation_rate = 0.1
)


# JC69, strict, Yule
generative_experiment <- create_experiment(
  inference_conditions = create_inference_conditions(
    model_type = "generative",
    run_if = "always",
    do_measure_evidence = TRUE
  ),
  inference_model = create_inference_model(
    site_model = create_jc69_site_model(),
    clock_model = create_strict_clock_model(),
    tree_prior = create_yule_tree_prior(),
    mcmc = create_mcmc(chain_length = 10e+7, store_every = 1000)
  )
)
generative_experiment <- create_gen_experiment()
check_experiment(generative_experiment)

# All non-Yule tree priors
candidate_experiments <- create_all_experiments(
  exclude_model = generative_experiment$inference_model
)
check_experiments(candidate_experiments)

experiments <- c(list(generative_experiment), candidate_experiments)
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

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments,
  twinning_params = create_twinning_params()
)

print("#######################################################################")
print("Settings to run on Peregrine cluster")
print("#######################################################################")
pir_params$alignment_params$fasta_filename <- file.path(example_folder, "true.fasta")
for (i in seq_along(pir_params$experiments)) {
  pir_params$experiments[[i]]$beast2_options$input_filename <- file.path(example_folder, "beast2_input_best.xml")
  pir_params$experiments[[i]]$beast2_options$output_log_filename <- file.path(example_folder, "beast2_output_best.log")
  pir_params$experiments[[i]]$beast2_options$output_trees_filenames <- file.path(example_folder, "beast2_output_best.trees")
  pir_params$experiments[[i]]$beast2_options$output_state_filename <- file.path(example_folder, "beast2_output_best.xml.state")
  pir_params$experiments[[i]]$beast2_options$beast2_working_dir <- example_folder
  pir_params$experiments[[i]]$errors_filename <- file.path(example_folder, "error_best.csv")
}
pir_params$experiments[[1]]$beast2_options$input_filename <- file.path(example_folder, "beast2_input_gen.xml")
pir_params$experiments[[1]]$beast2_options$output_log_filename <- file.path(example_folder, "beast2_output_gen.log")
pir_params$experiments[[1]]$beast2_options$output_trees_filenames <- file.path(example_folder, "beast2_output_gen.trees")
pir_params$experiments[[1]]$beast2_options$output_state_filename <- file.path(example_folder, "beast2_output_gen.xml.state")
pir_params$experiments[[1]]$errors_filename <- file.path(example_folder, "error_gen.csv")
pir_params$evidence_filename <- file.path(example_folder, "evidence_true.csv")
if (!is_one_na(pir_params$twinning_params)) {
  pir_params$twinning_params$twin_tree_filename <- file.path(example_folder, "twin.tree")
  pir_params$twinning_params$twin_alignment_filename <- file.path(example_folder, "twin.fasta")
  pir_params$twinning_params$twin_evidence_filename <- file.path(example_folder, "evidence_twin.csv")
}
rm_pir_param_files(pir_params)
print("#######################################################################")

Sys.time()
errors <- pir_run(
  phylogeny,
  pir_params = pir_params
)
Sys.time()

if (1 == 2) {
  errors <- utils::read.csv(
    file = file.path(example_folder, "errors.csv")
  )
  check_pir_out(errors)
  pir_plot(pir_out = errors)
}

utils::write.csv(
  x = errors,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plot(errors) +
  ggsave(file.path(example_folder, "errors.png"))

print("#######################################################################")
print("Evidence")
print("#######################################################################")

# Evidence, true
df_evidences <- utils::read.csv(pir_params$evidence_filename)[, c(-1, -6)]
df_evidences$site_model_name <- plyr::revalue(df_evidences$site_model_name, c("JC69" = "JC", "TN93" = "TN"))
df_evidences$clock_model_name <- plyr::revalue(
  df_evidences$clock_model_name,
  c("strict" = "Strict", "relaxed_log_normal" = "RLN")
)
df_evidences$tree_prior_name <- plyr::revalue(
  df_evidences$tree_prior_name,
  c(
    "yule" = "Yule",
    "birth_death" = "BD",
    "coalescent_bayesian_skyline" = "CBS",
    "coalescent_constant_population" = "CCP",
    "coalescent_exp_population" = "CEP"
  )
)
names(df_evidences) <- c("Site model", "Clock model", "Tree prior", "log(evidence)", "Weight")

sink(file.path(example_folder, "evidence_true.latex"))
xtable::print.xtable(
  xtable::xtable(
    df_evidences,
    caption = "Evidences of example 5", label = "tab:evidences_example_5", digits = 3
  ),
  include.rownames = FALSE
)
sink()

# Evidence, twin
df_evidences <- utils::read.csv(pir_params$twinning_params$twin_evidence_filename)[, c(-1, -6)]
df_evidences$site_model_name <- plyr::revalue(df_evidences$site_model_name, c("JC69" = "JC", "TN93" = "TN"))
df_evidences$clock_model_name <- plyr::revalue(
  df_evidences$clock_model_name,
  c("strict" = "Strict", "relaxed_log_normal" = "RLN")
)
df_evidences$tree_prior_name <- plyr::revalue(
  df_evidences$tree_prior_name,
  c(
    "yule" = "Yule",
    "birth_death" = "BD",
    "coalescent_bayesian_skyline" = "CBS",
    "coalescent_constant_population" = "CCP",
    "coalescent_exp_population" = "CEP"
  )
)
names(df_evidences) <- c("Site model", "Clock model", "Tree prior", "log(evidence)", "Weight")

sink(file.path(example_folder, "evidence_twin.latex"))
xtable::print.xtable(
  xtable::xtable(
    df_evidences,
    caption = "Evidences of example 5, twin tree", label = "tab:evidences_example_5_twin", digits = 3
  ),
  include.rownames = FALSE
)
sink()

print("#######################################################################")
print("ESSes")
print("#######################################################################")
testit::assert(pir_params$experiments[[1]]$inference_model$mcmc$store_every != -1)
esses_gen <- tracerer::calc_esses(
  traces = tracerer::parse_beast_log(pir_params$experiments[[1]]$beast2_options$output_log_filename),
  sample_interval = pir_params$experiments[[1]]$inference_model$mcmc$store_every
)
esses_best <- tracerer::calc_esses(
  traces = tracerer::parse_beast_log(pir_params$experiments[[2]]$beast2_options$output_log_filename),
  sample_interval = pir_params$experiments[[1]]$inference_model$mcmc$store_every
)
esses_twin_gen <- tracerer::calc_esses(
  traces = tracerer::parse_beast_log(to_twin_filename(pir_params$experiments[[1]]$beast2_options$output_log_filename)),
  sample_interval = pir_params$experiments[[1]]$inference_model$mcmc$store_every
)
esses_twin_best <- tracerer::calc_esses(
  traces = tracerer::parse_beast_log(to_twin_filename(pir_params$experiments[[2]]$beast2_options$output_log_filename)),
  sample_interval = pir_params$experiments[[1]]$inference_model$mcmc$store_every
)

df_esses_gen <- data.frame(parameter = colnames(esses_gen), ESS = as.character(esses_gen))
df_esses_best <- data.frame(parameter = colnames(esses_best), ESS = as.character(esses_best))
df_esses_twin_gen <- data.frame(parameter = colnames(esses_twin_gen), ESS = as.character(esses_twin_gen))
df_esses_twin_best <- data.frame(parameter = colnames(esses_twin_best), ESS = as.character(esses_twin_best))

sink(file.path(example_folder, "esses_gen.latex"))
xtable::print.xtable(
  xtable::xtable(
    df_esses_gen,
    caption = paste0("ESSes of example ", example_no, " for generative model"),
    label = paste0("tab:esses_example_", example_no, "_gen"),
    digits = 0
  ),
  include.rownames = FALSE
)
sink()

sink(file.path(example_folder, "esses_best.latex"))
xtable::print.xtable(
  xtable::xtable(
    df_esses_best,
    caption = paste0("ESSes of example ", example_no, " for best candidate model"),
    label = paste0("tab:esses_example_", example_no, "_best"),
    digits = 0
  ),
  include.rownames = FALSE
)
sink()

sink(file.path(example_folder, "esses_twin_gen.latex"))
xtable::print.xtable(
  xtable::xtable(
    df_esses_twin_gen,
    caption = paste0("ESSes of example ", example_no, " for generative model, twin tree"),
    label = paste0("tab:esses_example_", example_no, "_twin_gen"),
    digits = 0
  ),
  include.rownames = FALSE
)
sink()

sink(file.path(example_folder, "esses_twin_best.latex"))
xtable::print.xtable(
  xtable::xtable(
    df_esses_twin_best,
    caption = paste0("ESSes of example ", example_no, " for best candidate model, twin tree"),
    label = paste0("tab:esses_example_", example_no, "_twin__best"),
    digits = 0
  ),
  include.rownames = FALSE
)
sink()

print("#######################################################################")
print("Appendix")
print("#######################################################################")
pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)
