# Rebuild the summaries and rerun the additional studies.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && !identical(args, "--verify")) {
  stop("Use --verify to replay 44 original datasets before the additional studies.",
       call. = FALSE)
}

steps <- list(
  c("reproduce.R", args),
  "etude_complementaire.R", "etude_toroidale.R", "controle_calibration.R",
  "application_vent.R", "bandes_cv.R", "comparaison_vents_20261102.R",
  "etude_origine_tirages.R", "controle_independant_revision.R",
  "contre_verification_application.R", "contre_verification_modele.R",
  "verifier_rotation_20261102.R", "etude_rotation_20261102.R",
  "etude_sensibilite_20261102.R", "calibration_20261102.R",
  "valider_revision.R", "valider_renforcement.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
for (step in steps) {
  cat("\nRunning ", paste(step, collapse = " "), "\n", sep = "")
  status <- system2(rscript, c("--vanilla", shQuote(step)))
  if (status != 0L) {
    stop(paste("Failed:", paste(step, collapse = " ")), call. = FALSE)
  }
}
cat("\nReproduction checks passed.\n")
