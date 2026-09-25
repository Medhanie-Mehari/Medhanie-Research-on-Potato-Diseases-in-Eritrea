# ==============================================================================
# Disease prevalence, spatial distribution and environmental correlates
# — analysis of plant-level field survey data
# ==============================================================================
# Script  : 01_prevalence_analysis.R
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
# Tested  : R 4.4.1 (see SessionInfo.txt written at the end of the run)
#
# Output  : up to 8 PNG figures (400 dpi), 15 CSV tables, 1 Excel workbook
#
# All results, tables and figures are computed from the input data supplied
# in Section 0. No data are distributed with this script.
#
# Required input columns (one row per plant-level observation):
#   Disease_Name, Disease_Type (Fungal/Viral/Bacterial/Other/Unidentified),
#   Severity (ordinal score), Is_Named, Is_Specific_Named, Is_Other (TRUE/FALSE),
#   Region, Subzone, Subzone_Display, Name (farm identifier), Farm_ID,
#   Longitude, Latitude, Temperature_C, RH_pct, Soil_pH, Altitude, PlantAge,
#   Density_Category (Low/Medium/High), Density_Numeric, Plants_Assessed,
#   Area_Ha_farm, Farm_Incidence_pct, Subzone_N_Total, Subzone_Incidence_pct,
#   Subzone_Mean_DSI, Subzone_Disease_Richness,
#   Inc_<Disease> (farm-level 0/1 incidence, one column per disease)
#
# Implementation notes:
#   1. MASS is loaded before dplyr so that dplyr::select() is not masked.
#   2. Boolean columns exported as "True"/"False" text are converted to logical.
#   3. All dplyr verbs are namespaced (dplyr::) throughout.
#   4. Label offsets for Figure 3 are pre-computed (dsi_offset) outside aes().
#   5. Moran's I uses Longitude/Latitude vectors directly.
#   6. POMs are fitted for every disease with n >= MIN_N_POM observations.
# ==============================================================================


# ── 0. CONFIGURATION — replace the placeholders before running ───────────────

DATA_FILE <- "<PATH_TO_INPUT_DATA>/<INPUT_DATA_FILE>.csv"
OUT_DIR   <- "<PATH_TO_OUTPUT_FOLDER>"

# Environmental predictors used in both regression models
PREDICTORS <- c("Temperature_C", "RH_pct", "Soil_pH", "Altitude", "PlantAge", "Density_Numeric")

# Minimum number of observations for a disease to be modelled by the POM
MIN_N_POM <- 30

# Number of most frequent diseases shown in the frequency figures
N_TOP <- 10

# Optional: a disease tested for occurring at lower altitude than all others
# (Section 3, Figure 7). Leave the placeholder to skip these tests.
FOCAL_DISEASE     <- "<FOCAL_DISEASE_NAME>"
FOCAL_ALTERNATIVE <- "less"   # "less", "greater" or "two.sided"

# Optional: historical data for Figure 6. Leave the placeholders to skip.
#   Incidence file columns : Year, Group, Incidence_pct
#   Taxa file columns      : Year, Disease_Type, Cumulative_n
HIST_INCIDENCE_FILE <- "<PATH_TO_HISTORICAL_INCIDENCE_FILE>.csv"
HIST_TAXA_FILE      <- "<PATH_TO_CUMULATIVE_TAXA_FILE>.csv"

if (!file.exists(DATA_FILE))
  stop("Input data not found. Set DATA_FILE in Section 0 to your data file.")

FIG_DIR <- file.path(OUT_DIR, "Figures")
TAB_DIR <- file.path(OUT_DIR, "Tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)


# ── 0. PACKAGES ──────────────────────────────────────────────────────────────

need <- c("MASS", "DescTools", "brant", "spdep", "dplyr", "tidyr",
          "ggplot2", "patchwork", "scales", "openxlsx", "dunn.test")
for (p in need)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

# Load MASS first, then dplyr, so dplyr::select() is not masked by MASS::select()
suppressPackageStartupMessages({
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(DescTools)
  library(brant)
  library(spdep)
  library(dunn.test)
  library(openxlsx)
})

cat("=== Analysis started:", format(Sys.time()), "===\n")


# ── Helpers ──────────────────────────────────────────────────────────────────

save_fig <- function(p, name, w = 16, h = 8, dpi = 400) {
  ggsave(file.path(FIG_DIR, name), plot = p,
         width = w, height = h, dpi = dpi, units = "in", bg = "white")
  cat("  Figure:", name, "\n")
}

save_tab <- function(d, name) {
  write.csv(d, file.path(TAB_DIR, name), row.names = FALSE, fileEncoding = "UTF-8")
  cat("  Table:", name, "\n")
}

fmt_p <- function(p) ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

# Readable labels for predictors in Figure 5 (unlisted predictors keep their name)
PRED_LABELS <- c(Temperature_C   = "Temperature (per 1°C)",
                 RH_pct          = "Relative Humidity (per 1%)",
                 Soil_pH         = "Soil pH (per unit)",
                 Altitude        = "Altitude (per 1 m)",
                 PlantAge        = "Plant Age (per day)",
                 Density_Numeric = "Plant Density (per unit)")
lab_pred <- function(x) ifelse(x %in% names(PRED_LABELS), PRED_LABELS[x], x)


# ── Colour palettes ──────────────────────────────────────────────────────────

CF <- "#2D8653"; CV <- "#B5342A"; CB <- "#2166AC"
CA <- "#E07B39"; CG <- "#7B4F9E"
PAL <- c(CB, CF, CA, CG, CV, "#4D4D4D", "#1B9E77", "#E7298A", "#A6761D", "#66A61E")

DCOL <- c(Fungal = CF, Viral = CV, Bacterial = CB)


# ── Publication theme ────────────────────────────────────────────────────────

PUB <- theme_classic(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold", size = 13, margin = margin(b = 8)),
    axis.title        = element_text(size = 12),
    axis.text         = element_text(size = 10.5, colour = "black"),
    legend.position   = "bottom",
    legend.text       = element_text(size = 10.5),
    legend.title      = element_text(size = 11, face = "bold"),
    legend.key.size   = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", colour = "#CCCCCC", linewidth = 0.4),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 10.5),
    plot.caption      = element_text(size = 9, colour = "grey45", hjust = 0, margin = margin(t = 6)),
    plot.margin       = margin(10, 14, 10, 10)
  )


# ==============================================================================
# SECTION 1: DATA LOADING AND PREPARATION
# ==============================================================================
cat("\n=== SECTION 1: Loading data ===\n")

raw <- read.csv(DATA_FILE, stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM",
                na.strings = c("NA", "", "N/A"))
cat(sprintf("  Loaded: %d rows x %d columns\n", nrow(raw), ncol(raw)))

required <- c("Disease_Name", "Disease_Type", "Severity", "Is_Named", "Is_Specific_Named",
              "Is_Other", "Region", "Subzone", "Subzone_Display", "Name", "Farm_ID",
              "Longitude", "Latitude", PREDICTORS, "Density_Category", "Plants_Assessed",
              "Area_Ha_farm", "Farm_Incidence_pct", "Subzone_N_Total",
              "Subzone_Incidence_pct", "Subzone_Mean_DSI", "Subzone_Disease_Richness")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols) > 0)
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

# Boolean columns exported as "True"/"False" text are read as character
bool_cols <- c("Is_Named", "Is_Specific_Named", "Is_Other",
               "Is_Unidentified", "Is_Symptomatic", "For_Incidence_Calc")
for (bc in bool_cols) {
  if (bc %in% names(raw)) {
    v <- raw[[bc]]
    if (is.character(v)) v <- v %in% c("True", "TRUE", "true", "1", "T")
    raw[[bc]] <- as.logical(v)
  }
}

# Factor and numeric conversions
raw <- dplyr::mutate(raw,
  Disease_Type     = factor(Disease_Type,
                            levels = c("Fungal", "Viral", "Bacterial", "Other", "Unidentified")),
  Density_Category = factor(Density_Category, levels = c("Low", "Medium", "High")),
  Severity         = as.integer(Severity),
  Altitude         = as.numeric(Altitude),
  Temperature_C    = as.numeric(Temperature_C),
  RH_pct           = as.numeric(RH_pct),
  Soil_pH          = as.numeric(Soil_pH),
  PlantAge         = as.integer(PlantAge),
  Density_Numeric  = as.numeric(Density_Numeric),
  Longitude        = as.numeric(Longitude),
  Latitude         = as.numeric(Latitude),
  dplyr::across(dplyr::starts_with("Inc_"), as.integer)
)

# Sub-datasets
df     <- raw
named  <- dplyr::filter(raw, Is_Named)           # records with a named disease
df_vis <- dplyr::filter(raw, Is_Specific_Named)  # records with a specific disease name
symp   <- dplyr::filter(raw, !Is_Other)          # excluding "Other"

# Farm-level dataset (one row per farm)
df_farm <- raw %>%
  dplyr::group_by(Name, Subzone, Region) %>%
  dplyr::summarise(
    Farm_ID            = dplyr::first(Farm_ID),
    Altitude           = mean(Altitude, na.rm = TRUE),
    Longitude          = mean(Longitude, na.rm = TRUE),
    Latitude           = mean(Latitude, na.rm = TRUE),
    Temperature_C      = mean(Temperature_C, na.rm = TRUE),
    RH_pct             = mean(RH_pct, na.rm = TRUE),
    Soil_pH            = mean(Soil_pH, na.rm = TRUE),
    PlantAge           = mean(PlantAge, na.rm = TRUE),
    Density_Category   = dplyr::first(Density_Category),
    Density_Numeric    = dplyr::first(Density_Numeric),
    Plants_Assessed    = dplyr::first(Plants_Assessed),
    Area_Ha_farm       = dplyr::first(Area_Ha_farm),
    Farm_Incidence_pct = dplyr::first(Farm_Incidence_pct),
    dplyr::across(dplyr::starts_with("Inc_"), dplyr::first),
    .groups = "drop"
  )

# Study-level counts used as denominators
N_NAMED  <- nrow(named)
N_FARMS  <- nrow(df_farm)
N_PLANTS <- raw %>%
  dplyr::distinct(Subzone, Subzone_N_Total) %>%
  dplyr::pull(Subzone_N_Total) %>%
  sum(na.rm = TRUE)

cat(sprintf("  df=%d | named=%d | df_vis=%d | df_farm=%d | plants assessed=%d\n",
            nrow(df), N_NAMED, nrow(df_vis), N_FARMS, N_PLANTS))


# ==============================================================================
# SECTION 2: DESCRIPTIVE STATISTICS
# ==============================================================================
cat("\n=== SECTION 2: Descriptive statistics ===\n")

dis_freq <- named %>%
  dplyr::count(Disease_Name, Disease_Type) %>%
  dplyr::mutate(Pct_of_named  = round(n / N_NAMED * 100, 1),
                Pct_of_plants = round(n / N_PLANTS * 100, 2)) %>%
  dplyr::arrange(dplyr::desc(n))

type_totals <- named %>%
  dplyr::group_by(Disease_Type) %>%
  dplyr::summarise(N_diseases   = dplyr::n_distinct(Disease_Name),
                   N_records    = dplyr::n(),
                   Pct_of_named = round(dplyr::n() / N_NAMED * 100, 1),
                   .groups = "drop")
cat("  Type totals:\n"); print(as.data.frame(type_totals))

dsi_disease <- symp %>%
  dplyr::filter(!is.na(Disease_Name)) %>%
  dplyr::group_by(Disease_Name, Disease_Type) %>%
  dplyr::summarise(n        = dplyr::n(),
                   DSI_mean = round(mean(Severity, na.rm = TRUE), 3),
                   DSI_SD   = round(sd(Severity, na.rm = TRUE), 3),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(DSI_mean))

sub_summary <- raw %>%
  dplyr::group_by(Region, Subzone, Subzone_Display) %>%
  dplyr::summarise(
    Farms            = dplyr::n_distinct(Name),
    Plants_Assessed  = dplyr::first(Subzone_N_Total),
    Pct_of_plants    = round(dplyr::first(Subzone_N_Total) / N_PLANTS * 100, 1),
    W_Inc_pct        = dplyr::first(Subzone_Incidence_pct),
    Mean_DSI_0to5    = dplyr::first(Subzone_Mean_DSI),
    Disease_Richness = dplyr::first(Subzone_Disease_Richness),
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(W_Inc_pct))

env_vars <- c("Altitude", "Temperature_C", "RH_pct", "Soil_pH", "PlantAge")
env_summary <- df_farm %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(env_vars),
    list(Min  = ~round(min(.x, na.rm = TRUE), 2),
         Max  = ~round(max(.x, na.rm = TRUE), 2),
         Mean = ~round(mean(.x, na.rm = TRUE), 3),
         SD   = ~round(sd(.x, na.rm = TRUE), 3)))) %>%
  tidyr::pivot_longer(dplyr::everything(),
                      names_to = c("Variable", ".value"),
                      names_sep = "_(?=[^_]+$)")

save_tab(dis_freq,    "T01_Disease_Frequency.csv")
save_tab(type_totals, "T02_Disease_Type_Totals.csv")
save_tab(dsi_disease, "T03_DSI_per_Disease.csv")
save_tab(sub_summary, "T04_Subregion_Summary.csv")
save_tab(env_summary, "T05_Environmental_Summary.csv")


# ==============================================================================
# SECTION 3: NON-PARAMETRIC TESTS
# ==============================================================================
cat("\n=== SECTION 3: Non-parametric tests ===\n")

np_rows <- list()
add_np <- function(test, stat, p)
  np_rows[[length(np_rows) + 1]] <<- data.frame(Test = test, Statistic = stat,
                                                p_num = p, stringsAsFactors = FALSE)

kw_sev <- kruskal.test(Severity ~ Disease_Type, data = named)
cat(sprintf("  KW Severity x Type: H=%.3f, p=%.4g\n",
            kw_sev$statistic, kw_sev$p.value))
suppressMessages(dunn.test::dunn.test(named$Severity, named$Disease_Type,
                                      method = "bonferroni", kw = FALSE, list = FALSE))
add_np("KW Severity x Type", round(kw_sev$statistic, 3), kw_sev$p.value)

kw_alt <- kruskal.test(Altitude ~ Disease_Type, data = named)
cat(sprintf("  KW Altitude x Type: H=%.3f, p=%.4f\n",
            kw_alt$statistic, kw_alt$p.value))
add_np("KW Altitude x Type", round(kw_alt$statistic, 3), kw_alt$p.value)

has_focal <- FOCAL_DISEASE %in% named$Disease_Name
if (has_focal) {
  foc_alt <- named$Altitude[named$Disease_Name == FOCAL_DISEASE]
  oth_alt <- named$Altitude[named$Disease_Name != FOCAL_DISEASE]
  mw_foc  <- wilcox.test(foc_alt, oth_alt, alternative = FOCAL_ALTERNATIVE, exact = FALSE)
  cat(sprintf("  MW %s vs others (%s): U=%.0f, p=%.4g\n",
              FOCAL_DISEASE, FOCAL_ALTERNATIVE, mw_foc$statistic, mw_foc$p.value))
  cat(sprintf("  %s mean alt=%.1f m | Others mean alt=%.1f m\n", FOCAL_DISEASE,
              mean(foc_alt, na.rm = TRUE), mean(oth_alt, na.rm = TRUE)))
  add_np(sprintf("MW %s altitude vs others (%s)", FOCAL_DISEASE, FOCAL_ALTERNATIVE),
         round(mw_foc$statistic, 0), mw_foc$p.value)

  sp_sub <- named %>%
    dplyr::filter(Disease_Name == FOCAL_DISEASE) %>%
    dplyr::group_by(Subzone) %>%
    dplyr::summarise(N_foc    = dplyr::n(),
                     Mean_Alt = mean(Altitude, na.rm = TRUE), .groups = "drop")
  if (nrow(sp_sub) >= 3) {
    sp_foc <- cor.test(sp_sub$Mean_Alt, sp_sub$N_foc, method = "spearman", exact = FALSE)
    cat(sprintf("  Spearman rho=%.3f, p=%.4f, n=%d sub-regions\n",
                sp_foc$estimate, sp_foc$p.value, nrow(sp_sub)))
    add_np(sprintf("Spearman rho (altitude vs %s frequency)", FOCAL_DISEASE),
           round(sp_foc$estimate, 3), sp_foc$p.value)
  }
} else {
  cat("  Focal-disease tests skipped (FOCAL_DISEASE not set or not found)\n")
}

sub_corr <- raw %>%
  dplyr::group_by(Subzone) %>%
  dplyr::summarise(W_Inc = dplyr::first(Subzone_Incidence_pct),
                   DSI   = dplyr::first(Subzone_Mean_DSI), .groups = "drop")
cor_id <- cor.test(sub_corr$W_Inc, sub_corr$DSI, method = "pearson")
cat(sprintf("  Pearson r=%.3f, p=%.4f (n=%d)\n", cor_id$estimate, cor_id$p.value, nrow(sub_corr)))
add_np("Pearson r (weighted incidence vs DSI)", round(cor_id$estimate, 3), cor_id$p.value)

nonpar_tab <- do.call(rbind, np_rows)
nonpar_tab$p_value <- formatC(nonpar_tab$p_num, digits = 4, format = "g")
nonpar_tab$Sig     <- ifelse(nonpar_tab$p_num < 0.05, "Yes", "No")
nonpar_tab$p_num   <- NULL
save_tab(nonpar_tab, "T06_NonParametric_Tests.csv")


# ==============================================================================
# SECTION 4: PROPORTIONAL ODDS MODEL — DISEASE SEVERITY
# ==============================================================================
cat("\n=== SECTION 4: POM ===\n")

pom_formula <- as.formula(paste("Sev_f ~", paste(PREDICTORS, collapse = " + ")))

run_pom <- function(disease_name) {
  sub <- dplyr::filter(named, Disease_Name == disease_name) %>%
    dplyr::mutate(Sev_f = factor(Severity, ordered = TRUE))
  n <- nrow(sub)
  cat(sprintf("\n  POM: %s (n=%d)\n", disease_name, n))
  if (n < MIN_N_POM) { cat(sprintf("  SKIP: n < %d\n", MIN_N_POM)); return(NULL) }

  mod <- tryCatch(
    MASS::polr(pom_formula, data = sub, Hess = TRUE, method = "logistic"),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL })
  if (is.null(mod)) return(NULL)

  ctab <- coef(summary(mod))
  beta <- ctab[, "Value"]; se <- ctab[, "Std. Error"]
  pv   <- pnorm(abs(beta / se), lower.tail = FALSE) * 2
  ci   <- tryCatch(exp(confint(mod)),
                   error = function(e) exp(cbind(beta - 1.96 * se, beta + 1.96 * se)))
  OR   <- exp(beta)
  idx  <- !grepl("\\|", rownames(ctab))
  # confint() returns predictor rows only; match CI rows to predictors by name
  ci   <- ci[rownames(ctab)[idx], , drop = FALSE]

  null_ll <- as.numeric(logLik(MASS::polr(Sev_f ~ 1, data = sub, Hess = TRUE)))
  full_ll <- as.numeric(logLik(mod))
  mcfR2   <- round(1 - full_ll / null_ll, 3)
  cat(sprintf("  AIC=%.2f | McFadden R2=%.3f\n", AIC(mod), mcfR2))

  cat("  Brant test:\n")
  brant_p <- tryCatch({
    br <- brant::brant(mod)
    round(as.numeric(br["Omnibus", "probability"]), 3)
  }, error = function(e) NA_real_)

  list(model   = mod,
       disease = disease_name,
       results = data.frame(
         Disease     = disease_name, n = n,
         Predictor   = rownames(ctab)[idx],
         OR          = round(OR[idx], 3),
         CI_L        = round(ci[, 1], 3),
         CI_H        = round(ci[, 2], 3),
         p_value     = round(pv[idx], 4),
         Sig         = ifelse(pv[idx] < 0.05, "Yes", "No"),
         AIC         = round(AIC(mod), 2),
         McFadden_R2 = mcfR2,
         stringsAsFactors = FALSE),
       fit = data.frame(Disease = disease_name, n = n,
                        AIC = round(AIC(mod), 2), McFadden_R2 = mcfR2,
                        Brant_p_omnibus = brant_p, stringsAsFactors = FALSE))
}

# Fit a POM for every named disease, most frequent first
pom_diseases <- dis_freq$Disease_Name[dis_freq$n >= MIN_N_POM]
pom_list <- Filter(Negate(is.null), lapply(pom_diseases, run_pom))

pom_all <- if (length(pom_list) > 0)
  do.call(rbind, lapply(pom_list, function(x) x$results)) else data.frame()
pom_fit <- if (length(pom_list) > 0)
  do.call(rbind, lapply(pom_list, function(x) x$fit)) else data.frame()

pom_sig <- if (nrow(pom_all) > 0) dplyr::filter(pom_all, Sig == "Yes") else pom_all
cat("\n  POM significant results:\n")
if (nrow(pom_sig) > 0)
  print(pom_sig[, c("Disease", "Predictor", "OR", "CI_L", "CI_H", "p_value")], row.names = FALSE)

save_tab(pom_all, "T07_POM_Results_Full.csv")
save_tab(pom_sig, "T08_POM_Results_Significant.csv")
save_tab(pom_fit, "T08b_POM_Model_Fit.csv")


# ==============================================================================
# SECTION 5: BINOMIAL LOGISTIC REGRESSION — DISEASE INCIDENCE
# ==============================================================================
cat(sprintf("\n=== SECTION 5: Logistic regression (n=%d farms) ===\n", N_FARMS))

run_logit <- function(outcome, label) {
  y <- df_farm[[outcome]]
  if (length(unique(stats::na.omit(y))) < 2) {
    cat("  SKIP:", label, "- outcome has no variation\n"); return(NULL)
  }
  fml <- as.formula(paste(outcome, "~", paste(PREDICTORS, collapse = " + ")))
  mod <- tryCatch(
    glm(fml, data = df_farm, family = binomial(link = "logit")),
    error = function(e) { cat("  ERROR:", label, "-", e$message, "\n"); NULL })
  if (is.null(mod)) return(NULL)

  cf  <- coef(summary(mod))[-1, , drop = FALSE]
  cis <- confint.default(mod)[-1, , drop = FALSE]
  nag <- tryCatch(round(DescTools::PseudoR2(mod, which = "Nagelkerke"), 3),
                  error = function(e) NA)
  cat(sprintf("  %s: AIC=%.2f | Nagelkerke R2=%.3f\n", label, AIC(mod), nag))

  list(model   = mod,
       disease = label,
       results = data.frame(
         Disease       = label, n_farms = N_FARMS,
         Predictor     = rownames(cf),
         OR            = round(exp(cf[, "Estimate"]), 3),
         CI_L          = round(exp(cis[, 1]), 3),
         CI_H          = round(exp(cis[, 2]), 3),
         p_value       = round(cf[, "Pr(>|z|)"], 4),
         Sig           = ifelse(cf[, "Pr(>|z|)"] < 0.05, "Yes", "No"),
         Nagelkerke_R2 = nag,
         stringsAsFactors = FALSE),
       fit = data.frame(Disease = label, n_farms = N_FARMS,
                        AIC = round(AIC(mod), 2), Nagelkerke_R2 = nag,
                        stringsAsFactors = FALSE))
}

# One model per farm-level incidence column (Inc_<Disease>)
inc_cols   <- grep("^Inc_", names(df_farm), value = TRUE)
logit_list <- Filter(Negate(is.null),
                     lapply(inc_cols, function(v) run_logit(v, sub("^Inc_", "", v))))

logit_all <- if (length(logit_list) > 0)
  do.call(rbind, lapply(logit_list, function(x) x$results)) else data.frame()
logit_fit <- if (length(logit_list) > 0)
  do.call(rbind, lapply(logit_list, function(x) x$fit)) else data.frame()

logit_sig <- if (nrow(logit_all) > 0) dplyr::filter(logit_all, Sig == "Yes") else logit_all
cat("\n  Logistic significant results:\n")
if (nrow(logit_sig) > 0)
  print(logit_sig[, c("Disease", "Predictor", "OR", "CI_L", "CI_H", "p_value")], row.names = FALSE)

save_tab(logit_all, "T09_Logistic_Results_Full.csv")
save_tab(logit_sig, "T10_Logistic_Results_Significant.csv")
save_tab(logit_fit, "T10b_Logistic_Model_Fit.csv")


# ==============================================================================
# SECTION 6: MORAN'S I — SPATIAL AUTOCORRELATION OF MODEL RESIDUALS
# ==============================================================================
cat("\n=== SECTION 6: Moran's I ===\n")

build_weights <- function(lon, lat, k = 5) {
  coords <- cbind(lon, lat)
  ok     <- complete.cases(coords)
  coords <- coords[ok, , drop = FALSE]
  if (nrow(coords) < k + 1) return(NULL)
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  spdep::nb2listw(nb, style = "W")
}

run_mi <- function(model, label, lon_vec, lat_vec) {
  wts <- tryCatch(build_weights(lon_vec, lat_vec), error = function(e) NULL)
  if (is.null(wts)) { cat("  SKIP:", label, "\n"); return(NULL) }
  res <- tryCatch(residuals(model), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  mi <- tryCatch(spdep::moran.test(res, wts, randomisation = TRUE, alternative = "two.sided"),
                 error = function(e) NULL)
  if (is.null(mi)) return(NULL)
  cat(sprintf("  %s: I=%.3f, p=%.4f %s\n", label, mi$statistic, mi$p.value,
              ifelse(mi$p.value > 0.05, "[OK - independent]", "[!]")))
  data.frame(Model = label, Morans_I = round(mi$statistic, 3),
             p_value = round(mi$p.value, 4), OK = mi$p.value > 0.05,
             stringsAsFactors = FALSE)
}

mi_rows <- c(
  lapply(pom_list, function(x) {
    d <- dplyr::filter(named, Disease_Name == x$disease)
    run_mi(x$model, paste("POM", x$disease), d$Longitude, d$Latitude)
  }),
  lapply(logit_list, function(x)
    run_mi(x$model, paste("Logistic", x$disease), df_farm$Longitude, df_farm$Latitude))
)
morans_all <- do.call(rbind, Filter(Negate(is.null), mi_rows))
if (!is.null(morans_all) && nrow(morans_all) > 0)
  save_tab(morans_all, "T11_Morans_I.csv")


# ==============================================================================
# SECTION 7: FIGURES
# ==============================================================================
cat("\n=== SECTION 7: Figures ===\n")

freq_rank <- dis_freq %>%
  dplyr::group_by(Disease_Name) %>%
  dplyr::summarise(n = sum(n), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::pull(Disease_Name)
TOP10 <- head(freq_rank, N_TOP)
TOP8  <- head(freq_rank, 8)
TOP6  <- head(freq_rank, 6)

# ─── FIGURE 2 — Disease frequency ────────────────────────────────────────────
f2_dat <- named %>%
  dplyr::filter(Disease_Name %in% TOP10) %>%
  dplyr::count(Disease_Name, Disease_Type) %>%
  dplyr::mutate(
    Disease_Name = factor(Disease_Name, levels = rev(TOP10)),
    pct   = round(n / N_NAMED * 100, 1),
    label = paste0(n, " (", pct, "%)"))

p2 <- ggplot(f2_dat, aes(x = n, y = Disease_Name, fill = Disease_Type)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.8,
            fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = DCOL, name = "Category") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.33))) +
  labs(title = sprintf("Frequency of the Top %d Diseases by Plant-Level Observation Count", length(TOP10)),
       x = "Number of Plant Observations", y = NULL,
       caption = sprintf("n = %d named records", N_NAMED)) +
  PUB +
  theme(axis.text.y = element_text(size = 11.5), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
save_fig(p2, "Figure_2_Disease_Frequency.png", 13, 7)

# ─── FIGURE 3 — Incidence and DSI by sub-region ──────────────────────────────
# dsi_offset is pre-computed here because dplyr::n() cannot be used inside aes()
sub_f3 <- raw %>%
  dplyr::group_by(Subzone, Subzone_Display, Region) %>%
  dplyr::summarise(W_Inc = dplyr::first(Subzone_Incidence_pct),
                   DSI   = dplyr::first(Subzone_Mean_DSI), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(W_Inc)) %>%
  dplyr::mutate(
    Subzone_Display = factor(Subzone_Display, levels = unique(Subzone_Display)),
    dsi_offset      = ifelse(dplyr::row_number() %% 2 == 0, 1.5, -1.5))

regions <- sort(unique(as.character(sub_f3$Region)))
RCOL    <- setNames(rep(PAL, length.out = length(regions)), regions)
y3_max  <- max(c(sub_f3$W_Inc + 2, sub_f3$DSI * 10 + 3), na.rm = TRUE)

p3 <- ggplot(sub_f3, aes(x = Subzone_Display)) +
  geom_col(aes(y = W_Inc, fill = Region), width = 0.65,
           colour = "white", linewidth = 0.3, alpha = 0.85) +
  geom_line(aes(y = DSI * 10, group = 1),
            colour = CA, linewidth = 2.0, lineend = "round") +
  geom_point(aes(y = DSI * 10), colour = CA, size = 4.5,
             shape = 21, fill = "white", stroke = 2.0) +
  geom_text(aes(y = W_Inc + 0.7, label = sprintf("%.1f", W_Inc)),
            size = 2.9, fontface = "bold", colour = "grey20") +
  geom_text(aes(y = DSI * 10 + dsi_offset, label = sprintf("%.2f", DSI)),
            colour = CA, size = 2.9, fontface = "bold") +
  scale_fill_manual(values = RCOL, name = "Region") +
  scale_y_continuous(
    name = "Weighted Disease Incidence (%)", limits = c(0, y3_max),
    sec.axis = sec_axis(~ . / 10, name = "Mean DSI (0–5 raw ordinal score)")) +
  labs(title = "Area-Weighted Disease Incidence (%) and Mean Disease Severity Index (0–5) by Sub-region",
       x = NULL) +
  PUB +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.5),
        axis.title.y.right = element_text(colour = CA, size = 12))
save_fig(p3, "Figure_3_Incidence_DSI_Subregion.png", 18, 8)

# ─── FIGURE 4 — Disease proportion heatmap ───────────────────────────────────
row_ord <- sub_f3 %>%
  dplyr::arrange(dplyr::desc(W_Inc)) %>%
  dplyr::pull(Subzone_Display) %>% as.character() %>% unique()

hm4 <- df_vis %>%
  dplyr::filter(Disease_Name %in% TOP10) %>%
  dplyr::count(Subzone_Display, Disease_Name) %>%
  dplyr::group_by(Subzone_Display) %>%
  dplyr::mutate(pct = round(n / sum(n) * 100, 1)) %>%
  dplyr::ungroup() %>%
  tidyr::complete(Subzone_Display, Disease_Name = TOP10, fill = list(n = 0, pct = NA)) %>%
  dplyr::mutate(Disease_Name    = factor(Disease_Name, levels = TOP10),
                Subzone_Display = factor(Subzone_Display, levels = rev(row_ord)))

p4 <- ggplot(hm4, aes(x = Disease_Name, y = Subzone_Display, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = ifelse(is.na(pct), "", sprintf("%.1f", pct))),
            size = 3.6, colour = "grey10", na.rm = TRUE) +
  scale_fill_gradientn(
    colours = c("#FFFDE7", "#FFF176", "#FFCA28", "#FB8C00", "#E64A19", "#BF360C"),
    na.value = "#EFEFEF", name = "Proportion (%)",
    limits = c(0, max(hm4$pct, na.rm = TRUE)),
    guide = guide_colorbar(barwidth = 18, barheight = 0.8, title.position = "top")) +
  labs(title = sprintf("Proportion (%%) of Top %d Disease Observations per Sub-region", length(TOP10)),
       x = "Disease", y = NULL,
       caption = sprintf("Grey = not recorded | n = %d specific named records", nrow(df_vis))) +
  theme_minimal(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 13, margin = margin(b = 10)),
        plot.caption = element_text(size = 9, colour = "grey45", hjust = 0),
        axis.text.x  = element_text(angle = 38, hjust = 1, size = 11.5),
        axis.text.y  = element_text(size = 11),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", plot.margin = margin(10, 12, 10, 10))
save_fig(p4, "Figure_4_Disease_Heatmap.png", 18, 8.5)

# ─── FIGURE 5 — Forest plot of significant predictors (computed above) ───────
forest_dat <- function(sig_tab) {
  if (nrow(sig_tab) == 0) return(sig_tab)
  d <- dplyr::filter(sig_tab, is.finite(CI_L), is.finite(CI_H), CI_L > 0, OR > 0)
  d$Predictor <- lab_pred(d$Predictor)
  d$pcol <- ifelse(d$OR < 1, "#2166AC", "#B5342A")
  d$plab <- paste0(fmt_p(d$p_value), stars(d$p_value))
  d
}
no_sig <- function(all_tab, sig_tab)
  setdiff(unique(all_tab$Disease), unique(sig_tab$Disease))

POM_F <- forest_dat(pom_sig)
LOG_F <- forest_dat(logit_sig)
panels <- list()

if (nrow(POM_F) > 0) {
  POM_F$Disease <- factor(POM_F$Disease, levels = unique(POM_F$Disease))
  pom_lab <- with(pom_fit, setNames(
    sprintf("%s (n=%d; AIC=%.2f; McFadden R²=%.2f)", Disease, n, AIC, McFadden_R2), Disease))
  ns_pom  <- no_sig(pom_all, pom_sig)

  panels$A <- ggplot(POM_F, aes(x = OR, y = reorder(Predictor, -OR), xmin = CI_L, xmax = CI_H)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.8) +
    geom_errorbarh(aes(colour = pcol), height = 0.22, linewidth = 1.8, alpha = 0.9) +
    geom_point(aes(colour = pcol), size = 5, shape = 23, fill = "white", stroke = 2) +
    geom_text(aes(x = CI_H, label = paste0(" ", plab)), hjust = 0, size = 3.5, colour = "grey25") +
    scale_colour_identity() +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.55))) +
    facet_wrap(~Disease, scales = "free_y", ncol = min(3, nlevels(POM_F$Disease)),
               labeller = labeller(Disease = pom_lab)) +
    labs(title = "(A) Disease Severity — Proportional Odds Model | Plant-level observations",
         x = "Odds Ratio (log scale) with 95% profile-likelihood CI", y = NULL,
         caption = paste0("Blue = OR<1 p<0.05 (protective) | Red = OR>1 p<0.05 (positive)",
                          if (length(ns_pom) > 0)
                            paste0(" | No significant predictors: ", paste(ns_pom, collapse = ", ")))) +
    PUB + theme(panel.grid.major.y = element_blank(), strip.text = element_text(size = 9.5))
}

if (nrow(LOG_F) > 0) {
  LOG_F$Predictor <- sub(" \\(", "\n(", LOG_F$Predictor)
  LOG_F$Disease   <- factor(LOG_F$Disease, levels = unique(LOG_F$Disease))
  log_lab <- with(logit_fit, setNames(
    sprintf("%s\n(Nag.R²=%.2f)", Disease, Nagelkerke_R2), Disease))

  panels$B <- ggplot(LOG_F, aes(x = OR, y = reorder(Predictor, -OR), xmin = CI_L, xmax = CI_H)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.8) +
    geom_errorbarh(aes(colour = pcol), height = 0.22, linewidth = 1.8, alpha = 0.9) +
    geom_point(aes(colour = pcol), size = 5, shape = 23, fill = "white", stroke = 2) +
    geom_text(aes(x = pmax(CI_H, OR * 1.05), label = paste0(" ", plab)),
              hjust = 0, size = 3.4, colour = "grey25") +
    scale_colour_identity() +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.60))) +
    facet_wrap(~Disease, scales = "free_y", ncol = min(4, nlevels(LOG_F$Disease)),
               labeller = labeller(Disease = log_lab)) +
    labs(title = sprintf("(B) Disease Incidence — Binomial Logistic Regression | Farm-level binary (n=%d farms)", N_FARMS),
         x = "Odds Ratio (log scale) with 95% Wald CI", y = NULL,
         caption = "OR = exp(β); β = log-odds change per unit — NOT a probability change (Hosmer et al., 2013)") +
    PUB + theme(panel.grid.major.y = element_blank(), strip.text = element_text(size = 9.5))
}

if (length(panels) > 0) {
  p5 <- patchwork::wrap_plots(panels, ncol = 1)
  save_fig(p5, "Figure_5_ForestPlot_POM_Logistic.png", 19, 7 * length(panels))
} else {
  cat("  Figure 5 skipped: no significant predictors in either model\n")
}

# ─── FIGURE 6 — Historical trends (optional external data) ───────────────────
p6_parts <- list()

if (file.exists(HIST_INCIDENCE_FILE)) {
  av <- read.csv(HIST_INCIDENCE_FILE, stringsAsFactors = FALSE)
  av$Group <- factor(av$Group, levels = unique(av$Group))
  grp  <- levels(av$Group)
  p6_parts$a <- ggplot(av, aes(Year, Incidence_pct, colour = Group, shape = Group)) +
    geom_line(linewidth = 2.0) + geom_point(size = 3.5, stroke = 1.5, fill = "white") +
    scale_colour_manual(values = setNames(rep(PAL, length.out = length(grp)), grp)) +
    scale_shape_manual(values = setNames(rep(c(21, 22, 24, 23, 25), length.out = length(grp)), grp)) +
    labs(title = "(A) Historical Disease Incidence", x = "Year", y = "Incidence (%)", colour = NULL, shape = NULL) + PUB
}

if (file.exists(HIST_TAXA_FILE)) {
  cl <- read.csv(HIST_TAXA_FILE, stringsAsFactors = FALSE)
  cl$Disease_Type <- factor(cl$Disease_Type, levels = intersect(names(DCOL), unique(cl$Disease_Type)))
  p6_parts$b <- ggplot(cl, aes(Year, Cumulative_n, colour = Disease_Type, shape = Disease_Type)) +
    geom_line(linewidth = 2.0) + geom_point(size = 3.5, stroke = 1.5, fill = "white") +
    scale_colour_manual(values = DCOL) + scale_shape_manual(values = c(Fungal = 21, Viral = 22, Bacterial = 24)) +
    labs(title = "(B) Cumulative Disease Taxa Documented", x = "Year", y = "Cumulative taxa (n)",
         colour = "Type", shape = "Type") + PUB
}

if (length(p6_parts) > 0) {
  p6 <- patchwork::wrap_plots(p6_parts, nrow = 1) +
    plot_annotation(title = "Historical Trends in Disease Documentation",
                    theme = theme(plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6))))
  save_fig(p6, "Figure_6_Historical_Trends.png", 8.5 * length(p6_parts), 7)
} else {
  cat("  Figure 6 skipped: historical data files not provided\n")
}

# ─── FIGURE 7 — Altitude distribution ────────────────────────────────────────
alt8 <- named %>%
  dplyr::filter(Disease_Name %in% TOP8) %>%
  dplyr::mutate(Disease_Name = factor(Disease_Name,
    levels = names(sort(tapply(Altitude, Disease_Name, median, na.rm = TRUE), decreasing = TRUE))))

set.seed(42)
p7 <- ggplot(alt8, aes(x = Disease_Name, y = Altitude, fill = Disease_Type)) +
  geom_boxplot(alpha = 0.60, width = 0.55, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(colour = Disease_Type), width = 0.14, alpha = 0.20, size = 1.6, shape = 16) +
  scale_fill_manual(values = DCOL, name = "Category") +
  scale_colour_manual(values = DCOL, guide = "none") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = sprintf("Altitudinal Distribution of Disease Occurrences Across the Top %d Most Prevalent Diseases",
                       length(TOP8)),
       x = NULL, y = "Altitude (m a.s.l.)") +
  PUB + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 11.5))

if (has_focal) {
  foc_mean <- mean(named$Altitude[named$Disease_Name == FOCAL_DISEASE], na.rm = TRUE)
  p7 <- p7 +
    geom_hline(yintercept = foc_mean, colour = "#B5342A", linetype = "dashed", linewidth = 1.5, alpha = 0.85) +
    labs(caption = sprintf("Dashed red line = mean altitude of %s observations (%.0f m a.s.l.)",
                           FOCAL_DISEASE, foc_mean))
}
save_fig(p7, "Figure_7_Altitude_Distribution.png", 15, 7.5)

# ─── FIGURE 8 — DSI heatmap ──────────────────────────────────────────────────
dsi8_src <- dplyr::filter(df_vis, Disease_Name %in% TOP10, !is.na(Severity))
dsi8 <- dsi8_src %>%
  dplyr::group_by(Subzone_Display, Disease_Name) %>%
  dplyr::summarise(Mean_DSI = round(mean(Severity, na.rm = TRUE), 2), .groups = "drop") %>%
  tidyr::complete(Subzone_Display, Disease_Name = TOP10) %>%
  dplyr::mutate(Disease_Name    = factor(Disease_Name, levels = TOP10),
                Subzone_Display = factor(Subzone_Display, levels = rev(row_ord)))

p8 <- ggplot(dsi8, aes(x = Disease_Name, y = Subzone_Display, fill = Mean_DSI)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(Mean_DSI), "", sprintf("%.2f", Mean_DSI))),
            size = 3.5, colour = "grey10", na.rm = TRUE) +
  scale_fill_gradientn(colours = c("#FFFFFF", "#B7E1C3", "#FFEDA0", "#FC8D59", "#D73027", "#7B0000"),
                       na.value = "#F2F2F2", name = "Mean DSI (0–5)", limits = c(0.5, 5.5), breaks = 1:5,
                       guide = guide_colorbar(barwidth = 20, barheight = 0.8, title.position = "top")) +
  labs(title = "Mean Disease Severity Index (0–5) by Disease and Sub-region",
       x = "Disease", y = NULL,
       caption = sprintf("Grey = not recorded | DSI = mean raw ordinal score (0–5) | n=%d observations",
                         nrow(dsi8_src))) +
  theme_minimal(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 13, margin = margin(b = 10)),
        plot.caption = element_text(size = 9, colour = "grey45", hjust = 0),
        axis.text.x  = element_text(angle = 38, hjust = 1, size = 11.5),
        axis.text.y  = element_text(size = 11), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "bottom",
        plot.margin = margin(10, 12, 10, 10))
save_fig(p8, "Figure_8_DSI_Heatmap.png", 18, 8.5)

# ─── FIGURE 9 — Severity distribution ────────────────────────────────────────
sev6 <- named %>%
  dplyr::filter(Disease_Name %in% TOP6) %>%
  dplyr::mutate(Disease_Name = factor(Disease_Name, levels = TOP6))
sm6 <- sev6 %>%
  dplyr::group_by(Disease_Name, Disease_Type) %>%
  dplyr::summarise(Mean = round(mean(Severity, na.rm = TRUE), 2), n = dplyr::n(), .groups = "drop")

set.seed(42)
p9 <- ggplot(sev6, aes(x = Disease_Name, y = Severity, fill = Disease_Type)) +
  geom_boxplot(alpha = 0.60, width = 0.55, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(colour = Disease_Type), width = 0.13, alpha = 0.20, size = 1.6, shape = 16) +
  geom_point(data = sm6, aes(x = Disease_Name, y = Mean, colour = Disease_Type),
             size = 5, shape = 18, inherit.aes = FALSE) +
  geom_text(data = sm6, aes(x = Disease_Name, y = Mean + 0.14,
                            label = sprintf("μ=%.2f", Mean), colour = Disease_Type),
            size = 3.7, fontface = "bold", inherit.aes = FALSE, vjust = 0) +
  geom_text(data = sm6, aes(x = Disease_Name, y = 5.85, label = paste0("n=", n)),
            size = 3.2, colour = "grey45", inherit.aes = FALSE) +
  scale_fill_manual(values = DCOL, name = "Category") +
  scale_colour_manual(values = DCOL, guide = "none") +
  scale_y_continuous(breaks = 1:5,
                     labels = c("1\n(very mild)", "2\n(mild)", "3\n(moderate)", "4\n(severe)", "5\n(very severe)"),
                     limits = c(0.5, 6.5), expand = expansion(mult = c(0, 0.02))) +
  labs(title = sprintf("Severity Score Distribution for the %d Most Frequent Diseases", length(TOP6)),
       x = NULL, y = "Severity Score (0–5 ordinal scale)",
       caption = sprintf("Diamond = mean severity per disease | Kruskal-Wallis (severity by disease type) H=%.2f, %s",
                         kw_sev$statistic, fmt_p(kw_sev$p.value))) +
  PUB + theme(axis.text.x = element_text(angle = 16, hjust = 1, size = 11.5))
save_fig(p9, "Figure_9_Severity_Distribution.png", 14, 8)


# ==============================================================================
# SECTION 8: ADDITIONAL TABLES AND EXCEL WORKBOOK
# ==============================================================================
cat("\n=== SECTION 8: Additional tables + Excel workbook ===\n")

consolidated <- rbind(
  data.frame(Model = "Non-parametric",
             Test = nonpar_tab$Test, Statistic = as.character(nonpar_tab$Statistic),
             p_value = nonpar_tab$p_value, Sig = nonpar_tab$Sig,
             stringsAsFactors = FALSE),
  if (nrow(pom_all) > 0)
    data.frame(Model = "POM",
               Test = paste(pom_all$Disease, pom_all$Predictor, sep = " — "),
               Statistic = paste0("OR=", pom_all$OR),
               p_value = as.character(pom_all$p_value),
               Sig = pom_all$Sig, stringsAsFactors = FALSE),
  if (nrow(logit_all) > 0)
    data.frame(Model = "Logistic",
               Test = paste(logit_all$Disease, logit_all$Predictor, sep = " — "),
               Statistic = paste0("OR=", logit_all$OR),
               p_value = as.character(logit_all$p_value),
               Sig = logit_all$Sig, stringsAsFactors = FALSE)
)
save_tab(consolidated, "T12_Consolidated_All_Results.csv")

# Farm-level table: farm name and GPS coordinates are dropped from the export
# so that individual farms cannot be identified if outputs are shared.
save_tab(as.data.frame(dplyr::select(df_farm, -Name, -Longitude, -Latitude)),
         "T13_Farm_Level_Data.csv")

# Excel workbook with all sheets
wb <- createWorkbook()
tabs <- list(
  T01_DiseaseFreq   = dis_freq,
  T02_TypeTotals    = type_totals,
  T03_DSI_Disease   = dsi_disease,
  T04_Subregion     = sub_summary,
  T05_Environment   = env_summary,
  T06_NonParametric = nonpar_tab,
  T07_POM_Full      = pom_all,
  T08_POM_Sig       = pom_sig,
  T08b_POM_Fit      = pom_fit,
  T09_Logistic_Full = logit_all,
  T10_Logistic_Sig  = logit_sig,
  T10b_Logistic_Fit = logit_fit,
  T11_MoransI       = if (!is.null(morans_all)) morans_all else data.frame(),
  T12_Consolidated  = consolidated
)
hs <- createStyle(textDecoration = "bold", fgFill = "#1A3A5C",
                  fontColour = "#FFFFFF", border = "Bottom")
for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  if (nrow(tabs[[nm]]) > 0) {
    writeData(wb, nm, tabs[[nm]])
    addStyle(wb, nm, hs, rows = 1, cols = seq_len(ncol(tabs[[nm]])))
  }
}
saveWorkbook(wb, file.path(TAB_DIR, "All_Tables.xlsx"), overwrite = TRUE)
cat("  Excel workbook: All_Tables.xlsx\n")

writeLines(capture.output(sessionInfo()), file.path(TAB_DIR, "SessionInfo.txt"))

# ==============================================================================
cat(sprintf("\n=== ANALYSIS COMPLETE: %s ===\n", format(Sys.time())))
cat(sprintf("  Figures: %d PNG files written to the Figures folder\n",
            length(list.files(FIG_DIR, "\\.png$"))))
cat(sprintf("  Tables : %d CSV files + All_Tables.xlsx written to the Tables folder\n",
            length(list.files(TAB_DIR, "\\.csv$"))))
