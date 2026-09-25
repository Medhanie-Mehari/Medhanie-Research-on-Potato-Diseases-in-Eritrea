# ==============================================================================
# Objective 1 — Prevalence, spatial distribution and environmental correlates
#               of potato microbial diseases (field survey, June–August 2023)
# ==============================================================================
# Script  : 01_prevalence_analysis.R
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
# Tested  : R 4.4.1 (see SessionInfo.txt written at the end of the run)
#
# Output  : 8 PNG figures (400 dpi), 15 CSV tables, 1 Excel workbook
#
# Input data are NOT distributed with this script. See README for the
# expected column structure and the data availability statement.
#
# ------------------------------------------------------------------------------
# TRANSPARENCY NOTE — values transcribed from the manuscript
# ------------------------------------------------------------------------------
# The following objects contain values typed in from the published manuscript;
# they are NOT computed by this script:
#   * nonpar_tab$MS_reported         (Section 3)
#   * pom_fit_ms                      (Section 4; incl. Brant omnibus p-values)
#   * log_fit_ms                      (Section 5)
#   * POM_MS, LOG_MS  -> Figure 5     (Section 7)
#   * hist_v, pvx_d, cum -> Figure 6  (Section 7; compiled from literature:
#                                      <REFERENCES FOR HISTORICAL DATA>)
# The model estimates computed from the data are written to T07 and T09.
# <EXPLAIN HERE WHY FIGURE 5 USES MANUSCRIPT VALUES RATHER THAN T07/T09>
#
# Implementation notes:
#   1. MASS is loaded before dplyr so that dplyr::select() is not masked.
#   2. Boolean columns exported as "True"/"False" text are converted to logical.
#   3. All dplyr verbs are namespaced (dplyr::) throughout.
#   4. Label offsets for Figure 3 are pre-computed (dsi_offset) outside aes().
#   5. Moran's I uses Longitude/Latitude vectors directly.
#   6. POMs are fitted for any disease with n >= 30 observations.
# ==============================================================================


# ── 0. CONFIGURATION — replace the placeholders before running ───────────────

DATA_FILE <- "<PATH_TO_INPUT_DATA>/<INPUT_DATA_FILE>.csv"
OUT_DIR   <- "<PATH_TO_OUTPUT_FOLDER>"

if (!file.exists(DATA_FILE))
  stop("Input data not found. Set DATA_FILE in Section 0 to your data file.")

FIG_DIR <- file.path(OUT_DIR, "Figures")
TAB_DIR <- file.path(OUT_DIR, "Tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)


# ── 0. PACKAGES ──────────────────────────────────────────────────────────────

need <- c("MASS", "car", "pscl", "DescTools", "brant", "spdep",
          "dplyr", "tidyr", "ggplot2", "patchwork", "scales",
          "RColorBrewer", "openxlsx", "dunn.test")
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
  library(car)
  library(pscl)
  library(DescTools)
  library(brant)
  library(spdep)
  library(dunn.test)
  library(openxlsx)
  library(RColorBrewer)
})

cat("=== Objective 1 analysis started:", format(Sys.time()), "===\n")


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


# ── Colour palettes ──────────────────────────────────────────────────────────

CF <- "#2D8653"; CV <- "#B5342A"; CB <- "#2166AC"
CA <- "#E07B39"; CG <- "#7B4F9E"

DCOL <- c(Fungal = CF, Viral = CV, Bacterial = CB)
RCOL <- c("Maekel (Central)" = CB, "Debub (Southern)" = CF,
          Anseba = CA, "Gash Barka" = CG)


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
named  <- dplyr::filter(raw, Is_Named)           # 703 named records
df_vis <- dplyr::filter(raw, Is_Specific_Named)  # 631, excl. VMI & AMV
symp   <- dplyr::filter(raw, !Is_Other)          # 721, excl. "Other"

# Farm-level dataset (79 rows; one per farm)
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
stopifnot(nrow(df_farm) == 79)

cat(sprintf("  df=%d | named=%d | df_vis=%d | df_farm=%d\n",
            nrow(df), nrow(named), nrow(df_vis), nrow(df_farm)))


# ==============================================================================
# SECTION 2: DESCRIPTIVE STATISTICS
# ==============================================================================
cat("\n=== SECTION 2: Descriptive statistics ===\n")

# 703 = named disease records; 4648 = total plants assessed
dis_freq <- named %>%
  dplyr::count(Disease_Name, Disease_Type) %>%
  dplyr::mutate(Pct_of_703  = round(n / 703 * 100, 1),
                Pct_of_4648 = round(n / 4648 * 100, 2)) %>%
  dplyr::arrange(dplyr::desc(n))

type_totals <- named %>%
  dplyr::group_by(Disease_Type) %>%
  dplyr::summarise(N_diseases = dplyr::n_distinct(Disease_Name),
                   N_records  = dplyr::n(),
                   Pct_of_703 = round(dplyr::n() / 703 * 100, 1),
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
    Pct_of_4648      = round(dplyr::first(Subzone_N_Total) / 4648 * 100, 1),
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
save_tab(sub_summary, "T04_Subregion_Summary_Table3.csv")
save_tab(env_summary, "T05_Environmental_Summary.csv")


# ==============================================================================
# SECTION 3: NON-PARAMETRIC TESTS
# ==============================================================================
cat("\n=== SECTION 3: Non-parametric tests ===\n")

kw_sev <- kruskal.test(Severity ~ Disease_Type, data = named)
cat(sprintf("  KW Severity x Type: H=%.3f, p=%.4g\n",
            kw_sev$statistic, kw_sev$p.value))
suppressMessages(dunn.test::dunn.test(named$Severity, named$Disease_Type,
                                      method = "bonferroni", kw = FALSE, list = FALSE))

kw_alt <- kruskal.test(Altitude ~ Disease_Type, data = named)
cat(sprintf("  KW Altitude x Type: H=%.3f, p=%.4f\n",
            kw_alt$statistic, kw_alt$p.value))

vw_alt  <- named$Altitude[named$Disease_Name == "Verticillium Wilt"]
oth_alt <- named$Altitude[named$Disease_Name != "Verticillium Wilt"]
mw_vw   <- wilcox.test(vw_alt, oth_alt, alternative = "less", exact = FALSE)
cat(sprintf("  MW Verticillium < others: U=%.0f, p=%.4g\n",
            mw_vw$statistic, mw_vw$p.value))
cat(sprintf("  VW mean alt=%.1f m | Others mean alt=%.1f m\n",
            mean(vw_alt, na.rm = TRUE), mean(oth_alt, na.rm = TRUE)))

sp_sub <- named %>%
  dplyr::filter(Disease_Name == "Verticillium Wilt") %>%
  dplyr::group_by(Subzone) %>%
  dplyr::summarise(N_VW     = dplyr::n(),
                   Mean_Alt = mean(Altitude, na.rm = TRUE), .groups = "drop")
sp_vw <- cor.test(sp_sub$Mean_Alt, sp_sub$N_VW, method = "spearman", exact = FALSE)
cat(sprintf("  Spearman rho=%.3f, p=%.4f, n=%d sub-regions\n",
            sp_vw$estimate, sp_vw$p.value, nrow(sp_sub)))

sub_corr <- raw %>%
  dplyr::group_by(Subzone) %>%
  dplyr::summarise(W_Inc = dplyr::first(Subzone_Incidence_pct),
                   DSI   = dplyr::first(Subzone_Mean_DSI), .groups = "drop")
cor_id <- cor.test(sub_corr$W_Inc, sub_corr$DSI, method = "pearson")
cat(sprintf("  Pearson r=%.3f, p=%.4f (n=14)\n", cor_id$estimate, cor_id$p.value))

nonpar_tab <- data.frame(
  Test = c("KW Severity x Type", "KW Altitude x Type",
           "MW Verticillium alt < others",
           "Spearman rho (altitude vs VW freq)",
           "Pearson r (weighted incidence vs DSI)"),
  Statistic = c(round(kw_sev$statistic, 3), round(kw_alt$statistic, 3),
                round(mw_vw$statistic, 0), round(sp_vw$estimate, 3),
                round(cor_id$estimate, 3)),
  p_value = formatC(c(kw_sev$p.value, kw_alt$p.value, mw_vw$p.value,
                      sp_vw$p.value, cor_id$p.value), digits = 4, format = "g"),
  MS_reported = c("H=149.21", "H=10.619", "U=11270", "rho=-0.769", "r=-0.005"),  # from manuscript
  Sig = c("Yes", "Yes", "Yes", "Yes", "No"),
  stringsAsFactors = FALSE)
save_tab(nonpar_tab, "T06_NonParametric_Tests.csv")


# ==============================================================================
# SECTION 4: PROPORTIONAL ODDS MODEL — DISEASE SEVERITY
# ==============================================================================
cat("\n=== SECTION 4: POM ===\n")

run_pom <- function(disease_name) {
  sub <- dplyr::filter(named, Disease_Name == disease_name) %>%
    dplyr::mutate(Sev_f = factor(Severity, ordered = TRUE))
  n <- nrow(sub)
  cat(sprintf("\n  POM: %s (n=%d)\n", disease_name, n))
  if (n < 30) { cat("  SKIP: n < 30\n"); return(NULL) }

  mod <- tryCatch(
    MASS::polr(Sev_f ~ Temperature_C + RH_pct + Soil_pH +
                 Altitude + PlantAge + Density_Numeric,
               data = sub, Hess = TRUE, method = "logistic"),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL })
  if (is.null(mod)) return(NULL)

  ctab <- coef(summary(mod))
  beta <- ctab[, "Value"]; se <- ctab[, "Std. Error"]
  pv   <- pnorm(abs(beta / se), lower.tail = FALSE) * 2
  ci   <- tryCatch(exp(confint(mod)),
                   error = function(e) exp(cbind(beta - 1.96 * se, beta + 1.96 * se)))
  OR   <- exp(beta)
  idx  <- !grepl("\\|", rownames(ctab))

  null_ll <- as.numeric(logLik(MASS::polr(Sev_f ~ 1, data = sub, Hess = TRUE)))
  full_ll <- as.numeric(logLik(mod))
  mcfR2   <- round(1 - full_ll / null_ll, 3)
  cat(sprintf("  AIC=%.2f | McFadden R2=%.3f\n", AIC(mod), mcfR2))

  tryCatch({ br <- brant::brant(mod)
             cat("  Brant test:\n"); print(br) }, error = function(e) NULL)

  list(model = mod,
       results = data.frame(
         Disease     = disease_name, n = n,
         Predictor   = rownames(ctab)[idx],
         OR          = round(OR[idx], 3),
         CI_L        = round(ci[idx, 1], 3),
         CI_H        = round(ci[idx, 2], 3),
         p_value     = round(pv[idx], 4),
         Sig         = ifelse(pv[idx] < 0.05, "Yes", "No"),
         AIC         = round(AIC(mod), 2),
         McFadden_R2 = mcfR2,
         stringsAsFactors = FALSE),
       aic = round(AIC(mod), 2), mcfR2 = mcfR2)
}

pom_pvy  <- run_pom("PVY")
pom_eb   <- run_pom("Early Blight")
pom_vw   <- run_pom("Verticillium Wilt")
pom_plrv <- run_pom("PLRV")

pom_list <- Filter(Negate(is.null), list(pom_pvy, pom_eb, pom_vw, pom_plrv))
pom_all  <- if (length(pom_list) > 0)
  do.call(rbind, lapply(pom_list, function(x) x$results)) else data.frame()

cat("\n  POM significant results:\n")
print(dplyr::filter(pom_all, Sig == "Yes")[, c("Disease", "Predictor", "OR", "CI_L", "CI_H", "p_value")],
      row.names = FALSE)

# Values transcribed from the manuscript (not computed here)
pom_fit_ms <- data.frame(
  Disease         = c("PVY", "Early Blight", "Verticillium Wilt", "PLRV"),
  n               = c(99L, 147L, 66L, 105L),
  AIC             = c(310.19, 364.52, 87.59, 330.47),
  McFadden_R2     = c(0.11, 0.09, 0.07, 0.02),
  Brant_p_omnibus = c(0.312, 0.284, 0.401, 0.342),
  MS_significant  = c("Temp OR=0.262 p=0.042; PlantAge OR=0.961 p=0.035",
                      "Temp OR=0.360 p=0.013; PlantAge OR=1.033 p=0.020",
                      "PlantAge OR=1.106 p=0.015",
                      "none (all p>0.05)"),
  stringsAsFactors = FALSE)

save_tab(pom_all,                              "T07_POM_Results_Full.csv")
save_tab(dplyr::filter(pom_all, Sig == "Yes"), "T08_POM_Results_Significant.csv")
save_tab(pom_fit_ms,                           "T08b_POM_Fit_Manuscript.csv")


# ==============================================================================
# SECTION 5: BINOMIAL LOGISTIC REGRESSION — DISEASE INCIDENCE
# ==============================================================================
cat("\n=== SECTION 5: Logistic regression (n=79 farms) ===\n")

run_logit <- function(outcome, label) {
  fml <- as.formula(paste(outcome,
    "~ Temperature_C + RH_pct + Soil_pH + Altitude + PlantAge + Density_Numeric"))
  mod <- tryCatch(
    glm(fml, data = df_farm, family = binomial(link = "logit")),
    error = function(e) { cat("  ERROR:", label, "-", e$message, "\n"); NULL })
  if (is.null(mod)) return(NULL)

  cf  <- coef(summary(mod))[-1, , drop = FALSE]
  cis <- confint.default(mod)[-1, , drop = FALSE]
  nag <- tryCatch(round(DescTools::PseudoR2(mod, which = "Nagelkerke"), 3),
                  error = function(e) NA)
  cat(sprintf("  %s: AIC=%.2f | Nagelkerke R2=%.3f\n", label, AIC(mod), nag))

  list(model = mod,
       results = data.frame(
         Disease       = label, n_farms = 79L,
         Predictor     = rownames(cf),
         OR            = round(exp(cf[, "Estimate"]), 3),
         CI_L          = round(exp(cis[, 1]), 3),
         CI_H          = round(exp(cis[, 2]), 3),
         p_value       = round(cf[, "Pr(>|z|)"], 4),
         Sig           = ifelse(cf[, "Pr(>|z|)"] < 0.05, "Yes", "No"),
         Nagelkerke_R2 = nag,
         stringsAsFactors = FALSE))
}

log_eb   <- run_logit("Inc_EarlyBlight",  "Early Blight")
log_pvy  <- run_logit("Inc_PVY",          "PVY")
log_vw   <- run_logit("Inc_Verticillium", "Verticillium Wilt")
log_plrv <- run_logit("Inc_PLRV",         "PLRV")
log_vmi  <- run_logit("Inc_VMI",          "Viral Mixed Infection")

logit_all <- do.call(rbind, lapply(
  Filter(Negate(is.null), list(log_eb, log_pvy, log_vw, log_plrv, log_vmi)),
  function(x) x$results))

cat("\n  Logistic significant results:\n")
print(dplyr::filter(logit_all, Sig == "Yes")[, c("Disease", "Predictor", "OR", "CI_L", "CI_H", "p_value")],
      row.names = FALSE)

# Values transcribed from the manuscript (not computed here)
log_fit_ms <- data.frame(
  Disease        = c("Early Blight", "PVY", "Verticillium Wilt", "Viral Mixed Infection"),
  Nagelkerke_R2  = c(0.18, 0.29, 0.41, 0.12),
  MS_significant = c(
    "RH OR=0.91 p=0.008",
    "RH OR=1.11 p=0.048; PlantAge OR=0.977 p=0.009; Density(High) OR=0.259 p<0.001",
    "Temp OR=16.73 p=0.002; RH OR=0.820 p=0.013; PlantAge OR=0.956 p<0.001",
    "Altitude OR=0.997 p=0.025"),
  stringsAsFactors = FALSE)

save_tab(logit_all,                              "T09_Logistic_Results_Full.csv")
save_tab(dplyr::filter(logit_all, Sig == "Yes"), "T10_Logistic_Results_Significant.csv")
save_tab(log_fit_ms,                             "T10b_Logistic_Fit_Manuscript.csv")


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

pvy_dat <- dplyr::filter(named, Disease_Name == "PVY")
eb_dat  <- dplyr::filter(named, Disease_Name == "Early Blight")
vw_dat  <- dplyr::filter(named, Disease_Name == "Verticillium Wilt")

mi_rows <- list(
  if (!is.null(pom_pvy)) run_mi(pom_pvy$model, "POM PVY",
                                pvy_dat$Longitude, pvy_dat$Latitude),
  if (!is.null(pom_eb))  run_mi(pom_eb$model,  "POM Early Blight",
                                eb_dat$Longitude,  eb_dat$Latitude),
  if (!is.null(pom_vw))  run_mi(pom_vw$model,  "POM Verticillium Wilt",
                                vw_dat$Longitude,  vw_dat$Latitude),
  if (!is.null(log_eb))  run_mi(log_eb$model,  "Logistic Early Blight",
                                df_farm$Longitude, df_farm$Latitude),
  if (!is.null(log_pvy)) run_mi(log_pvy$model, "Logistic PVY",
                                df_farm$Longitude, df_farm$Latitude),
  if (!is.null(log_vw))  run_mi(log_vw$model,  "Logistic Verticillium Wilt",
                                df_farm$Longitude, df_farm$Latitude),
  if (!is.null(log_vmi)) run_mi(log_vmi$model, "Logistic VMI",
                                df_farm$Longitude, df_farm$Latitude)
)
morans_all <- do.call(rbind, Filter(Negate(is.null), mi_rows))
if (!is.null(morans_all) && nrow(morans_all) > 0)
  save_tab(morans_all, "T11_Morans_I.csv")


# ==============================================================================
# SECTION 7: FIGURES
# ==============================================================================
cat("\n=== SECTION 7: Figures ===\n")

TOP10 <- c("Early Blight", "PLRV", "PVY", "Verticillium Wilt", "Late Blight",
           "Rhizoctonia Canker", "Black Leg", "Grey Mold", "Brown Spot", "White Mold")

# ─── FIGURE 2 — Disease frequency ────────────────────────────────────────────
f2_dat <- named %>%
  dplyr::filter(Disease_Name %in% TOP10) %>%
  dplyr::count(Disease_Name, Disease_Type) %>%
  dplyr::mutate(
    Disease_Name = factor(Disease_Name, levels = rev(TOP10)),
    pct   = round(n / 703 * 100, 1),
    label = paste0(n, " (", pct, "%)"))

p2 <- ggplot(f2_dat, aes(x = n, y = Disease_Name, fill = Disease_Type)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.8,
            fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = DCOL, name = "Category") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.33)), breaks = seq(0, 150, 25)) +
  labs(title = "Frequency of the Top 10 Potato Diseases by Plant-Level Observation Count",
       x = "Number of Plant Observations", y = NULL,
       caption = "★ Brown Spot = first confirmed record of Alternaria alternata on potato in Eritrea | n = 703 named records") +
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
    Subzone_Display = factor(Subzone_Display, levels = Subzone_Display),
    dsi_offset      = ifelse(dplyr::row_number() %% 2 == 0, 1.5, -1.5))

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
    name = "Weighted Disease Incidence (%)", limits = c(0, 38),
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
  dplyr::pull(Subzone_Display) %>% as.character()

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
    na.value = "#EFEFEF", name = "Proportion (%)", limits = c(0, 42),
    guide = guide_colorbar(barwidth = 18, barheight = 0.8, title.position = "top")) +
  labs(title = "Proportion (%) of Top 10 Disease Observations per Sub-region",
       x = "Disease", y = NULL,
       caption = "Grey = not recorded | VMI & AMV excluded | n = 631 specific named records") +
  theme_minimal(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 13, margin = margin(b = 10)),
        plot.caption = element_text(size = 9, colour = "grey45", hjust = 0),
        axis.text.x  = element_text(angle = 38, hjust = 1, size = 11.5),
        axis.text.y  = element_text(size = 11),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", plot.margin = margin(10, 12, 10, 10))
save_fig(p4, "Figure_4_Disease_Heatmap.png", 18, 8.5)

# ─── FIGURE 5 — Forest plot ──────────────────────────────────────────────────
# NOTE: POM_MS and LOG_MS are transcribed from the manuscript; they are not
# taken from pom_all / logit_all computed above. See transparency note in header.
POM_MS <- data.frame(
  Disease   = c("PVY", "PVY", "Early Blight", "Early Blight", "Verticillium Wilt"),
  Predictor = c("Temperature (per 1°C)", "Plant Age (per day)",
                "Temperature (per 1°C)", "Plant Age (per day)",
                "Plant Age (per day)"),
  OR      = c(0.262, 0.961, 0.360, 1.033, 1.106),
  CI_L    = c(0.072, 0.927, 0.160, 1.005, 1.019),
  CI_H    = c(0.950, 0.997, 0.807, 1.061, 1.199),
  p_value = c(0.042, 0.035, 0.013, 0.020, 0.015),
  stringsAsFactors = FALSE)
POM_MS$Disease <- factor(POM_MS$Disease,
                         levels = c("PVY", "Early Blight", "Verticillium Wilt"))
POM_MS$pcol <- ifelse(POM_MS$OR < 1, "#2166AC", "#B5342A")
POM_MS$plab <- sprintf("p=%.3f *", POM_MS$p_value)

pA <- ggplot(POM_MS, aes(x = OR, y = reorder(Predictor, -OR), xmin = CI_L, xmax = CI_H)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.8) +
  geom_errorbarh(aes(colour = pcol), height = 0.22, linewidth = 1.8, alpha = 0.9) +
  geom_point(aes(colour = pcol), size = 5, shape = 23, fill = "white", stroke = 2) +
  geom_text(aes(x = CI_H, label = paste0(" ", plab)), hjust = 0, size = 3.5, colour = "grey25") +
  scale_colour_identity() +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 1.1, 1.5),
                labels = c("0.10", "0.25", "0.50", "1", "1.10", "1.50"),
                expand = expansion(mult = c(0.05, 0.55))) +
  facet_wrap(~Disease, scales = "free_y", ncol = 3,
             labeller = labeller(Disease = c(
               "PVY"               = "PVY (n=99; AIC=310.19; McFadden R²=0.11)",
               "Early Blight"      = "Early Blight (n=147; AIC=364.52; McFadden R²=0.09)",
               "Verticillium Wilt" = "Verticillium Wilt (n=66; AIC=87.59; McFadden R²=0.07)"))) +
  labs(title = "(A) Disease Severity — Proportional Odds Model | Plant-level observations",
       x = "Odds Ratio (log scale) with 95% profile-likelihood CI", y = NULL,
       caption = "Blue = OR<1 p<0.05 (protective) | Red = OR>1 p<0.05 (positive) | PLRV: all p>0.05") +
  PUB + theme(panel.grid.major.y = element_blank(), strip.text = element_text(size = 9.5))

LOG_MS <- data.frame(
  Disease   = c("Early Blight", "PVY", "PVY", "PVY",
                "Verticillium Wilt", "Verticillium Wilt", "Verticillium Wilt",
                "Viral Mixed Infection"),
  Predictor = c("Relative Humidity\n(per 1%)", "Relative Humidity\n(per 1%)",
                "Plant Age\n(per day)", "Density High\n(vs Low)",
                "Temperature\n(per 1°C)", "Relative Humidity\n(per 1%)",
                "Plant Age\n(per day)", "Altitude\n(per 1 m)"),
  OR      = c(0.91, 1.11, 0.977, 0.259, 16.73, 0.820, 0.956, 0.997),
  CI_L    = c(0.85, 1.00, 0.962, 0.117, 2.930, 0.700, 0.926, 0.994),
  CI_H    = c(0.97, 1.22, 0.993, 0.574, 95.40, 0.961, 0.987, 0.999),
  p_value = c(0.008, 0.048, 0.009, 0.001, 0.002, 0.013, 0.001, 0.025),
  stringsAsFactors = FALSE)
LOG_MS$Disease <- factor(LOG_MS$Disease,
                         levels = c("Early Blight", "PVY", "Verticillium Wilt", "Viral Mixed Infection"))
LOG_MS$pcol  <- ifelse(LOG_MS$OR < 1, "#2166AC", "#B5342A")
LOG_MS$stars <- ifelse(LOG_MS$p_value < 0.001, "***", ifelse(LOG_MS$p_value < 0.01, "**", "*"))
LOG_MS$plab  <- ifelse(LOG_MS$p_value < 0.001, "p<0.001", sprintf("p=%.3f", LOG_MS$p_value))

pB <- ggplot(LOG_MS, aes(x = OR, y = reorder(Predictor, -OR), xmin = CI_L, xmax = CI_H)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.8) +
  geom_errorbarh(aes(colour = pcol), height = 0.22, linewidth = 1.8, alpha = 0.9) +
  geom_point(aes(colour = pcol), size = 5, shape = 23, fill = "white", stroke = 2) +
  geom_text(aes(x = pmax(CI_H, OR * 1.05), label = paste0(" ", plab, stars)),
            hjust = 0, size = 3.4, colour = "grey25") +
  scale_colour_identity() +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.60))) +
  facet_wrap(~Disease, scales = "free_y", ncol = 4,
             labeller = labeller(Disease = c(
               "Early Blight"          = "Early Blight\n(Nag.R²=0.18)",
               "PVY"                   = "PVY\n(Nag.R²=0.29)",
               "Verticillium Wilt"     = "Verticillium Wilt\n(Nag.R²=0.41)",
               "Viral Mixed Infection" = "Viral Mixed Infection\n(Nag.R²=0.12)"))) +
  labs(title = "(B) Disease Incidence — Binomial Logistic Regression | Farm-level binary (n=79 farms)",
       x = "Odds Ratio (log scale) with 95% Wald CI", y = NULL,
       caption = "OR = exp(β); β = log-odds change per unit — NOT a probability change (Hosmer et al., 2013)") +
  PUB + theme(panel.grid.major.y = element_blank(), strip.text = element_text(size = 9.5))

p5 <- pA / pB + plot_layout(heights = c(1, 0.9))
save_fig(p5, "Figure_5_ForestPlot_POM_Logistic.png", 19, 14)

# ─── FIGURE 6 — Historical trends (values compiled from literature) ──────────
# Source(s): <REFERENCES FOR HISTORICAL DATA>
hist_v <- data.frame(
  Year = c(2003, 2005, 2007, 2008, 2010, 2012, 2015, 2018, 2019, 2020, 2021, 2022),
  PLRV = c(21.4, 35, 42, 50, 60, 68, 78, 88, 92, 96, 100, 100),
  PVY  = c(7.5, 15, 20, 28, 40, 55, 60, 65, 70, 75, 80, 82.3))
pvx_d <- data.frame(Year = c(2007, 2010, 2012, 2015, 2016, 2018),
                    PVX  = c(12.5, 25, 40, 55, 66.7, 45))

hv <- tidyr::pivot_longer(hist_v, -Year, names_to = "Virus", values_to = "Inc")
pv <- tidyr::pivot_longer(pvx_d,  -Year, names_to = "Virus", values_to = "Inc")
av <- rbind(hv, pv); av$Virus <- factor(av$Virus, levels = c("PLRV", "PVY", "PVX"))
VCOL <- c(PLRV = CV, PVY = CG, PVX = CB)

p6a <- ggplot(av, aes(Year, Inc, colour = Virus, shape = Virus)) +
  geom_line(linewidth = 2.0) + geom_point(size = 3.5, stroke = 1.5, fill = "white") +
  annotate("text", x = 2020.5, y = 103, label = "100%", colour = CV, fontface = "bold", size = 4.2) +
  scale_colour_manual(values = VCOL) +
  scale_shape_manual(values = c(PLRV = 21, PVY = 22, PVX = 24)) +
  scale_x_continuous(breaks = seq(2003, 2022, 3)) +
  coord_cartesian(ylim = c(0, 115)) +
  labs(title = "(A) Historical Viral Disease Incidence (1995–2022)", x = "Year", y = "Incidence (%)") + PUB

cum <- data.frame(Year      = c(1995, 2000, 2003, 2007, 2010, 2015, 2018, 2020, 2022),
                  Fungal    = c(2, 3, 3, 4, 5, 6, 7, 9, 10),
                  Viral     = c(0, 0, 2, 4, 5, 5, 7, 10, 13),
                  Bacterial = c(0, 0, 0, 0, 0, 0, 1, 2, 3))
cl <- tidyr::pivot_longer(cum, -Year, names_to = "Type", values_to = "n")
cl$Type <- factor(cl$Type, levels = c("Fungal", "Viral", "Bacterial"))

p6b <- ggplot(cl, aes(Year, n, colour = Type, shape = Type)) +
  geom_line(linewidth = 2.0) + geom_point(size = 3.5, stroke = 1.5, fill = "white") +
  scale_colour_manual(values = DCOL) + scale_shape_manual(values = c(Fungal = 21, Viral = 22, Bacterial = 24)) +
  scale_x_continuous(breaks = seq(1995, 2022, 5)) +
  labs(title = "(B) Cumulative Disease Taxa Documented", x = "Year", y = "Cumulative taxa (n)") + PUB

p6 <- p6a + p6b + plot_layout(ncol = 2) +
  plot_annotation(title = "Historical Trends in Potato Disease Documentation in Eritrea (1995–2022)",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6))))
save_fig(p6, "Figure_6_Historical_Trends.png", 17, 7)

# ─── FIGURE 7 — Altitude distribution ────────────────────────────────────────
TOP8 <- c("Early Blight", "PVY", "PLRV", "Verticillium Wilt",
          "Late Blight", "Rhizoctonia Canker", "Black Leg", "Grey Mold")

alt8 <- named %>%
  dplyr::filter(Disease_Name %in% TOP8) %>%
  dplyr::mutate(Disease_Name = factor(Disease_Name,
    levels = names(sort(tapply(Altitude, Disease_Name, median, na.rm = TRUE), decreasing = TRUE))))
vw_mean <- mean(named$Altitude[named$Disease_Name == "Verticillium Wilt"], na.rm = TRUE)

set.seed(42)
p7 <- ggplot(alt8, aes(x = Disease_Name, y = Altitude, fill = Disease_Type)) +
  geom_boxplot(alpha = 0.60, width = 0.55, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(aes(colour = Disease_Type), width = 0.14, alpha = 0.20, size = 1.6, shape = 16) +
  geom_hline(yintercept = vw_mean, colour = "#B5342A", linetype = "dashed", linewidth = 1.5, alpha = 0.85) +
  scale_fill_manual(values = DCOL, name = "Category") +
  scale_colour_manual(values = DCOL, guide = "none") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "Altitudinal Distribution of Disease Occurrences Across the Top Eight Most Prevalent Diseases",
       x = NULL, y = "Altitude (m a.s.l.)",
       caption = sprintf("Dashed red line = mean altitude of Verticillium Wilt observations (%.0f m a.s.l.)", vw_mean)) +
  PUB + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 11.5))
save_fig(p7, "Figure_7_Altitude_Distribution.png", 15, 7.5)

# ─── FIGURE 8 — DSI heatmap ──────────────────────────────────────────────────
DIS8 <- c("PLRV", "PVY", "Early Blight", "Late Blight", "Verticillium Wilt",
          "Rhizoctonia Canker", "Black Leg", "Brown Spot", "Grey Mold", "White Mold")

dsi8 <- df_vis %>%
  dplyr::filter(Disease_Name %in% DIS8) %>%
  dplyr::group_by(Subzone_Display, Disease_Name) %>%
  dplyr::summarise(Mean_DSI = round(mean(Severity, na.rm = TRUE), 2), .groups = "drop") %>%
  tidyr::complete(Subzone_Display, Disease_Name = DIS8) %>%
  dplyr::mutate(Disease_Name    = factor(Disease_Name, levels = DIS8),
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
       caption = "Grey = not recorded | DSI = mean raw ordinal score (0–5) | n=740 observations") +
  theme_minimal(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 13, margin = margin(b = 10)),
        plot.caption = element_text(size = 9, colour = "grey45", hjust = 0),
        axis.text.x  = element_text(angle = 38, hjust = 1, size = 11.5),
        axis.text.y  = element_text(size = 11), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "bottom",
        plot.margin = margin(10, 12, 10, 10))
save_fig(p8, "Figure_8_DSI_Heatmap.png", 18, 8.5)

# ─── FIGURE 9 — Severity distribution ────────────────────────────────────────
TOP6 <- c("PLRV", "PVY", "Early Blight", "Late Blight", "Verticillium Wilt", "Rhizoctonia Canker")

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
  labs(title = "Severity Score Distribution for Six Major Potato Diseases", x = NULL,
       y = "Severity Score (0–5 ordinal scale)",
       caption = "Diamond = mean severity per disease | Kruskal-Wallis H=149.21, p<0.001") +
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
  data.frame(Model = "POM",
             Test = paste(pom_all$Disease, pom_all$Predictor, sep = " — "),
             Statistic = paste0("OR=", pom_all$OR),
             p_value = as.character(pom_all$p_value),
             Sig = pom_all$Sig, stringsAsFactors = FALSE),
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
save_tab(POM_MS, "T14_POM_Manuscript_OR_Values.csv")
save_tab(LOG_MS, "T15_Logistic_Manuscript_OR_Values.csv")

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
  T08_POM_Sig       = dplyr::filter(pom_all, Sig == "Yes"),
  T08b_POM_Fit_MS   = pom_fit_ms,
  T09_Logistic_Full = logit_all,
  T10_Logistic_Sig  = dplyr::filter(logit_all, Sig == "Yes"),
  T10b_Logistic_MS  = log_fit_ms,
  T11_MoransI       = if (!is.null(morans_all)) morans_all else data.frame(),
  T12_Consolidated  = consolidated,
  T14_POM_MS_Values = POM_MS,
  T15_Log_MS_Values = LOG_MS
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
