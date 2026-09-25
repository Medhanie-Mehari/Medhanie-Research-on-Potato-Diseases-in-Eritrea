# =============================================================================
# POTATO VARIETY EVALUATION FOR FOLIAR DISEASE RESISTANCE, TOLERANCE AND YIELD
# Multi-season RCBD variety trial
#
# COMPLETE END-TO-END ANALYSIS SCRIPT
#
# Script  : 07_variety_evaluation.R
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
# Outputs : tables/ (CSV + one Excel workbook), figures/ (300 dpi PNG),
#           sessionInfo.txt
#
# ---------------------------------------------------------------------------
# INPUT FILES (set paths in Section 0)
# ---------------------------------------------------------------------------
# PLOT_FILE           one row per plot: Season, Block, Variety, Plant_Stand,
#                     Yield_t_ha, AUDPC_TOT, sAUDPC_TOT, sAUDPC_<CODE> per disease
# TIMESERIES_FILE     one row per plot x assessment: Season, Block, Variety,
#                     DAP, Total_Severity
# VARIETY_MEANS_FILE  one row per variety: Variety, sAUDPC_<CODE>, sAUDPC_TOT,
#                     DRI_percent, Peak_Total_Severity, Mean_Total_Severity,
#                     Yield_t_ha, Yield_adj_t_ha, Disease_Reaction, Plant_Stand
# WEATHER_FILE        daily weather: Date (yyyy-mm-dd), Season_window (season
#                     name, empty outside seasons), Temp_C_daily_avg,
#                     Humidity_pct_daily_avg, Rainfall_mm
# VARIETY_INFO_FILE   optional: any descriptive columns per variety (e.g.
#                     Variety, Source, Breeder, Origin); copied to table T01
#
# Design: randomised complete block design repeated over seasons; blocks are
# nested within seasons. Design sizes (seasons, blocks, varieties, error df)
# are taken from the data.
#
# Note: bootstrap confidence intervals use seed = 1. Parametric bootstrapping
# of near-zero variance components is unstable, so upper limits may shift
# slightly between machines or package versions.
# =============================================================================

## ---- 0. SETTINGS (edit this section) ---------------------------------------

PLOT_FILE          <- "<PATH_TO_PLOT_LEVEL_DATA>.csv"
TIMESERIES_FILE    <- "<PATH_TO_PLOT_TIMESERIES_DATA>.csv"
VARIETY_MEANS_FILE <- "<PATH_TO_VARIETY_MEANS_DATA>.csv"
WEATHER_FILE       <- "<PATH_TO_DAILY_WEATHER_DATA>.csv"
VARIETY_INFO_FILE  <- NULL                         # optional, e.g. "<PATH>.csv"
OUTPUT_DIR         <- "<PATH_TO_OUTPUT_FOLDER>"

# Diseases: names = column suffix (sAUDPC_<CODE>), values = labels
DISEASES     <- c(EB = "Early blight", LB = "Late blight", BS = "Brown spot")
SEVERITY_MAX <- 5          # per-disease severity scale; total = SEVERITY_MAX x diseases

# Planting date of each season (names must match Season in the data), used to
# place weather and disease on a common days-after-planting axis. Leave the
# placeholders to skip the weather figure.
SEASON_START_DATES <- c("<SEASON_1_NAME>" = "<YYYY-MM-DD>",
                        "<SEASON_2_NAME>" = "<YYYY-MM-DD>")

# Variety classification thresholds
YIELD_HIGH   <- 11         # t/ha: >= high yield
YIELD_MEDIUM <- 7          # t/ha: >= medium yield
DRI_GOOD     <- 25         # %: disease resistance index for "good" resistance

ALPHA      <- 0.05         # significance level (Tukey HSD)
POWER      <- 0.80         # power for the minimum detectable difference
N_BOOT     <- 500          # parametric bootstrap samples for heritability CIs
N_CLUSTERS <- 3            # clusters outlined on the UPGMA dendrogram


## ---- 0b. SETUP -------------------------------------------------------------

for (v in c("PLOT_FILE", "TIMESERIES_FILE", "VARIETY_MEANS_FILE", "WEATHER_FILE", "OUTPUT_DIR"))
  if (grepl("<", get(v), fixed = TRUE)) stop("Set ", v, " in Section 0 before running.")
for (f in c(PLOT_FILE, TIMESERIES_FILE, VARIETY_MEANS_FILE, WEATHER_FILE, VARIETY_INFO_FILE))
  if (!file.exists(f)) stop("Input file not found: ", f)

fig_dir <- file.path(OUTPUT_DIR, "figures")
tab_dir <- file.path(OUTPUT_DIR, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

need <- c("ggplot2", "dplyr", "tidyr", "lme4", "lmerTest", "car", "emmeans",
          "openxlsx", "RColorBrewer", "scales", "multcompView")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(lme4); library(lmerTest)
  library(car); library(emmeans); library(openxlsx); library(RColorBrewer); library(scales)
})
set.seed(1)   # reproducibility (bootstrap)
options(stringsAsFactors = FALSE)

# publication theme -----------------------------------------------------------
theme_pub <- function(base = 14) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
          axis.title  = element_text(face = "bold", size = base + 2),
          axis.text   = element_text(colour = "grey15", size = base),
          plot.title  = element_text(face = "bold", size = base + 4, hjust = 0.5),
          legend.title = element_text(face = "bold", size = base),
          legend.text = element_text(size = base - 1),
          strip.background = element_rect(fill = "grey92", colour = "grey60"),
          strip.text = element_text(face = "bold", size = base))
}
save_png <- function(plot, file, w = 12, h = 7) {
  ggsave(file.path(fig_dir, file), plot, width = w, height = h, dpi = 300, units = "in", bg = "white")
  cat("  saved:", file, "\n")
}
save_csv <- function(d, file) write.csv(d, file.path(tab_dir, file), row.names = FALSE)

CODES     <- names(DISEASES)
SAUDPC    <- paste0("sAUDPC_", CODES)
SEV_TOTAL <- SEVERITY_MAX * length(DISEASES)
PALETTE   <- c("#E6862A", "#3B6EA5", "#2E7D32", "#7570B3", "#D95F02", "#1B9E77", "#C0392B")


## ---- 1. IMPORT DATA ---------------------------------------------------------
cat("\n=== 1. IMPORTING DATA ===\n")
dat   <- read.csv(PLOT_FILE)             # plot level  <-- PRIMARY
ts    <- read.csv(TIMESERIES_FILE)       # plot x assessment
vmean <- read.csv(VARIETY_MEANS_FILE)    # variety level
wx    <- read.csv(WEATHER_FILE)          # daily weather

check_cols <- function(d, cols, what) {
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(what, " lacks column(s): ", paste(miss, collapse = ", "))
}
check_cols(dat, c("Season", "Block", "Variety", "Plant_Stand", "Yield_t_ha", "AUDPC_TOT",
                  "sAUDPC_TOT", SAUDPC), "Plot file")
check_cols(ts, c("Season", "Block", "Variety", "DAP", "Total_Severity"), "Time-series file")
check_cols(vmean, c("Variety", SAUDPC, "sAUDPC_TOT", "DRI_percent", "Peak_Total_Severity",
                    "Mean_Total_Severity", "Yield_t_ha", "Yield_adj_t_ha", "Disease_Reaction",
                    "Plant_Stand"), "Variety means file")

fac <- function(d) { d$Season <- factor(d$Season); d$Block <- factor(d$Block)
                     d$Variety <- factor(d$Variety); d }
dat <- fac(dat); ts <- fac(ts)
dat$Plot <- factor(paste(dat$Season, dat$Block, dat$Variety))

# design sizes, taken from the data
N_SEASON  <- nlevels(dat$Season)
N_BLOCK   <- dat %>% group_by(Season) %>% summarise(b = n_distinct(Block)) %>% pull(b) %>% min()
N_VARIETY <- nlevels(dat$Variety)
cat("  plots:", nrow(dat), " varieties:", N_VARIETY, " seasons:", N_SEASON,
    " blocks per season:", N_BLOCK, "\n")

SEASONS    <- levels(dat$Season)
COL_SEASON <- setNames(rep(PALETTE, length.out = N_SEASON), SEASONS)
COL_DISEASE <- setNames(c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E")[seq_along(DISEASES)],
                        DISEASES)

# ordering used throughout
ord_disease <- vmean$Variety[order(vmean$sAUDPC_TOT)]   # least -> most diseased
ord_yield   <- vmean$Variety[order(-vmean$Yield_t_ha)]  # high -> low yield


## ---- 2. HELPER: Tukey HSD with compact letter display ------------------------
# HSD = q(1 - alpha; k, df_error) x sqrt(MSE / r). Two varieties differ when
# |mean_i - mean_j| > HSD; letters are assigned with multcompView so that
# varieties sharing a letter are never significantly different.
tukey_cld <- function(model, dfr, means, alpha = ALPHA) {
  mse <- summary(model)$sigma^2; dfe <- df.residual(model)
  hsd <- qtukey(1 - alpha, length(means), dfe) * sqrt(mse / dfr)
  m <- sort(means, decreasing = TRUE)
  sig <- outer(m, m, function(a, b) abs(a - b) > hsd)
  dimnames(sig) <- list(names(m), names(m))
  lt <- multcompView::multcompLetters(sig)$Letters
  data.frame(Variety = names(m), Mean = as.numeric(m),
             Group = unname(lt[names(m)]), HSD = hsd, row.names = NULL)
}


## ---- 3. VARIETIES EVALUATED (optional descriptive table) -----------------------
if (!is.null(VARIETY_INFO_FILE)) {
  cat("\n=== 3. Varieties ===\n")
  T1 <- read.csv(VARIETY_INFO_FILE)
  save_csv(T1, "T01_varieties.csv")
}


## ---- 4. VARIETY MEANS (sAUDPC, DRI, severity, yield) -------------------------
cat("=== 4. Variety means ===\n")
T2 <- vmean %>%
  select(Variety, all_of(SAUDPC), sAUDPC_TOT,
         DRI_percent, Peak_Total_Severity, Mean_Total_Severity,
         Yield_t_ha, Yield_adj_t_ha, Disease_Reaction) %>%
  arrange(desc(DRI_percent)) %>% mutate(across(where(is.numeric), ~round(.x, 2)))
save_csv(T2, "T02_variety_means.csv")
print(T2, row.names = FALSE)


## ---- 5. COMBINED OVER-SEASON ANOVA (mean squares) ---------------------------
cat("\n=== 5. Combined ANOVA ===\n")
combined_anova <- function(resp, d = dat) {
  f   <- as.formula(paste(resp, "~ Season + Season:Block + Variety + Variety:Season"))
  a   <- anova(lm(f, data = d))
  ms  <- setNames(a[, "Mean Sq"], rownames(a)); dfs <- setNames(a[, "Df"], rownames(a))
  mse <- ms[["Residuals"]]; dfe <- dfs[["Residuals"]]
  Fv  <- ms[["Variety"]] / mse; Fi <- ms[["Season:Variety"]] / mse
  Fs  <- ms[["Season"]] / ms[["Season:Block"]]
  data.frame(Trait = resp,
             MS_Season = ms[["Season"]], MS_RepSeason = ms[["Season:Block"]],
             MS_Variety = ms[["Variety"]], MS_VxS = ms[["Season:Variety"]], MS_Error = mse,
             DF_Error = dfe,
             F_Variety = Fv, P_Variety = pf(Fv, dfs[["Variety"]], dfe, lower.tail = FALSE),
             F_VxS = Fi, P_VxS = pf(Fi, dfs[["Season:Variety"]], dfe, lower.tail = FALSE),
             F_Season = Fs, P_Season = pf(Fs, dfs[["Season"]], dfs[["Season:Block"]], lower.tail = FALSE),
             CV_percent = 100 * sqrt(mse) / mean(d[[resp]], na.rm = TRUE))
}
T3 <- bind_rows(lapply(c(SAUDPC, "sAUDPC_TOT", "Yield_t_ha"), combined_anova))
T3$Signif_Variety <- ifelse(T3$P_Variety < 0.01, "**", ifelse(T3$P_Variety < 0.05, "*", "ns"))
T3n <- T3 %>% mutate(across(where(is.numeric), ~round(.x, 4)))
save_csv(T3n, "T03_combined_ANOVA.csv")
print(T3n[, c("Trait", "MS_Variety", "MS_Error", "F_Variety", "P_Variety", "Signif_Variety", "CV_percent")],
      row.names = FALSE)


## ---- 6. SENSITIVITY TO AUDPC FORMULATION --------------------------------------
cat("\n=== 6. AUDPC sensitivity ===\n")
sens <- lapply(c("AUDPC_TOT", "sAUDPC_TOT"), function(v) {
  out <- lapply(levels(dat$Season), function(s) {
    ds <- droplevels(subset(dat, Season == s))
    a  <- anova(lm(as.formula(paste(v, "~ Variety + Block")), data = ds))
    data.frame(Metric = v, Season = s, F = a["Variety", "F value"], P = a["Variety", "Pr(>F)"])
  }); bind_rows(out)
})
# relative AUDPC (per-season maximum) - third variant
dat <- dat %>% group_by(Season) %>% mutate(rAUDPC_TOT = AUDPC_TOT / max(AUDPC_TOT)) %>% ungroup()
dat <- fac(as.data.frame(dat)); dat$Plot <- factor(paste(dat$Season, dat$Block, dat$Variety))
rel <- lapply(levels(dat$Season), function(s) {
  ds <- droplevels(subset(dat, Season == s))
  a  <- anova(lm(rAUDPC_TOT ~ Variety + Block, data = ds))
  data.frame(Metric = "rAUDPC_TOT", Season = s, F = a["Variety", "F value"], P = a["Variety", "Pr(>F)"])
})
T4 <- bind_rows(bind_rows(sens), bind_rows(rel)) %>% mutate(across(where(is.numeric), ~round(.x, 4)))
save_csv(T4, "T04_AUDPC_sensitivity.csv")
print(T4, row.names = FALSE)


## ---- 7. YIELD BY SEASON WITH TUKEY GROUPS -------------------------------------
cat("\n=== 7. Yield + Tukey HSD ===\n")
tuk <- lapply(levels(dat$Season), function(s) {
  ds  <- droplevels(subset(dat, Season == s))
  m   <- lm(Yield_t_ha ~ Variety + Block, data = ds)
  mu  <- tapply(ds$Yield_t_ha, ds$Variety, mean)
  reps <- mean(table(ds$Variety))                    # replicates per variety
  cld <- tukey_cld(m, dfr = reps, means = mu)
  cld$Season <- s; cld
})
T5 <- bind_rows(tuk) %>% select(Season, Variety, Mean, Group) %>%
  pivot_wider(names_from = Season, values_from = c(Mean, Group)) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
save_csv(T5, "T05_yield_Tukey.csv")
print(as.data.frame(T5), row.names = FALSE)


## ---- 8. ANCOVA (plant stand covariate) ----------------------------------------
cat("\n=== 8. ANCOVA ===\n")
# (a) homogeneity of slopes
slope_test <- lapply(levels(dat$Season), function(s) {
  ds <- droplevels(subset(dat, Season == s))
  a  <- anova(lm(Yield_t_ha ~ Plant_Stand + Variety, data = ds),
              lm(Yield_t_ha ~ Plant_Stand * Variety, data = ds))
  data.frame(Season = s, F = a$F[2], P = a$`Pr(>F)`[2])
}) %>% bind_rows()
cat("  slope homogeneity:\n"); print(slope_test, row.names = FALSE)

m_anc <- lm(Yield_t_ha ~ Plant_Stand + Variety * Season, data = dat)
A  <- car::Anova(m_anc, type = 2)
T6 <- data.frame(Source = rownames(A), df = A$Df, SS = A$`Sum Sq`,
                 MS = A$`Sum Sq` / A$Df, F = A$`F value`, P = A$`Pr(>F)`) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))
save_csv(T6, "T06_ANCOVA.csv")
print(T6, row.names = FALSE)

# adjusted (least-squares) means
adj <- as.data.frame(emmeans(m_anc, ~ Variety))[, c("Variety", "emmean", "SE")]
names(adj) <- c("Variety", "Yield_adjusted", "SE")
unadj <- dat %>% group_by(Variety) %>% summarise(Yield_unadjusted = mean(Yield_t_ha),
                                                 Plant_Stand = mean(Plant_Stand), .groups = "drop")
T6b <- left_join(unadj, adj, by = "Variety") %>% mutate(across(where(is.numeric), ~round(.x, 2)))
save_csv(T6b, "T06b_adjusted_means.csv")
cv_before <- 100 * sqrt(summary(lm(Yield_t_ha ~ Season + Season:Block + Variety + Variety:Season, dat))$sigma^2) / mean(dat$Yield_t_ha)
cv_after  <- 100 * summary(m_anc)$sigma / mean(dat$Yield_t_ha)
cat(sprintf("  residual CV: %.1f%% -> %.1f%% after covariance adjustment\n", cv_before, cv_after))


## ---- 9. REML VARIANCE COMPONENTS + HERITABILITY (bootstrap CI) --------------
cat("\n=== 9. REML variance components & heritability ===\n")
# broad-sense heritability on a variety-mean basis:
#   H2 = s2_g / (s2_g + s2_gs / seasons + s2_e / (seasons x blocks))
h2_fun <- function(fit) {
  v  <- as.data.frame(VarCorr(fit))
  g  <- v$vcov[v$grp == "Variety"]; gs <- v$vcov[v$grp == "Variety:Season"]
  e  <- attr(VarCorr(fit), "sc")^2
  g / (g + gs / N_SEASON + e / (N_SEASON * N_BLOCK))
}
vc_row <- function(resp, label, boot = TRUE, nsim = N_BOOT) {
  f <- as.formula(paste(resp, "~ (1|Variety) + (1|Season) + (1|Variety:Season) + (1|Season:Block)"))
  m <- suppressWarnings(suppressMessages(lmer(f, data = dat)))
  v <- as.data.frame(VarCorr(m))
  g <- v$vcov[v$grp == "Variety"]; gs <- v$vcov[v$grp == "Variety:Season"]
  e <- attr(VarCorr(m), "sc")^2; h2 <- 100 * h2_fun(m)
  lo <- NA; hi <- NA
  if (boot) {
    set.seed(1)
    bb <- suppressWarnings(suppressMessages(
      bootMer(m, h2_fun, nsim = nsim, type = "parametric", use.u = FALSE)))
    q <- quantile(bb$t, c(.025, .975), na.rm = TRUE); lo <- 100 * q[1]; hi <- 100 * q[2]
  }
  data.frame(Trait = label, sigma2_g = g, sigma2_gs = gs, sigma2_e = e,
             H2_percent = h2, CI_lower = lo, CI_upper = hi, row.names = NULL)
}
T7 <- bind_rows(
  bind_rows(lapply(CODES, function(k) vc_row(paste0("sAUDPC_", k), paste(DISEASES[[k]], "sAUDPC"), boot = FALSE))),
  vc_row("sAUDPC_TOT",  "Total disease sAUDPC",  boot = TRUE),
  vc_row("Yield_t_ha",  "Tuber yield (t/ha)",    boot = TRUE),
  vc_row("Plant_Stand", "Establishment (stand)", boot = TRUE)) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
save_csv(T7, "T07_variance_components.csv")
print(T7, row.names = FALSE)


## ---- 10. CORRELATIONS AMONG VARIETY MEANS -------------------------------------
cat("\n=== 10. Correlations ===\n")
cm <- vmean[, c(SAUDPC, "sAUDPC_TOT", "Yield_t_ha")]
names(cm) <- c(unname(DISEASES), "Total Disease", "Yield")
R  <- cor(cm); Pm <- outer(1:ncol(cm), 1:ncol(cm),
                           Vectorize(function(i, j) cor.test(cm[[i]], cm[[j]])$p.value))
dimnames(Pm) <- dimnames(R)
T8 <- as.data.frame(round(R, 2)); T8 <- cbind(Trait = rownames(T8), T8)
save_csv(T8, "T08_correlations.csv")
save_csv(cbind(Trait = rownames(Pm), as.data.frame(round(Pm, 4))), "T08b_correlation_pvalues.csv")
print(T8, row.names = FALSE)


## ---- 11. VARIETY CLASSIFICATION ----------------------------------------------
cat("\n=== 11. Classification ===\n")
T9 <- vmean %>%
  mutate(Yield_class = ifelse(Yield_t_ha >= YIELD_HIGH, "High",
                         ifelse(Yield_t_ha >= YIELD_MEDIUM, "Medium", "Low")),
         Overall = ifelse(Yield_t_ha >= YIELD_HIGH & DRI_percent >= DRI_GOOD, "Excellent",
                     ifelse(Yield_t_ha >= YIELD_HIGH | DRI_percent >= DRI_GOOD, "Good",
                       ifelse(Yield_t_ha >= YIELD_MEDIUM, "Fair", "Poor")))) %>%
  select(Variety, DRI_percent, Mean_Total_Severity, Disease_Reaction,
         Yield_t_ha, Yield_adj_t_ha, Yield_class, Overall) %>%
  arrange(desc(DRI_percent)) %>% mutate(across(where(is.numeric), ~round(.x, 2)))
save_csv(T9, "T09_classification.csv")
print(T9, row.names = FALSE)


## ---- 12. MIXED MODELS, POWER AND REPEATED MEASURES ----------------------------
cat("\n=== 12. Mixed models / power / repeated measures ===\n")
mm_y <- suppressWarnings(lmer(Yield_t_ha ~ Variety * Season + (1|Season:Block), data = dat))
mm_d <- suppressWarnings(lmer(sAUDPC_TOT ~ Variety * Season + (1|Season:Block), data = dat))
cat("  Mixed model yield  - Variety F =", round(anova(mm_y)["Variety", "F value"], 2), "\n")
cat("  Mixed model disease- Variety F =", round(anova(mm_d)["Variety", "F value"], 2), "\n")

# minimum detectable difference between two variety means (two-sided alpha, given power)
mse_y <- T3$MS_Error[T3$Trait == "Yield_t_ha"]; mse_d <- T3$MS_Error[T3$Trait == "sAUDPC_TOT"]
dfe   <- T3$DF_Error[T3$Trait == "Yield_t_ha"]
mdd <- function(mse, dfe, r = N_BLOCK, s = N_SEASON)
  (qt(1 - ALPHA / 2, dfe) + qt(POWER, dfe)) * sqrt(2 * mse / (r * s))
T10 <- data.frame(
  Trait = c("Tuber yield (t/ha)", "Total disease sAUDPC"),
  Grand_mean = c(mean(dat$Yield_t_ha), mean(dat$sAUDPC_TOT)),
  MDD = c(mdd(mse_y, dfe), mdd(mse_d, dfe))) %>%
  mutate(Percent_of_mean = 100 * MDD / Grand_mean) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
names(T10)[names(T10) == "MDD"] <- sprintf("MDD_%dpct_power", round(100 * POWER))
save_csv(T10, "T10_power_MDD.csv")
print(T10, row.names = FALSE)

# repeated measures: does disease progress RATE differ by variety?
ts$Plot <- factor(paste(ts$Season, ts$Block, ts$Variety))
ts$cDAP <- ts$DAP - mean(ts$DAP)
m_ri <- suppressWarnings(lmer(Total_Severity ~ Variety * cDAP + Season + (1|Plot), data = ts))
m_rs <- suppressWarnings(lmer(Total_Severity ~ Variety * cDAP + Season + (1 + cDAP|Plot), data = ts,
                              control = lmerControl(optimizer = "bobyqa")))
slopes <- ts %>% group_by(Season, Block, Variety) %>%
  summarise(slope = coef(lm(Total_Severity ~ DAP))[2], .groups = "drop")
a_sl <- anova(lm(slope ~ Variety + Season + Season:Block, data = slopes))
T11 <- data.frame(
  Model = c("Random intercept only (anti-conservative)", "Random slopes per plot",
            "Per-plot slope ANOVA"),
  F = c(anova(m_ri)["Variety:cDAP", "F value"], anova(m_rs)["Variety:cDAP", "F value"], a_sl["Variety", "F value"]),
  P = c(anova(m_ri)["Variety:cDAP", "Pr(>F)"],  anova(m_rs)["Variety:cDAP", "Pr(>F)"],  a_sl["Variety", "Pr(>F)"])) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))
save_csv(T11, "T11_repeated_measures.csv")
print(T11, row.names = FALSE)


## ---- 13. EXCEL WORKBOOK OF ALL TABLES -----------------------------------------
cat("\n=== 13. Writing Excel workbook of all tables ===\n")
wb <- createWorkbook()
add_sheet <- function(nm, d) { addWorksheet(wb, nm); writeData(wb, nm, d)
  addStyle(wb, nm, createStyle(textDecoration = "bold", fgFill = "#2E5496", fontColour = "white"),
           rows = 1, cols = 1:ncol(d), gridExpand = TRUE); setColWidths(wb, nm, 1:ncol(d), "auto") }
if (exists("T1")) add_sheet("T01_Varieties", T1)
add_sheet("T02_VarietyMeans", T2)
add_sheet("T03_CombinedANOVA", T3n);  add_sheet("T04_AUDPCsensitivity", T4)
add_sheet("T05_Yield_Tukey", as.data.frame(T5)); add_sheet("T06_ANCOVA", T6)
add_sheet("T06b_AdjustedMeans", T6b); add_sheet("T07_VarianceComponents", T7)
add_sheet("T08_Correlations", T8);    add_sheet("T09_Classification", T9)
add_sheet("T10_Power_MDD", T10);      add_sheet("T11_RepeatedMeasures", T11)
add_sheet("SlopeHomogeneity", slope_test)
saveWorkbook(wb, file.path(tab_dir, "All_Tables.xlsx"), overwrite = TRUE)
cat("  saved: tables/All_Tables.xlsx\n")


## =============================================================================
## FIGURES
## =============================================================================
cat("\n=== 14. FIGURES ===\n")

# --- Weather and disease severity on a common time axis ----------------------
has_dates <- !any(grepl("<", c(names(SEASON_START_DATES), SEASON_START_DATES), fixed = TRUE))
if (has_dates) {
  check_cols(wx, c("Date", "Season_window", "Temp_C_daily_avg", "Humidity_pct_daily_avg",
                   "Rainfall_mm"), "Weather file")
  start <- as.Date(SEASON_START_DATES)
  names(start) <- names(SEASON_START_DATES)
  wx$Date <- as.Date(wx$Date)
  wx2 <- wx %>% filter(Season_window %in% names(start)) %>%
    mutate(DAP = as.numeric(Date - start[Season_window]))
  sev_scale <- 10 * 15 / SEV_TOTAL          # puts total severity on the weather axis
  sev_dap <- ts %>% group_by(Season, DAP) %>% summarise(Sev = mean(Total_Severity), .groups = "drop")
  wlong <- wx2 %>% select(Season = Season_window, DAP, Temp_C_daily_avg,
                          Humidity_pct_daily_avg, Rainfall_mm) %>%
    pivot_longer(-c(Season, DAP), names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = dplyr::recode(Variable, Temp_C_daily_avg = "Temperature (°C)",
                                    Humidity_pct_daily_avg = "Humidity (%)", Rainfall_mm = "Rainfall (mm)"))
  f1 <- ggplot() +
    geom_line(data = wlong, aes(DAP, Value, colour = Variable), linewidth = 0.8) +
    geom_point(data = sev_dap, aes(DAP, Sev * sev_scale), colour = "black", size = 2.6) +
    geom_line(data = sev_dap, aes(DAP, Sev * sev_scale), colour = "black", linewidth = 1) +
    facet_wrap(~Season, scales = "free_x") +
    scale_y_continuous(name = "Weather variable",
                       sec.axis = sec_axis(~ . / sev_scale,
                                           name = sprintf("Mean total disease severity (0-%g)", SEV_TOTAL))) +
    scale_colour_brewer(palette = "Set2", name = NULL) +
    labs(x = "Days after planting", title = "Weather and disease progress on a common time axis") +
    theme_pub() + theme(legend.position = "bottom")
  save_png(f1, "F01_weather_and_disease.png", 13, 6.5)
} else {
  cat("  weather figure skipped (SEASON_START_DATES not set)\n")
}

# --- Severity heat map (variety x crop age) ------------------------------------
hm <- ts %>% group_by(Season, Variety, DAP) %>%
  summarise(Sev = mean(Total_Severity), .groups = "drop") %>%
  mutate(Variety = factor(Variety, levels = rev(ord_disease)))
f2 <- ggplot(hm, aes(factor(DAP), Variety, fill = Sev)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", Sev)), size = 3.4) +
  facet_wrap(~Season, scales = "free_x") +
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd"), name = sprintf("Severity\n(0-%g)", SEV_TOTAL)) +
  labs(x = "Days after planting", y = NULL,
       title = "Mean total disease severity by variety and crop age") +
  theme_pub(13)
save_png(f2, "F02_severity_heatmap.png", 14, 7.5)

# --- Standardised AUDPC by variety and disease ----------------------------------
sa <- vmean %>% select(Variety, all_of(SAUDPC)) %>%
  pivot_longer(-Variety, names_to = "Disease", values_to = "sAUDPC") %>%
  mutate(Disease = factor(unname(DISEASES[sub("^sAUDPC_", "", Disease)]), levels = unname(DISEASES)),
         Variety = factor(Variety, levels = ord_disease))
f3 <- ggplot(sa, aes(Variety, sAUDPC, fill = Disease)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = COL_DISEASE, name = NULL) +
  labs(x = NULL, y = "Standardised AUDPC (severity units/day)",
       title = "Disease pressure by variety (ordered least to most diseased)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 40, hjust = 1),
                      legend.position = "top")
save_png(f3, "F03_standardised_AUDPC.png", 13, 7)

# --- Tuber yield by variety with Tukey letters ------------------------------------
ybar <- dat %>% group_by(Season, Variety) %>%
  summarise(Mean = mean(Yield_t_ha), SE = sd(Yield_t_ha) / sqrt(n()), .groups = "drop") %>%
  left_join(bind_rows(tuk)[, c("Season", "Variety", "Group")], by = c("Season", "Variety")) %>%
  mutate(Variety = factor(Variety, levels = ord_yield))
f4 <- ggplot(ybar, aes(Variety, Mean, fill = Season)) +
  geom_col(position = position_dodge(0.85), width = 0.78) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(0.85), width = 0.22) +
  geom_text(aes(label = Group, y = Mean + SE + 0.04 * max(Mean + SE)),
            position = position_dodge(0.85), size = 4.2) +
  scale_fill_manual(values = COL_SEASON, name = NULL) +
  labs(x = NULL, y = expression(bold("Tuber yield (t ha"^-1*")")),
       title = "Tuber yield by variety and season (Tukey HSD groups)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 40, hjust = 1),
                      legend.position = "top")
save_png(f4, "F04_yield_Tukey.png", 13.5, 7.5)

# --- Variety x season interaction ------------------------------------------------
inter <- dat %>% group_by(Season, Variety) %>%
  summarise(Mean = mean(Yield_t_ha), .groups = "drop")
f5 <- ggplot(inter, aes(Season, Mean, group = Variety, colour = Variety)) +
  geom_line(linewidth = 1.05) + geom_point(size = 3) +
  geom_text(data = subset(inter, Season == tail(SEASONS, 1)),
            aes(label = Variety), hjust = -0.12, size = 4, show.legend = FALSE) +
  scale_x_discrete(expand = expansion(mult = c(0.08, 0.42))) +
  labs(x = NULL, y = expression(bold("Tuber yield (t ha"^-1*")")),
       title = "Variety × season interaction for tuber yield") +
  theme_pub() + theme(legend.position = "none")
save_png(f5, "F05_variety_season_interaction.png", 11, 8)

# --- ANCOVA: unadjusted vs stand-adjusted yield ----------------------------------
ancp <- T6b %>% mutate(Variety = factor(Variety, levels = Variety[order(Yield_unadjusted)]))
f6 <- ggplot(ancp) +
  geom_segment(aes(x = Yield_unadjusted, xend = Yield_adjusted,
                   y = Variety, yend = Variety), colour = "grey65", linewidth = 1.4) +
  geom_point(aes(Yield_unadjusted, Variety, colour = "Unadjusted (realised)"), size = 4.4) +
  geom_point(aes(Yield_adjusted, Variety, colour = "Stand-adjusted (intrinsic)"), size = 4.4) +
  scale_colour_manual(values = c("Unadjusted (realised)" = "#E6862A",
                                 "Stand-adjusted (intrinsic)" = "#3B6EA5"), name = NULL) +
  labs(x = expression(bold("Tuber yield (t ha"^-1*")")), y = NULL,
       title = "ANCOVA: realised versus stand-adjusted tuber yield") +
  theme_pub() + theme(legend.position = "top")
save_png(f6, "F06_ANCOVA_adjusted_yield.png", 12, 7.5)

# --- Crop establishment ----------------------------------------------------------
est <- dat %>% group_by(Variety) %>%
  summarise(Stand = mean(Plant_Stand), SE = sd(Plant_Stand) / sqrt(n()), .groups = "drop") %>%
  arrange(desc(Stand)) %>% mutate(Variety = factor(Variety, levels = Variety))
f7 <- ggplot(est, aes(Variety, Stand, fill = Stand)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = Stand - SE, ymax = Stand + SE), width = 0.22) +
  geom_hline(yintercept = mean(dat$Plant_Stand), linetype = "dashed", colour = "grey35") +
  scale_fill_gradient(low = "#C0392B", high = "#2E7D32", guide = "none") +
  labs(x = NULL, y = "Plants established per plot",
       title = "Crop establishment by variety (dashed line = trial mean)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_png(f7, "F07_establishment.png", 13, 7)

# --- Resistance vs yield classification ------------------------------------------
cls <- vmean %>% mutate(Reaction = Disease_Reaction)
reactions <- unique(cls$Reaction)
f8 <- ggplot(cls, aes(DRI_percent, Yield_t_ha)) +
  annotate("rect", xmin = median(cls$DRI_percent), xmax = Inf,
           ymin = median(cls$Yield_t_ha), ymax = Inf, fill = "#E8F5E9") +
  geom_vline(xintercept = median(cls$DRI_percent), linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = median(cls$Yield_t_ha), linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = Reaction), size = 6, alpha = 0.95) +
  geom_text(aes(label = Variety), hjust = -0.14, size = 4.3) +
  scale_colour_manual(values = setNames(rep(c("#2E7D32", "#E65100", "#3B6EA5", "#7570B3"),
                                            length.out = length(reactions)), reactions), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.24))) +
  labs(x = "Disease Resistance Index (%) → more resistant",
       y = expression(bold("Tuber yield (t ha"^-1*")")),
       title = sprintf("Resistance and yield classification of the %d varieties", nrow(cls))) +
  theme_pub() + theme(legend.position = "bottom")
save_png(f8, "F08_resistance_vs_yield.png", 12, 8.5)

# --- Correlation heat map --------------------------------------------------------
Rl <- as.data.frame(as.table(R)); names(Rl) <- c("V1", "V2", "r")
f9 <- ggplot(Rl, aes(V1, V2, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 5) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(x = NULL, y = NULL, title = sprintf("Correlations among variety means (n = %d)", nrow(cm))) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_png(f9, "F09_correlation_heatmap.png", 10, 8)

# --- UPGMA dendrogram ------------------------------------------------------------
cl_in <- vmean %>% select(Variety, all_of(SAUDPC), Yield_t_ha, Plant_Stand)
rownames(cl_in) <- cl_in$Variety
Z    <- hclust(dist(scale(cl_in[, -1])), method = "average")
coph <- cor(cophenetic(Z), dist(scale(cl_in[, -1])))
png(file.path(fig_dir, "F10_UPGMA_dendrogram.png"), width = 12, height = 7.5,
    units = "in", res = 300, bg = "white")
par(mar = c(8, 5, 4, 2), cex.lab = 1.3, cex.axis = 1.15, font.lab = 2)
plot(Z, main = "", xlab = "", sub = "", ylab = "Euclidean distance", cex = 1.2, hang = -1)
title(main = sprintf("UPGMA clustering of the %d varieties (cophenetic r = %.2f)", nrow(cl_in), coph),
      cex.main = 1.5, font.main = 2)
k <- min(N_CLUSTERS, nrow(cl_in) - 1)
rect.hclust(Z, k = k, border = rep(c("#2E7D32", "#E6862A", "#3B6EA5", "#7570B3"), length.out = k))
invisible(dev.off()); cat("  saved: F10_UPGMA_dendrogram.png\n")

# --- Disease progress curves -----------------------------------------------------
f11 <- ggplot(ts, aes(DAP, Total_Severity, colour = Variety)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~Season, scales = "free_x") +
  labs(x = "Days after planting", y = sprintf("Total disease severity (0-%g)", SEV_TOTAL),
       title = "Disease progress curves by variety") +
  theme_pub(13) + theme(legend.position = "right", legend.key.size = unit(0.5, "cm"))
save_png(f11, "F11_disease_progress_curves.png", 14, 7)

# --- ANCOVA slope-homogeneity diagnostic ------------------------------------------
f12 <- ggplot(dat, aes(Plant_Stand, Yield_t_ha, colour = Season)) +
  geom_point(size = 2.6, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9) +
  facet_wrap(~Season, scales = "free_x") +
  scale_colour_manual(values = COL_SEASON, guide = "none") +
  labs(x = "Plant stand (plants per plot)", y = expression(bold("Tuber yield (t ha"^-1*")")),
       title = "ANCOVA slope-homogeneity check") +
  theme_pub(13)
save_png(f12, "F12_ANCOVA_slope_diagnostic.png", 13, 6)

# --- Yield distribution boxplot ----------------------------------------------------
f13 <- ggplot(dat, aes(factor(Variety, levels = ord_yield), Yield_t_ha, fill = Season)) +
  geom_boxplot(outlier.shape = 21, width = 0.72) +
  scale_fill_manual(values = COL_SEASON, name = NULL) +
  labs(x = NULL, y = expression(bold("Tuber yield (t ha"^-1*")")),
       title = "Distribution of plot yields by variety and season") +
  theme_pub() + theme(axis.text.x = element_text(angle = 40, hjust = 1),
                      legend.position = "top")
save_png(f13, "F13_yield_boxplot.png", 13.5, 7)

# --- Disease-yield relationship (tolerance) ---------------------------------------
f14 <- ggplot(vmean, aes(sAUDPC_TOT, Yield_t_ha)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey40", fill = "grey85") +
  geom_point(size = 5.5, colour = "#3B6EA5") +
  geom_text(aes(label = Variety), hjust = -0.13, size = 4.2) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.26))) +
  labs(x = "Total disease pressure (standardised AUDPC)",
       y = expression(bold("Tuber yield (t ha"^-1*")")),
       title = "Disease-yield relationship (flat slope indicates tolerance/escape)") +
  theme_pub()
save_png(f14, "F14_disease_yield_relationship.png", 12, 7.5)


## ---- 15. SESSION SUMMARY ----------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
cat("\n=============================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("  Tables  ->", tab_dir, sprintf("(%d CSV + All_Tables.xlsx)\n", length(list.files(tab_dir, "\\.csv$"))))
cat("  Figures ->", fig_dir, sprintf("(%d PNG, 300 dpi)\n", length(list.files(fig_dir, "\\.png$"))))
cat("=============================================================\n")
