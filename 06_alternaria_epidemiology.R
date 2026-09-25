################################################################################
#
# EPIDEMIOLOGY OF FOLIAR DISEASES OF POTATO ACROSS FIELD SITES
# Climate, disease progress, apparent infection rate, AUDPC and yield
#
# Clean, end-to-end reproducible analysis script
# (raw data -> all tables, figures and statistics)
#
# Script  : 06_alternaria_epidemiology.R
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
# Tested  : R 4.3.3 (versions of every run are written to session_info.txt)
#
# ---------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
# ---------------------------------------------------------------------------
# Starting from THREE RAW DATA FILES only, it runs, in order:
#   Step 1  Data loading and quality checks
#   Step 2  Site climate characterisation and favourable/high-risk days
#   Step 3  Disease severity and incidence summaries
#   Step 4  Climate-disease relationship (site-level correlations)
#   Step 5  Disease progression rates and progress curves
#   Step 6  Apparent infection rate (Van der Plank) and site comparison
#   Step 7  AUDPC and treatment ANOVA (RCBD within site)
#   Step 8  Yield, yield loss, economic loss and yield-disease regressions
#   Step 9  Session information (R and package versions)
#
# Each step writes CSV table(s) to <OUTPUT_DIR>/tables/, each figure is saved
# to <OUTPUT_DIR>/figures/ as a PNG and an LZW-compressed TIFF, the whole
# console output is saved to <OUTPUT_DIR>/analysis_log.txt, and the package
# versions to <OUTPUT_DIR>/session_info.txt.
#
# ---------------------------------------------------------------------------
# SOFTWARE
# ---------------------------------------------------------------------------
# R 4.1 or newer (https://cran.r-project.org); RStudio is optional.
# The five R packages used (ggplot2, dplyr, tidyr, zoo, multcompView) are
# installed automatically by Section B on the first run. No compiler needed.
# Optional: freeze package versions with renv:
#   install.packages("renv"); renv::init(); renv::snapshot()
#
# ---------------------------------------------------------------------------
# HOW TO RUN
# ---------------------------------------------------------------------------
# 1. Set the file paths and study settings in Section A.
# 2. RStudio: open this file and click "Source", or
#    Terminal: Rscript 06_alternaria_epidemiology.R
#
# ---------------------------------------------------------------------------
# INPUT DATA STRUCTURE
# ---------------------------------------------------------------------------
# Disease file : Site, Block, Treatment, Sample_ID, Time_Point, and one
#                severity column per disease named <CODE>_Sev (codes set in
#                DISEASES, Section A). One row = one tagged plant at one
#                assessment; severity on the 0-SEVERITY_MAX scale;
#                Time_Point 1..n matches ASSESSMENT_DAYS.
# Climate file : Date, Temp, Humidity, Site (optional: Rain, Dew)
#                Daily means; date as dd/mm/yyyy, yyyy-mm-dd or "01 May 2025".
# Yield file   : Site, Treatment, Replication (or Block), Yield (t/ha)
# Site and treatment codes must be identical in the three files.
#
# Experimental unit: the PLOT (Site x Block x Treatment). Tagged plants are
# sub-samples and are averaged to plot means before any ANOVA (RCBD analysed
# within each site, treatment fixed, block random). Incidence is derived from
# severity (a plant is diseased when severity > 0). AUDPC is always computed
# from the raw severity readings.
################################################################################


################################################################################
# SECTION A. USER SETTINGS (the ONLY section you normally need to edit)
################################################################################

# A1. Input files and output folder -------------------------------------------
DISEASE_FILE <- "<PATH_TO_DISEASE_DATA>.csv"
CLIMATE_FILE <- "<PATH_TO_CLIMATE_DATA>.csv"
YIELD_FILE   <- "<PATH_TO_YIELD_DATA>.csv"
OUTPUT_DIR   <- "<PATH_TO_OUTPUT_FOLDER>"

# A2. Study design --------------------------------------------------------------
# Diseases: names = column prefix in the disease file (<CODE>_Sev), values = labels
DISEASES        <- c(EB = "Early blight", BS = "Brown spot")
ASSESSMENT_DAYS <- c(15, 30, 45, 60, 75)   # crop age (days) at Time_Point 1, 2, ...
SEVERITY_MAX    <- 5                       # 0-5 ordinal severity scale

# Order of sites in tables and figures. NULL = order of first appearance in the
# climate file; or give the codes, e.g. c("SiteA", "SiteB", "SiteC").
SITES <- NULL

# Optional display names for sites in tables and figures, as code = label, e.g.
#   c(SiteA = "Site A (highland)"). NULL = codes with "_" shown as spaces.
SITE_DISPLAY <- NULL

# Optional: report each site's infection rate as % difference from this site
# (a site code). NULL = not reported.
REFERENCE_SITE <- NULL

# Optional legend labels for treatments, e.g.
#   c("1" = "T1: inoculated at planting", "2" = "T2: ...", "4" = "T4: control")
# NULL = "Treatment 1", "Treatment 2", ...
TREATMENT_LABELS <- NULL

# A3. Thresholds -----------------------------------------------------------------
FAV_TEMP_MIN <- 15      # favourable daily mean temperature range (deg C)
FAV_TEMP_MAX <- 25
FAV_RH_MIN   <- 60      # favourable relative humidity (> %)
RATE_LOW     <- 0.03    # apparent infection rate classes: r < RATE_LOW = low
RATE_HIGH    <- 0.10    # RATE_LOW <= r < RATE_HIGH = moderate; >= RATE_HIGH = high
ALPHA        <- 0.05    # significance level

# A4. Economics (optional) --------------------------------------------------------
# Leave PRICE_PER_TONNE as NA to omit the economic-loss columns.
PRICE_PER_TONNE <- NA_real_              # farm-gate price per tonne, local currency
CURRENCY        <- "<CURRENCY_CODE>"     # e.g. "USD"
LOCAL_PER_USD   <- NA_real_              # exchange rate; NA = no USD columns

# A5. Figure options -------------------------------------------------------------
FIG_DPI           <- 600       # journal resolution
FIG_FONT          <- "serif"   # maps to Times New Roman on Windows
ADD_DAY0_BASELINE <- FALSE     # TRUE adds an assumed severity of 0 at day 0 to the
                               # progress curves. It is NOT an observed value; if
                               # switched on, say so in the figure caption.
ROLLING_WINDOW    <- 7         # days, for the climate moving averages


################################################################################
# SECTION B. INSTALL AND LOAD THE REQUIRED R PACKAGES
################################################################################

cat("\n==================================================================\n")
cat("  DISEASE EPIDEMIOLOGY - REPRODUCIBLE ANALYSIS\n")
cat("==================================================================\n\n")

if (getRversion() < "4.1.0") {
  stop("R 4.1.0 or newer is required (4.3+ recommended). ",
       "Please update R from https://cran.r-project.org")
}

options(repos = c(CRAN = "https://cloud.r-project.org"),
        stringsAsFactors = FALSE,
        scipen = 999,
        dplyr.summarise.inform = FALSE)

required_packages <- c("ggplot2", "dplyr", "tidyr", "zoo", "multcompView")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package '%s' from CRAN ...\n", pkg))
    install.packages(pkg, dependencies = c("Depends", "Imports"))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("B. Packages loaded:",
    paste(sprintf("%s %s", required_packages,
                  sapply(required_packages, function(p) as.character(packageVersion(p)))),
          collapse = ", "), "\n\n")

set.seed(2025)   # no random numbers are used, but fixed for good practice


################################################################################
# SECTION C. FOLDERS, LOG FILE AND HELPER FUNCTIONS
################################################################################

# C1. Check settings and create output folders --------------------------------
for (v in c("DISEASE_FILE", "CLIMATE_FILE", "YIELD_FILE", "OUTPUT_DIR")) {
  if (grepl("<", get(v), fixed = TRUE)) stop("Set ", v, " in Section A before running.")
}
for (f in c(DISEASE_FILE, CLIMATE_FILE, YIELD_FILE)) {
  if (!file.exists(f)) stop("Input file not found: ", f, "\nPlease check Section A.")
}
HAS_PRICE <- !is.na(PRICE_PER_TONNE)
HAS_USD   <- HAS_PRICE && !is.na(LOCAL_PER_USD)
if (HAS_PRICE && grepl("<", CURRENCY, fixed = TRUE)) stop("Set CURRENCY in Section A.")

TABLE_DIR  <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# C2. Copy everything printed to the console into a log file -----------------
while (sink.number() > 0) sink()   # close any log left open by an interrupted run
log_file <- file.path(OUTPUT_DIR, "analysis_log.txt")
log_con  <- file(log_file, open = "wt")
sink(log_con, split = TRUE)

cat(sprintf("Run started : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output      : %s\n\n", OUTPUT_DIR))

# C3. Small helper functions ---------------------------------------------------
banner <- function(txt) {
  cat("\n------------------------------------------------------------------\n")
  cat(" ", txt, "\n")
  cat("------------------------------------------------------------------\n")
}

# Standard error of the mean (ignores missing values)
se <- function(x) { x <- x[!is.na(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }

save_table <- function(df, name) {
  path <- file.path(TABLE_DIR, paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("  -> table saved : tables/%s.csv (%d rows)\n", name, nrow(df)))
  invisible(path)
}

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(FIGURE_DIR, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = FIG_DPI, bg = "white")
  ggsave(file.path(FIGURE_DIR, paste0(name, ".tiff")), plot,
         width = width, height = height, units = "in", dpi = FIG_DPI, bg = "white",
         device = "tiff", compression = "lzw")
  cat(sprintf("  -> figure saved: figures/%s.png / .tiff\n", name))
}

fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# Significance stars computed from the P value (never typed by hand)
p_stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",
                         ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

# Parse dates robustly: tries day-first formats first
parse_dates <- function(x) {
  for (fmt in c("%d/%m/%Y", "%d-%m-%Y", "%Y-%m-%d", "%d %B %Y", "%d %b %Y")) {
    d <- as.Date(as.character(x), format = fmt)
    if (all(!is.na(d[!is.na(x) & x != ""]))) return(d)
  }
  stop("Could not recognise the date format in the climate file.")
}

# Colour-blind-safe palette; colour is always paired with shape/line type so
# figures also read in grey-scale print.
PALETTE    <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8f5bd6", "#d6457a", "#5a5a5a")
SHAPES     <- c(16, 17, 15, 18, 8, 3, 4)
FILL_SHAPES <- c(21, 24, 22, 23, 25)
LINETYPES  <- c("solid", "dashed", "dotdash", "dotted", "longdash", "twodash")
recycle    <- function(v, n) rep(v, length.out = n)

# Common publication theme for all figures
theme_pub <- theme_bw(base_size = 12, base_family = FIG_FONT) +
  theme(plot.title        = element_blank(),     # titles go in the caption
        strip.text        = element_text(face = "bold", size = 11),
        strip.background  = element_rect(fill = "grey95", colour = "grey70"),
        axis.title        = element_text(face = "bold"),
        panel.grid.major  = element_line(colour = "grey88", linewidth = 0.3),
        panel.grid.minor  = element_blank(),
        legend.position   = "top",
        legend.title      = element_text(face = "bold"),
        legend.key.width  = unit(1.6, "lines"))


################################################################################
# STEP 1. LOAD THE RAW DATA AND CHECK ITS QUALITY
################################################################################
banner("STEP 1. Loading raw data and quality checks")

disease_raw <- read.csv(DISEASE_FILE, check.names = FALSE)
climate_raw <- read.csv(CLIMATE_FILE, check.names = FALSE)
yield_raw   <- read.csv(YIELD_FILE,   check.names = FALSE)

# Remove stray quotes/spaces that spreadsheet software sometimes adds to names
clean_names <- function(df) { names(df) <- trimws(gsub('"', "", names(df))); df }
disease_raw <- clean_names(disease_raw)
climate_raw <- clean_names(climate_raw)
yield_raw   <- clean_names(yield_raw)

# 1.1 Check the required columns exist ----------------------------------------
SEV_COLS <- paste0(names(DISEASES), "_Sev")
need <- list(disease = c("Site", "Block", "Treatment", "Sample_ID", "Time_Point", SEV_COLS),
             climate = c("Date", "Temp", "Humidity", "Site"),
             yield   = c("Treatment", "Site", "Yield"))   # plus Replication or Block
for (nm in names(need)) {
  df <- get(paste0(nm, "_raw"))
  miss <- setdiff(need[[nm]], names(df))
  if (length(miss)) stop(sprintf("The %s file lacks column(s): %s", nm, paste(miss, collapse = ", ")))
}
if (max(disease_raw$Time_Point, na.rm = TRUE) > length(ASSESSMENT_DAYS))
  stop("Time_Point exceeds the number of ASSESSMENT_DAYS set in Section A.")

# 1.2 Sites and treatments ------------------------------------------------------
SITE_CODES  <- if (is.null(SITES)) unique(as.character(climate_raw$Site)) else SITES
SITE_LABELS <- if (is.null(SITE_DISPLAY)) gsub("_", " ", SITE_CODES) else
  unname(ifelse(SITE_CODES %in% names(SITE_DISPLAY), SITE_DISPLAY[SITE_CODES], gsub("_", " ", SITE_CODES)))
site_factor <- function(x) factor(SITE_LABELS[match(x, SITE_CODES)], levels = SITE_LABELS)
SITE_COLOURS <- setNames(recycle(PALETTE, length(SITE_LABELS)), SITE_LABELS)

TRT_LEVELS <- as.character(sort(unique(disease_raw$Treatment)))
TRT_LABELS <- if (is.null(TREATMENT_LABELS)) setNames(paste("Treatment", TRT_LEVELS), TRT_LEVELS) else
  TREATMENT_LABELS[TRT_LEVELS]
TRT_COLOURS <- setNames(recycle(PALETTE, length(TRT_LEVELS)), TRT_LEVELS)
TRT_SHAPES  <- setNames(recycle(SHAPES, length(TRT_LEVELS)), TRT_LEVELS)
TRT_LINES   <- setNames(recycle(LINETYPES, length(TRT_LEVELS)), TRT_LEVELS)

# 1.3 Put every variable in the right type -------------------------------------
disease <- disease_raw %>%
  mutate(Site       = site_factor(Site),
         Block      = factor(Block),
         Treatment  = factor(Treatment, levels = TRT_LEVELS),
         Sample_ID  = as.integer(Sample_ID),
         Time_Point = as.integer(Time_Point),
         Day        = ASSESSMENT_DAYS[Time_Point],
         across(all_of(SEV_COLS), as.numeric),
         Plot       = interaction(Site, Block, Treatment, drop = TRUE))

climate <- climate_raw %>%
  mutate(Site     = site_factor(Site),
         Date     = parse_dates(Date),
         Temp     = as.numeric(Temp),
         Humidity = as.numeric(Humidity)) %>%
  arrange(Site, Date)

# the block column may be called "Replication" or "Block"
if (!"Replication" %in% names(yield_raw)) {
  if ("Block" %in% names(yield_raw)) yield_raw$Replication <- yield_raw$Block
  else stop("Yield file needs a 'Replication' or 'Block' column.")
}
yield <- yield_raw %>%
  mutate(Site      = site_factor(Site),
         Treatment = factor(Treatment, levels = TRT_LEVELS),
         Block     = factor(Replication),   # replication = block of the RCBD
         Yield     = as.numeric(Yield))

if (any(is.na(disease$Site)) || any(is.na(climate$Site)) || any(is.na(yield$Site))) {
  stop("Site name not found in all files. Expected: ", paste(SITE_CODES, collapse = ", "))
}

# 1.4 Data-quality table --------------------------------------------------------
plants_per_plot <- disease %>% distinct(Site, Block, Treatment, Sample_ID) %>%
  count(Site, Block, Treatment, name = "plants")

quality <- bind_rows(
  disease %>% group_by(Site) %>% summarise(
    `Disease: records (plant x assessment)` = n(),
    `Disease: plots`                        = n_distinct(Plot),
    `Disease: tagged plants`                = n_distinct(paste(Block, Treatment, Sample_ID)),
    `Disease: assessments`                  = n_distinct(Time_Point),
    across(all_of(SEV_COLS), ~ sum(is.na(.x)), .names = "Disease: missing {.col}")) %>%
    pivot_longer(-Site, names_to = "Check", values_to = "Value") %>% mutate(Value = as.character(Value)),
  climate %>% group_by(Site) %>% summarise(
    `Climate: daily records`       = as.character(n()),
    `Climate: first date`          = as.character(min(Date)),
    `Climate: last date`           = as.character(max(Date)),
    `Climate: missing temperature` = as.character(sum(is.na(Temp))),
    `Climate: missing humidity`    = as.character(sum(is.na(Humidity)))) %>%
    pivot_longer(-Site, names_to = "Check", values_to = "Value"),
  yield %>% group_by(Site) %>% summarise(
    `Yield: plots`          = as.character(n()),
    `Yield: missing values` = as.character(sum(is.na(Yield)))) %>%
    pivot_longer(-Site, names_to = "Check", values_to = "Value")
) %>% pivot_wider(names_from = Site, values_from = Value)

print(as.data.frame(quality))
cat(sprintf("\nTagged plants per plot: min %d, max %d, mean %.1f\n",
            min(plants_per_plot$plants), max(plants_per_plot$plants),
            mean(plants_per_plot$plants)))
save_table(quality,         "T01_data_quality")
save_table(plants_per_plot, "T01b_plants_per_plot")


################################################################################
# STEP 2. SITE CLIMATE CHARACTERISATION
################################################################################
banner("STEP 2. Site climate characterisation")

# 2.1 Flag favourable days ---------------------------------------------------
climate <- climate %>%
  mutate(fav_temp  = Temp >= FAV_TEMP_MIN & Temp <= FAV_TEMP_MAX,
         fav_rh    = Humidity > FAV_RH_MIN,
         high_risk = fav_temp & fav_rh)          # both conditions on the same day

# 2.2 Climate and disease-risk days per site ---------------------------------
climate_tab <- climate %>% group_by(Site) %>% summarise(
  Days                     = n(),
  Temp_mean                = round(mean(Temp, na.rm = TRUE), 1),
  Temp_min_daily           = round(min(Temp, na.rm = TRUE), 1),
  Temp_max_daily           = round(max(Temp, na.rm = TRUE), 1),
  Humidity_mean            = round(mean(Humidity, na.rm = TRUE), 1),
  Favourable_Temp_Days     = sum(fav_temp, na.rm = TRUE),
  Favourable_Humidity_Days = sum(fav_rh, na.rm = TRUE),
  High_Risk_Days           = sum(high_risk, na.rm = TRUE),
  Percent_High_Risk_Days   = round(100 * High_Risk_Days / Days, 1),
  Total_Rain_mm            = if ("Rain" %in% names(climate)) round(sum(as.numeric(Rain), na.rm = TRUE), 1) else NA_real_)
print(as.data.frame(climate_tab))
save_table(climate_tab, "T02_site_climate_and_risk_days")

# 2.3 Moving averages of temperature and RH ----------------------------------
TEMP_LAB <- "Temperature (°C)"   # unicode escape keeps the script ASCII-safe
climate_fig <- climate %>% group_by(Site) %>% arrange(Date) %>%
  mutate(!!TEMP_LAB := zoo::rollmeanr(Temp, k = ROLLING_WINDOW, fill = NA),
         `Relative humidity (%)` = zoo::rollmeanr(Humidity, k = ROLLING_WINDOW, fill = NA)) %>%
  ungroup() %>%
  pivot_longer(c(all_of(TEMP_LAB), `Relative humidity (%)`),
               names_to = "Variable", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Variable = factor(Variable, levels = c(TEMP_LAB, "Relative humidity (%)")))

# Reference bands showing the favourable range
fav_bands <- data.frame(
  Variable = factor(c(TEMP_LAB, "Relative humidity (%)"), levels = levels(climate_fig$Variable)),
  ymin = c(FAV_TEMP_MIN, FAV_RH_MIN), ymax = c(FAV_TEMP_MAX, Inf))

fig_climate <- ggplot(climate_fig, aes(Date, Value, colour = Site, linetype = Site)) +
  geom_rect(data = fav_bands, inherit.aes = FALSE,
            aes(xmin = min(climate_fig$Date), xmax = max(climate_fig$Date),
                ymin = ymin, ymax = ymax),
            fill = "grey92", colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_colour_manual(values = SITE_COLOURS) +
  scale_linetype_manual(values = recycle(LINETYPES, length(SITE_LABELS))) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  labs(x = "Date", y = NULL, colour = "Site", linetype = "Site") +
  theme_pub + theme(strip.placement = "outside", strip.background = element_blank(),
                    plot.margin = margin(5.5, 18, 5.5, 5.5))
save_figure(fig_climate, "F01_climate_moving_average", width = 7.5, height = 6)
cat(sprintf(paste0("  Caption: %d-day moving averages of daily mean temperature (A) and relative\n",
                   "  humidity (B); shaded bands = favourable range (%g-%g deg C; RH > %g %%).\n"),
            ROLLING_WINDOW, FAV_TEMP_MIN, FAV_TEMP_MAX, FAV_RH_MIN))


################################################################################
# STEP 3. DISEASE SEVERITY AND INCIDENCE SUMMARIES
################################################################################
banner("STEP 3. Disease severity and incidence summaries")

# 3.1 Long format: one row per plant x assessment x disease -------------------
code_of <- setNames(names(DISEASES), DISEASES)          # label -> code
disease_long <- disease %>%
  select(Site, Block, Treatment, Plot, Sample_ID, Time_Point, Day, all_of(SEV_COLS)) %>%
  pivot_longer(all_of(SEV_COLS), names_to = "Disease", values_to = "Severity") %>%
  mutate(Disease  = factor(unname(DISEASES[sub("_Sev$", "", Disease)]), levels = unname(DISEASES)),
         Infected = ifelse(is.na(Severity), NA, Severity > 0))

# 3.2 Plot means at each assessment (the experimental unit) -------------------
plot_time <- disease_long %>%
  group_by(Site, Block, Treatment, Plot, Disease, Time_Point, Day) %>%
  summarise(Severity  = mean(Severity, na.rm = TRUE),
            Incidence = 100 * mean(Infected, na.rm = TRUE),
            n_plants  = sum(!is.na(Infected)), .groups = "drop")

# 3.3 Site x treatment x assessment (mean over plots) --------------------------
sev_trt <- plot_time %>% group_by(Site, Treatment, Disease, Time_Point, Day) %>%
  summarise(Severity_mean = mean(Severity), Severity_SE = se(Severity),
            Incidence_mean = mean(Incidence), n_plots = n(), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
save_table(sev_trt, "T03_severity_by_site_treatment_time")

# 3.4 Site x assessment, pooled over plants -------------------------------------
sev_site <- disease_long %>% group_by(Site, Disease, Time_Point, Day) %>%
  summarise(Severity_mean = mean(Severity, na.rm = TRUE),
            Severity_SE   = se(Severity),
            Incidence_pct = 100 * mean(Infected, na.rm = TRUE),
            n_plants      = sum(!is.na(Severity)), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 2)))
print(as.data.frame(sev_site))
save_table(sev_site, "T04_severity_incidence_by_site_time")

# Largest increase between consecutive assessments (vulnerability window)
window <- sev_site %>% group_by(Site, Disease) %>% arrange(Day) %>%
  mutate(increase = Severity_mean - lag(Severity_mean),
         interval = paste0(lag(Day), "-", Day, " DAP")) %>%
  filter(!is.na(increase)) %>% slice_max(increase, n = 1, with_ties = FALSE) %>%
  select(Site, Disease, interval, increase)
cat("\nInterval with the largest severity increase (vulnerability window):\n")
print(as.data.frame(window))
save_table(window, "T04b_largest_increase_interval")


################################################################################
# STEP 4. CLIMATE-DISEASE RELATIONSHIP
################################################################################
banner("STEP 4. Climate-disease relationship")

FINAL_DAY <- max(disease$Day)

# 4.1 Climate means next to final-assessment disease values -------------------
final_wide <- sev_site %>% filter(Time_Point == max(Time_Point)) %>%
  transmute(Site, Code = code_of[as.character(Disease)],
            Sev = Severity_mean, Inc = Incidence_pct) %>%
  pivot_longer(c(Sev, Inc), names_to = "metric", values_to = "value") %>%
  mutate(column = paste0(Code, "_", metric)) %>%
  select(Site, column, value) %>%
  pivot_wider(names_from = column, values_from = value)

SEV_FINAL <- paste0(names(DISEASES), "_Sev")
INC_FINAL <- paste0(names(DISEASES), "_Inc")
climate_disease <- climate_tab %>% select(Site, Temp_mean, Humidity_mean) %>%
  left_join(final_wide, by = "Site") %>%
  select(Site, Temp_mean, Humidity_mean, all_of(SEV_FINAL), all_of(INC_FINAL))
print(as.data.frame(climate_disease))
save_table(climate_disease, "T05_climate_disease_final_assessment")
cat(sprintf("  Note: severity and incidence are from the final assessment (%d DAP).\n", FINAL_DAY))

# 4.2 Pearson correlations between site climate means and disease -------------
# Climate is recorded once per site, so n = number of sites (df = n - 2).
# With few sites these coefficients are descriptive only: P values have very
# little power and a single site can drive the sign.
site_audpc_tmp <- disease_long %>% arrange(Time_Point) %>%
  group_by(Site, Block, Treatment, Sample_ID, Disease) %>%
  summarise(AUDPC = sum((head(Severity, -1) + tail(Severity, -1)) / 2), .groups = "drop") %>%
  group_by(Site, Disease) %>% summarise(AUDPC = mean(AUDPC, na.rm = TRUE), .groups = "drop")

corr_input <- climate_disease %>%
  left_join(site_audpc_tmp %>% pivot_wider(names_from = Disease, values_from = AUDPC,
                                           names_prefix = "AUDPC_"), by = "Site")
climate_vars <- c("Temp_mean", "Humidity_mean")
AUDPC_COLS   <- paste0("AUDPC_", DISEASES)
disease_vars <- c(SEV_FINAL, INC_FINAL, AUDPC_COLS)

corr_tab <- expand.grid(Climate = climate_vars, Disease_metric = disease_vars,
                        stringsAsFactors = FALSE) %>% rowwise() %>%
  mutate(n = nrow(corr_input),
         r = cor(corr_input[[Climate]], corr_input[[Disease_metric]]),
         P = { t <- r * sqrt((n - 2) / (1 - r^2)); 2 * pt(-abs(t), df = n - 2) }) %>%
  ungroup() %>% mutate(r = round(r, 3), P = round(P, 3), Sig = p_stars(P))
print(as.data.frame(corr_tab))
save_table(corr_tab, "T06_climate_disease_correlations")

# 4.3 Correlation heat map (diverging, grey midpoint) ---------------------------
metric_labels <- setNames(
  c(sprintf("%s severity (%d DAP)", names(DISEASES), FINAL_DAY),
    sprintf("%s incidence (%d DAP)", names(DISEASES), FINAL_DAY),
    sprintf("%s AUDPC", names(DISEASES))),
  disease_vars)
fig_corr <- ggplot(corr_tab, aes(Climate, Disease_metric, fill = r)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f%s", r, ifelse(Sig %in% c("", "ns"), "", Sig))),
            family = FIG_FONT, size = 3.8) +
  scale_fill_gradient2(low = "#2a78d6", mid = "grey95", high = "#eb6834",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  scale_x_discrete(labels = c(Temp_mean = "Mean temperature", Humidity_mean = "Mean RH")) +
  scale_y_discrete(limits = rev(disease_vars), labels = metric_labels) +
  labs(x = NULL, y = NULL) + theme_pub + theme(panel.grid = element_blank())
save_figure(fig_corr, "F02_climate_disease_correlation_heatmap", width = 6, height = 5)


################################################################################
# STEP 5. DISEASE PROGRESSION RATES AND PROGRESS CURVES
################################################################################
banner("STEP 5. Disease progression rates and progress curves")

# 5.1 Progression rate = slope of plot-mean severity on assessment number -----
#     (severity units per assessment interval), one slope per plot.
slope <- function(y, x) { ok <- !is.na(y); if (sum(ok) < 2) NA_real_ else unname(coef(lm(y[ok] ~ x[ok]))[2]) }

prog_plot <- plot_time %>% group_by(Site, Block, Treatment, Disease) %>%
  summarise(rate_per_assessment = slope(Severity, Time_Point),
            R2 = { ok <- !is.na(Severity); if (sum(ok) > 2 && var(Severity[ok]) > 0)
                     summary(lm(Severity[ok] ~ Time_Point[ok]))$r.squared else NA_real_ },
            .groups = "drop") %>%
  mutate(rate_per_day = rate_per_assessment / diff(ASSESSMENT_DAYS)[1])
save_table(prog_plot %>% mutate(across(where(is.double), ~ round(.x, 3))),
           "T07_progression_rate_per_plot")

prog_site <- prog_plot %>% group_by(Site, Disease) %>%
  summarise(mean_rate = mean(rate_per_assessment), SE = se(rate_per_assessment),
            min_rate = min(rate_per_assessment), max_rate = max(rate_per_assessment),
            mean_R2 = mean(R2, na.rm = TRUE), n_plots = n(), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 2)))
print(as.data.frame(prog_site))
save_table(prog_site, "T07b_progression_rate_by_site")

# 5.2 Disease progress curves (mean of plots +/- SE) ----------------------------
curve_data <- sev_trt %>% mutate(Treatment = as.character(Treatment))
if (ADD_DAY0_BASELINE) {
  curve_data <- bind_rows(
    curve_data %>% distinct(Site, Treatment, Disease) %>%
      mutate(Day = 0, Severity_mean = 0, Severity_SE = NA), curve_data)
}
dodge <- position_dodge(width = diff(ASSESSMENT_DAYS)[1] / 5)

fig_progress <- ggplot(curve_data, aes(Day, Severity_mean, colour = Treatment,
                                       shape = Treatment, linetype = Treatment, group = Treatment)) +
  geom_errorbar(aes(ymin = pmax(0, Severity_mean - Severity_SE),
                    ymax = Severity_mean + Severity_SE),
                width = 2.5, linewidth = 0.4, linetype = "solid", position = dodge) +
  geom_line(linewidth = 0.7, position = dodge) +
  geom_point(size = 2.4, position = dodge) +
  facet_grid(Disease ~ Site) +
  scale_colour_manual(values = TRT_COLOURS, labels = TRT_LABELS) +
  scale_shape_manual(values = TRT_SHAPES, labels = TRT_LABELS) +
  scale_linetype_manual(values = TRT_LINES, labels = TRT_LABELS) +
  scale_x_continuous(breaks = c(if (ADD_DAY0_BASELINE) 0, ASSESSMENT_DAYS)) +
  scale_y_continuous(limits = c(0, SEVERITY_MAX), breaks = 0:SEVERITY_MAX) +
  labs(x = "Crop age (days after planting)",
       y = sprintf("Disease severity (0-%g scale)", SEVERITY_MAX),
       colour = NULL, shape = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 2)) + theme_pub
save_figure(fig_progress, "F03_disease_progress_curves", width = 9, height = 6.5)


################################################################################
# STEP 6. APPARENT INFECTION RATE
################################################################################
banner("STEP 6. Apparent infection rate")

# 6.1 Van der Plank exponential model:
#       r = [ln(x2) - ln(x1)] / (t2 - t1)
#     x = plot-mean severity / SEVERITY_MAX (proportion), t = days after planting.
#     r is computed for every consecutive pair of assessments in which the plot
#     already shows disease (x1 > 0 and x2 > 0, because ln(0) is undefined), then
#     averaged per plot. Plots that never show disease have no defined r. The
#     logistic rate [logit(x2) - logit(x1)]/(t2 - t1) is given as a sensitivity check.
logit <- function(p) log(p / (1 - p))

rates_interval <- plot_time %>%
  mutate(x = pmin(Severity / SEVERITY_MAX, 0.999)) %>%
  arrange(Site, Block, Treatment, Disease, Day) %>%
  group_by(Site, Block, Treatment, Plot, Disease) %>%
  mutate(x1 = lag(x), t1 = lag(Day)) %>%
  filter(!is.na(x1), x1 > 0, x > 0) %>%
  mutate(r_exponential = (log(x) - log(x1)) / (Day - t1),
         r_logistic    = (logit(x) - logit(x1)) / (Day - t1)) %>% ungroup()

rates_plot <- rates_interval %>% group_by(Site, Block, Treatment, Disease) %>%
  summarise(r_exponential = mean(r_exponential), r_logistic = mean(r_logistic),
            n_intervals = n(), .groups = "drop")
save_table(rates_plot %>% mutate(across(where(is.double), ~ round(.x, 4))),
           "T08_apparent_infection_rate_per_plot")

# 6.2 Site summary and classification (low / moderate / high) -----------------
classify_rate <- function(r) cut(r, c(-Inf, RATE_LOW, RATE_HIGH, Inf), right = FALSE,
                                 labels = c("Low", "Moderate", "High"))
rates_site <- rates_plot %>% group_by(Site, Disease) %>%
  summarise(r_mean = mean(r_exponential), r_SE = se(r_exponential),
            r_logistic_mean = mean(r_logistic), n_plots = n(), .groups = "drop") %>%
  mutate(doubling_time_days = log(2) / r_mean,
         class = classify_rate(r_mean))

if (!is.null(REFERENCE_SITE)) {
  if (!REFERENCE_SITE %in% SITE_CODES) stop("REFERENCE_SITE '", REFERENCE_SITE, "' not found.")
  ref_label <- SITE_LABELS[match(REFERENCE_SITE, SITE_CODES)]
  ref <- rates_site %>% filter(Site == ref_label) %>% select(Disease, r_ref = r_mean)
  rates_site <- rates_site %>% left_join(ref, by = "Disease") %>%
    mutate(pct_difference_from_reference = 100 * (r_mean - r_ref) / r_ref) %>% select(-r_ref)
}

# 6.3 Do sites differ? One-way ANOVA on plot rates + Tukey HSD, per disease ---
rate_tests <- list(); rate_letters <- list()
for (dz in levels(rates_plot$Disease)) {
  d <- droplevels(filter(rates_plot, Disease == dz))
  if (nlevels(d$Site) < 2) { cat(sprintf("\n%s - fewer than 2 sites with a defined r; test skipped\n", dz)); next }
  m  <- aov(r_exponential ~ Site, data = d)
  s  <- summary(m)[[1]]
  tk <- TukeyHSD(m, "Site")$Site
  lt <- multcompView::multcompLetters(setNames(tk[, "p adj"], gsub(" ", "_", rownames(tk))),
                                      threshold = ALPHA)$Letters
  rate_tests[[dz]]   <- data.frame(Disease = dz, F = s[1, "F value"], df1 = s[1, "Df"],
                                   df2 = s[2, "Df"], P = s[1, "Pr(>F)"])
  rate_letters[[dz]] <- data.frame(Disease = dz, Site = gsub("_", " ", names(lt)), letter = unname(lt))
  cat(sprintf("\n%s - site effect on r: F(%d,%d) = %.2f, P = %s\n", dz,
              s[1, "Df"], s[2, "Df"], s[1, "F value"], fmt_p(s[1, "Pr(>F)"])))
  print(round(tk, 4))
}
rate_tests   <- bind_rows(rate_tests)
rate_letters <- bind_rows(rate_letters)
if (nrow(rate_letters) > 0) {
  rate_letters <- rate_letters %>%
    mutate(Site = factor(Site, levels = SITE_LABELS),
           Disease = factor(Disease, levels = levels(rates_plot$Disease)))
  rates_site <- rates_site %>% left_join(rate_letters, by = c("Site", "Disease"))
} else {
  rates_site$letter <- ""
}

print(as.data.frame(rates_site %>% mutate(across(where(is.double), ~ round(.x, 3)))))
save_table(rates_site %>% mutate(across(where(is.double), ~ round(.x, 4))),
           "T08b_apparent_infection_rate_by_site")
if (nrow(rate_tests) > 0)
  save_table(rate_tests %>% mutate(P = round(P, 4)), "T08c_infection_rate_site_ANOVA")

# 6.4 Apparent infection rate by site, letters from Tukey HSD -----------------
fig_rate <- ggplot(rates_site, aes(Site, r_mean, fill = Site)) +
  geom_hline(yintercept = c(RATE_LOW, RATE_HIGH), linetype = "dashed",
             colour = "grey45", linewidth = 0.4) +
  geom_col(width = 0.6, colour = "white") +
  geom_errorbar(aes(ymin = r_mean - r_SE, ymax = r_mean + r_SE), width = 0.15, linewidth = 0.4) +
  geom_text(aes(y = r_mean + ifelse(is.na(r_SE), 0, r_SE), label = letter),
            vjust = -0.6, family = FIG_FONT, size = 4.5) +
  geom_text(data = data.frame(y = c(RATE_LOW, RATE_HIGH),
                              label = c(sprintf("low / moderate (%g)", RATE_LOW),
                                        sprintf("moderate / high (%g)", RATE_HIGH))),
            aes(x = 0.45, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = -0.4, size = 3, colour = "grey30", family = FIG_FONT) +
  facet_wrap(~ Disease) +
  scale_fill_manual(values = SITE_COLOURS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = expression(bold("Apparent infection rate, r (day"^-1*")"))) +
  theme_pub
save_figure(fig_rate, "F04_apparent_infection_rate", width = 7.5, height = 4.5)
cat(sprintf(paste0("  Caption: bars = mean of plots with a defined rate (n = %s) +/- SE;\n",
                   "  different letters = significant difference between sites (Tukey HSD,\n",
                   "  P < %g) within each disease.\n"),
            paste(range(rates_site$n_plots), collapse = "-"), ALPHA))


################################################################################
# STEP 7. AUDPC AND TREATMENT EFFECTS
################################################################################
banner("STEP 7. AUDPC and treatment ANOVA")

# 7.1 AUDPC per plant with the trapezoidal rule:
#       AUDPC = sum[(y_i + y_i+1)/2 x (t_i+1 - t_i)]
#     Two equivalent scales are reported - state in the Methods which is used:
#       AUDPC_assess : t = assessment number (interval counted as 1)
#       AUDPC_days   : t = days after planting
#     A plant with a missing assessment gets AUDPC = NA (no interpolation).
audpc_plant <- disease_long %>% arrange(Time_Point) %>%
  group_by(Site, Block, Treatment, Plot, Sample_ID, Disease) %>%
  summarise(AUDPC_assess = sum((head(Severity, -1) + tail(Severity, -1)) / 2 * diff(Time_Point)),
            AUDPC_days   = sum((head(Severity, -1) + tail(Severity, -1)) / 2 * diff(Day)),
            .groups = "drop")
save_table(audpc_plant, "T09_AUDPC_per_plant")

# 7.2 Plot means (experimental unit) -------------------------------------------
audpc_plot <- audpc_plant %>% group_by(Site, Block, Treatment, Disease) %>%
  summarise(AUDPC_assess = mean(AUDPC_assess, na.rm = TRUE),
            AUDPC_days   = mean(AUDPC_days, na.rm = TRUE),
            n_plants     = sum(!is.na(AUDPC_assess)), .groups = "drop")
save_table(audpc_plot %>% mutate(across(where(is.double), ~ round(.x, 3))), "T09b_AUDPC_per_plot")

# 7.3 Site x treatment means (plant-level mean and plot-level SE) --------------
audpc_trt <- audpc_plant %>% group_by(Site, Treatment, Disease) %>%
  summarise(AUDPC_mean = mean(AUDPC_assess, na.rm = TRUE),
            AUDPC_SD   = sd(AUDPC_assess, na.rm = TRUE),
            n_plants   = sum(!is.na(AUDPC_assess)), .groups = "drop") %>%
  left_join(audpc_plot %>% group_by(Site, Treatment, Disease) %>%
              summarise(AUDPC_SE_plots = se(AUDPC_assess), .groups = "drop"),
            by = c("Site", "Treatment", "Disease")) %>%
  mutate(AUDPC_days_mean = AUDPC_mean * diff(ASSESSMENT_DAYS)[1]) %>%
  mutate(across(where(is.double), ~ round(.x, 2)))
print(as.data.frame(audpc_trt))
save_table(audpc_trt, "T09c_AUDPC_by_site_treatment")

audpc_site <- audpc_plant %>% group_by(Site, Disease) %>%
  summarise(AUDPC_mean = mean(AUDPC_assess, na.rm = TRUE), .groups = "drop") %>%
  left_join(audpc_plot %>% group_by(Site, Disease) %>%
              summarise(AUDPC_SE = se(AUDPC_assess), .groups = "drop"), by = c("Site", "Disease"))

# 7.4 ANOVA of treatment within each site and disease --------------------------
#     RCBD on PLOT means: AUDPC ~ Treatment + Block. With a balanced RCBD,
#     treating Block as random gives exactly the same F-test for Treatment.
aov_term <- function(fit, term) {   # extract F, df and P for one model term
  tab <- summary(fit)[[1]]; rownames(tab) <- trimws(rownames(tab))
  list(F = tab[term, "F value"], df = paste0(tab[term, "Df"], ",", tab["Residuals", "Df"]),
       P = tab[term, "Pr(>F)"])
}

anova_rows <- list()
for (s in levels(audpc_plot$Site)) for (dz in levels(audpc_plot$Disease)) {
  pd <- filter(audpc_plot, Site == s, Disease == dz)
  a1 <- aov_term(aov(AUDPC_assess ~ Treatment + Block, data = pd), "Treatment")
  anova_rows[[paste(s, dz)]] <- data.frame(Site = s, Disease = dz,
                                           F = a1$F, df = a1$df, P = a1$P)
}
anova_audpc <- bind_rows(anova_rows) %>%
  mutate(Sig = p_stars(P), across(where(is.double), ~ round(.x, 3)))
print(anova_audpc, row.names = FALSE)
save_table(anova_audpc, "T10_AUDPC_treatment_ANOVA")


################################################################################
# STEP 8. YIELD, YIELD LOSS AND ECONOMIC LOSS
################################################################################
banner("STEP 8. Yield, yield loss and economic loss")

# 8.1 Yield by site and by site x treatment --------------------------------------
yield_trt <- yield %>% group_by(Site, Treatment) %>%
  summarise(Yield_mean = mean(Yield), Yield_SE = se(Yield), n = n(), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 2)))
save_table(yield_trt, "T11_yield_by_site_treatment")

# 8.2 Treatment effect on yield within each site (RCBD) ------------------------
yield_anova <- bind_rows(lapply(levels(yield$Site), function(s) {
  a <- aov_term(aov(Yield ~ Treatment + Block, data = filter(yield, Site == s)), "Treatment")
  data.frame(Site = s, F = a$F, df = a$df, P = a$P)
})) %>% mutate(Sig = p_stars(P), across(where(is.double), ~ round(.x, 3)))
print(yield_anova, row.names = FALSE)
save_table(yield_anova, "T11b_yield_treatment_ANOVA")

# 8.3 Yield loss relative to the highest-yielding site ---------------------------
#     Yield loss (%)  = (Y_max - Y_site) / Y_max x 100
#     Economic loss   = (Y_max - Y_site) x price per tonne
yield_loss <- yield %>% group_by(Site) %>%
  summarise(Yield_t_ha = mean(Yield), Yield_SE = se(Yield), n = n(), .groups = "drop") %>%
  mutate(Yield_gap_t_ha = max(Yield_t_ha) - Yield_t_ha,
         Pct_loss       = 100 * Yield_gap_t_ha / max(Yield_t_ha)) %>%
  arrange(desc(Yield_t_ha))

if (HAS_PRICE) {
  yield_loss <- yield_loss %>%
    mutate(Price_per_t         = PRICE_PER_TONNE,
           Value_local_ha      = Yield_t_ha * PRICE_PER_TONNE,
           Economic_loss_local_ha = Yield_gap_t_ha * PRICE_PER_TONNE)
  if (HAS_USD) yield_loss <- yield_loss %>%
    mutate(Value_USD_ha         = Value_local_ha / LOCAL_PER_USD,
           Economic_loss_USD_ha = Economic_loss_local_ha / LOCAL_PER_USD)
}

yield_loss_print <- yield_loss %>% mutate(
  across(c(Yield_t_ha, Yield_SE, Yield_gap_t_ha), ~ round(.x, 2)),
  Pct_loss = round(Pct_loss, 1),
  across(any_of(c("Value_local_ha", "Economic_loss_local_ha", "Value_USD_ha", "Economic_loss_USD_ha")),
         ~ round(.x, 0)))
if (HAS_PRICE) names(yield_loss_print) <- sub("_local_", paste0("_", CURRENCY, "_"), names(yield_loss_print))
print(as.data.frame(yield_loss_print))
save_table(yield_loss_print, "T12_yield_and_economic_loss")

# 8.4 Plot-level data set linking disease, infection rate and yield -------------
plot_link <- audpc_plot %>% select(Site, Block, Treatment, Disease, AUDPC_assess) %>%
  left_join(rates_plot %>% select(Site, Block, Treatment, Disease, r_exponential),
            by = c("Site", "Block", "Treatment", "Disease")) %>%
  left_join(yield %>% select(Site, Block, Treatment, Yield), by = c("Site", "Block", "Treatment"))
save_table(plot_link %>% mutate(across(where(is.double), ~ round(.x, 4))),
           "T13_plot_level_disease_rate_yield")

# 8.5 Linear regressions: yield ~ AUDPC and yield ~ r -----------------------------
reg_one <- function(d, x, label) {
  d <- d[!is.na(d[[x]]) & !is.na(d$Yield), ]
  if (nrow(d) < 3) return(data.frame(Predictor = label, n = nrow(d), intercept = NA,
                                     slope = NA, R2 = NA, r = NA, P = NA))
  m <- lm(reformulate(x, "Yield"), data = d); s <- summary(m)
  ct <- cor.test(d[[x]], d$Yield)
  data.frame(Predictor = label, n = nrow(d),
             intercept = coef(m)[1], slope = coef(m)[2], R2 = s$r.squared,
             r = unname(ct$estimate), P = ct$p.value)
}
regs <- bind_rows(lapply(levels(plot_link$Disease), function(dz) {
  d <- filter(plot_link, Disease == dz)
  bind_rows(reg_one(d, "AUDPC_assess",  paste(dz, "AUDPC")),
            reg_one(d, "r_exponential", paste(dz, "apparent infection rate")))
})) %>% mutate(Sig = p_stars(P), across(where(is.double), ~ signif(.x, 4)))
print(regs, row.names = FALSE)
save_table(regs, "T13b_yield_regressions")
cat("  Note: infection-rate regressions use only plots with a defined r\n",
    "  (plots that stayed disease-free have no rate), see column n.\n")

# 8.6 AUDPC vs yield. Small points = plots; large points = site means +/- SE;
#     line = least-squares fit to the plot-level data.
site_means <- audpc_site %>%
  left_join(yield_loss %>% select(Site, Yield_t_ha, Yield_SE), by = "Site")
reg_lab <- regs %>% filter(grepl("AUDPC", Predictor)) %>%
  mutate(Disease = factor(sub(" AUDPC", "", Predictor), levels = levels(plot_link$Disease)),
         label = sprintf("r = %.2f, P %s\nn = %d plots", r,
                         ifelse(P < 0.001, "< 0.001", paste("=", sprintf("%.3f", P))), n))

fig_yield <- ggplot(plot_link, aes(AUDPC_assess, Yield)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey30",
              fill = "grey85", linewidth = 0.6) +
  geom_point(aes(colour = Site, shape = Site), size = 1.8, alpha = 0.55) +
  geom_errorbar(data = site_means, inherit.aes = FALSE, width = 0,
                aes(x = AUDPC_mean, ymin = Yield_t_ha - Yield_SE, ymax = Yield_t_ha + Yield_SE)) +
  geom_errorbarh(data = site_means, inherit.aes = FALSE, height = 0,
                 aes(y = Yield_t_ha, xmin = AUDPC_mean - AUDPC_SE, xmax = AUDPC_mean + AUDPC_SE)) +
  geom_point(data = site_means, aes(AUDPC_mean, Yield_t_ha, fill = Site, shape = Site),
             size = 4, colour = "black", stroke = 0.6) +
  geom_text(data = reg_lab, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.3,
            size = 3.4, family = FIG_FONT, inherit.aes = FALSE) +
  facet_wrap(~ Disease, scales = "free_x") +
  scale_colour_manual(values = SITE_COLOURS) +
  scale_fill_manual(values = SITE_COLOURS) +
  scale_shape_manual(values = recycle(FILL_SHAPES, length(SITE_LABELS))) +
  labs(x = "AUDPC (severity x assessment interval)", y = "Tuber yield (t/ha)",
       colour = "Site", fill = "Site", shape = "Site") + theme_pub
save_figure(fig_yield, "F05_AUDPC_vs_yield", width = 8, height = 4.8)

# 8.7 Yield loss relative to the highest-yielding site -------------------------
loss_data <- yield_loss %>% filter(Pct_loss > 0)
if (nrow(loss_data) > 0) {
  loss_data <- loss_data %>%
    mutate(label = if (HAS_PRICE)
      sprintf("%.1f %%\n%s %s/ha", Pct_loss, CURRENCY,
              format(round(Economic_loss_local_ha), big.mark = ",", trim = TRUE))
      else sprintf("%.1f %%", Pct_loss))
  fig_loss <- ggplot(loss_data, aes(Site, Pct_loss, fill = Site)) +
    geom_col(width = 0.55, colour = "white") +
    geom_text(aes(label = label), vjust = -0.3, family = FIG_FONT, size = 3.6) +
    scale_fill_manual(values = SITE_COLOURS, guide = "none") +
    scale_y_continuous(limits = c(0, max(loss_data$Pct_loss) * 1.3),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = sprintf("Yield loss relative to %s (%%)", yield_loss$Site[1])) + theme_pub
  save_figure(fig_loss, "F06_yield_loss", width = 5.5, height = 4.5)
}


################################################################################
# STEP 9. SESSION INFORMATION
################################################################################
banner("STEP 9. Session information")
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
cat("  -> session_info.txt written (R and package versions used)\n")

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Tables : %d CSV files in %s\n", length(list.files(TABLE_DIR, "\\.csv$")), TABLE_DIR))
cat(sprintf("Figures: %d PNG + %d TIFF in %s\n", length(list.files(FIGURE_DIR, "\\.png$")),
            length(list.files(FIGURE_DIR, "\\.tiff$")), FIGURE_DIR))
sink(); close(log_con)
cat("Log saved to:", log_file, "\n")

################################################################################
# END OF SCRIPT
################################################################################
