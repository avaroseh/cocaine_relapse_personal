# COCAINE BRAIN CONCENTRATION PHARMACOKINETIC MODEL
  # Author: Ava Rose Holmes

# This script simulates brain cocaine concentration over time
# during self-administration behavioral sessions using a
# two-compartment pharmacokinetic (PK) model.

# The model treats each cocaine infusion as an instantaneous
# IV bolus into a central compartment (blood), with drug
# exchange between the central and brain compartments and
# first-order elimination from the central compartment.

# Behavioral timestamps are used to map infusion
# timing for each animal and session, and the model generates
# a continuous brain concentration time course across sessions.

# From these simulations, the script calculates PK
# exposure metrics including:
#   *Cmax (maximum brain concentration)
#   *Tmax (time of maximum concentration)
#   *AUC (area under the concentration curve)

# What this project demonstrates:
#   *Data cleaning and preprocessing
#   *Event-based ODE simulation
#   *Pharmacokinetic modeling
#   *Time-series analysis
#   *Data visualization

# --------------------------
# LOAD PACKAGES
library(tidyverse)
library(deSolve)
library(pracma)
library(ggplot2)
library(viridis)
library(stringr)


# --------------------------
# SET PARAMETERS
# pk_params are taken from Pan et al. (1991) and D'Ottavio et al. (2025)
infusion_label <- "infusion"
sim_dt_min <- 0.5
session_length_min <- 180
MW <- 339.81        # cocaine HCl MW; not critical for RELATIVE units
dose_mgkg_inf <- 0.75
exclude_sessions <- c("S21")
pk_params <- list(
  Kel = 0.090295,   # elimination from central (min^-1)
  K12 = 0.067599,   # central -> brain (min^-1)
  K21 = 0.047336,   # brain -> central (min^-1)
  V1  = 0.122182    # central volume (L) - scaling for relative concentration
)


# --------------------------
# INPUT NECESSARY FILES
# df should include rat ID, measure names (infusions, active lever presses, etc.), 
# timestamp of behavior
# If not included in df, include a table with rat weights (Our weights are taken 
# from S01 of acquisition training)
mRFP_timestamps <- "/Users/holmesar/Desktop/rcode_012226/mRFP_timestamps.RDS"
df <- readRDS(mRFP_timestamps)

rat_info <- read.csv("/Users/holmesar/Desktop/rat_info.csv",
                     stringsAsFactors = FALSE)

if (!all(c("rat", "weight_kg") %in% names(rat_info))) {
  stop("rat_info.csv must have columns: rat, weight_kg")
}


# --------------------------
# NORMALIZE DF AND ADD WEIGHTS
df <- readRDS(mRFP_timestamps) %>%
  mutate(
    rat = str_squish(as.character(rat)),
    session = str_to_upper(str_squish(as.character(session))),
    measure_norm = str_to_lower(str_squish(measure_name)),
    time_min = timestamp / 60   # change to timestamp if your timestamps are already minutes
  ) %>%
  left_join(rat_info, by = "rat") %>%
  mutate(
    dose_mg = dose_mgkg_inf * weight_kg,
    dose_umol = ((dose_mg * 1e-3) / MW) * 1e6
  )

# Observed sessions (includes zero-infusions)
observed_sessions <- df %>%
  distinct(rat, session, weight_kg) %>%
  filter(!session %in% exclude_sessions)

# Infusion rows (used to create dosing events)
infusions <- df %>%
  filter(!session %in% exclude_sessions) %>%
  filter(measure_norm == infusion_label) %>%
  arrange(rat, session, time_min)

# Infusion counts per rat-session (pad zeros for sessions w/o infusions)
inf_counts <- infusions %>%
  count(rat, session, name = "n_infusions") %>%
  right_join(observed_sessions %>% distinct(rat, session), by = c("rat", "session")) %>%
  mutate(n_infusions = replace_na(n_infusions, 0L),
         has_infusions = n_infusions > 0)

# Build rat_sessions from observed_sessions (keeps weight_kg)
rat_sessions <- observed_sessions %>%
  mutate(
    session_num = as.integer(str_remove(session, "^S")),
    day = ceiling(session_num / 2),
    tod = if_else(session_num %% 2 == 1, "AM", "PM")
  ) %>%
  arrange(rat, session_num)

# --------------------------
# ODE: two-compartment IV bolus (CENTRAL <-> BRAIN; elimination from CENTRAL)
# States are AMOUNTS (same units as dose_col)
# Compartments:
#   CENTRAL : well-mixed blood/central compartment
#   BRAIN   : brain extracellular fluid compartment
#
# Assumptions:
#   Drug is delivered as instantaneous IV boluses into CENTRAL
#   Drug exchanges between CENTRAL and BRAIN with first-order kinetics
#   Drug is eliminated from CENTRAL with first-order kinetics
#
# State variables are AMOUNTS (not concentrations)
# Units of amount are the same as the dose input (dose_col)
pk_ode_iv <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    # Rate of change of CENTRAL amount:
    #   − Kel * CENTRAL     : elimination from central compartment
    #   − K12 * CENTRAL     : transfer from central → brain
    #   + K21 * BRAIN       : transfer from brain → central
    dCENTRAL <- -(Kel + K12) * CENTRAL + K21 * BRAIN
    # Rate of change of BRAIN amount:
    #   + K12 * CENTRAL     : transfer from central → brain
    #   − K21 * BRAIN       : transfer from brain → central
    dBRAIN   <-  K12 * CENTRAL - K21 * BRAIN
    list(c(dCENTRAL, dBRAIN))
  })
}


# --------------------------
# SIMULATE SESSION: add doses directly to CENTRAL at infusion times
# At time time_min, add dose into the CENTRAL compartment
# akes infusion timestamps for a rat/session, treats each infusion 
# as an instantaneous dose added to CENTRAL, solves the PK ODEs across 180 minutes, 
# and returns a brain-level time course.
simulate_session <- function(inf_sess, params,
                             session_len = 180, dt = 0.5,
                             dose_col = c("dose_umol", "dose_mg")) {
  
  dose_col <- match.arg(dose_col)
  tvec <- seq(0, session_len, by = dt)
  
  # If no infusions, return zero trace (we will NA it later if desired)
  if (nrow(inf_sess) == 0) {
    return(tibble(time = tvec, CENTRAL = 0, BRAIN = 0) %>%
             mutate(time_panel = time,
                    CP_BRAIN = 0))
  }
  
  # Build events: add dose amount directly to CENTRAL
  ev <- inf_sess %>%
    transmute(
      var = "CENTRAL",
      time = pmin(time_min, session_len),
      value = .data[[dose_col]],
      method = "add"
    ) %>%
    as.data.frame()
  
  out <- ode(
    y = c(CENTRAL = 0, BRAIN = 0),
    times = tvec,
    func = pk_ode_iv,
    parms = c(Kel = params$Kel, K12 = params$K12, K21 = params$K21, V1 = params$V1),
    events = list(data = ev),
    method = "lsoda"
  ) %>%
    as.data.frame() %>%
    as_tibble()
  
  out %>%
    mutate(time_panel = time,
           CP_BRAIN = BRAIN / params$V1)  # relative brain level (a.u.)
}

# --------------------------
# Run simulations for every observed rat-session (including zero-infusions)
# Note: select only the columns your pmap function expects to avoid unused-argument errors
# Block repeats that for every rat/session and combines the results
sim_df_day <- purrr::pmap_dfr(
  rat_sessions %>% select(rat, session, weight_kg, session_num, day, tod),
  function(rat, session, weight_kg, session_num, day, tod) {
    
    inf_sess <- infusions %>% filter(rat == !!rat, session == !!session)
    
    simulate_session(
      inf_sess = inf_sess,
      params = pk_params,
      session_len = session_length_min,
      dt = sim_dt_min,
      dose_col = "dose_umol"   # or "dose_mg" for relative-only units
    ) %>%
      mutate(rat = rat, session = session, session_num = session_num, day = day, tod = tod)
  }
) %>%
  mutate(tod = factor(tod, levels = c("AM", "PM"))) %>%
  left_join(inf_counts %>% select(rat, session, has_infusions), by = c("rat", "session")) %>%
  mutate(
    has_infusions = replace_na(has_infusions, FALSE),
    CP_BRAIN = if_else(has_infusions, CP_BRAIN, NA_real_)
  )


# --------------------------
# QUICK CHECKS
inf_counts %>% count(has_infusions)
sim_df_day %>% distinct(rat, session, has_infusions) %>% arrange(rat, session) %>% head(50)


# --------------------------
# SUMMARY
summary_df_day <- sim_df_day %>%
  group_by(rat, day, tod, session) %>%
  summarise(
    Cmax = max(CP_BRAIN, na.rm = TRUE),
    Tmax_min = time_panel[which.max(CP_BRAIN)],
    AUC = pracma::trapz(time_panel, CP_BRAIN),
    .groups = "drop"
  )

cat("Simulation finished. Created: sim_df_day, summary_df_day\n",
    "Rows - sim_df_day:", nrow(sim_df_day),
    " | summary_df_day:", nrow(summary_df_day), "\n")


# --------------------------
# HEATMAP STACKED BY DAY (AM/BREAK/PM columns)
plot_heatmap_by_day <- function(df_day, rat_id, dt = 0.5,
                                days = NULL,
                                x_breaks = c(0, 60, 120, 180)) {
  plot_df <- df_day %>%
    filter(rat == rat_id) %>%
    filter(is.finite(time_panel), is.finite(CP_BRAIN)) %>%
    mutate(
      time_panel = round(time_panel / dt) * dt,
      day = factor(day)
    )
  
  if (!is.null(days)) plot_df <- plot_df %>% filter(as.integer(as.character(day)) %in% days)
  
  ggplot(plot_df, aes(time_panel, y = 1, fill = CP_BRAIN)) +
    geom_raster() +
    facet_grid(rows = vars(day), cols = vars(tod), scales = "free_x") +
    scale_fill_viridis_c(option = "magma", na.value = "white", name = "CP_BRAIN (µM)") +
    scale_x_continuous(breaks = x_breaks) +
    labs(title = paste0(rat_id, " — heatmap stacked by day (AM / PM)"),
         x = "Time in segment (min)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

# --------------------------
# EXAMPLE PLOTS
plot_heatmap_by_day(sim_df_day, rat_id = "mRFP05R01")
plot_heatmap_by_day(sim_df_day, rat_id = "mRFP05R01", days = 1:10)

