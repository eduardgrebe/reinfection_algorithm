# SARS-CoV-2 reinfection detection algorithm
# Copyright (C) 2025-2026 Eduard Grebe Consulting and South African National
# Blood Service. Author: Eduard Grebe <eduard@grebe.consulting>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Large language model usage disclosure: The original analytic code, from
# which this routine is derived, was written by Eduard Grebe without AI
# assistance. AI models were used to prepare a more readable script for
# public release, to identify and correct defects in that script, and to
# verify that it reproduces the original classifications and incidence
# estimates exactly when run against the analytic dataset used for the
# manuscript. Models used: DeepSeek-V4.1-Flash and Anthropic Claude Opus 5
# (claude-opus-5). All model-generated changes were reviewed by Eduard Grebe.

# ============================================================================
# Custom routine used to identify first infections and reinfections in the
# SANBS repeat blood donor cohort (South Africa, 2022-2023), and to estimate
# infection incidence with exact Poisson confidence intervals.
#
# Companion code for:
#   Vermeulen M, Swanevelder R, Mmenu C, Brits T, Mitchel J, Grebe E. Rates
#   of SARS-CoV-2 infection in a highly exposed South African blood donor
#   cohort, 2022-2023. Vox Sang. Forthcoming 2026.
#
# How the algorithm works
# -----------------------
# The routine works on longitudinal anti-nucleocapsid (anti-N) S/CO results
# from repeat donors, in which infection is detected as a rise in anti-N
# reactivity between consecutive donations (a reactive sample is one with an
# anti-N S/CO >= 1.0):
#
#   * A first infection is called at a seroconversion: a reactive donation
#     preceded by a non-reactive donation. Donors whose first tested donation
#     is already reactive ("entered positive") are assumed to have been
#     infected before enrolment; they are not at risk of a first infection.
#   * A reinfection is called when a reactive donation in a donor who is
#     already infected (i.e. not a seroconversion) has an anti-N S/CO at least
#     1.5-fold the previous anti-N S/CO. The denominator of the ratio is
#     floored at the reactivity cut-off, so that a rise from a near-zero
#     non-reactive result cannot by itself exceed the threshold. Because
#     continuing antibody maturation after a single infection commonly produces
#     a series of rising results, consecutive boosting donations are treated as
#     one infection event: a reinfection is only called once reactivity has
#     stabilised or declined since the previous boosting donation.
#   * A donor who seroreverts remains "previously infected": a later reactive
#     donation is not called as a second first infection.
#   * The date of infection (first infection or reinfection) is estimated as
#     the midpoint of the interdonation interval in which it was detected.
#
# Incidence is then estimated as infection events / time at risk, separately
# for first infections and for reinfections, and for all infections combined.
# Donors are considered at risk of a first infection if they were never
# previously infected; first infections are censoring events. Donors are
# considered at risk of reinfection from their estimated first infection date
# (or from their first tested donation, if they entered the study already
# infected); reinfections are not censoring events, since a donor remains at
# risk after an infection. Donor series are right-censored by the last tested
# donation sample. Donors contribute whole interdonation intervals during which
# no infection occurred and half of the interdonation interval during which an
# infection was detected, since infection time is presumed to be at the
# midpoint of that interval. Only donors with at least two tested donations
# (so that time at risk is observed) contribute time at risk. Confidence
# intervals are exact Poisson intervals, using the relationship between the
# Poisson and chi-squared distributions.
#
# Requirements: R (>= 4.1) with dplyr, tidyr, lubridate and tibble. The
# published analysis was run under R version 4.5.1-3.
#
# Input
# -----
# `detect_infections()` takes a data frame with one row per donation:
#   donor_urn      donor identifier
#   donation_date  date of the donation (Date)
#   anti_n_sco     anti-N S/CO result of the tested sample (numeric, NA if the
#                  sample was not tested for anti-N, in which case the
#                  donation is ignored by the routine)
#
# There must be exactly one row per donation: repeated test records for the
# same donation would be read as separate donations.
#
# Output
# ------
# detect_infections() returns a list with
#   $donations     donation-level data with the classification
#   $infections    one row per first infection / reinfection detected
# estimate_incidence() returns a data frame with incidence estimates
# ============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(tibble)

# ----------------------------------------------------------------------------
# 1. Identify first infections and reinfections
# ----------------------------------------------------------------------------
#
# `reactive_cutoff` is the anti-N S/CO at or above which a sample is reactive
# and `boost_ratio_cutoff` is the fold-rise in anti-N S/CO over the previous
# sample required to call a reinfection.
detect_infections <- function(donations,
                              reactive_cutoff = 1.0,
                              boost_ratio_cutoff = 1.5) {

  if (!is.data.frame(donations)) {
    stop("`donations` must be a data frame")
  }
  required <- c("donor_urn", "donation_date", "anti_n_sco")
  if (length(setdiff(required, names(donations))) > 0) {
    stop("`donations` must contain the columns: ", paste(required, collapse = ", "))
  }
  if (!inherits(donations$donation_date, "Date")) {
    stop("`donation_date` must be a Date (see ?as.Date)")
  }

  # Anti-N tested donations, in donation order within each donor -------------
  dat <- donations |>
    filter(!is.na(anti_n_sco)) |>
    mutate(
      anti_n_interp = if_else(anti_n_sco >= reactive_cutoff, "reactive", "non-reactive")
    ) |>
    group_by(donor_urn) |>
    arrange(donation_date) |>
    mutate(
      prev_n_interp = lag(anti_n_interp, 1),
      prev_n_sco = lag(anti_n_sco, 1),
      entered_pos = first(anti_n_interp) == "reactive",
      # Serological status of each donation relative to the previous donation:
      #   not_infected   non-reactive, no previous reactive result
      #   seroconversion first reactive result (a first infection)
      #   infected       reactive, previous result also reactive (potentially
      #                  a reinfection, see the boosting rule below)
      #   seroreversion  non-reactive after a reactive result
      sero_status = case_when(
        is.na(prev_n_interp) & anti_n_interp == "reactive" ~ "infected",
        is.na(prev_n_interp) ~ "not_infected",
        anti_n_interp == "reactive" & prev_n_interp == "non-reactive" ~ "seroconversion",
        anti_n_interp == "reactive" & prev_n_interp == "reactive" ~ "infected",
        anti_n_interp == "non-reactive" & prev_n_interp == "reactive" ~ "seroreversion",
        TRUE ~ "not_infected"
      ),
      # Rise in anti-N S/CO relative to the previous donation, with the
      # denominator floored at the reactivity cut-off.
      n_ratio_2 = anti_n_sco / ifelse(prev_n_sco < reactive_cutoff, reactive_cutoff, prev_n_sco)
    ) |>
    ungroup()

  # After a seroreversion a donor is still treated as previously infected, so
  # that a later reactive result is not called as a second first infection.
  seroreversion_dates <- dat |>
    filter(sero_status == "seroreversion") |>
    group_by(donor_urn) |>
    arrange(donation_date) |>
    slice(1) |>
    ungroup() |>
    select(donor_urn, seroreversion_date = donation_date)

  dat <- dat |>
    left_join(seroreversion_dates, by = "donor_urn") |>
    mutate(
      sero_status = case_when(
        !is.na(seroreversion_date) & donation_date > seroreversion_date ~ "infected",
        TRUE ~ sero_status
      )
    )

  # Reinfection: a >= boost_ratio_cutoff-fold rise in anti-N S/CO in a
  # previously infected donor. Only the first boosting donation of a run of
  # consecutive boosting donations is counted, as a single infection event.
  dat <- dat |>
    mutate(
      n_boost = case_when(
        sero_status != "infected" ~ FALSE,
        n_ratio_2 < boost_ratio_cutoff ~ FALSE,
        n_ratio_2 >= boost_ratio_cutoff ~ TRUE
      )
    ) |>
    group_by(donor_urn) |>
    arrange(donation_date) |>
    mutate(
      prev_boost = lag(n_boost, 1),
      reinfection = case_when(
        n_boost == FALSE ~ FALSE,
        n_boost == TRUE & prev_boost == TRUE ~ FALSE,
        n_boost == TRUE ~ TRUE
      )
    ) |>
    ungroup()

  # Date of infection: midpoint of the interdonation interval in which the
  # infection was detected.
  dat <- dat |>
    group_by(donor_urn) |>
    arrange(donation_date) |>
    mutate(
      idi = interval(lag(donation_date, 1), donation_date) / ddays(1),
      infection_date = case_when(
        sero_status == "seroconversion" ~ lag(donation_date, 1) + idi / 2,
        reinfection == TRUE ~ lag(donation_date, 1) + idi / 2,
        TRUE ~ as.Date(NA)
      )
    ) |>
    ungroup()

  # Infection events, numbered in order of estimated infection date. Donors who
  # were already infected at enrolment start at 2, since their first infection
  # occurred before the study period.
  infections <- dat |>
    filter(sero_status == "seroconversion" | reinfection == TRUE) |>
    mutate(
      type = if_else(sero_status == "seroconversion", "First infection", "Reinfection")
    ) |>
    group_by(donor_urn) |>
    arrange(infection_date) |>
    mutate(
      infection_number = seq_len(n()) + as.integer(entered_pos)
    ) |>
    ungroup() |>
    select(donor_urn, donation_date, infection_date, infection_number, type)

  # Infection history in wide format (inf1, inf2, ...), used to define time at
  # risk. Up to four infection events per donor (the maximum observed in the
  # cohort); further events would not be used for time at risk.
  if (nrow(infections) > 0 && max(infections$infection_number) > 4) {
    warning("Some donors have more than 4 infection events; only the first four are used.")
  }

  infections_wide <- infections |>
    mutate(inf = paste0("inf", infection_number)) |>
    select(donor_urn, inf, infection_date) |>
    pivot_wider(names_from = inf, values_from = infection_date)

  for (k in 1:4) {
    column <- paste0("inf", k)
    if (!column %in% names(infections_wide)) {
      infections_wide[[column]] <- as.Date(NA)
    }
  }

  dat <- dat |>
    left_join(select(infections_wide, donor_urn, inf1, inf2, inf3, inf4), by = "donor_urn")

  list(donations = dat, infections = infections)
}

# ----------------------------------------------------------------------------
# 2. Time at risk, and number of events, for first infections and reinfections
# ----------------------------------------------------------------------------
#
# Returns a list of two donor-level data frames (first_infection,
# reinfection), each with one row per donor at risk and the number of events
# and days at risk contributed by that donor.
time_at_risk <- function(classified) {

  dat <- classified$donations |>
    # Time at risk must have been observed: at least two tested donations.
    # (A donor with a single tested donation contributes neither time at risk
    # nor events, so this does not change the estimates.)
    group_by(donor_urn) |>
    filter(n() >= 2) |>
    ungroup()

  # First infections: donors who were not previously infected contribute time
  # from their first tested donation until their estimated first infection date
  # (first infections are censoring events), or until their last tested
  # donation if they were not infected during follow-up.
  first_infection <- dat |>
    filter(!entered_pos & (is.na(inf1) | donation_date < inf1)) |>
    group_by(donor_urn) |>
    arrange(donation_date) |>
    summarise(
      events = if_else(!is.na(first(inf1)), 1, 0),
      time_at_risk = if_else(
        is.na(first(inf1)),
        interval(first(donation_date), last(donation_date)) / ddays(1),
        interval(first(donation_date), first(inf1)) / ddays(1)
      )
    ) |>
    ungroup()

  # Reinfections: donors contribute time at risk from their first tested
  # donation if they entered the study already infected, or from their
  # estimated first infection date if they were infected during follow-up,
  # until their last tested donation. Reinfections are not censoring events.
  reinfection <- dat |>
    filter(entered_pos | (!is.na(inf1) & donation_date > inf1)) |>
    group_by(donor_urn) |>
    arrange(donation_date) |>
    summarise(
      events = case_when(
        !is.na(first(inf4)) ~ 3,
        !is.na(first(inf3)) ~ 2,
        !is.na(first(inf2)) ~ 1,
        TRUE ~ 0
      ),
      time_at_risk = if_else(
        first(entered_pos),
        interval(first(donation_date), last(donation_date)) / ddays(1),
        interval(first(inf1), last(donation_date)) / ddays(1)
      )
    ) |>
    ungroup()

  list(first_infection = first_infection, reinfection = reinfection)
}

# ----------------------------------------------------------------------------
# 3. Incidence and exact Poisson confidence interval
# ----------------------------------------------------------------------------
#
# Incidence for a count of events Y and a total time at risk Texp ('exposure',
# in days by default). The confidence interval is derived from the properties
# of the Poisson and chi-squared distributions:
#   lower = qchisq(alpha / 2, 2 * Y) / 2
#   upper = qchisq(1 - alpha / 2, 2 * (Y + 1)) / 2
incidence_n <- function(Y, Texp, unit_conversion = 365.25, per = 100, alpha = 0.05) {
  Texp <- Texp / unit_conversion
  Yl <- qchisq(alpha / 2, 2 * Y) / 2
  Yu <- qchisq(1 - alpha / 2, 2 * (Y + 1)) / 2
  c(
    I = Y / Texp * per,
    CI_lower = Yl / Texp * per,
    CI_upper = Yu / Texp * per
  )
}

# ----------------------------------------------------------------------------
# 4. Incidence for first infections, reinfections and all infections combined
# ----------------------------------------------------------------------------
estimate_incidence <- function(classified, unit_conversion = 365.25,
                               per = 100, alpha = 0.05) {

  at_risk <- time_at_risk(classified)

  one_row <- function(type, events, time_at_risk) {
    estimate <- incidence_n(events, time_at_risk,
                            unit_conversion = unit_conversion,
                            per = per, alpha = alpha)
    tibble(
      infection_type = type,
      events = events,
      person_time = time_at_risk / unit_conversion,
      incidence = estimate[["I"]],
      ci_lower = estimate[["CI_lower"]],
      ci_upper = estimate[["CI_upper"]]
    )
  }

  first <- at_risk$first_infection
  reinf <- at_risk$reinfection

  bind_rows(
    one_row("First infection", sum(first$events), sum(first$time_at_risk)),
    one_row("Reinfection", sum(reinf$events), sum(reinf$time_at_risk)),
    one_row("Total", sum(first$events) + sum(reinf$events),
            sum(first$time_at_risk) + sum(reinf$time_at_risk))
  )
}

# ----------------------------------------------------------------------------
# Example
# ----------------------------------------------------------------------------
# donations <- read.csv("donations.csv")
# donations$donation_date <- as.Date(donations$donation_date)
#
# classified <- detect_infections(donations)
# classified$infections                     # first infections and reinfections
# estimate_incidence(classified)            # incidence per 100 person-years
