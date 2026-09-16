# SARS-CoV-2 reinfection detection algorithm

Companion code for:

> Vermeulen M, Swanevelder R, Mmenu C, Brits T, Mitchel J, Grebe E. Rates of
> SARS-CoV-2 infection in a highly exposed South African blood donor cohort,
> 2022-2023. Vox Sang. Forthcoming 2026.

`reinfection_algorithm.R` is the custom routine used to identify first
SARS-CoV-2 infections and reinfections from longitudinal anti-nucleocapsid
(anti-N) S/CO results in repeat blood donors, and to estimate first-infection,
reinfection and total infection incidence with exact Poisson 95% confidence
intervals.

The file is self-contained (R version 4.1 or later, with dplyr, tidyr,
lubridate and tibble); the published analysis was run under R version
4.5.1-3. It takes a data frame with one row per tested donation
(`donor_urn`, `donation_date`, `anti_n_sco`) and returns the donation-level
classification, the detected infection events, and the incidence estimates.
The input must contain exactly one row per donation: repeated test records for
the same donation would be read as separate donations. See the comments at the
top of the file for the algorithm and the expected input.

Licensed under GPL-3.0-or-later (see `LICENSE`).
