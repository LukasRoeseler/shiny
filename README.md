# Shiny Apps (Static)

Static HTML versions of Lukas Röseler's interactive R Shiny apps, hosted on GitHub Pages. These were originally hosted on `shinyapps.io` and have been migrated to static sites that require no server.

## Migration: old Shiny pages → new static pages

| App | Old (shinyapps.io) | New (GitHub Pages) |
|-----|--------------------|--------------------|
| OpAQ — Open Anchoring Quest | <https://metaanalyses.shinyapps.io/OpAQ/> | <https://lukasroeseler.github.io/shiny/opaq/> |
| Toolbox — Trustworthiness of Published Findings | <https://metaanalyses.shinyapps.io/toolbox/> | <https://lukasroeseler.github.io/shiny/toolbox/> |
| Dynamic Meta-Analysis of Body Positions | <https://metaanalyses.shinyapps.io/bodypositions/> | <https://lukasroeseler.github.io/shiny/bodypositions/> |
| Studienfeedback (RHO-Daten) — *not converted* | <https://l-air.shinyapps.io/feedback/> | — (archived, see note) |

## Notes about function

Each app below is a static reimplementation of the original interactive Shiny app. Analyses that required the R server backend are precomputed and embedded; straightforward computations are done in client-side JavaScript.

| App | Function (what it does) | Migration notes |
|-----|-------------------------|-----------------|
| **OpAQ** | Interactive meta-analysis of anchoring effects: dataset overview, effect-size distribution, moderators (e.g., DV type, design, culture), reliabilities, anchor-extremeness, and file-drawer analyses. | Full static dashboard preserved (data precomputed into `data/*.json`). Original analysis logic and figures intact. |
| **Toolbox** | Determines the trustworthiness of a set of published findings. Upload a dataset of test statistics and it computes effect sizes/z-values and runs a **p-curve**, **z-curve**, **caliper test**, and **relative proximity** (significance tests + mediation CIs) analysis. | Effect-size conversion, z-values, caliper test and relative proximity run in client-side JS (verified against R). **p-curve and z-curve are documented only** — they require the R backend and are not reimplemented in JS. |
| **Body Positions** | Dynamic meta-analysis of experimentally induced body position effects on **behavior**, **self-report**, and **physiology** outcomes, with forest/funnel plots, p-curve, z-curve, and Egger's test. | Static results dashboard with the precomputed REML meta-analysis (overall, by DV type, by preregistration) and a forest-plot visualization. The interactive filtering/adding of studies from the original app is not available. |
| **Studienfeedback** | German-language participant feedback app ("RHO-Daten: Ihr Studienfeedback"). Participants enter a personal code to view their individual results vs. other participants (thinking styles/narcissism, and opinion formation on YouTube). | **Not converted** to static HTML (per migration request). Source code was not located in public repos/OSF; archived in the parent `shinymigration/feedback` folder. |

## Hosting on GitHub Pages

This repository root is published with GitHub Pages (deploy from the `main` branch, root `/` directory). Each app lives in its own subdirectory so it gets its own URL:

- `https://lukasroeseler.github.io/shiny/`
- `https://lukasroeseler.github.io/shiny/opaq/`
- `https://lukasroeseler.github.io/shiny/toolbox/`
- `https://lukasroeseler.github.io/shiny/bodypositions/`

### To enable
1. Push this folder's contents to a repository named **`shiny`** under the `LukasRoeseler` account.
2. In the repo **Settings → Pages**, set the source to *Deploy from a branch*, branch `main`, directory `/ (root)`.
3. The apps will be served at the URLs above.

## About the migration

The original Shiny apps (and their source code, data, and OSF projects) are archived in the parent `shinymigration` folder. Each static app here reproduces the original app's core functionality with precomputed analysis results (generated in R) and client-side JavaScript where the analyses are straightforward (e.g., caliper test, relative proximity).

> **Note:** The `feedback` (Studienfeedback) app on `l-air.shinyapps.io` was **not** converted to static HTML, per the migration request. Its source and content are noted in the parent folder.
