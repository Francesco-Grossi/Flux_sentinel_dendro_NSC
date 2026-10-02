# Dendrometer data at flux sites - what is on disk and what is not

Downloaded on 2026-10-02 into `data_raw/dendro/` (not tracked by git; re-download from the sources
below). Nothing here is used by the pipeline yet.

## Downloaded

| Folder | Source and licence | What it is | Flux sites it covers | In the pipeline? |
|---|---|---|---|---|
| `dendrought2018/` | Salomón et al. 2022, Zenodo record 5711706, CC-BY 4.0 | Daily stem growth (`gr_dmax`) and tree water deficit from point and band dendrometers; 77 sites, 429 trees, 2015-2018 | CH-Dav (9 trees, 2015-2018), CZ-BK1 (8, 2016-2018), CZ-RAJ (7, 2017-2018), CZ-Stn (8, 2017-2018), FR-Fon (12, 2015-2018) | CH-Dav, CZ-BK1, CZ-RAJ, CZ-Stn yes; FR-Fon fails the flux QC |
| `harvard_forest/hf069-*` | Harvard Forest Data Archive HF069, CC0 | Dendrometer bands on all trees of the EMS tower plots: DBH and biomass per tree, 1993-2024, 2-7 readings per year | US-Ha1 | yes (1991-2025) |
| `harvard_forest/hf149-*` | Harvard Forest Data Archive HF149, CC0 | Dendrometer bands at the Hemlock and Little Prospect Hill towers, since 2001 | US-Ha2, US-LPH | no (not QC-passing) |
| `zoebelboden_lter/` | LTER Zöbelboden, Zenodo record 22940658, CC-BY 4.0 | Monthly stem, branch and foliage growth of spruce and beech from weekly dendrometer readings, 3 plots, 1996-2019 | AT-Zoe | yes (2014-2024) |
| `hyytiala_liu2023/` | Liu et al., Zenodo record 10037224, CC-BY 4.0 | 30-minute stem size and sap flow of two Scots pines, 2015-2019 | FI-Hyy | yes (1997-2025) |
| `automated_bands_four_forests/` | Zenodo record 4944203, CC0 | 15-minute automated dendrometer bands, 12-40 trees per site, from 2014 | SCBI, SERC, Wind River (US-Wrc), BCI | no (none of these is in the pipeline) |

`harvard_forest/` also holds the other HF069 tables (litter, soil respiration, soil water) that came
with the same record.

## Usable now: overlap with the pipeline

| Flux site | Dataset | Years | Time step |
|---|---|---|---|
| US-Ha1 | HF069 | 1993-2024 | 2-7 readings per year |
| AT-Zoe | Zöbelboden LTER | 2014-2019 (overlap) | monthly |
| CH-Dav | DenDrought2018 | 2015-2018 | daily |
| FI-Hyy | Liu et al. | 2015-2019 | 30 minutes, 2 trees |
| CZ-BK1 | DenDrought2018 | 2016-2018 | daily |
| CZ-RAJ, CZ-Stn | DenDrought2018 | 2017-2018 | daily |

That is 7 sites and roughly 55 site-years, of which about 20 have daily or finer resolution. The
within-site models of the pipeline need at least 3 years per site, so only US-Ha1, AT-Zoe, CH-Dav,
FI-Hyy and CZ-BK1 qualify.

## Not downloaded

| Dataset | Why not | How to get it |
|---|---|---|
| Rao et al. 2026, oaks incl. US-Ton (Dryad, doi:10.5061/dryad.m63xsj4h4, 384 MB, CC0) | The Dryad download API refused an anonymous request (HTTP 401) | Download in a browser from the Dryad page and put the zip in `data_raw/dendro/rao_oaks_dryad/` |
| TreeNet (WSL): CH-Dav and CH-Lae, 10-minute data since 2011 | On request only | Ask the TreeNet team at WSL |
| Paired flux-tower and dendrometer network (M. Rao; includes SE-Svb) | No public release found | Contact the author |
| Davos Seehornwald on EnviDat | The public resource is a picture, the data are linked to ICOS / TreeNet | Covered by DenDrought2018 and TreeNet |
