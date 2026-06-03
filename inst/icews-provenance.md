# ICEWS bundled example provenance

This package bundles a reduced illustrative ICEWS example as `data/icews.rda`.
It is not the published estimation sample for Minhas and Hoff (2025).

Source paper:

- Minhas, S. and Hoff, P. D. (2025). "Decomposing Network Influence: Social Influence Regression." Political Analysis.
- Cambridge DOI: https://doi.org/10.1017/pan.2025.10013

Replication archive:

- Harvard Dataverse DOI: https://doi.org/10.7910/DVN/VTFDX6
- Repository mirror referenced by the paper: https://github.com/s7minhas/sir_paper
- Raw file expected by `data-raw/icews.R`: `replArchive/data/socRegData.rda`

Bundled object:

- File: `data/icews.rda`
- SHA256: `8c717ea32ff7ceb430434d1f36e076cce07bc7b5b201a7842bf3e5dc4bc819b4`
- Size: 267,944 bytes
- Dimensions: `Y` and `X` are `50 x 50 x 95`; `W` is `50 x 50 x 4`; `Z` is `50 x 50 x 5 x 95`
- Metadata: the bundled list includes a `metadata` entry with source DOI, replication DOI, raw-file path, access date, country-selection rule, dates, and transforms.
- Dates: 2005-02-01 through 2012-12-01

Derivation:

- Script in the source repository: `data-raw/icews.R` (excluded from the installed package by `.Rbuildignore`)
- Script SHA256 at this package snapshot: `dc845efa339576410b1a1a63afbbff74c3d21d0ebab49680c230e24f33bb6acc`
- Country rule: keep the 50 countries with the highest total sent plus received material-conflict event counts over the source archive, then sort retained indices into source order.
- Transforms: `W[, , "verbCoop"]` and `Z[, , "verbCoop", ]` are transformed with `log(x + 1)`.
- The first bundled `X` period uses the retained January 2005 lag from the source archive and cannot be reconstructed from bundled `Y` alone.

To regenerate from a source checkout, set `SIR_ICEWS_RAW` to the local path of
`socRegData.rda`, run `data-raw/icews.R`, then verify the source archive version
and raw SHA256 before publishing a rebuilt data object.
