# Heritable variance in fitness and neutral diversity in wild vertebrates

Code and data for:

> Mawass, W. *Heritable variance in fitness in wild vertebrates accounts for a small share of their
> missing neutral diversity and exceeds what recurrent mutation is estimated to maintain.* Preprint
> (bioRxiv, 2026).

The paper takes the additive genetic variance in relative fitness, $V_A$, estimated by Bonnet et al.
(2022, *Science* 376:1012) in 19 wild vertebrate pedigrees, and runs it through two equilibrium tests:
how much neutral diversity that variance removes (relative to the diversity each population is
missing), and how much of it recurrent deleterious mutation can maintain at mutation–selection–drift
balance.

## Contents

```
scripts/R/          analysis, figure and table scripts (run from the repository root)
data/               curated inputs, with provenance in data/README.md
data/processed/     the tables every analysis reads (built by 00_prepare_data.R)
main/, supp/        main-text and supplementary figures (PDF), as printed in the paper
```

## Two inputs are not included

### Posterior draws of $V_A$ (Bonnet et al. 2022)

Every analysis starts from the 4,000 MCMC draws of $V_A$ per population from the animal-model fits of
Bonnet et al. (2022). They were shared by T. Bonnet and are **not redistributed here**; they
are available on request from the authors of that study. (Their phenotypic records and relatedness
matrices, from which the models can be refitted, are archived as Data S1 of Bonnet et al. 2022.)

To reproduce the analysis, save the draws at the repository root as `alldistzip`: an R save file
holding a list named `alldistzip`, one element per population code (below), whose `$id` component is
the vector of draws of $V_A$ on the observed scale of relative fitness. Then run
`scripts/R/00_prepare_data.R`, which writes `data/processed/va_posterior.csv` (the draws, also not
included) alongside the other processed tables.

Without the draws, every other processed table is included (posterior modes and 95% HPD intervals of
$V_A$ are in `data/processed/populations.csv`), and every script runs from them except the five that
use the draws themselves: `fig1_unlinked.R`, `fig4_selection.R`, `figS_gap.R`,
`sens_maintenance_joint.R` and `audit_numbers.R`.

### Lewin & Eyre-Walker (2026) supplementary workbook

`00_prepare_data.R` also reads the nucleotide diversity, mutation rates, density and range data of
Lewin & Eyre-Walker (2026, *Molecular Ecology* 35:e70265). Download their Supporting Information
workbook (Supplementary Tables 1–7) from the article page, https://doi.org/10.1111/mec.70265, and save
it as `data/lewin_eyre_walker_2026_data.xlsx`. The values the analysis takes from it are already in
`data/processed/`, so it is needed only to rebuild those tables.

## Reproducing the paper

Requirements: R 4.6.0 with dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, readxl 1.4.5, ggplot2 4.0.3,
patchwork 1.3.2, scales 1.4.0, MCMCglmm 2.36, coda 0.19.4.1 and ragg. All commands run from the
repository root.

```sh
Rscript scripts/R/00_prepare_data.R        # ~10 s; rebuilds data/processed/ (needs both inputs above)
Rscript scripts/R/fig1_unlinked.R          # then any figure or table script, in any order
```

| script | writes | in the paper |
|---|---|---|
| `fig1_unlinked.R` | `main/Fig1_unlinked` | Figure 1 |
| `fig3_ceiling.R` | `main/Fig3_ceiling`, `supp/FigS_required`, `supp/FigS_ceiling` | Figure 2, Figures S1, S2 |
| `fig2_linked.R` | `main/Fig2_linked`, `supp/FigS_omega`, `supp/FigS_linked_fgene` | Figure 3, Figures S4, S5 |
| `fig4_selection.R` | `main/Fig4_selection` | Figure 4 |
| `fig5_rescaling.R` | `main/Fig5_rescaling` | Figure 5 |
| `tab_reduction.R` | `manuscript/tables/tab_reduction.tex` (LaTeX) | Table 3 |
| `figS_kappa.R` | `supp/FigS_kappa` | Figure S3 |
| `figS_noncoding.R` | `supp/FigS_noncoding`, `supp/FigS_noncoding_linked` | Figures S6, S7 |
| `figS_gap.R` | `supp/FigS_gap` | Figure S8 |
| `sens_maintenance_joint.R` | `supp/FigS_maintenance_joint`, `data/processed/maintenance_joint_sensitivity*.csv` | Figure S9 |
| `figS_rescaling_linked.R` | `supp/FigS_rescaling_linked` | Figure S10 |
| `figS_maps.R` | `supp/FigS_maps` | Figure S11 |
| `shares.R` | `data/processed/shares.csv` | the two shares quoted in Results |
| `audit_numbers.R` | (prints) | other numbers quoted in the text |
| `check_q2bar_quadrature.R` | (prints) | check of the closed form in Appendix D |
| `check_tail_realization.R` | (prints) | simulation in Appendix B |

The figure file names do not follow the printed numbering (`Fig3_ceiling` prints as Figure 2 and
`Fig2_linked` as Figure 3). Figures are written as PDF and TIFF; only the PDFs are kept here.

Shared code: `kappa.R` (mating-system scalar κ, $Q$, $N_e/N$), `linked.R` (DFE assignment, Haldane
averaging of $Q^2$, the linked term), `load_data.R` (the only route by which scripts read data) and
`theme_genetics.R` (figure style).

### Upstream steps (already run; outputs included)

These rebuild inputs that are committed under `data/` and need network access or local copies of
NCBI genome annotations. See `data/README.md` for what each one does.

| script | writes |
|---|---|
| `fetch_extract_assemblies.R` | `data/assembly_stats_added.csv` |
| `lg_orthology.R` | `data/lg_orthology.csv` |
| `census_buffalo.R` | `data/census/census_buffalo_method.csv` (a reimplementation of Buffalo 2021, not used in the paper) |
| `functional_target.R` (+ `functional_target_gff.py`, Python 3 standard library) | `data/functional_target_lengths.csv` |

NCBI assemblies (accessions in Table S3 of the paper) are expected as
`genomic_data/<species>/<accession>.zip`, as downloaded by the NCBI `datasets` tool.

## Population codes

| code | population | species |
|---|---|---|
| RD | Red deer (Isle of Rum, UK) | *Cervus elaphus* |
| SV | Snow voles (Graubünden, Switzerland) | *Chionomys nivalis* |
| NH | Spotted hyenas (Ngorongoro Crater, Tanzania) | *Crocuta crocuta* |
| BTR | Blue tits (La Rouvière, France) | *Cyanistes caeruleus* |
| BTP | Blue tits (Pirio, Corsica, France) | *Cyanistes caeruleus* |
| BTM | Blue tits (Muro, Corsica, France) | *Cyanistes caeruleus* |
| CF | Collared flycatchers (Gotland, Sweden) | *Ficedula albicollis* |
| CSM | Rhesus macaques (Cayo Santiago, Puerto Rico) | *Macaca mulatta* |
| SFW | Superb fairy-wrens (Canberra, Australia) | *Malurus cyaneus* |
| SSP | Song sparrows (Mandarte Island, Canada) | *Melospiza melodia* |
| HT | Hihi (Tiritiri Matangi Island, New Zealand) | *Notiomystis cincta* |
| HK | Hihi (Karori, New Zealand) | *Notiomystis cincta* |
| SSH | Soay sheep (St Kilda, UK) | *Ovis aries* |
| BHS | Bighorn sheep (Ram Mountain, Canada) | *Ovis canadensis* |
| AB | Yellow baboons (Amboseli National Park, Kenya) | *Papio cynocephalus* |
| GTU | Great tits (Wytham Woods, UK) | *Parus major* |
| GTN | Great tits (Hoge Veluwe, the Netherlands) | *Parus major* |
| KS | Red squirrels (Kluane, Canada) | *Tamiasciurus hudsonicus* |
| MK | Meerkats (Kalahari, South Africa) | *Suricata suricatta* |

## Third-party data

Several inputs are compiled from published sources, cited in the paper and in `data/README.md`:
genome assemblies and annotations (NCBI), mutation rates (Wang & Obbard 2023; Lewin & Eyre-Walker
2026), nucleotide diversity (Lewin & Eyre-Walker 2026 and the studies in `data/pi_literature.csv`),
census sizes (Buffalo 2021) and distributions of fitness effects (Lin et al. 2026). Cite the
original publications if you use these values.

## Citation

Cite the paper above and this repository (`CITATION.cff`; GitHub's "Cite this repository" button
gives the formatted citation).

## License

Code: MIT License (see `LICENSE`). Data compiled for this study: CC BY 4.0. Third-party data remain
under the terms of their original publications.
