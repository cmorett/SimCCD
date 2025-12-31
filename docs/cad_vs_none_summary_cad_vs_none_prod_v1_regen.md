# CAD vs none summary (cad_vs_none_prod_v1_regen)

## Run config
- tag: cad_vs_none_prod_v1_regen
- thrown_none: 1000000
- thrown_cad: 500000
- compare_dir: paper_outputs\cad_vs_none_prod_v1_regen\compare

## Cutflow (counts)
| stage | none | cad |
| --- | --- | --- |
| hits | 23100 | 11504 |
| quality | 23099 | 10688 |
| throughgoing | 22643 | 8460 |
| thrown | 1000000 | 500000 |

## Headline effects
| metric | none | cad | ratio |
| --- | --- | --- | --- |
| hits_per_thrown | 0.0231 +/- 0.00015 | 0.023008 +/- 0.00021 | 0.996017 +/- 0.011 |
| through_per_thrown | 0.022643 +/- 0.00015 | 0.01692 +/- 0.00018 | 0.747251 +/- 0.0094 |
| through_fraction_of_hits | 0.980216 +/- 0.00092 | 0.735396 +/- 0.0041 | 0.750239 +/- 0.0043 |
| charge_p99_hits | 207715 +/- 0 | 437323 +/- 0 | 2.1054 +/- 0 |
| charge_p99_through | 210163 +/- 0 | 444435 +/- 0 | 2.11471 +/- 0 |

## Key plots
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_charge_hits_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_charge_hits_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_charge_through_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_charge_through_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_dedx_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_dedx_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_edep_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_edep_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_hit_efficiency_vs_coszen.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_lcos_distribution.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_through_fraction_of_hits_vs_coszen.pdf`
- `paper_outputs\cad_vs_none_prod_v1_regen\compare\compare_through_fraction_per_thrown_vs_coszen.pdf`

