# CAD vs none summary (cad_vs_none_prod_v1)

## Run config
- tag: cad_vs_none_prod_v1
- thrown_none: 1000000
- thrown_cad: 500000
- compare_dir: paper_outputs\cad_vs_none_prod_v1\compare

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
| hits_per_thrown | 0.023099 +/- 0.00015 | 0.021376 +/- 0.0002 | 0.925408 +/- 0.011 |
| through_per_thrown | 0.022643 +/- 0.00015 | 0.01692 +/- 0.00018 | 0.747251 +/- 0.0094 |
| through_fraction_of_hits | 0.980259 +/- 0.00092 | 0.791542 +/- 0.0039 | 0.807483 +/- 0.0041 |
| charge_p99_hits | 207715 +/- 0 | 437323 +/- 0 | 2.1054 +/- 0 |
| charge_p99_through | 208081 +/- 0 | 439548 +/- 0 | 2.11239 +/- 0 |

## Key plots
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_charge_hits_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_charge_hits_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_charge_through_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_charge_through_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_dedx_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_dedx_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_edep_core.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_edep_tail.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_hit_efficiency_vs_coszen.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_lcos_distribution.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_through_fraction_of_hits_vs_coszen.pdf`
- `paper_outputs\cad_vs_none_prod_v1\compare\compare_through_fraction_per_thrown_vs_coszen.pdf`

