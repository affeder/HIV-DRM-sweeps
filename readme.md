# How to generate the analysis/figures from “More efficient drugs lead to harder sweeps of drug resistance in HIV-1”

### To be run once: 
```
align.r
```
This script cleans the data and generates various parsed and clean versions of the data that are stored in `/tmp` and will be used by other scripts. After `align.r` is run once, everything else requires only `read.in.data.r` besides other dependencies listed explicitly.

### To do the validations for the ambiguous reads, run in order:
```
ambread-validation-prep.r
ambread-validation-graph.r
```
Note: the graph from this figure isn’t actually in the paper (only the correlation coefficient).

### Producing figures:
To produce the following figures that don’t require the fitting of GLMMs:
- Figure 2 (F2, F2-S1, F2-S2) 
- Figure 5 (F5-S1, F5-S2)
```
make-nonmodel-figs.r
```
The rest of the figures in Figure 5 (F5-S3, F5-S4, F5-S5) are produced in the top of `subsample-run.r`. `subsample-run.r` also fits and saves all of the models that are eventually plotted in F3 and F4. It fits subsampled GLMs and GLMMs for three different iterations of the data:

1. All data, p-thinning to the 1995 level
2. Data truncated to 4 DRMs, p-thinning to the 1989 level
3. Data truncated to 4 DRMs, p-thinning to the 1995 level

Most of `subsample-run.r` runs slowly (i.e., days when parallelized to four cores on a 2016 macbook pro). I’ve therefore lowered the iterations from the 1000 shown in the text to 20, which should give you an idea without running for days.

Versions of Figure 3 and Figure 4 for each of data iterations 1, 2 and 3 above can be plotted with `subsample-plot-etc.r`.