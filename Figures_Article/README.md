# Figures for the manuscript  

***

### Scripts and data used to generate the figures of the manuscript [*Theulot et al, 2024*](https://doi.org/10.1038/s41467-024-55520-3)

The large file containing the telomeric sequences and the metric for filtering used in the figures 4 and S15 to S19 is on the Zenodo repository under the DOI [*10.5281/zenodo.12668295*](https://doi.org/10.5281/zenodo.12668295). These data were generated using the [*telo_tools scripts*](https://github.com/touala/telo_tools)

To generate session_info file  

```R
library(sessioninfo)  
session_info() %>% capture.output(file="nanotiming_figures_session_info.txt")
```  

***
