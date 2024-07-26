# Figures for the manuscript  

***

### Scripts and data used to generate the figures of the manuscript [*Theulot et al, 2024*](https://doi.org/10.1101/2024.07.05.602252)

The large file containing the telomeric sequences and the metric for filtering used in the figures 4 and S19 to s23 is on the Zenodo repository under the DOI [*10.5281/zenodo.12668295*](https://doi.org/10.5281/zenodo.12668295)

To generate session_info file  

```R
library(sessioninfo)  
session_info() %>% capture.output(file="nanotiming_figures_session_info.txt")
```  

***
