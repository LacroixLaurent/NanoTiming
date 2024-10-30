# NanoTiming
### Laurent Lacroix (laurent.lacroix@inserm.fr)
***
### BrdU Base calling

Raw data are available from ENA repository under accession number PRJEB76824.  
Expected outputs using the fast5 from WT_rep3 is available from the Zenodo repository under the DOI [10.5281/zenodo.12668295](https://doi.org/10.5281/zenodo.12668295).  

In order to test this procedure, the WT_rep3.tar.gz file should be downloaded and expanded into a WT_rep3/fast5 folder.  

This procedure is for ONT R9.4.1 data type acquired in a fast5 format.

*basecalling_sample.sh* contains a example of the base calling procedure going from RawData to the megalodon mod_mappings.bam file.  

The mod_mappings.bam data should then be processed using the script in the *NanoTiming_Data* folder.  

For this project, we have built a reference genome from ONT and Ilumina reads the the BT1 yeast strain. This strain was reported in Theulot et all 2022. The procedure used to build this reference genome is explained in the [manuscript](https://doi.org/10.1101/2024.07.05.602252) and the corresponding annotation is detailed in the *Reference_Genome* folder.  

***
