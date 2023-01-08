
# Title of Dataset: Diet can alter the cost of resistance to a natural parasite in Caenorhabditis elegans
---

This dataset contains the data and analysis files for the paper "Diet can alter the cost of resistance to a natural parasite in Caenorhabditis elegans." 

fecundityFinal.csv is the dataset for the fecundity assay, which compared the fecundity of susceptible and disease-resistant individuals.

pg.csv is the dataset for the population expansion assay, which compared population expansion between resistant and susceptible inidividuals over a standardized period.

cost_final.R is the R file containing the analyses of these datasets. glmm_funs.R is an additional file with functions that are used to test for overdispersion in the data. 


## Description of the Data and file structure

fecundityFinal.csv contains data in 8 columns. "ID" is a blinded identity marker (the genotype was unknown to people scoring the assay). "Fecundity" refers to the number of offspring produced on that day of the assay. It is left blank if there is no data from that day, indicating that the worm was censored or died. "Genotype" refers to whether the individual was susceptible (N2) or resistant (ERT250) to the parasite. "Food"" is the bacterial diet the individual was fed on throughout the assay. "Day" is the day of the assay. "Censored" refers to whether the individual was included (0) or removed (1) from the assay due to non-focal damage. "Counter" is the person performing the assay. "Notes" includes miscellaneous observations and further information on other columns.

pg.csv contains data in 12 columns. "ID" is a blinded identity marker (the genotype and food type was unknown to people scoring the assay). "Replicate" indicates which replicate of the 6 count replicates. "Count" indicates the number of individuals per 20 uL aliquot, a subset of the population contained the in the 14.5 mL total volume. An average of these counts was multiplied by 725 to calculate the estimated population size. "Counter" is the assay scorer. "Food" is the individual's diet. "Genotype" refers to whether the individual was susceptible (N2) or resistant (ERT250) to the parasite. "Died" marks whether the focal individual survived to reproduction (0) or died (1). "Fungus" indicates that there was (1) or was not (0) a fungal contaminant. "Bacteria" indicates that there was (1) or was not (0) a bacterial contaminant. "Censored" refers to whether the individual was included (0) or removed (1) from the assay due to non-focal damage or recording errors. "Contaminated" indicates presence (1) or absence (0) of either bacterial or fungal contaminants. "Notes" includes miscellaneous observations and further information on other columns.


## Sharing/access Information

Links to other publicly accessible locations of the data: doi:10.5061/dryad.2v6wwpzsv![image](https://user-images.githubusercontent.com/41449473/211208015-16701e10-5dff-4824-8e0f-c12be12a51c2.png)
