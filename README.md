# pumba-datamining

This project contains the R code used to generate datamining results from Pumba.
It requests data from [Pumba API](https://pumba.dcsr.unil.ch/backend) and [UniProt](https://www.uniprot.org).

## Run the R scripts

You will have to adapt `res_path` in `well_behavioud_proteins.R`.

First you have to run `well_behavioud_proteins.R` to create the main table, and then `add_uniprot_data.R` to add the additional information parsed from UniProt.
