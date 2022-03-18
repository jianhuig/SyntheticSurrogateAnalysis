# Short description of the dataset

Link to folder: https://utoronto-my.sharepoint.com/:f:/g/personal/jianhui_gao_mail_utoronto_ca/EghIzbuI_blGjBlNuHYdwQ0BZt2EMZ69PJ7jEmklC8zl9g?e=IvMhPt

Unfortunately, Onedrive has very strict sharing policy. This link will expire in one month.

## Text Format Data

* plink2.PHENO1.glm.linear

Plink linear regression result adjusted for age,sex,and 10 genetic PC

* covariate.txt

Contains family ID, age, sex, and 10 genetic PC used in plink asscoation test

* phenotype.txt

Contains family ID and inverse normal transfornmed height used in plink association test

* id.txt


List of patient ids that passed QC controls and have no missing covariate/phenotype data

* Weight.txt

Contains family id and weight extracted from UKBB

* Phenotype.tab

Contains family id, age, sex, ethinicty, and height extracted from UKBB

* 50_irnt.gwas.imputed_v3.both_sexes.tsv

GWAS summary statistics of height downloaded from Nealab

## Processed Data in RData format

* height_nlab.RData

Cleaned GWAS summary statistics of height by Nealab

* height_prediction.RData

Contains id, height, age, sex, and weight to build prediction model.
