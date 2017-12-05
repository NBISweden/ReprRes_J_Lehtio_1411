# j_lehtio_1411

Repository for NBIS LT project j_lehtio_1411 with Janne Lehtiö group 2015-17

This repository provides the data-files and R-scripts necessary to
perform analyses and create the tables and plots for the manuscript
*Johansson et al. Breast cancer quantitative proteome and
proteogenomic landscape*.

More specifically, the script __FinalFigures.Rmd__ Main text figure
5a-b, and Supplementary figure S5 and Suppementary Table S3.


## Usage
The easiest way of 



## Repository content:

* __README.md__
This file.

##### Folders

* __Data/__ Directory holding the data provided by the Lehtiö group.
See further README.Rmd in that folder.

##### Main R-files
The two .RMD documents perform the main analyses in this repository.
They double as the main documentation tool:

* __FinalFigures.Rmd__
R markdown document for manuscript figures relating to analysis of
correlation between protein expression and mRNA expression.

##### Helper R-files

* __Helper.R__
Helper R-script defining useful convenience fuctions
* __formatData.R__
R-script that reads the data files obtained from Henrik Johansson,
Lehtio's group (usually in Excel format), reformat data to format 
useful in the present analyses, and save them as a R-object, for 
later use in .Rmd-files.

* __updating_gene_symbols.R__
Fixing some gene symbol errors in mRNA data

##### Miscellaneous files

* __GRCh37_chromLengths_NCBI.txt__
Chromosome lengths copied from http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/?build=37