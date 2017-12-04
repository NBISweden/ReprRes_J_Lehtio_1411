# j_lehtio_1411

Repository for NBIS LT project j_lehtio_1411 with Janne Lehtiö 2015-17


## Repository content:

* __README.md__
This file.

##### Folders

* Data/ Directory holding the data provided by the Lehtiö group

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

##### Miscellaneous files

* __GRCh37_chromLengths_NCBI.txt__
Chromosome lengths copied from http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/?build=37
