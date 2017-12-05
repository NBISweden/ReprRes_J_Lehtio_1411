# Breast cancer quantitative proteome and proteogenomic landscape -- supplementary data and software

This is a collaborative project between the [Janne Lehtiö
group](http://ki.se/en/onkpat/janne-lehtios-group), Department of
Oncology-Pathology, Karolinska Institutet and the
[NBIS](http://nbis.se) long-term support (a.k.a. WABI), performed
during 2015-17 (NBIS LT project name: j_lehtio_1411)

As a part of the NBIS reproducible research policy, this repository
provides the data-files and R-scripts necessary to perform the
underlying analyses and create the tables and plots present in the
manuscript *Johansson et al. Breast cancer quantitative proteome and
proteogenomic landscape*.

More specifically, the script __FinalFigures.Rmd__ reproduces parts of
Main text figure 5, and Supplementary figure S5 and Suppementary Table
S3 (as well as some additional figures not used in the manuscript).


## Usage

#### Download

The content of this repository can be obtained in two ways:

* Download the current release as a zip-archive by clicking [download zip]()

* Clone the repository from the github page (reached by clicking [View
  the Project on GitHub](https://github.com/NBISweden/j_lehtio_1411))

#### Running the scripts

The easiest way of running the script is to open __FinalFigures.Rmd__ in Rstudio and
press the knitr button.

Alternatively, the following command can be issued inside a
R-environment (command-line or script):


`require(rmarkdown)`
 
`render('FinalFigures.Rmd', output_format='html_document', output_file='FinalFigures.html')`

*Note: it seems like this alternative requires (at least) pandoc to
 be installed (e.g., using [macports](https://www.macports.org)).*

The script are designed to install required R-packages as needed;
however, this feature has not been tested extensively on different
operative systems (scripts were developed on MacOSX).


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
