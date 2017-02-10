# R Package for Community Matrix Analysis

In development. 

##Use my time to save your time

The package contains functions to read, write, manipulate, visualize, analyse community data and its meta data, 
such as taxonomy, environment variables, phylogenetic trees. The analyses aim to assessing diversities,
finding primary effective metadata, discover the pattern of high related samples, etc.   

This package also import package [gg1L](https://github.com/walterxie/gg1L) providing one-line plotting functions.

To see all exported functions:
```R
help(package = "ComMA")
```

##Installation

You can use the **devtools** *install\_github()* function to install the lastest development version directly from the GitHub.

```R
library("devtools")
devtools::install_github("walterxie/gg1L")
devtools::install_github("walterxie/ComMA")
library("gg1L")
library("ComMA")
```

##Tutorials

It may take time to load the page from GitHub.

1. [One-line plotting] (https://github.com/walterxie/gg1L)

2. [RDP + Greengenes vs. BLAST + MEGAN] (https://github.com/walterxie/ComMA/blob/master/tutorials/GreengenesVsBLAST.ipynb)
