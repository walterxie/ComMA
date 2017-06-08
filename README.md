# R Package for Community Matrix Analysis

In development. 

## Tools to anlyse community matrix especially from eDNA 

The package contains functions to read, write, manipulate, visualize, analyse community data and its meta data, 
such as taxonomy, environment variables, phylogenetic trees. The analyses aim to assessing diversities,
finding primary effective metadata, discover the pattern of high related samples, etc.   

Some visualizations have been moved to the [gg1L](https://github.com/walterxie/gg1L) package which provides one-line plotting functions. Recommend to load gg1L as well when using a [pipeline](http://walterxie.github.io/eDNA-pipeline/). 

## Installation

You can use the **devtools::install\_github()** function to install the latest release directly from the GitHub:
```R
library(devtools)
devtools::install_github("walterxie/ComMA@*release")
library(ComMA)
library(gg1L)
```

Or a perticular release version:
```R
devtools::install_github("walterxie/ComMA@v0.1.0")
```

Or the latest development version:
```R
devtools::install_github("walterxie/ComMA")
```

Read the package introduction:
```R
package?ComMA
```

To see all exported functions:
```R
help(package = "ComMA")
```

## Tutorials

It may take time to load the page from GitHub.

1. [One-line plotting](https://github.com/walterxie/gg1L)

2. [RDP + Greengenes vs. BLAST + MEGAN](https://cdn.rawgit.com/walterxie/ComMA/master/tutorials/GreengenesVsBLAST.html)
