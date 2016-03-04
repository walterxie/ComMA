# Developer Note

The following commands create __*.Rd__ files to the __man__ directory, 
and update __NAMESPACE__ file in the package directory.

```
setwd("???/ComMA")

library("devtools")
library(roxygen2)

document()
```

The following commands install package from local.
```
setwd("???/ComMA")

library("devtools")
install("../ComMA")
```
