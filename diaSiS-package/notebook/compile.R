# usethis::use_mit_license()

library(devtools);
load_all(".");


library(roxygen2); # Read in the roxygen2 R package
roxygenise(clean = T);      # Builds the help files


install.packages("../diaSiS-package", repos = NULL, type = 'source')

library(diaSiS)

