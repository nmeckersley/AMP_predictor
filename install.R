# install.r 

# Packages required for the R shiny app
packages <- c(
  "shiny",
  "Biostrings",
  "randomForest",
  "Peptides",
  "dplyr",
  "stringr",
  "ggplot2",
  "DT",
  "here",
  "httr"
)

# Install any missing packages
installed <- rownames(installed.packages())
for(pkg in packages){
  if(!pkg %in% installed) install.packages(pkg)
}
message("All required packages are installed!")
