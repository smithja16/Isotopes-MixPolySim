########################################################################
##  Mixing Polygon Simulation (2 isotopes)                            ##
##  v1.3, March 2026, R version 4.4.1                                 ##
##  From Smith et al 2013, Methods in Ecology and Evolution           ##
##  https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12048
##  Old version was found at http://www.famer.unsw.edu.au/downloads.html
##  New in v1.1: added Fig. 4: black and white plot of mixing region  ##
##  New in v1.2: added Fig. 5: standard bi-plot with 95% contour      ##
##  New in v1.3: updates to file loading, saving, plotting            ##
##  contact: james.smith@unsw.edu.au; james.a.smith@dpird.nsw.gov.au  ##
########################################################################

## This is the top script, and is used to call the main function and
## explore results.
## The code is focused on two isotope systems. The C and N syntax is used for
## a common two isotope system, but this code can be used for any two
## isotopes.


## Housekeeping
# install.packages(c("sp","ggplot2","progress"))
library(sp)
library(ggplot2)
library(progress)

data_dir <- "C:/.../MixPolySim"  #where your data and R scripts are saved
data_dir <- "C:/Users/smithj08/OneDrive - DPIE/Other/MixPolySim"
out_dir  <- file.path(data_dir, "output")  #folder for saving results and figures

source(paste0(data_dir,"/MixPolySim_2iso_1.3_Function.r"))


## Load data
# Match the format of the example data
sources <- read.csv(file.path(data_dir, "Sources_example.csv"), header=TRUE) #always put 13C(x) before 15N(y)
mixture <- read.csv(file.path(data_dir, "Mixture_example.csv"), header=TRUE) #some error is required for every value
TEF     <- read.csv(file.path(data_dir, "TEF_example.csv"),     header=TRUE)


## Run the function
results <- mix_poly_sim(
  sources = sources,
  mixture = mixture,
  TEF = TEF,
  its = 1500,  # number of iterations (1500 ususally enough, but check area plot)
  min_C = -50, # specify the dimensions and resolution for the mixing region figures
  max_C = -20, # choose values outside the 95% mixing region
  min_N = -2,
  max_N = 10,
  res = 250,  # resolution of the mixing region figure; reducing this improves performance
  out_dir = out_dir,   # for saving results
  seed = 117 ) 


## Examine results
results$probabilities
print(results$area_plot)
print(results$mix_plot)
print(results$biplot_95)

