################################################################################

## 5-0_SDM-wallace.R

### Author: Emily Beckman Bruns
### Funding:
# Cooperative agreement between the United States Botanic Garden and
#   San Diego Botanic Garden (subcontracted to The Morton Arboretum), to work
#   on a conservation gap analysis for North American fruit and nut tree
#   crop wild relatives

### Creation date: 07 September 2022
### Last updated:

### R version 4.1.3

### DESCRIPTION:
  #
  # This script
  #

### DATA IN:
  #

### DATA OUT:
  #

################################################################################
# Load libraries
################################################################################

#rm(list=ls())
my.packages <- c('wallace')
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

run_wallace()

## alternativley, you can open a termainal window and run this (so you can
#   still use the R consol while wallace is running)
R
library(wallace)
run_wallace()

################################################################################
# Set working directory
################################################################################


################################################################################
# Load functions
################################################################################


################################################################################
################################################################################
#
################################################################################
