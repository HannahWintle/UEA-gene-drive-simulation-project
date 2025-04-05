###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   MEREA
#   Hannah Wintle, Philip Leftwich
#   March 2025
#   cnm20syu@uea.ac.uk
#
#   Inheritance Cube: MEREA (Maternal Effect RECESSIVE Embryonic Arrest)
#
#   This drive contains two different transgenes on the MEREA allele. Exu-Cas9 
#   (acting as the toxin) remains the same, but the rescue genes are split:
#    * Ma: Contains rescuegenepromoter-tTAV
#    * Mb: Contains teto-recoded.rescuegene
#
#   Rescue Conditions:
#   - Males must inherit both copies of each transgene to survive (MaMb genotype).
#   - Instead of a simple homozygous rescue system (as in MEREA), males must be 
#     trans-homozygous.
#   - Males with MaMa or MbMb **are not viable**.
#   - Only MaMb males **survive**.
#
#   Male Genotype Key:
#    * ZZ      : Wild-type male
#    * MaZ     : Male with only Ma transgene 
#    * MbZ     : Male with only Mb transgene 
#    * MaMa    : Male with two Ma alleles
#    * MbMb    : Male with two Mb alleles
#    * MaMb    : Male with both Ma and Mb
#    * RZ      : Male with resistance allele
#    * MaR     : Male with resistance + Ma
#    * MbR     : Male with resistance + Mb
#    * RR      : Male with two resistance alleles
#   
#   Female Genotype Key:
#    * ZW   : Wild-type female
#    * MaW  : Female with only Ma transgene
#    * MbW  : Female with only Mb transgene
#    * RW   : Female with resistance allele
#
#' @param rM    Breakdown of MEREA allele, no homing/toxin/antidote, M -> R conversion
#' @param Teff  Efficacy of the toxin
#' @param eta   Genotype-specific mating fitness
#' @param phi   Genotype-specific sex ratio at emergence
#' @param omega Genotype-specific multiplicative modifier of adult mortality
#' @param xiF   Genotype-specific female pupatory success
#' @param xiM   Genotype-specific male pupatory success
#' @param s     Genotype-specific fractional reduction (increase) in fertility
#'
#' @return Named list containing the inheritance cube, transition matrix, genotypes, 
#' wild-type allele, and all genotype-specific parameters.
#' @export
#
###############################################################################

cubeMEREA_2L <- function(rM = 0, Teff = 1.0, eta = NULL, phi = NULL,
                         omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  ## safety checks
  if(any(c(rM, Teff)<0) || any(c(rM, Teff)>1)){
    stop("Parameters are rates.
         0 <= x <= 1")
  }
  
  ## Define genotype list
  gtype <- c('ZW', 'ZZ', 'MaZ', 'MbZ', 'MaMb', 'MaMa', 'MbMb', 'RZ', 'MaR', 'MbR', 'RR', 'MaW', 'MbW', 'RW')
  size <- length(gtype)
  
  ## Initialize transition matrix
  tMatrix <- array(data = 0, dim = c(size, size, size), dimnames = list(gtype, gtype, gtype))
  
  ## Fill tMatrix with inheritance probabilities
  tMatrix['ZW', 'ZZ', c('ZW', 'ZZ')] <- c(1, 1) / 2  # Wild-type cross
  
  tMatrix['ZW', 'MaZ', c('MaZ', 'ZZ', 'MaW', 'ZW', 'RZ', 'RW')] <- c((1 - rM), 1, (1 - rM), 1, rM, rM) / 4 # checked
  tMatrix['ZW', 'MbZ', c('MbZ', 'ZZ', 'MbW', 'ZW', 'RZ', 'RW')] <- c((1 - rM), 1, (1 - rM), 1, rM, rM) / 4 # checked
  tMatrix['ZW', 'MaMa', c('MaZ', 'MaW', 'RZ', 'RW')] <- c((1-rM)^2, (1-rM)^2, 1-((1-rM)^2), 1-((1-rM)^2))/2 # checked
  tMatrix['ZW', 'MbMb', c('MbZ', 'MbW', 'RZ', 'RW')] <- c((1-rM)^2, (1-rM)^2, 1-((1-rM)^2), 1-((1-rM)^2))/2 # checked
  tMatrix['ZW', 'MaMb', c('MaZ', 'MbZ', 'MaW', 'MbW', 'RZ', 'RW')] <- c((1 - rM), (1 - rM), (1 - rM), (1 - rM), 2*rM, 2*rM) / 4 # checked
  tMatrix['ZW', 'RZ', c('ZZ', 'RZ', 'ZW', 'RW')] <- c(1, 1, 1, 1) / 4 # checked
  tMatrix['ZW', 'MaR', c('MaZ', 'RZ', 'MaW', 'RW')] <- c(1-rM, (1-rM)+rM*2, 1-rM, (1-rM)+rM*2)/4 # checked
  tMatrix['ZW', 'MbR', c('MbZ', 'RZ', 'MbW', 'RW')] <- c(1-rM, (1-rM)+rM*2, 1-rM, (1-rM)+rM*2)/4 # checked
  tMatrix['ZW', 'RR', c('RZ', 'RW')] <- c(1, 1) / 2 # checked
  
  tMatrix['MaW','ZZ', c('MaZ', 'ZW', 'RZ')] <- c(1-rM, 1, rM)/2 # checked
  tMatrix['MaW','MaZ', c('MaW', 'MaMa', 'MaZ', 'ZW', 'RZ', 'RW', 'RR', 'MaR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4 # checked
  tMatrix['MaW','MbZ', c('MbW','MaMb','MaZ','ZW','RZ','RW','RR','MaR','MbR')] <- c((1 - rM), (1 - rM)^2, (1 - rM), 1, rM, rM, rM^2, rM*(1 - rM), rM*(1 - rM)) / 4 # checked
  tMatrix['MaW','MaMa', c('MaMa', 'MaW', 'RR', 'RW', 'MaR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2 # checked
  tMatrix['MaW','MbMb', c('MaMb','MbW','RR','MaR','MbR','RW')] <- c((1 - rM)^2, (1 - rM), rM^2, rM*(1 - rM), rM*(1 - rM), rM) / 2 # checked
  tMatrix['MaW','MaMb', c('MaMa', 'MaMb', 'MaW', 'MbW', 'RR', 'RW', 'MaR', 'MbR')] <- c((1 - rM)^2, (1 - rM)^2, (1 - rM), (1 - rM), 2*rM-(2*rM*(1-rM)), 2*rM, 3*rM*(1 - rM), rM*(1 - rM))/4 # checked
  tMatrix['MaW','RZ', c('MaZ', 'MaR', 'ZW', 'RW', 'RZ', 'RR')] <- c((1-rM), (1-rM), 1, 1, rM, rM)/4 # checked
  tMatrix['MaW','MaR', c('MaMa', 'MaR', 'MaW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4 # checked
  tMatrix['MaW','MbR', c('MaR', 'MaMb', 'MbW', 'RW', 'RR', 'MbR')] <- c((1 - rM)*(1 + rM), (1 - rM)^2, (1 - rM), (1 + rM), rM*(1 + rM), rM*(1 - rM)) / 4 # checked
  tMatrix['MaW','RR', c('MaR', 'RW', 'RR')] <- c((1-rM), 1, rM)/2 # checked
  
  tMatrix['MbW','ZZ', c('MbZ', 'ZW', 'RZ')] <- c(1-rM, 1, rM)/2 # checked
  tMatrix['MbW','MaZ', c('MaW','MaMb','MbZ','ZW','RZ','RW','RR','MaR','MbR')] <- c((1 - rM), (1 - rM)^2, (1 - rM), 1, rM, rM, rM^2, rM*(1 - rM), rM*(1 - rM)) / 4 # checked
  tMatrix['MbW','MbZ', c('MbW', 'MbMb', 'MbZ', 'ZW', 'RZ', 'RW', 'RR', 'MbR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4 # checked
  tMatrix['MbW','MaMa', c('MaMb', 'MaW', 'RR', 'MaR', 'MbR', 'RW')] <- c((1 - rM)^2, (1 - rM), rM^2, rM*(1 - rM), rM*(1 - rM), rM) / 2 # checked
  tMatrix['MbW','MbMb', c('MbMb', 'MbW', 'RR', 'RW', 'MbR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2 # checked
  tMatrix['MbW','MaMb', c('MbMb', 'MaMb', 'MaW', 'MbW', 'RR', 'RW', 'MaR', 'MbR')] <- c((1 - rM)^2, (1 - rM)^2, (1 - rM), (1 - rM), 2*rM-(2*rM*(1-rM)), 2*rM, 3*rM*(1 - rM), rM*(1 - rM))/4 # checked
  tMatrix['MbW','RZ', c('MbZ', 'MbR', 'ZW', 'RW', 'RZ', 'RR')] <- c((1-rM), (1-rM), 1, 1, rM, rM)/4 # checked
  tMatrix['MbW','MaR', c('MbR', 'MaMb', 'MaW', 'RW', 'RR', 'MaR')] <- c((1 - rM)*(1 + rM), (1 - rM)^2, (1 - rM), (1 + rM), rM*(1 + rM), rM*(1 - rM)) / 4 # checked
  tMatrix['MbW','MbR', c('MbMb', 'MbR', 'MbW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4 # checked
  tMatrix['MbW','RR', c('MbR', 'RW', 'RR')] <- c((1-rM), 1, rM)/2 # checked
  
  tMatrix['RW','ZZ', c('RZ', 'ZW')] <- c(1, 1)/2 # checked
  tMatrix['RW','MaZ', c('MaR', 'RZ', 'MaW', 'ZW', 'RR', 'RW')] <- c((1-rM), 1, (1-rM), 1, rM, rM)/4 # checked
  tMatrix['RW','MbZ', c('MbR', 'RZ', 'MbW', 'ZW', 'RR', 'RW')] <- c((1-rM), 1, (1-rM), 1, rM, rM)/4 # checked
  tMatrix['RW','MaMa', c('MaR', 'MaW', 'RR', 'RW')] <- c((1-rM), (1-rM), rM, rM)/2 # checked
  tMatrix['RW','MbMb', c('MbR', 'MbW', 'RR', 'RW')] <- c((1-rM), (1-rM), rM, rM)/2 # checked
  tMatrix['RW', 'MaMb', c('MaR', 'MbR', 'RR', 'MaW', 'MbW', 'RW')] <- c((1 - rM), (1 - rM), 2*rM, (1 - rM), (1 - rM), 2*rM) / 4 # checked
  tMatrix['RW','RZ', c('RR', 'RZ', 'RW', 'ZW')] <- c(1, 1, 1, 1)/4 # checked
  tMatrix['RW','MaR', c('RW', 'MaR', 'RR', 'MaW')] <- c((1+rM), (1-rM), (1+rM), (1-rM))/4 # checked
  tMatrix['RW','MbR', c('RW', 'MbR', 'RR', 'MbW')] <- c((1+rM), (1-rM), (1+rM), (1-rM))/4 # checked
  tMatrix['RW','RR', c('RW', 'RR')] <- c(1, 1)/2 # checked
  
  ## Define viability mask
  viabilityMask <- array(data = 1, dim = c(size, size, size), dimnames = list(gtype, gtype, gtype))
  viabilityMask[c('MaMa', 'MbMb', 'RZ', 'MaR', 'MbR'), , ] <- 0 
  
  ## Genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)
  
  ## Define sex modifiers
  # 0 is male and 1 is female
  sexes <- c(
    "ZW" = 1, 
    "ZZ" = 0,
    "MaZ" = 0, 
    "MbZ" = 0, 
    "MaMb" = 0,
    "MaMa" = 0, 
    "MbMb" = 0, 
    "RZ" = 0, 
    "MaR" = 0, 
    "MbR" = 0, 
    "RR" = 0,
    "MaW" = 1, 
    "MbW" = 1, 
    "RW" = 1
  )
  phi <- sexes
  
  ## Return cube structure
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c("ZZ", "ZW"),
    eta = modifiers$eta,
    phi = phi,  # Assign sex ratios
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = c("MaMb", "MaMa", "MbMb", "MaW", "MbW")
  ))
}
