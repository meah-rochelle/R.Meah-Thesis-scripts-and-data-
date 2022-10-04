library(plyr)
library(tidyverse)
library(sf)
library(raster)
#library(spData)
#library(spDataLarge)
#library(rnaturalearth)

library(reshape2)
library(cluster)
#library(ggfortify)
#library(FactoMineR)
#library(factoextra)
library(corrplot)
library(tidyr)
library(viridis)
#library(hexbin)
#library(broom)


library(dplyr)
library(ggplot2)
library(reshape)
library(tiff)
library(fields)
library(imager)
library(magrittr)
library(mmand)
library(jpeg)
library(raster)
library(RNiftyReg)
library(RColorBrewer)


DoP <- function(z,f,n,o) {
  S.1 <- (z - n)/(z+n)
  S.2 <- (f-o)/(f+o)
  d <- sqrt(S.1^2+S.2^2)
  return(d)
}

AoP <- function(z,f,n,o) {
  S.1 <- (z - n)/(z+n)
  S.2 <- (f-o)/(f+o)
  da <- 0.5*atan(S.2/S.1)
  return(da)
}

rang <- function(x) {
  f <- median(x)
  sq <- sd(x)
  a <- c(f-sq, f+sq)
  return(a)
}

ecovis_colors <- c(
  `red`        = "#FF0000",
  `endred`        = "#550000",
  `darkgreen`   = "#003300",
  `brightgreen` = "#80FF00",
  `brightblue`   = "#00FFFF",
  `darkblue`   = "#002EFF",
  `endblue` = "#002151",
  `orange`     = "#FF8000",
  `yellow`     = "#FFF300",
  `light grey` = "#cccccc",
  `dark grey`  = "#8c8c8c",
  `zero` = "#f6e653",
  `fifteen` = "#98d25e",
  `thrity` = "#5ab080",
  `fortyfive` = "#458b8c",
  `sixty` = "#3e6289",
  `seventyfive` = "#42347a",
  `ninety` = "#3e6289",
  `oneofive` = "#7e2f9f",
  `onetwoone` = "#a73c88",
  `onethreefive` = "#c65c6f",
  `onefifty` = "#e68b57",
  `onesixfive` = "#f5c14e",
  `oneeighty` = "#f6e653",
  `neg` = "#091be2",
  `less` = "#020a66",
  `ze` = "#1b1b1c",
  `more` = "#7a7517",
  `pos` = "#efe62f")


#' Function to extract ecovis colors as hex codes
#'
#' @param ... Character names of ecovisvis_colors 
#'
ecovis_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (ecovis_colors)
  
  ecovis_colors[cols]
}


ecovis_palettes <- list(
  `main`  = ecovis_cols("blue", "brightgreen", "red"),
  
  `Cont`  = ecovis_cols("darkgreen", "brightgreen", "neg", "ze", "ze", "pos", "red", "endred"),
  
  `AoP`   = ecovis_cols("zero", "fifteen", "thirty", "fortyfive",
                        "sixty", "seventyfive", "ninety", "oneofive",
                        "onetwoone", "onethreefive", "onefifty", "onesixfive",
                        "oneeighty"),
  
  `DoP` = ecovis_cols("endblue", "darkblue", "brightblue", "yellow", "orange", "red", "endred"),
  
  `DoPmin` = ecovis_cols("endblue", "yellow", "orange", "red", "endred", "endred", "endred", "endred", "endred", "endred"),
  
  
  `Light` = ecovis_cols("darkblue","brightblue", "brightgreen","yellow", 
                        "orange", "red", 
                        "endred"),
  
  `grey`  = ecovis_cols("light grey", "dark grey")
)

#' Return function to interpolate a ecovis color palette
#'
#' @param palette Character name of palette in ecovis_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette()
#'
ecovis_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- ecovis_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}


#' Color scale constructor for ecovis colors
#'
#' @param palette Character name of palette in ecovis_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_ecovis <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- ecovis_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("ecovis_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for ecovis colors
#'
#' @param palette Character name of palette in ecovis_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_ecovis <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- ecovis_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("ecovis_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}


#Constants
Planck = 6.626070040e-34
speed = 299792458

# create two functions to give the number of photons absorbed in the photoreceptor 
sen.par <- function(a_par, a_perp,lambda, lambda.max, P0, P1) {
  x.alpha <- log10(lambda/lambda.max)
  x.beta <- log10(lambda/350)
  a.alpha <- 380
  b.alpha <- 6.09
  a.beta <- 247
  b.beta <- 3.59
  alpha <- exp(-a.alpha*x.alpha*x.alpha*(1+b.alpha*x.alpha+(3/8)*(b.alpha*x.alpha)*(b.alpha*x.alpha)))
  beta<- 0.29*exp(-a.beta*x.beta*x.beta*(1+b.beta*x.beta+(3/8)*(b.beta*x.beta)*(b.beta*x.beta)))
  q <- 10^(-a_par*(alpha+beta))
  r <- 10^(-a_perp*(alpha+beta))
  result <- P0-(0.5*(P0*(q[row(P0)]+r[row(P0)])+P1*(q[row(P1)]-r[row(P1)])))
  return(result)
}

sen.perp <- function(a_par, a_perp, lambda, lambda.max, P0, P1) {
  x.alpha <- log10(lambda/lambda.max)
  x.beta <- log10(lambda/350)
  a.alpha <- 380
  b.alpha <- 6.09
  a.beta <- 247
  b.beta <- 3.59
  alpha <- exp(-a.alpha*x.alpha*x.alpha*(1+b.alpha*x.alpha+(3/8)*(b.alpha*x.alpha)*(b.alpha*x.alpha)))
  beta<- 0.29*exp(-a.beta*x.beta*x.beta*(1+b.beta*x.beta+(3/8)*(b.beta*x.beta)*(b.beta*x.beta)))
  q <- 10^(-a_perp*(alpha+beta))
  r <- 10^(-a_par*(alpha+beta))
  result <- P0-(0.5*(P0*(q[row(P0)]+r[row(P0)])+P1*(q[row(P1)]-r[row(P1)])))
  return(result)
}


fop <- function(filename) {
  da <- read.table(filename, header=F, sep="\t")
  da$V1 <- abs(da$V1-mean(da$V1[1:100]))
  ou <- cbind(da, rep(filename, nrow(da)))
  return(ou)
}

nr <- 1044 # rows in each file

