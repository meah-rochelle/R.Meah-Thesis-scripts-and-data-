
### signal
#### set working directory to be FINALVIIRSdata folder with csv files
setwd("D:\\Oxford Sky Pol Data 2018\\FINALVIIRSdata\\test")
library(broom)
library(plyr)
library(ggplot2)

path <- getwd()

#### for moth DRA
filenames <- list.files(path)
filenames <- filenames[-c(1,5)]
DRAc <- c(0.407, 0.603, 0.692, 0.239, 0.159, 0.093, 0.448, 0.206, 0.378, 0.539,
          0.620, 0.629, 0.142, 0.211, 0.071, 0.109, 0.216)

# these are the 14 sites in the thesis (except Gavarnie which 
# PRC coudln't be modelled for due to acute agnle of the pol pattern
# in the amiages)
filenames <- filenames[-c(6,11,13,15)]
DRAc <- c(0.407, 0.603, 0.692, 0.239, 0.159, 0.448, 0.206, 0.378, 0.539,
          0.629, 0.211,0.109, 0.216)
Sig.data <- data.frame(name = filenames, signal = DRAc) # DRA PRC moth


ppr = laply(filenames, function(filename){
  dap <- read.csv(filename)
  dar <- na.omit(dap)
  t.data <- ddply(dar, .(), summarise,
                  me = mean(value),
                  std = sd(value) 
  )
  dan <- c(t.data$me, t.data$std)
  return(dan)
}
)


final.data.moth <- cbind(Sig.data, ppr)
colnames(final.data.moth) <- c("name, total", "signal", "rad", "sd")

####  plots with the x-axis with error bars 
#p<- ggplot(final.data, aes(x=(rad), y=signal)) + geom_point() +
#  geom_errorbarh(aes(xmin=rad-sd, xmax=rad+sd))
#p

o.moth <- order(final.data.moth$rad)
df_o.moth <- final.data.moth[o.moth, ] # order it by increasing x


fit_moth <- nls(signal ~ SSasymp(rad, yf, y0, log_alpha), data = df_o.moth)

qplot(rad, signal, data = augment(fit_moth)) + 
  geom_line(aes(y = .fitted)) + 
  xlim(c(0,70))+ylim(c(0,0.7))

summary(fit_moth)
modelr::rsquare(fit_moth, df_o.moth)

############## for spiders

filenames <- list.files(path)
filenames <- filenames[-1]

sigg <- c(0.274, 0.439, 0.419, 0.470, 0.165, 0.193, 0.149, 0.321, 0.192, 0.262,
          0.276, 0.455, 0.423, 0.157, 0.122, 0.123, 0.144, 0.236)

# these are the 14 sites in the thesis (except Gavarnie which 
# PRC coudln't be modelled for due to acute agnle of the pol pattern
# in the amiages)
filenames <- filenames[-c(7,12,14,16)]
sigg <- c(0.274, 0.439, 0.419, 0.470, 0.165, 0.193, 0.321, 0.192, 0.262,
          0.276, 0.423, 0.122, 0.144, 0.236)
Sig.data.spid <- data.frame(name = filenames, signal = sigg)  # spider

ppr = laply(filenames, function(filename){
  dap <- read.csv(filename)
  dar <- na.omit(dap)
  t.data <- ddply(dar, .(), summarise,
                  me = mean(value),
                  std = sd(value) 
  )
  dan <- c(t.data$me, t.data$std)
  return(dan)
}
)


final.data.spid <- cbind(Sig.data.spid, ppr)
colnames(final.data.spid) <- c("name, total", "signal", "rad", "sd")


#p<- ggplot(final.data, aes(x=(rad), y=signal)) + geom_point() +
#  geom_errorbarh(aes(xmin=rad-sd, xmax=rad+sd))
#p

o.spid <- order(final.data.spid$rad)
df_o.spid <- final.data.spid[o.spid, ] # order it by increasing x


fit_spid <- nls(signal ~ SSasymp(rad, yf, y0, log_alpha), data = df_o.spid)

qplot(rad, signal, data = augment(fit_spid)) + 
  geom_line(aes(y = .fitted)) + 
  xlim(c(0,70)) + ylim(c(0,0.7))

summary(fit_spid)
modelr::rsquare(fit_spid, df_o.spid)

anova(fit_moth, fit_spid)

#ggplot(final.data, aes(x=(rad), y=signal)) + geom_point() +
 ## geom_errorbarh(aes(xmin=rad-sd, xmax=rad+sd)) + 
  #xlim(c(0,70)) 

#### plot both species together ####


qplot() + 
  geom_point(size=2, shape=17,aes(rad, signal, colour="spider"), data = augment(fit_spid)) +
  geom_line(size=1.5, aes(x=rad, y = .fitted, colour="spider"), data = augment(fit_spid)) + 
  geom_point(size=2, shape=16, aes(rad, signal, colour="moth"), data = augment(fit_moth)) +
  geom_line(size=1.5, aes(x=rad, y = .fitted, colour="moth"), data = augment(fit_moth)) +
  scale_x_continuous("Light pollution radiance (nW/cm2/sr)", limits=c(0,50)) +
  scale_y_continuous("Photoreceptor contrast", limits=c(0,0.7)) +
  scale_colour_manual(name="Species", labels= c("Moth", "Spider"), values = c("thistle", "slategray3")) +
  theme_classic() +
  theme(text = element_text(size=22))
