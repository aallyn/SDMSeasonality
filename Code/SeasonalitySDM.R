#####
## Accounting for Seasonality in Species Distribution Models
######
# Libraries and preliminaries ---------------------------------------------
# Libraries
library(mgcv)
library(tidyverse)
library(ncdf4)
library(raster)
library(zoo)
library(maptools)
library(ggpmisc)
library(rgeos)
library(forecast)
library(sdm)
library(ROCR)
library(SimDesign)
library(akima)
library(sf)
library(viridis)
library(cowplot)
library(SDMTools)
library(lme4)
library(bestNormalize)
library(adehabitatHR)
library(geosphere)

# Used functions
temp.scale<- function(new.temp, base.temp.mean, base.temp.sd){
  if(is.na(new.temp)){
    temp.scaled<- NA
    return(temp.scaled)
  } else {
    temp.scaled<- (new.temp - base.temp.mean)/base.temp.sd
    return(temp.scaled)
  }
}

gam_fit_full_func<- function(df, response){
  if(response == "Presence"){
    gam.mod0<- gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = binomial(link = logit), select = TRUE)
    return(gam.mod0)
  } 
  
  if(response == "Biomass"){
    gam.mod0<- gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = gaussian, select = TRUE)
    return(gam.mod0)
  }
}

gam_fit_red_func<- function(df, response){
  if(response == "Presence"){
    gam.mod0<- gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = binomial(link = logit), select = TRUE)
    return(gam.mod0)
  } 
  
  if(response == "Biomass"){
    gam.mod0<- gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = df, family = gaussian, select = TRUE)
    return(gam.mod0)
  }
}

resids_map_func<- function(test.data, response, predicted, type) {
  # Preliminary spatial stuff
  # Spatial projections
  proj.wgs84<- "+init=epsg:4326" #WGS84
  proj.utm<- "+init=epsg:2960" #UTM 19
  
  # NELME domaine
  nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
  st_crs(nelme)<- proj.wgs84
  
  #Bounds
  xlim.use<- c(-77, -65)
  ylim.use<- c(35, 45)
  
  states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
  provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")
  
  us <- raster::getData("GADM",country="USA",level=1)
  us.states <- us[us$NAME_1 %in% states,]
  us.states <- gSimplify(us.states, tol = 0.025, topologyPreserve = TRUE)
  canada <- raster::getData("GADM",country="CAN",level=1)
  ca.provinces <- canada[canada$NAME_1 %in% provinces,]
  ca.provinces <- gSimplify(ca.provinces, tol = 0.025, topologyPreserve = TRUE)
  
  us.states.f<- fortify(us.states, NAME_1)
  ca.provinces.f<- fortify(ca.provinces, NAME_1)
  
  if(response == "Presence" & type == "response"){
    resids<- predicted - test.data$PRESENCE.BIOMASS
    dat.resids<- data.frame("year" = test.data$EST_YEAR, "x" = test.data$DECDEG_BEGLON, "y" = test.data$DECDEG_BEGLAT, "resid" = resids)
    
    plots.out<- vector("list", length(unique(dat.resids$year)))
    
    zlim.use<- c(min(dat.resids$resid, na.rm = T), max(dat.resids$resid, na.rm = T))
    
    for(k in 1:length(unique(dat.resids$year))){
      data.use<- dat.resids[dat.resids$year == unique(dat.resids$year)[k],]
      pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$resid))
      pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                              xo=seq(-87.99457, -57.4307, length = 115),
                              yo=seq(22.27352, 48.11657, length = 133))
      pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
      pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
      
      # Clip to nelme
      pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
      coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
      row.names(coords.keep)<- NULL
      pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
      names(pred.df.use)<- c("X", "Y", "z")
      
      # Plot
      plots.out[[k]]<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(name = paste(unique(dat.resids$year)[k], " pred - obs", sep = ""), low = "blue", high = "red", mid = "white", midpoint = 0.0, na.value = "white", limits = c(-1,1)) +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("Lat") +
        scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) + 
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "in"))
    }
    
    out<- plot_grid(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], plots.out[[5]], nrow = 2, ncol = 3, scale = 1)
  }
  
  if(response == "Biomass" & type == "response"){
    resids<- predicted - test.data$BIOMASS
    dat.resids<- data.frame("year" = test.data$EST_YEAR, "x" = test.data$DECDEG_BEGLON, "y" = test.data$DECDEG_BEGLAT, "resid" = resids)
    
    plots.out<- vector("list", length(unique(dat.resids$year)))
    
    zlim.use<- c(min(dat.resids$resid, na.rm = T), max(dat.resids$resid, na.rm = T))
    
    for(k in 1:length(unique(dat.resids$year))){
      data.use<- dat.resids[dat.resids$year == unique(dat.resids$year)[k],]
      pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$resid))
      pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                              xo=seq(-87.99457, -57.4307, length = 115),
                              yo=seq(22.27352, 48.11657, length = 133))
      pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
      pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
      
      # Clip to nelme
      pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
      coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
      row.names(coords.keep)<- NULL
      pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
      names(pred.df.use)<- c("X", "Y", "z")
      
      # Plot
      plots.out[[k]]<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(name = paste(unique(dat.resids$year)[k], " pred - obs", sep = ""), low = "blue", high = "red", mid = "white", midpoint = 0.0, na.value = "white") +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("Lat") +
        scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) + 
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "in"))
    }
    
    out<- plot_grid(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], plots.out[[5]], nrow = 2, ncol = 3, scale = 1)
  }
  
  return(out)
}

predict_func<- function(mod.fitted.p, mod.fitted.b, response, percentile, test.data) {
  temp<- dplyr::select(test.data, one_of(c("DEPTH.Scale", percentile)))
  test.data<- data.frame(na.omit(temp))
  out.p<- round(as.numeric(predict.gam(mod.fitted.p, newdata = test.data, type = "response", se.fit = TRUE)$fit), 3)
  
  if(response == "Biomass"){
    out.b<- exp(round(as.numeric(predict.gam(mod.fitted.b, newdata = test.data, type = "response", se.fit = TRUE)$fit), 3))
    out.c<- out.p * out.b
    return(out.c)
  } else {
    return(out.p)
  }
}

pred_ranges_func<- function(predicted){
  pred.ranges<- data.frame("Min.Pred" = min(predicted, na.rm = T), "Max.Pred" = max(predicted, na.rm = T), "Mean.Pred" = mean(predicted, na.rm = T))
  return(pred.ranges)
}

auc_func<- function(test.data, predicted) {
  library(ROCR)
  temp<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS")))
  test.data<- data.frame(na.omit(temp))
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
    return(performance(dat, measure = "auc")@y.values[[1]])
  }
}

rmse_func<- function(test.data, predicted, response) {
  if(response == "Presence"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
    if(all(test.data$PRESENCE.BIOMASS == 0)){
      return(NA)
    } else {
      dat<- prediction(predictions = predicted, labels = test.data[,col.ind])
      return(performance(dat, measure = "rmse")@y.values[[1]])
    }
  }
  
  if(response == "Biomass"){
    test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) 
    test.data<- data.frame(na.omit(test.data))
    col.ind<- which(colnames(test.data) == "BIOMASS")
    if(all(test.data$BIOMASS == 0)){
      return(NA)
    } else {
      return(nrmse(sim = as.numeric(predicted), obs = test.data$BIOMASS))
    }
  }
}

calib_stat_func<- function(test.data, predicted){
  test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
  test.data<- data.frame(na.omit(test.data))
  col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    calib.stat<- round(calibration(x = test.data[,col.ind], p = predicted)@statistic, 3)
    return(calib.stat)
  }
}

calib_plot_func<- function(test.data, predicted){
  test.data<- dplyr::select(test.data, one_of(c("DEPTH.Scale", "SEASONALMU.OISST.Scale", "PRESENCE.BIOMASS"))) 
  test.data<- data.frame(na.omit(test.data))
  col.ind<- which(colnames(test.data) == "PRESENCE.BIOMASS")
  if(all(test.data$PRESENCE.BIOMASS == 0)){
    return(NA)
  } else {
    calib.all<- calibration(x = test.data[,col.ind], p = predicted)
    calib.df<- calib.all@calibration
    calib.plot<- ggplot() +
      geom_point(data = calib.df, aes(x = pedicted_pobability_bin, y = observed_poportion)) +
      xlim(c(0, 1)) + 
      ylim(c(0, 1)) +
      ylab("Proportion of Observed Occurrences") +
      xlab("Predicted Probability of Occurrence") +
      geom_abline(intercept = 0, linetype = "dashed")
    return(calib.plot)
  }
}

# Temp data exploration ---------------------------------------------------------------
## Get oisst data
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19
oisst.dat<- raster::stack("~/GitHub/SDMSeasonality/Data/Projected_OISST.grd")
proj4string(oisst.dat)<- proj.utm
oisst.dat<- projectRaster(oisst.dat, crs = proj.wgs84)
names.use<- names(oisst.dat)

# Need to get climatology from the OISST data -- set up OISST stack as time series
oisst.min<- gsub("X", "", min(names(oisst.dat)))
oisst.min.date<- as.Date(gsub("[.]", "-", oisst.min))
oisst.max<- gsub("X", "", max(names(oisst.dat)))
oisst.max.date<- as.Date(gsub("[.]", "-", oisst.max))

oisst.dates<- seq.Date(from = oisst.min.date, to = oisst.max.date, by = "day")
oisst.dat<- setZ(oisst.dat, oisst.dates)

# Mask 
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
nelme<- readShapePoly("~/GitHub/SDMSeasonality/Data/nelme.shp")
proj4string(nelme)<- proj.wgs84

oisst.m<- mask(oisst.dat, nelme)

# Filter full years only
oisst.m<- oisst.m[[which(getZ(oisst.m) >= "1982-01-01" & getZ(oisst.m) <= "2016-12-31")]]

## Baseline
## Need baseline daily average temp across years between 1982-2011 for every day
oisst.m.base<- oisst.m[[which(getZ(oisst.m) >= "1982-01-01" & getZ(oisst.m) <= "2011-12-31")]]
oisst.m.base<- setZ(oisst.m.base, seq.Date(from = as.Date("1982-01-01"), to = as.Date("2011-12-31"), by = "day"))

dates.unique<- unique(format(getZ(oisst.m), "%m-%d"))
daily.means<- stack(lapply(seq(length(dates.unique)), function(x) calc(oisst.m.base[[which(format(getZ(oisst.m.base), "%m-%d") == dates.unique[x])]], fun = mean)))
names(daily.means)<- dates.unique
daily.sd<- stack(lapply(seq(length(dates.unique)), function(x) calc(oisst.m.base[[which(format(getZ(oisst.m.base), "%m-%d") == dates.unique[x])]], fun = sd)))
names(daily.sd)<- dates.unique

# Alright, now substract each daily OISST from the daily climatology to get the anomaly
anom.type<- "Non.standardized"
daily.anoms<- switch(anom.type, 
                     Standardized = stack(lapply(seq(1:nlayers(oisst.m)), function(x) (oisst.m[[x]] - daily.means[[match(format(getZ(oisst.m)[x], "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.means))))]])/daily.sd[[match(format(getZ(oisst.m)[x], "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.sd))))]])),
                     Non.standardized = stack(lapply(seq(1:nlayers(oisst.m)), function(x) (oisst.m[[x]] - daily.means[[match(format(getZ(oisst.m)[x], "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.means))))]]))))
names(daily.anoms)<- getZ(oisst.m)

## Get the data ready for plotting
# Getting the anomalies
ts.wide.daily<- do.call("cbind", lapply(seq(1:nlayers(daily.anoms)), function(x) as.data.frame(daily.anoms[[x]], xy = TRUE)))
ts.df.daily<- ts.wide.daily %>%
  subset(., select=which(!duplicated(names(.)))) %>%
  gather(., Year, SST, -x, -y) 
names(ts.df.daily)[3]<- "Date"
ts.df.daily$Date<- gsub("X", "", gsub("[.]", "-", ts.df.daily$Date))
ts.df.daily<- ts.df.daily %>%
  separate(., Date, c("YYYY", "MM", "DD")) %>%
  group_by(., YYYY, MM, DD) %>%
  summarize_at(., "SST", mean, na.rm = T)

# Keep going for individual region...
ts.df.dailymu<- ts.df.daily %>%
  mutate(., "Plot.Date" = as.Date(paste(YYYY, MM, DD, sep = "-"))) %>%
  data.frame

# Smooth daily values
# Creat zoo series
zoo.dates<- as.Date(ts.df.dailymu$Plot.Date)
zoo.vals<- zoo(ts.df.dailymu$SST, zoo.dates)
zoo.ts<- as.ts(zoo.vals)
ts.df.dailymu$smoothed<- as.numeric(ma(zoo.vals, order = 15, centre = FALSE))

ts.df.monthlymu<- ts.df.dailymu %>%
  group_by(YYYY, MM) %>%
  dplyr::summarize(., Mean.SST = mean(SST, na.rm =T)) %>%
  mutate(., Plot.Date = as.Date(paste(YYYY, MM, "15", sep = "-"))) %>%
  data.frame

ts.df.yearlymu<- ts.df.dailymu %>%
  group_by(YYYY) %>%
  dplyr::summarize(., Mean.SST = mean(SST, na.rm = T)) %>%
  mutate(., Plot.Date = as.Date(paste(YYYY, "06", "15", sep = "-"))) %>%
  filter(., YYYY != format(Sys.time(), "%Y")) %>%
  data.frame()

ts.df.yearlymu$Year.Model<- as.integer(ts.df.yearlymu$YYYY) - 1982

# Plots
# Full trend line
sst.anom.lm.full<- lm(Mean.SST ~ Year.Model, data = ts.df.yearlymu)
summary(sst.anom.lm.full)
adj.r2.full<- paste("Adj R2 = ", round(summary(sst.anom.lm.full)$adj.r.squared, 3), sep = "")

# Fitted model formula
my.formula.full<- paste("Mean.SST = ", signif(round(sst.anom.lm.full$coef[[1]], 3), 5), " + ", signif(round(sst.anom.lm.full$coef[[2]], 3), 5), "*Year", sep = "")
base.ts<- ggplot(data = subset(ts.df.dailymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date("2016-12-31", format = "%Y-%m-%d")), aes(x = Plot.Date, y = smoothed)) + 
  geom_line(col = "#d9d9d9", lwd = 0.15) +
  geom_point(data = subset(ts.df.yearlymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date("2016-12-31", format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), col = "black") +
  geom_line(data = subset(ts.df.yearlymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date("2016-12-31", format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), col = "black", lwd = 0.25) +
  geom_smooth(data = subset(ts.df.yearlymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date("2016-12-31", format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), method = "lm", formula = y ~ x, col = "#969696", size = 0.75, se = FALSE) +
  ylim(c(-3, 4)) +
  ylab("SST Anomaly (1982-2011 baseline)") + 
  xlab("Year") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(aes(x = as.Date("1987-06-15"), y = -2.5, label = my.formula.full), col = "#969696") +
  geom_text(aes(x = as.Date("1997-06-15"), y = -2.5, label = adj.r2.full), col = "#969696") 

# 2004-2016 trend line
sst.anom.lm.plot2<- lm(Mean.SST ~ Year.Model, data = subset(ts.df.yearlymu, YYYY >= 2004 & YYYY <= 2016))
summary(sst.anom.lm.plot2)
adj.r2.plot2<- paste("Adj R2 = ", round(summary(sst.anom.lm.plot2)$adj.r.squared, 3), sep = "")

# Fitted model formula
my.formula.plot2<- paste("Mean.SST = ", signif(round(sst.anom.lm.plot2$coef[[1]], 3), 5), " + ", signif(round(sst.anom.lm.plot2$coef[[2]], 3), 5), "*Year", sep = "")
base.ts.form<- base.ts +
  geom_smooth(data = subset(ts.df.yearlymu, Plot.Date >= as.Date("2004-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date("2016-12-31", format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), method = "lm", formula = y ~ x, col = "red", size = 0.75, se = FALSE) +
  geom_text(aes(x = as.Date("1987-06-15"), y = -2.8, label = my.formula.plot2), col = "red") +
  geom_text(aes(x = as.Date("1997-06-15"), y = -2.8, label = adj.r2.plot2), col = "red") 
base.ts.form
ggsave("~/GitHub/SDMSeasonality/Results/Figures/SSTTrend.jpg", width = 8, height = 6, units = "in")
dev.off()

# Data preparation --------------------------------------------------------
# Data path
dat.path<- "~/GitHub/SDMSeasonality/Data/model.dat.rds"

# Read it in, filter to species with at least 250 observations and observed in at least 10 years, do some quick formatting to fit the GAM
dat<- readRDS(dat.path)

dat.filter.a<- dat %>%
  group_by(., SVSPP) %>%
  summarize_at(vars(PRESENCE), sum, na.rm = T) %>%
  filter(., PRESENCE >= 250)

dat.filter.b<- dat %>%
  filter(., SVSPP %in% dat.filter.a$SVSPP) %>%
  group_by(., SEASON, SVSPP) %>% 
  filter(., PRESENCE > 0) %>%
  summarize_at(vars(EST_YEAR), n_distinct, na.rm = T) %>%
  filter(., EST_YEAR >= 10) 
dat.filter.b$Count<- ifelse(dat.filter.b$EST_YEAR >= 10, 1, 0)
dat.filter.b<- dat.filter.b %>%
  group_by(., SVSPP) %>%
  summarize_at(vars(Count), sum, na.rm = T) %>%
  filter(., Count == 2)
  
dat<- dat %>%
  filter(., SVSPP %in% dat.filter.b$SVSPP)

# Add species common name
species.names<- read_csv("~/GitHub/SDMSeasonality/Data/svspp_spp names.csv") %>%
  dplyr::select(., -X4) 
species.names$SVSPP<- as.numeric(species.names$SVSPP)
dat<- dat %>%
  left_join(., species.names) 
dat<- dat[-which(grepl("UNCL", dat$COMNAME)),]
dat<- dat[-which(grepl("UNKNOWN 01", dat$COMNAME)),]
dat<- dat[-which(is.na(dat$COMNAME)),]

# Training vs. testing
train.start<- "1982-01-01"
train.end<- "2003-12-31"
test.start<- "2004-01-01"
test.end<- "2016-01-01"

# Durations?
as.Date(train.end) - as.Date(train.start)
as.Date(test.end) - as.Date(test.start)

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither"))

# Bottom trawl strata
bstrat<- st_read("~/GitHub/SDMSeasonality/Data/BottomTrawlStrata/BTS_Strata.shp")

# Get names of strata
bstrat.names<- unique(bstrat$STRATA)

# Reduce dataset
dat<- dat[dat$STRATUM %in% bstrat.names,]

# Training data
dat.train.f<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST)))
dat.train.f$PRESENCE.ABUNDANCE<- ifelse(dat.train.f$ABUNDANCE > 0, 1, 0)
dat.train.f$WT.ABUNDANCE<- ifelse(dat.train.f$PRESENCE.ABUNDANCE == 0, 1, dat.train.f$ABUNDANCE)
dat.train.f$PRESENCE.BIOMASS<- ifelse(dat.train.f$BIOMASS > 0, 1, 0)

# Temps to rescale other time periods
base.depth.mean.f<- mean(abs(dat.train.f$DEPTH))
base.depth.sd.f<- sd(abs(dat.train.f$DEPTH))
base.shelf.mean.f<- mean(dat.train.f$SHELF_POS)
base.shelf.sd.f<- sd(dat.train.f$SHELF_POS)
base.temp.mean.f<- mean(dat.train.f$SEASONALMU.OISST, na.rm = T)
base.temp.sd.f<- sd(dat.train.f$SEASONALMU.OISST, na.rm = T)
fall.rescale.df<- data.frame(SEASON = "FALL", mean.t = base.temp.mean.f, sd.t = base.temp.sd.f, mean.depth = base.depth.mean.f, sd.depth = base.depth.sd.f, mean.shelf = base.shelf.mean.f, sd.shelf = base.shelf.sd.f)

dat.train.s<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST)))
dat.train.s$PRESENCE.ABUNDANCE<- ifelse(dat.train.s$ABUNDANCE > 0, 1, 0)
dat.train.s$WT.ABUNDANCE<- ifelse(dat.train.s$PRESENCE.ABUNDANCE == 0, 1, dat.train.s$ABUNDANCE)
dat.train.s$PRESENCE.BIOMASS<- ifelse(dat.train.s$BIOMASS > 0, 1, 0)

# Temps to rescale other variables
base.depth.mean.sp<- mean(abs(dat.train.s$DEPTH))
base.depth.sd.sp<- sd(abs(dat.train.s$DEPTH))
base.shelf.mean.sp<- mean(dat.train.s$SHELF_POS)
base.shelf.sd.sp<- sd(dat.train.s$SHELF_POS)
base.temp.mean.sp<- mean(dat.train.s$SEASONALMU.OISST, na.rm = T)
base.temp.sd.sp<- sd(dat.train.s$SEASONALMU.OISST, na.rm = T)
spring.rescale.df<- data.frame(SEASON = "SPRING", mean.t = base.temp.mean.sp, sd.t = base.temp.sd.sp, mean.depth = base.depth.mean.sp, sd.depth = base.depth.sd.sp, mean.shelf = base.shelf.mean.sp, sd.shelf = base.shelf.sd.sp)

## Testing dataframes
dat.test.f<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
dat.test.f$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
dat.test.f$SHELF_POS.Scale<- mapply(temp.scale, dat.test.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
dat.test.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.f$SEASONALMU.OISST, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
dat.test.f$SEASON<- factor(dat.test.f$SEASON, levels = c("FALL", "SPRING"))

dat.test.f<- dat.test.f 
dat.test.f$PRESENCE.ABUNDANCE<- ifelse(dat.test.f$ABUNDANCE > 0, 1, 0)
dat.test.f$WT.ABUNDANCE<- ifelse(dat.test.f$PRESENCE.ABUNDANCE == 0, 1, dat.test.f$ABUNDANCE)
dat.test.f$PRESENCE.BIOMASS<- ifelse(dat.test.f$BIOMASS > 0, 1, 0)

## Testing dataframes
dat.test.s<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
dat.test.s$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.s$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
dat.test.s$SHELF_POS.Scale<- mapply(temp.scale, dat.test.s$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
dat.test.s$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.s$SEASONALMU.OISST, spring.rescale.df$mean.t, spring.rescale.df$sd.t)
dat.test.s$SEASON<- factor(dat.test.s$SEASON, levels = c("FALL", "SPRING"))

dat.test.s<- dat.test.s 
dat.test.s$PRESENCE.ABUNDANCE<- ifelse(dat.test.s$ABUNDANCE > 0, 1, 0)
dat.test.s$WT.ABUNDANCE<- ifelse(dat.test.s$PRESENCE.ABUNDANCE == 0, 1, dat.test.s$ABUNDANCE)
dat.test.s$PRESENCE.BIOMASS<- ifelse(dat.test.s$BIOMASS > 0, 1, 0)

# Create nested dataframes, one for testing, one for training
# Training
dat.train<- dat.train.f %>%
  bind_rows(., dat.train.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TRAIN.DATA") %>%
  arrange(COMNAME)

# Testing
dat.test<- dat.test.f %>%
  bind_rows(., dat.test.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TEST.DATA") %>%
  arrange(COMNAME) 

# Need base.preds, fut.preds, nevaD and nevaV
spring.preds = "~/GitHub/SDMSeasonality/Data/spring.rast.seasonality.rds"
spring.preds<- readRDS(spring.preds)

fall.preds = "~/GitHub/SDMSeasonality/Data/fall.rast.seasonality.rds"
fall.preds<- readRDS(fall.preds)

base.preds.sp<- spring.preds %>%
  dplyr::select(., -Spring.2055, -TRI, -SED.SIZE, -SED.TYPE, -SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.)))

base.preds.sp<- base.preds.sp %>%
  left_join(., spring.rescale.df, by = "SEASON")
base.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)

sst.loop<- paste("Spring.", seq(2004, 2017, by = 1), sep = "")
res.ind<- ncol(base.preds.sp)

for(i in seq_along(sst.loop)){
  sst.use<- sst.loop[i]
  col.ind<- which(colnames(base.preds.sp) == sst.use)
  base.preds.sp[,res.ind+1]<- mapply(temp.scale, base.preds.sp[,col.ind], spring.rescale.df$mean.t, spring.rescale.df$sd.t)
  names(base.preds.sp)[res.ind+1]<- paste(sst.use, ".Scale", sep = "")
  res.ind<- res.ind+1
}

base.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.sp$Baseline, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

base.preds.sp<- base.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds.f<- fall.preds %>%
  dplyr::select(., -Fall.2055, -TRI, -SED.SIZE, -SED.TYPE, -SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) 
base.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)

sst.loop<- paste("Fall.", seq(2004, 2017, by = 1), sep = "")
res.ind<- ncol(base.preds.f)

for(i in seq_along(sst.loop)){
  sst.use<- sst.loop[i]
  col.ind<- which(colnames(base.preds.f) == sst.use)
  base.preds.f[,res.ind+1]<- mapply(temp.scale, base.preds.f[,col.ind], fall.rescale.df$mean.t, fall.rescale.df$sd.t)
  names(base.preds.f)[res.ind+1]<- paste(sst.use, ".Scale", sep = "")
  res.ind<- res.ind+1
}

base.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.f$Baseline, fall.rescale.df$mean.t, fall.rescale.df$sd.t)

base.preds.f<- base.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.f %>%
  bind_rows(., base.preds.sp)

dat.full<- dat.train %>%
  left_join(., dat.test, by = c("COMNAME", "SEASON"))

dat.full<- dat.full[complete.cases(dat.full$COMNAME, dat.full$SEASON),]

# Model fitting and prediction -----------------------------------------------------------
species<- unique(dat.full$COMNAME)
scenarios<- c("SeasonalSST", "SeasonalSST.SeasonFactor", "SeasonalSST.SeasonInteraction", "SeasonalSST.Independent.Fall", "SeasonalSST.Independent.Spring")

mod.results<- data.frame("Species" = species, 
                         "SeasonalSST.DevExp.P" = rep(NA, length(species)), 
                         "SeasonalSST.AUC.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.AUC.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.RMSE.P.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.RMSE.P.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Calib.SDM.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Calib.SDM.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.DevExp.B" = rep(NA, length(species)), 
                         "SeasonalSST.RMSE.B.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.RMSE.B.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.COG.Long.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.COG.Lat.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.COG.Long.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.COG.Lat.Spring" = rep(NA, length(species)),  
                         "SeasonalSST.SeasonFactor.DevExp.P" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.AUC.Fall" = rep(NA, length(species)),
                         "SeasonalSST.SeasonFactor.AUC.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.RMSE.P.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.RMSE.P.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.Calib.SDM.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.Calib.SDM.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.DevExp.B" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.RMSE.B.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.RMSE.B.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.COG.Long.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.COG.Lat.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.COG.Long.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonFactor.COG.Lat.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.DevExp.P" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.AUC.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.AUC.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.RMSE.P.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.RMSE.P.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.Calib.SDM.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.Calib.SDM.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.DevExp.B" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.RMSE.B.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.RMSE.B.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.COG.Long.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.COG.Lat.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.COG.Long.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.SeasonInteraction.COG.Lat.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.DevExp.P" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.AUC.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.AUC.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.RMSE.P.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.RMSE.P.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.Calib.SDM.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.Calib.SDM.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.DevExp.B" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.RMSE.B.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.RMSE.B.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.COG.Long.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.COG.Lat.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.COG.Long.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Fall.COG.Lat.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.DevExp.P" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.AUC.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.AUC.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.RMSE.P.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.RMSE.P.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.Calib.SDM.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.Calib.SDM.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.DevExp.B" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.RMSE.B.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.RMSE.B.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.COG.Long.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.COG.Lat.Spring" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.COG.Long.Fall" = rep(NA, length(species)), 
                         "SeasonalSST.Independent.Spring.COG.Lat.Fall" = rep(NA, length(species)))

# Prediction stuff
# Maps
# Spatial projections
proj.wgs84<- "+init=epsg:4326" #WGS84
proj.utm<- "+init=epsg:2960" #UTM 19

# NELME domain
nelme<- st_read("~/GitHub/SDMSeasonality/Data/NELME_clipped.shp")
st_crs(nelme)<- proj.wgs84
nelme.sp<- as(nelme, "Spatial")

#Bounds
xlim.use<- c(-77, -65)
ylim.use<- c(35, 45)

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us <- raster::getData("GADM",country="USA",level=1)
us.states <- us[us$NAME_1 %in% states,]
us.states <- gSimplify(us.states, tol = 0.025, topologyPreserve = TRUE)
canada <- raster::getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]
ca.provinces <- gSimplify(ca.provinces, tol = 0.025, topologyPreserve = TRUE)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

pred.cols.f<- c("SEASON", "DEPTH.Scale", paste("Fall.", seq(from = 2004, to = 2017, by = 1), ".Scale", sep = ""))
pred.cols.s<- c("SEASON", "DEPTH.Scale", paste("Spring.", seq(from = 2004, to = 2017, by = 1), ".Scale", sep = ""))
pred.years<- seq(from = 2004, to = 2017, by = 1)

set.seed(13)

for(i in seq_along(unique(dat.full$COMNAME))){
  # Species
  spp.use<- unique(dat.full$COMNAME)[i]
  
  # Get the data
  dat.temp<- dat.full[dat.full$COMNAME %in% spp.use,]
  
  res.ind<- i
  
  # Loop over species
  for(j in seq_along(scenarios)){
    
    scenario.use<- scenarios[j]
    
    # Get the data and fit the model
    if(scenario.use == "SeasonalSST"){
      dat.train<- bind_rows(dat.temp$TRAIN.DATA[[1]], dat.temp$TRAIN.DATA[[2]])
      dat.test<- bind_rows(dat.temp$TEST.DATA[[1]], dat.temp$TEST.DATA[[2]])
      dat.test$SEASON<- c(rep("FALL", nrow(dat.temp$TEST.DATA[[1]])), rep("SPRING", nrow(dat.temp$TEST.DATA[[2]])))
      
      # Fit the models
      mod.p<- try(gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = binomial(link = logit), select = TRUE), silent = TRUE)
      mod.b<- try(gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = gaussian, select = TRUE), silent = TRUE)
      
      if(any(class(mod.p) == "try-error") | any(class(mod.b) == "try-error") | is.infinite(summary(mod.p)$dev.expl)){
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
        next()
      } else {
        # Make predictions
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        pred.data.f<- data.frame(na.omit(temp.f))
        
        pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.f<- pred.p.f * pred.logb.f
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        pred.data.s<- data.frame(na.omit(temp.s))
        
        pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.s<- pred.p.s * pred.logb.s
        
        # Evaluate and validate model: DevianceExplained, AUC, RMSE, Calibration
        # Presence/absence component
        mod.results$SeasonalSST.DevExp.P[i]<- round(summary(mod.p)$dev.expl, 3)
        
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        test.data.f<- data.frame(na.omit(temp.f))
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        test.data.s<- data.frame(na.omit(temp.s))
        
        # AUC
        col.ind<- which(colnames(test.data.f) == "PRESENCE.BIOMASS")
        dat.f<- try(prediction(predictions = pred.p.f, labels = test.data.f[,col.ind]), silent = TRUE)
        dat.s<- try(prediction(predictions = pred.p.s, labels = test.data.s[,col.ind]), silent = TRUE)
        
        # Exit out if we can't at least get an AUC value
        if(class(dat.f) == "try-error" | class(dat.s) == "try-error") {
          next()
        } 
        
        mod.results$SeasonalSST.AUC.Fall[i]<- performance(dat.f, measure = "auc")@y.values[[1]]
        mod.results$SeasonalSST.AUC.Spring[i]<- performance(dat.s, measure = "auc")@y.values[[1]]
        
        # RMSE
        mod.results$SeasonalSST.RMSE.P.Fall[i]<- performance(dat.f, measure = "rmse")@y.values[[1]]
        mod.results$SeasonalSST.RMSE.P.Spring[i]<- performance(dat.s, measure = "rmse")@y.values[[1]]
        
        # Calibration
        mod.results$SeasonalSST.Calib.SDM.Fall[i]<- round(calibration(x = test.data.f[,col.ind], p = pred.p.f)@statistic, 3)
        mod.results$SeasonalSST.Calib.SDM.Spring[i]<- round(calibration(x = test.data.s[,col.ind], p = pred.p.s)@statistic, 3)
        
        # Biomass component
        mod.results$SeasonalSST.DevExp.B[i]<- round(summary(mod.b)$dev.expl, 3)
        
        # RMSE
        test.data.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "FALL")
        test.data.f<- data.frame(na.omit(test.data.f))
        
        test.data.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        test.data.s<- data.frame(na.omit(test.data.s))
        
        mod.results$SeasonalSST.RMSE.B.Fall[i]<- RMSE(as.numeric(pred.b.f), test.data.f$BIOMASS, "NRMSE")
        mod.results$SeasonalSST.RMSE.B.Spring[i]<- RMSE(as.numeric(pred.b.s), test.data.s$BIOMASS, "NRMSE")
        
        # Prediction maps...
        preds.f<- data.frame("x" = base.preds$Data[[1]]$x, "y" = base.preds$Data[[1]]$y)
        res.ind.f<- ncol(preds.f)
        preds.s<- data.frame("x" = base.preds$Data[[2]]$x, "y" = base.preds$Data[[2]]$y)
        res.ind.s<- ncol(preds.s)
        
        for(k in seq_along(pred.years)){
          fall.use<- c("x", "y", "DEPTH.Scale", pred.cols.f[grepl(pred.years[k], pred.cols.f)])
          spring.use<- c("x", "y", "DEPTH.Scale", pred.cols.s[grepl(pred.years[k], pred.cols.s)])
          
          pred.df.f<- data.frame(base.preds$Data[[1]])
          pred.df.f<- pred.df.f[,colnames(pred.df.f) %in% fall.use]
          names(pred.df.f)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.df.s<- data.frame(base.preds$Data[[2]])
          pred.df.s<- pred.df.s[,colnames(pred.df.s) %in% spring.use]
          names(pred.df.s)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3))
          preds.f[,res.ind.f+1]<- pred.p.f * pred.logb.f
          colnames(preds.f)[res.ind.f+1]<- paste("Pred.", pred.years[k], sep = "")
          
          pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3))
          preds.s[,res.ind.s+1]<- pred.p.s * pred.logb.s
          colnames(preds.s)[res.ind.s+1]<- paste("Pred.", pred.years[k], sep = "")
          
          res.ind.f<- res.ind.f+1
          res.ind.s<- res.ind.s+1
        }
        
        # Prediction maps
        SeasonalSST.pred.out.f<- data.frame("x" = preds.f$x, "y" = preds.f$y, "pred" = rowMeans(preds.f[,3:16]))
        SeasonalSST.pred.out.f$pred[!is.na(SeasonalSST.pred.out.f$pred) & SeasonalSST.pred.out.f$pred > 1000]<- NA
        SeasonalSST.pred.out.s<- data.frame("x" = preds.s$x, "y" = preds.s$y, "pred" = rowMeans(preds.s[,3:16]))
        SeasonalSST.pred.out.s$pred[!is.na(SeasonalSST.pred.out.s$pred) & SeasonalSST.pred.out.s$pred > 1000]<- NA
        
        # Center of Biomass
        cog.fall<- COGravity(x = SeasonalSST.pred.out.f$x, y = SeasonalSST.pred.out.f$y, wt = SeasonalSST.pred.out.f$pred)
        cog.spring<- COGravity(x = SeasonalSST.pred.out.s$x, y = SeasonalSST.pred.out.s$y, wt = SeasonalSST.pred.out.s$pred)
        mod.results$SeasonalSST.COG.Long.Fall[i]<- as.numeric(cog.fall[1])
        mod.results$SeasonalSST.COG.Lat.Fall[i]<- as.numeric(cog.fall[3])
        mod.results$SeasonalSST.COG.Long.Spring[i]<- as.numeric(cog.spring[1])
        mod.results$SeasonalSST.COG.Lat.Spring[i]<- as.numeric(cog.spring[3])
        
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
      }
    } 
    
    if(scenario.use == "SeasonalSST.SeasonFactor"){
      dat.train<- bind_rows(dat.temp$TRAIN.DATA[[1]], dat.temp$TRAIN.DATA[[2]])
      dat.train$SEASON<- factor(c(rep("FALL", nrow(dat.temp$TRAIN.DATA[[1]])), rep("SPRING", nrow(dat.temp$TRAIN.DATA[[2]]))), levels = c("FALL", "SPRING"))
      dat.test<- bind_rows(dat.temp$TEST.DATA[[1]], dat.temp$TEST.DATA[[2]])
      dat.test$SEASON<- factor(c(rep("FALL", nrow(dat.temp$TEST.DATA[[1]])), rep("SPRING", nrow(dat.temp$TEST.DATA[[2]]))), levels = c("FALL", "SPRING"))
      
      # Fit the models
      mod.p<- try(gam(PRESENCE.BIOMASS ~ SEASON + s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = binomial(link = logit), select = TRUE), silent = TRUE)
      mod.b<- try(gam(BIOMASS.MOD ~ SEASON + s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = gaussian, select = TRUE), silent = TRUE)
      
      if(any(class(mod.p) == "try-error") | any(class(mod.b) == "try-error") | is.infinite(summary(mod.p)$dev.expl)){
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
        next()
      } else {
        # Make predictions
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        pred.data.f<- data.frame(na.omit(temp.f))
        
        pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.f<- pred.p.f * pred.logb.f
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        pred.data.s<- data.frame(na.omit(temp.s))
        
        pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.s<- pred.p.s * pred.logb.s
        
        # Evaluate and validate model: DevianceExplained, AUC, RMSE, Calibration
        # Presence/absence component
        mod.results$SeasonalSST.SeasonFactor.DevExp.P[i]<- round(summary(mod.p)$dev.expl, 3)
        
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(temp.f))
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(temp.s))
        
        # AUC
        col.ind<- which(colnames(test.data.f) == "PRESENCE.BIOMASS")
        dat.f<- try(prediction(predictions = pred.p.f, labels = test.data.f[,col.ind]), silent = TRUE)
        dat.s<- try(prediction(predictions = pred.p.s, labels = test.data.s[,col.ind]), silent = TRUE)
        
        # Exit out if we can't at least get an AUC value
        if(class(dat.f) == "try-error" | class(dat.s) == "try-error") {
          next()
        } 
        mod.results$SeasonalSST.SeasonFactor.AUC.Fall[i]<- performance(dat.f, measure = "auc")@y.values[[1]]
        mod.results$SeasonalSST.SeasonFactor.AUC.Spring[i]<- performance(dat.s, measure = "auc")@y.values[[1]]
        
        # RMSE
        mod.results$SeasonalSST.SeasonFactor.RMSE.P.Fall[i]<- performance(dat.f, measure = "rmse")@y.values[[1]]
        mod.results$SeasonalSST.SeasonFactor.RMSE.P.Spring[i]<- performance(dat.s, measure = "rmse")@y.values[[1]]
        
        # Calibration
        mod.results$SeasonalSST.SeasonFactor.Calib.SDM.Fall[i]<- round(calibration(x = test.data.f[,col.ind], p = pred.p.f)@statistic, 3)
        mod.results$SeasonalSST.SeasonFactor.Calib.SDM.Spring[i]<- round(calibration(x = test.data.s[,col.ind], p = pred.p.s)@statistic, 3)
        
        # Biomass component
        mod.results$SeasonalSST.SeasonFactor.DevExp.B[i]<- round(summary(mod.b)$dev.expl, 3)
        
        # RMSE
        test.data.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "FALL")
        test.data.f$SEASON<- factor(test.data.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(test.data.f))
        
        test.data.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        test.data.s<- factor(test.data.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(test.data.s))
        
        mod.results$SeasonalSST.SeasonFactor.RMSE.B.Fall[i]<- RMSE(as.numeric(pred.b.f), test.data.f$BIOMASS, "NRMSE")
        mod.results$SeasonalSST.SeasonFactor.RMSE.B.Spring[i]<- RMSE(as.numeric(pred.b.s), test.data.s$BIOMASS, "NRMSE")
        
        # Prediction maps...
        preds.f<- data.frame("x" = base.preds$Data[[1]]$x, "y" = base.preds$Data[[1]]$y)
        res.ind.f<- ncol(preds.f)
        preds.s<- data.frame("x" = base.preds$Data[[2]]$x, "y" = base.preds$Data[[2]]$y)
        res.ind.s<- ncol(preds.s)
        
        for(k in seq_along(pred.years)){
          fall.use<- c("x", "y", "DEPTH.Scale", pred.cols.f[grepl(pred.years[k], pred.cols.f)])
          spring.use<- c("x", "y", "DEPTH.Scale", pred.cols.s[grepl(pred.years[k], pred.cols.s)])
          
          pred.df.f<- data.frame(base.preds$Data[[1]])
          pred.df.f<- pred.df.f[,colnames(pred.df.f) %in% fall.use]
          pred.df.f$SEASON<- factor(rep("FALL", nrow(pred.df.f)), levels = c("FALL", "SPRING"))
          names(pred.df.f)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.df.s<- data.frame(base.preds$Data[[2]])
          pred.df.s<- pred.df.s[,colnames(pred.df.s) %in% spring.use]
          pred.df.s$SEASON<- factor(rep("SPRING", nrow(pred.df.s)), levels = c("FALL", "SPRING"))
          names(pred.df.s)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3))
          preds.f[,res.ind.f+1]<- pred.p.f * pred.logb.f
          colnames(preds.f)[res.ind.f+1]<- paste("Pred.", pred.years[k], sep = "")
          
          pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3))
          preds.s[,res.ind.s+1]<- pred.p.s * pred.logb.s
          colnames(preds.s)[res.ind.s+1]<- paste("Pred.", pred.years[k], sep = "")
          
          res.ind.f<- res.ind.f+1
          res.ind.s<- res.ind.s+1
        }
        
        # Prediction maps
        SeasonalSST.SeasonFactor.pred.out.f<- data.frame("x" = preds.f$x, "y" = preds.f$y, "pred" = rowMeans(preds.f[,3:16]))
        SeasonalSST.SeasonFactor.pred.out.f$pred[!is.na(SeasonalSST.SeasonFactor.pred.out.f$pred) & SeasonalSST.SeasonFactor.pred.out.f$pred > 1000]<- NA
        SeasonalSST.SeasonFactor.pred.out.s<- data.frame("x" = preds.s$x, "y" = preds.s$y, "pred" = rowMeans(preds.s[,3:16]))
        SeasonalSST.SeasonFactor.pred.out.s$pred[!is.na(SeasonalSST.SeasonFactor.pred.out.s$pred) & SeasonalSST.SeasonFactor.pred.out.s$pred > 1000]<- NA
        
        # Center of Biomass
        cog.fall<- COGravity(x = SeasonalSST.SeasonFactor.pred.out.f$x, y = SeasonalSST.SeasonFactor.pred.out.f$y, wt = SeasonalSST.SeasonFactor.pred.out.f$pred)
        cog.spring<- COGravity(x = SeasonalSST.SeasonFactor.pred.out.s$x, y = SeasonalSST.SeasonFactor.pred.out.s$y, wt = SeasonalSST.SeasonFactor.pred.out.s$pred)
        mod.results$SeasonalSST.SeasonFactor.COG.Long.Fall[i]<- as.numeric(cog.fall[1])
        mod.results$SeasonalSST.SeasonFactor.COG.Lat.Fall[i]<- as.numeric(cog.fall[3])
        mod.results$SeasonalSST.SeasonFactor.COG.Long.Spring[i]<- as.numeric(cog.spring[1])
        mod.results$SeasonalSST.SeasonFactor.COG.Lat.Spring[i]<- as.numeric(cog.spring[3])
        
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
      }
    } 
    
    if(scenario.use == "SeasonalSST.SeasonInteraction"){
      dat.train<- bind_rows(dat.temp$TRAIN.DATA[[1]], dat.temp$TRAIN.DATA[[2]])
      dat.train$SEASON<- factor(c(rep("FALL", nrow(dat.temp$TRAIN.DATA[[1]])), rep("SPRING", nrow(dat.temp$TRAIN.DATA[[2]]))), levels = c("FALL", "SPRING"))
      dat.test<- bind_rows(dat.temp$TEST.DATA[[1]], dat.temp$TEST.DATA[[2]])
      dat.test$SEASON<- factor(c(rep("FALL", nrow(dat.temp$TEST.DATA[[1]])), rep("SPRING", nrow(dat.temp$TEST.DATA[[2]]))), levels = c("FALL", "SPRING"))
      
      # Fit the models
      mod.p<- try(gam(PRESENCE.BIOMASS ~ SEASON + s(DEPTH.Scale, fx = FALSE, bs = 'cs', by = SEASON) + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs', by = SEASON), drop.unused.levels = T, data = dat.train, family = binomial(link = logit), select = TRUE), silent = TRUE)
      mod.b<- try(gam(BIOMASS.MOD ~ SEASON + s(DEPTH.Scale, fx = FALSE, bs = 'cs', by = SEASON) + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs', by = SEASON), drop.unused.levels = T, data = dat.train, family = gaussian, select = TRUE), silent = TRUE)
      
      if(any(class(mod.p) == "try-error") | any(class(mod.b) == "try-error") | is.infinite(summary(mod.p)$dev.expl)){
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
        next()
      } else {
        # Make predictions
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        pred.data.f<- data.frame(na.omit(temp.f))
        
        pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.f<- pred.p.f * pred.logb.f
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        pred.data.s<- data.frame(na.omit(temp.s))
        
        pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.s<- pred.p.s * pred.logb.s
        
        # Evaluate and validate model: DevianceExplained, AUC, RMSE, Calibration
        # Presence/absence component
        mod.results$SeasonalSST.SeasonInteraction.DevExp.P[i]<- round(summary(mod.p)$dev.expl, 3)
        
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(temp.f))
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(temp.s))
        
        # AUC
        col.ind<- which(colnames(test.data.f) == "PRESENCE.BIOMASS")
        dat.f<- try(prediction(predictions = pred.p.f, labels = test.data.f[,col.ind]), silent = TRUE)
        dat.s<- try(prediction(predictions = pred.p.s, labels = test.data.s[,col.ind]), silent = TRUE)
        
        # Exit out if we can't at least get an AUC value
        if(class(dat.f) == "try-error" | class(dat.s) == "try-error") {
          next()
        }
        mod.results$SeasonalSST.SeasonInteraction.AUC.Fall[i]<- performance(dat.f, measure = "auc")@y.values[[1]]
        mod.results$SeasonalSST.SeasonInteraction.AUC.Spring[i]<- performance(dat.s, measure = "auc")@y.values[[1]]
        
        # RMSE
        mod.results$SeasonalSST.SeasonInteraction.RMSE.P.Fall[i]<- performance(dat.f, measure = "rmse")@y.values[[1]]
        mod.results$SeasonalSST.SeasonInteraction.RMSE.P.Spring[i]<- performance(dat.s, measure = "rmse")@y.values[[1]]
        
        # Calibration
        mod.results$SeasonalSST.SeasonInteraction.Calib.SDM.Fall[i]<- round(calibration(x = test.data.f[,col.ind], p = pred.p.f)@statistic, 3)
        mod.results$SeasonalSST.SeasonInteraction.Calib.SDM.Spring[i]<- round(calibration(x = test.data.s[,col.ind], p = pred.p.s)@statistic, 3)
        
        # Biomass component
        mod.results$SeasonalSST.SeasonInteraction.DevExp.B[i]<- round(summary(mod.b)$dev.expl, 3)
        
        # RMSE
        test.data.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "FALL")
        test.data.f$SEASON<- factor(test.data.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(test.data.f))
        
        test.data.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        test.data.s$SEASON<- factor(test.data.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(test.data.s))
        
        mod.results$SeasonalSST.SeasonInteraction.RMSE.B.Fall[i]<- RMSE(as.numeric(pred.b.f), test.data.f$BIOMASS, "NRMSE")
        mod.results$SeasonalSST.SeasonInteraction.RMSE.B.Spring[i]<- RMSE(as.numeric(pred.b.s), test.data.s$BIOMASS, "NRMSE")
        
        # Prediction maps...
        preds.f<- data.frame("x" = base.preds$Data[[1]]$x, "y" = base.preds$Data[[1]]$y)
        res.ind.f<- ncol(preds.f)
        preds.s<- data.frame("x" = base.preds$Data[[2]]$x, "y" = base.preds$Data[[2]]$y)
        res.ind.s<- ncol(preds.s)
        
        for(k in seq_along(pred.years)){
          fall.use<- c("x", "y", "DEPTH.Scale", pred.cols.f[grepl(pred.years[k], pred.cols.f)])
          spring.use<- c("x", "y", "DEPTH.Scale", pred.cols.s[grepl(pred.years[k], pred.cols.s)])
          
          pred.df.f<- data.frame(base.preds$Data[[1]])
          pred.df.f<- pred.df.f[,colnames(pred.df.f) %in% fall.use]
          pred.df.f$SEASON<- factor(rep("FALL", nrow(pred.df.f)), levels = c("FALL", "SPRING"))
          names(pred.df.f)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.df.s<- data.frame(base.preds$Data[[2]])
          pred.df.s<- pred.df.s[,colnames(pred.df.s) %in% spring.use]
          pred.df.s$SEASON<- factor(rep("SPRING", nrow(pred.df.s)), levels = c("FALL", "SPRING"))
          names(pred.df.s)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3))
          preds.f[,res.ind.f+1]<- pred.p.f * pred.logb.f
          colnames(preds.f)[res.ind.f+1]<- paste("Pred.", pred.years[k], sep = "")
          
          pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3))
          preds.s[,res.ind.s+1]<- pred.p.s * pred.logb.s
          colnames(preds.s)[res.ind.s+1]<- paste("Pred.", pred.years[k], sep = "")
          
          res.ind.f<- res.ind.f+1
          res.ind.s<- res.ind.s+1
        }
        
        # Prediction maps
        SeasonalSST.SeasonInteraction.pred.out.f<- data.frame("x" = preds.f$x, "y" = preds.f$y, "pred" = rowMeans(preds.f[,3:16]))
        SeasonalSST.SeasonInteraction.pred.out.f$pred[!is.na(SeasonalSST.SeasonInteraction.pred.out.f$pred) &  SeasonalSST.SeasonInteraction.pred.out.f$pred > 1000]<- NA
        SeasonalSST.SeasonInteraction.pred.out.s<- data.frame("x" = preds.s$x, "y" = preds.s$y, "pred" = rowMeans(preds.s[,3:16]))
        SeasonalSST.SeasonInteraction.pred.out.s$pred[!is.na(SeasonalSST.SeasonInteraction.pred.out.s$pred) &  SeasonalSST.SeasonInteraction.pred.out.s$pred > 1000]<- NA
        
        # Center of Biomass
        cog.fall<- COGravity(x = SeasonalSST.SeasonInteraction.pred.out.f$x, y = SeasonalSST.SeasonInteraction.pred.out.f$y, wt = SeasonalSST.SeasonInteraction.pred.out.f$pred)
        cog.spring<- COGravity(x = SeasonalSST.SeasonInteraction.pred.out.s$x, y = SeasonalSST.SeasonInteraction.pred.out.s$y, wt = SeasonalSST.SeasonInteraction.pred.out.s$pred)
        mod.results$SeasonalSST.SeasonInteraction.COG.Long.Fall[i]<- as.numeric(cog.fall[1])
        mod.results$SeasonalSST.SeasonInteraction.COG.Lat.Fall[i]<- as.numeric(cog.fall[3])
        mod.results$SeasonalSST.SeasonInteraction.COG.Long.Spring[i]<- as.numeric(cog.spring[1])
        mod.results$SeasonalSST.SeasonInteraction.COG.Lat.Spring[i]<- as.numeric(cog.spring[3])
        
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
      }
    } 
    
    if(scenario.use == "SeasonalSST.Independent.Fall"){
      dat.train<- bind_rows(dat.temp$TRAIN.DATA[[1]])
      dat.train$SEASON<- factor(rep("FALL", nrow(dat.train)), levels = c("FALL", "SPRING"))
      dat.test<- bind_rows(dat.temp$TEST.DATA[[1]], dat.temp$TEST.DATA[[2]])
      dat.test$SEASON<- factor(c(rep("FALL", nrow(dat.temp$TEST.DATA[[1]])), rep("SPRING", nrow(dat.temp$TEST.DATA[[2]]))), levels = c("FALL", "SPRING"))
      
      # Fit the models
      mod.p<- try(gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = binomial(link = logit), select = TRUE), silent = TRUE)
      mod.b<- try(gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = gaussian, select = TRUE), silent = TRUE)
      
      if(any(class(mod.p) == "try-error") | any(class(mod.b) == "try-error")| is.infinite(summary(mod.p)$dev.expl)){
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
        next()
      } else {
        # Make predictions
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        pred.data.f<- data.frame(na.omit(temp.f))
        
        pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.f<- pred.p.f * pred.logb.f
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        pred.data.s<- data.frame(na.omit(temp.s))
        
        pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.s<- pred.p.s * pred.logb.s
        
        # Evaluate and validate model: DevianceExplained, AUC, RMSE, Calibration
        # Presence/absence component
        mod.results$SeasonalSST.Independent.Fall.DevExp.P[i]<- round(summary(mod.p)$dev.expl, 3)
        
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(temp.f))
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(temp.s))
        
        # AUC
        col.ind<- which(colnames(test.data.f) == "PRESENCE.BIOMASS")
        dat.f<- try(prediction(predictions = pred.p.f, labels = test.data.f[,col.ind]), silent = TRUE)
        dat.s<- try(prediction(predictions = pred.p.s, labels = test.data.s[,col.ind]), silent = TRUE)
        
        # Exit out if we can't at least get an AUC value
        if(class(dat.f) == "try-error" | class(dat.s) == "try-error") {
          next()
        } 
        mod.results$SeasonalSST.Independent.Fall.AUC.Spring[i]<- performance(dat.s, measure = "auc")@y.values[[1]]
        mod.results$SeasonalSST.Independent.Fall.AUC.Fall[i]<- performance(dat.f, measure = "auc")@y.values[[1]]
        
        # RMSE
        mod.results$SeasonalSST.Independent.Fall.RMSE.P.Fall[i]<- performance(dat.f, measure = "rmse")@y.values[[1]]
        mod.results$SeasonalSST.Independent.Fall.RMSE.P.Spring[i]<- performance(dat.s, measure = "rmse")@y.values[[1]]
        
        # Calibration
        mod.results$SeasonalSST.Independent.Fall.Calib.SDM.Fall[i]<- round(calibration(x = test.data.f[,col.ind], p = pred.p.f)@statistic, 3)
        mod.results$SeasonalSST.Independent.Fall.Calib.SDM.Spring[i]<- round(calibration(x = test.data.s[,col.ind], p = pred.p.s)@statistic, 3)
        
        # Biomass component
        mod.results$SeasonalSST.Independent.Fall.DevExp.B[i]<- round(summary(mod.b)$dev.expl, 3)
        
        # RMSE
        test.data.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "FALL")
        test.data.f$SEASON<- factor(test.data.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(test.data.f))
        
        test.data.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        test.data.s$SEASON<- factor(test.data.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(test.data.s))
        
        mod.results$SeasonalSST.Independent.Fall.RMSE.B.Fall[i]<- RMSE(as.numeric(pred.b.f), test.data.f$BIOMASS, "NRMSE")
        mod.results$SeasonalSST.Independent.Fall.RMSE.B.Spring[i]<- RMSE(as.numeric(pred.b.s), test.data.s$BIOMASS, "NRMSE")
        
        # Prediction maps
        preds.f<- data.frame("x" = base.preds$Data[[1]]$x, "y" = base.preds$Data[[1]]$y)
        res.ind.f<- ncol(preds.f)
        preds.s<- data.frame("x" = base.preds$Data[[2]]$x, "y" = base.preds$Data[[2]]$y)
        res.ind.s<- ncol(preds.s)
        
        for(k in seq_along(pred.years)){
          fall.use<- c("x", "y", "DEPTH.Scale", pred.cols.f[grepl(pred.years[k], pred.cols.f)])
          spring.use<- c("x", "y", "DEPTH.Scale", pred.cols.s[grepl(pred.years[k], pred.cols.s)])
          
          pred.df.f<- data.frame(base.preds$Data[[1]])
          pred.df.f<- pred.df.f[,colnames(pred.df.f) %in% fall.use]
          pred.df.f$SEASON<- factor(rep("FALL", nrow(pred.df.f)), levels = c("FALL", "SPRING"))
          names(pred.df.f)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.df.s<- data.frame(base.preds$Data[[2]])
          pred.df.s<- pred.df.s[,colnames(pred.df.s) %in% spring.use]
          pred.df.s$SEASON<- factor(rep("SPRING", nrow(pred.df.s)), levels = c("FALL", "SPRING"))
          names(pred.df.s)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3))
          preds.f[,res.ind.f+1]<- pred.p.f * pred.logb.f
          colnames(preds.f)[res.ind.f+1]<- paste("Pred.", pred.years[k], sep = "")
          
          pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3))
          preds.s[,res.ind.s+1]<- pred.p.s * pred.logb.s
          colnames(preds.s)[res.ind.s+1]<- paste("Pred.", pred.years[k], sep = "")
          
          res.ind.f<- res.ind.f+1
          res.ind.s<- res.ind.s+1
        }
        
        SeasonalSST.Independent.Fall.pred.out.f<- data.frame("x" = preds.f$x, "y" = preds.f$y, "pred" = rowMeans(preds.f[,3:16]))
        SeasonalSST.Independent.Fall.pred.out.f$pred[!is.na(SeasonalSST.Independent.Fall.pred.out.f$pred) & SeasonalSST.Independent.Fall.pred.out.f$pred > 1000]<- NA
        SeasonalSST.Independent.Fall.pred.out.s<- data.frame("x" = preds.s$x, "y" = preds.s$y, "pred" = rowMeans(preds.s[,3:16]))
        SeasonalSST.Independent.Fall.pred.out.s$pred[!is.na(SeasonalSST.Independent.Fall.pred.out.s$pred) & SeasonalSST.Independent.Fall.pred.out.s$pred > 1000]<- NA
        
        # Center of Biomass
        cog.fall<- COGravity(x = SeasonalSST.Independent.Fall.pred.out.f$x, y = SeasonalSST.Independent.Fall.pred.out.f$y, wt = SeasonalSST.Independent.Fall.pred.out.f$pred)
        cog.spring<- COGravity(x = SeasonalSST.Independent.Fall.pred.out.s$x, y = SeasonalSST.Independent.Fall.pred.out.s$y, wt = SeasonalSST.Independent.Fall.pred.out.s$pred)
        mod.results$SeasonalSST.Independent.Fall.COG.Long.Fall[i]<- as.numeric(cog.fall[1])
        mod.results$SeasonalSST.Independent.Fall.COG.Lat.Fall[i]<- as.numeric(cog.fall[3])
        mod.results$SeasonalSST.Independent.Fall.COG.Long.Spring[i]<- as.numeric(cog.spring[1])
        mod.results$SeasonalSST.Independent.Fall.COG.Lat.Spring[i]<- as.numeric(cog.spring[3])
        
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
      }
    } 
    
    if(scenario.use == "SeasonalSST.Independent.Spring"){
      dat.train<- bind_rows(dat.temp$TRAIN.DATA[[2]])
      dat.train$SEASON<- factor(rep("SPRING", nrow(dat.train)), levels = c("FALL", "SPRING"))
      dat.test<- bind_rows(dat.temp$TEST.DATA[[1]], dat.temp$TEST.DATA[[2]])
      dat.test$SEASON<- factor(c(rep("FALL", nrow(dat.temp$TEST.DATA[[1]])), rep("SPRING", nrow(dat.temp$TEST.DATA[[2]]))), levels = c("FALL", "SPRING"))
      
      # Fit the models
      mod.p<- try(gam(PRESENCE.BIOMASS ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = binomial(link = logit), select = TRUE), silent = TRUE)
      mod.b<- try(gam(BIOMASS.MOD ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.train, family = gaussian, select = TRUE), silent = TRUE)
      
      if(any(class(mod.p) == "try-error") | any(class(mod.b) == "try-error")| is.infinite(summary(mod.p)$dev.expl)){
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
        next()
      } else {
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        pred.data.f<- data.frame(na.omit(temp.f))
        
        pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.f, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.f<- pred.p.f * pred.logb.f
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        pred.data.s<- data.frame(na.omit(temp.s))
        
        pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3)
        pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.data.s, type = "response", se.fit = TRUE)$fit), 3))
        pred.b.s<- pred.p.s * pred.logb.s
        
        # Evaluate and validate model: DevianceExplained, AUC, RMSE, Calibration
        # Presence/absence component
        mod.results$SeasonalSST.Independent.Spring.DevExp.P[i]<- round(summary(mod.p)$dev.expl, 3)
        
        temp.f<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "FALL")
        temp.f$SEASON<- factor(temp.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(temp.f))
        
        temp.s<- dplyr::select(dat.test, one_of(c("SEASON", "PRESENCE.BIOMASS", "DEPTH.Scale", "SEASONALMU.OISST.Scale"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        temp.s$SEASON<- factor(temp.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(temp.s))
        
        # AUC
        col.ind<- which(colnames(test.data.f) == "PRESENCE.BIOMASS")
        dat.f<- try(prediction(predictions = pred.p.f, labels = test.data.f[,col.ind]), silent = TRUE)
        dat.s<- try(prediction(predictions = pred.p.s, labels = test.data.s[,col.ind]), silent = TRUE)
        
        # Exit out if we can't at least get an AUC value
        if(class(dat.f) == "try-error" | class(dat.s) == "try-error") {
          next()
        } 
        
        mod.results$SeasonalSST.Independent.Spring.AUC.Fall[i]<- performance(dat.f, measure = "auc")@y.values[[1]]
        mod.results$SeasonalSST.Independent.Spring.AUC.Spring[i]<- performance(dat.s, measure = "auc")@y.values[[1]]
        
        # RMSE
        mod.results$SeasonalSST.Independent.Spring.RMSE.P.Fall[i]<- performance(dat.f, measure = "rmse")@y.values[[1]]
        mod.results$SeasonalSST.Independent.Spring.RMSE.P.Spring[i]<- performance(dat.s, measure = "rmse")@y.values[[1]]
        
        # Calibration
        mod.results$SeasonalSST.Independent.Spring.Calib.SDM.Fall[i]<- round(calibration(x = test.data.f[,col.ind], p = pred.p.f)@statistic, 3)
        mod.results$SeasonalSST.Independent.Spring.Calib.SDM.Spring[i]<- round(calibration(x = test.data.s[,col.ind], p = pred.p.s)@statistic, 3)
        
        # Biomass component
        mod.results$SeasonalSST.Independent.Spring.DevExp.B[i]<- round(summary(mod.b)$dev.expl, 3)
        
        # RMSE
        test.data.f<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "FALL")
        test.data.f$SEASON<- factor(test.data.f$SEASON, levels = c("FALL", "SPRING"))
        test.data.f<- data.frame(na.omit(test.data.f))
        
        test.data.s<- dplyr::select(dat.test, one_of(c("SEASON", "DEPTH.Scale", "SEASONALMU.OISST.Scale", "BIOMASS"))) %>%
          dplyr::filter(., SEASON == "SPRING")
        test.data.s$SEASON<- factor(test.data.s$SEASON, levels = c("FALL", "SPRING"))
        test.data.s<- data.frame(na.omit(test.data.s))
        
        mod.results$SeasonalSST.Independent.Spring.RMSE.B.Fall[i]<- RMSE(as.numeric(pred.b.f), test.data.f$BIOMASS, "NRMSE")
        mod.results$SeasonalSST.Independent.Spring.RMSE.B.Spring[i]<- RMSE(as.numeric(pred.b.s), test.data.s$BIOMASS, "NRMSE")
        
        # Prediction maps
        preds.f<- data.frame("x" = base.preds$Data[[1]]$x, "y" = base.preds$Data[[1]]$y)
        res.ind.f<- ncol(preds.f)
        preds.s<- data.frame("x" = base.preds$Data[[2]]$x, "y" = base.preds$Data[[2]]$y)
        res.ind.s<- ncol(preds.s)
        
        for(k in seq_along(pred.years)){
          fall.use<- c("x", "y", "DEPTH.Scale", pred.cols.f[grepl(pred.years[k], pred.cols.f)])
          spring.use<- c("x", "y", "DEPTH.Scale", pred.cols.s[grepl(pred.years[k], pred.cols.s)])
          
          pred.df.f<- data.frame(base.preds$Data[[1]])
          pred.df.f<- pred.df.f[,colnames(pred.df.f) %in% fall.use]
          pred.df.f$SEASON<- rep("FALL", nrow(pred.df.f))
          names(pred.df.f)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.df.s<- data.frame(base.preds$Data[[2]])
          pred.df.s<- pred.df.s[,colnames(pred.df.s) %in% spring.use]
          pred.df.s$SEASON<- rep("SPRING", nrow(pred.df.s))
          names(pred.df.s)[4]<- "SEASONALMU.OISST.Scale"
          
          pred.p.f<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.f<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.f, type = "response", se.fit = TRUE)$fit), 3))
          preds.f[,res.ind.f+1]<- pred.p.f * pred.logb.f
          colnames(preds.f)[res.ind.f+1]<- paste("Pred.", pred.years[k], sep = "")
          
          pred.p.s<- round(as.numeric(predict.gam(mod.p, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3)
          pred.logb.s<- exp(round(as.numeric(predict.gam(mod.b, newdata = pred.df.s, type = "response", se.fit = TRUE)$fit), 3))
          preds.s[,res.ind.s+1]<- pred.p.s * pred.logb.s
          colnames(preds.s)[res.ind.s+1]<- paste("Pred.", pred.years[k], sep = "")
          
          res.ind.f<- res.ind.f+1
          res.ind.s<- res.ind.s+1
        }
        
        SeasonalSST.Independent.Spring.pred.out.f<- data.frame("x" = preds.f$x, "y" = preds.f$y, "pred" = rowMeans(preds.f[,3:16]))
        SeasonalSST.Independent.Spring.pred.out.f$pred[!is.na(SeasonalSST.Independent.Spring.pred.out.f$pred) & SeasonalSST.Independent.Spring.pred.out.f$pred > 1000]<- NA
        SeasonalSST.Independent.Spring.pred.out.s<- data.frame("x" = preds.s$x, "y" = preds.s$y, "pred" = rowMeans(preds.s[,3:16]))
        SeasonalSST.Independent.Spring.pred.out.s$pred[!is.na(SeasonalSST.Independent.Spring.pred.out.s$pred) & SeasonalSST.Independent.Spring.pred.out.s$pred > 1000]<- NA
        
        # Center of Biomass
        cog.fall<- COGravity(x = SeasonalSST.Independent.Spring.pred.out.f$x, y = SeasonalSST.Independent.Spring.pred.out.f$y, wt = SeasonalSST.Independent.Spring.pred.out.f$pred)
        cog.spring<- COGravity(x = SeasonalSST.Independent.Spring.pred.out.s$x, y = SeasonalSST.Independent.Spring.pred.out.s$y, wt = SeasonalSST.Independent.Spring.pred.out.s$pred)
        mod.results$SeasonalSST.Independent.Spring.COG.Long.Fall[i]<- as.numeric(cog.fall[1])
        mod.results$SeasonalSST.Independent.Spring.COG.Lat.Fall[i]<- as.numeric(cog.fall[3])
        mod.results$SeasonalSST.Independent.Spring.COG.Long.Spring[i]<- as.numeric(cog.spring[1])
        mod.results$SeasonalSST.Independent.Spring.COG.Lat.Spring[i]<- as.numeric(cog.spring[3])
        
        print(paste(scenario.use, " for ", spp.use, " is done!", sep = ""))
      }
    } 
  }
  
  if(all(is.na(mod.results[i,]))){
    next
  }
  
  # Panel plot
  # Interpolation grid first
  coords.df<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y)
  pred.df<- na.omit(data.frame("x" = coords.df$x, "y" = coords.df$y, "layer" = rep(0, length(coords.df$x))))
  pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                          xo=seq(-87.99457, -57.4307, length = 115),
                          yo=seq(22.27352, 48.11657, length = 133))
  pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(pred.df.interp$z))
  pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
  
  # Get differences first...
  # Fall
  seas.seasfac.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.pred.out.f$pred - SeasonalSST.SeasonFactor.pred.out.f$pred)
  seas.seasint.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.pred.out.f$pred - SeasonalSST.SeasonInteraction.pred.out.f$pred)
  seas.seasindf.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.pred.out.f$pred - SeasonalSST.Independent.Fall.pred.out.f$pred)
  seas.seasinds.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.pred.out.f$pred - SeasonalSST.Independent.Spring.pred.out.f$pred)
  
  seasfac.seasint.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.SeasonFactor.pred.out.f$pred - SeasonalSST.SeasonInteraction.pred.out.f$pred)
  seasfac.seasindf.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.SeasonFactor.pred.out.f$pred - SeasonalSST.Independent.Fall.pred.out.f$pred)
  seasfac.seasinds.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.SeasonFactor.pred.out.f$pred - SeasonalSST.Independent.Spring.pred.out.f$pred)
  
  seasint.seasindf.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.SeasonInteraction.pred.out.f$pred - SeasonalSST.Independent.Fall.pred.out.f$pred)
  seasint.seasinds.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.SeasonInteraction.pred.out.f$pred - SeasonalSST.Independent.Spring.pred.out.f$pred)
  
  seasindf.seasinds.f<- data.frame("x" = SeasonalSST.pred.out.f$x, "y" = SeasonalSST.pred.out.f$y, "pred" = SeasonalSST.Independent.Fall.pred.out.f$pred - SeasonalSST.Independent.Spring.pred.out.f$pred)
  
  # Spring
  seas.seasfac.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.pred.out.s$pred - SeasonalSST.SeasonFactor.pred.out.s$pred)
  seas.seasint.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.pred.out.s$pred - SeasonalSST.SeasonInteraction.pred.out.s$pred)
  seas.seasindf.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.pred.out.s$pred - SeasonalSST.Independent.Fall.pred.out.s$pred)
  seas.seasinds.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.pred.out.s$pred - SeasonalSST.Independent.Spring.pred.out.s$pred)
  
  seasfac.seasint.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.SeasonFactor.pred.out.s$pred - SeasonalSST.SeasonInteraction.pred.out.s$pred)
  seasfac.seasindf.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.SeasonFactor.pred.out.s$pred - SeasonalSST.Independent.Fall.pred.out.s$pred)
  seasfac.seasinds.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.SeasonFactor.pred.out.s$pred - SeasonalSST.Independent.Spring.pred.out.s$pred)
  
  seasint.seasindf.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.SeasonInteraction.pred.out.s$pred - SeasonalSST.Independent.Fall.pred.out.s$pred)
  seasint.seasinds.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.SeasonInteraction.pred.out.s$pred - SeasonalSST.Independent.Spring.pred.out.s$pred)
  
  seasindf.seasinds.s<- data.frame("x" = SeasonalSST.pred.out.s$x, "y" = SeasonalSST.pred.out.s$y, "pred" = SeasonalSST.Independent.Fall.pred.out.s$pred - SeasonalSST.Independent.Spring.pred.out.s$pred)
  
  # Get range limits
  range.norm.f<- range(c(SeasonalSST.pred.out.f$pred, SeasonalSST.SeasonFactor.pred.out.f$pred, SeasonalSST.SeasonInteraction.pred.out.f$pred, SeasonalSST.Independent.Fall.pred.out.f$pred, SeasonalSST.Independent.Spring.pred.out.f$pred), na.rm = T)
  range.norm.f<- c(floor(range.norm.f[1]), ceiling(range.norm.f[2]))
  range.norm.s<- range(c(SeasonalSST.pred.out.s$pred, SeasonalSST.SeasonFactor.pred.out.s$pred, SeasonalSST.SeasonInteraction.pred.out.s$pred, SeasonalSST.Independent.Fall.pred.out.s$pred, SeasonalSST.Independent.Spring.pred.out.s$pred), na.rm = T)
  range.norm.s<- c(floor(range.norm.s[1]), ceiling(range.norm.s)[2])
  
  range.diff.f<- range(c(seas.seasfac.f$pred, seas.seasint.f$pred, seas.seasindf.f$pred, seas.seasinds.f$pred, seasfac.seasint.f$pred, seasfac.seasindf.f$pred, seasfac.seasinds.f$pred, seasint.seasindf.f$pred, seasint.seasinds.f$pred, seasindf.seasinds.f$pred), na.rm = T)
  range.diff.f<- c(floor(range.diff.f[1]), ceiling(range.diff.f[2]))
  range.diff.s<- ceiling(range(c(seas.seasfac.s$pred, seas.seasint.s$pred, seas.seasindf.s$pred, seas.seasinds.s$pred, seasfac.seasint.s$pred, seasfac.seasindf.s$pred, seasfac.seasinds.s$pred, seasint.seasindf.s$pred, seasint.seasinds.s$pred, seasindf.seasinds.s$pred), na.rm = T))
  range.diff.s<- c(floor(range.diff.s[1]), ceiling(range.diff.s[2]))
  
  # Make into list by season
  fall.list<- list(SeasonalSST.pred.out.f, SeasonalSST.SeasonFactor.pred.out.f, SeasonalSST.SeasonInteraction.pred.out.f, SeasonalSST.Independent.Fall.pred.out.f, SeasonalSST.Independent.Spring.pred.out.f, seas.seasfac.f, seas.seasint.f, seas.seasindf.f, seas.seasinds.f, seasfac.seasint.f, seasfac.seasindf.f,  seasfac.seasinds.f, seasint.seasindf.f,  seasint.seasinds.f, seasindf.seasinds.f)
  names(fall.list)<- c("SeasonalSST", "SeasonalSST.SeasonFactor", "SeasonalSST.SeasonInteraction", "SeasonalSST.Independent.Fall", "SeasonalSST.Independent.Spring", rep("", 10))
  fall.res<- vector("list", length = length(fall.list))
  
  for(l in seq_along(fall.list)){
    data.use<- fall.list[[l]]
    if(l<=5){
      range.use<- range.norm.f
      auc.name<- paste(names(fall.list)[l], ".AUC.Fall", sep = "")
      auc.use<- round(mod.results[i, match(auc.name, colnames(mod.results))], 3)
      rmse.p.name<- paste(names(fall.list)[l], ".RMSE.P.Fall", sep = "")
      rmse.p.use<- round(mod.results[i, match(rmse.p.name, colnames(mod.results))], 3)
      rmse.b.name<- paste(names(fall.list)[l], ".RMSE.B.Fall", sep = "")
      rmse.b.use<- round(mod.results[i, match(rmse.b.name, colnames(mod.results))], 3)
      label.use<- paste("AUC: ", auc.use, "\nRMSE.P: ", rmse.p.use, "\nRMSE.B: ", rmse.b.use, sep = "")
    } else {
      range.use<- range.diff.f
    }
    
    pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    pred.df.use$z[pred.df.use$z > 1000]<- NA
    pred.df.use$z[pred.df.use$z < -1000]<- NA
    
    if(l<= 5){
      fall.res[[l]]<- ggplot() +
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_viridis(option = "viridis", name = "", na.value = "white", limits = range.use, direction = 1) +
        annotate("text", label = label.use, x = -70, y = 37.5) +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("") +
        scale_x_continuous("", breaks = c(-75.0, -65.0), labels = c("-75.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) +
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5, 0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(size = 8))
    } else {
      fall.res[[l]]<- ggplot() +
      geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
      scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, name = "", na.value = "white", limits = range.use) +
      geom_map(data = us.states.f, map = us.states.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      geom_map(data = ca.provinces.f, map = ca.provinces.f,
               aes(map_id = id, group = group),
               fill = "gray65", color = "gray45", size = 0.15) +
      ylim(ylim.use) + ylab("") +
      scale_x_continuous("", breaks = c(-75.0, -65.0), labels = c("-75.0", "-65.0"), limits = xlim.use) +
      coord_fixed(1.3) +
      theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5, 0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(size = 8))
    }
  }

  # Make into list by season
  spring.list<- list(SeasonalSST.pred.out.s, SeasonalSST.SeasonFactor.pred.out.s, SeasonalSST.SeasonInteraction.pred.out.s, SeasonalSST.Independent.Fall.pred.out.s, SeasonalSST.Independent.Spring.pred.out.s, seas.seasfac.s, seas.seasint.s, seas.seasindf.s, seas.seasinds.s, seasfac.seasint.s, seasfac.seasindf.s,  seasfac.seasinds.s, seasint.seasindf.s,  seasint.seasinds.s, seasindf.seasinds.s)
  names(spring.list)<- c("SeasonalSST", "SeasonalSST.SeasonFactor", "SeasonalSST.SeasonInteraction", "SeasonalSST.Independent.Fall", "SeasonalSST.Independent.Spring", rep("", 10))
  spring.res<- vector("list", length = length(spring.list))
  
  for(l in seq_along(spring.list)){
    data.use<- spring.list[[l]]
    if(l<=5){
      range.use<- range.norm.s
      auc.name<- paste(names(spring.list)[l], ".AUC.Spring", sep = "")
      auc.use<- round(mod.results[i, match(auc.name, colnames(mod.results))], 3)
      rmse.p.name<- paste(names(spring.list)[l], ".RMSE.P.Spring", sep = "")
      rmse.p.use<- round(mod.results[i, match(rmse.p.name, colnames(mod.results))], 3)
      rmse.b.name<- paste(names(spring.list)[l], ".RMSE.B.Spring", sep = "")
      rmse.b.use<- round(mod.results[i, match(rmse.b.name, colnames(mod.results))], 3)
      label.use<- paste("AUC: ", auc.use, "\nRMSE.P: ", rmse.p.use, "\nRMSE.B: ", rmse.b.use, sep = "")
    } else {
      range.use<- range.diff.s
    }
    
    pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
    pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    if(l<= 5){
      spring.res[[l]]<- ggplot() +
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_viridis(option = "viridis", name = "Biomass\n", na.value = "white", limits = range.use, direction = 1) +
        annotate("text", label = label.use, x = -70, y = 37.5) +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("") +
        scale_x_continuous("", breaks = c(-75.0, -65.0), labels = c("-75.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) +
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5, 0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(size = 8))
    } else {
      spring.res[[l]]<- ggplot() +
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, name = "Biomass\n difference", na.value = "white", limits = range.use) +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("") +
        scale_x_continuous("", breaks = c(-75.0, -65.0), labels = c("-75.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) +
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5, 0.25), legend.text=element_text(size=10), legend.title=element_text(size=10), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(size = 8))
    }
  }
  
  # Combine all plots into one layout for each season
  margins_left<- c(-0.4, 0, 0, 0)
  margins_right<- c(0.4, 0, -0.2, 0)
  out.f<- plot_grid(fall.res[[1]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, NULL, NULL, NULL, fall.res[[6]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), fall.res[[2]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, NULL, NULL, fall.res[[7]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")),  fall.res[[10]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), fall.res[[3]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, NULL, fall.res[[8]]+ theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), fall.res[[11]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), fall.res[[13]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), fall.res[[4]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, fall.res[[9]] + theme(legend.position="none", plot.margin = unit(margins_right, "cm")), fall.res[[12]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), fall.res[[14]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), fall.res[[15]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), fall.res[[5]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), nrow = 5, labels = c("SeasonalSST", "Factor", "Interaction", "Fall", "Spring", rep("", 20)), hjust = -0.75, vjust = 0.75, scale = c(rep(0.9, 6), rep(0.85, 4), 0.9, rep(0.85, 4), 0.9, rep(0.85, 4), rep(1.15, 5)))
  legend.norm<- get_legend(fall.res[[1]] + theme(legend.position = c(0.0, 0.13), plot.margin = unit(c(0,0,0,0), "cm")))
  legend.diff<- get_legend(fall.res[[6]] + theme(legend.position = c(0.05, 0.13), plot.margin = unit(c(0,0,0,0), "cm")))
  legends<- plot_grid(legend.norm, legend.diff, nrow = 1, ncol = 2, align = "h")
  out.f2<- plot_grid(out.f, legends, nrow = 1, rel_widths = c(7, 0.9), rel_heights = c(7, 0.9), ncol = 2)
  ggsave(filename = paste("~/GitHub/SDMSeasonality/Results/Maps", spp.use, "fall.jpg", sep = ""), out.f2, width = 14, height = 10)
  
  out.s<- plot_grid(spring.res[[1]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, NULL, NULL, NULL, spring.res[[6]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), spring.res[[2]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, NULL, NULL, spring.res[[7]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")),  spring.res[[10]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), spring.res[[3]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, NULL, spring.res[[8]]+ theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(margins_left, "cm")), spring.res[[11]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), spring.res[[13]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), spring.res[[4]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_left, "cm")), NULL, spring.res[[9]] + theme(legend.position="none", plot.margin = unit(margins_right, "cm")), spring.res[[12]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), spring.res[[14]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), spring.res[[15]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), spring.res[[5]] + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(margins_right, "cm")), nrow = 5, labels = c("SeasonalSST", "Factor", "Interaction", "Fall", "Spring", rep("", 20)), hjust = -0.75, vjust = 0.75, scale = c(rep(0.9, 6), rep(0.85, 4), 0.9, rep(0.85, 4), 0.9, rep(0.85, 4), rep(1.15, 5)))
  legend.norm<- get_legend(spring.res[[1]] + theme(legend.position = c(0.0, 0.13), plot.margin = unit(c(0,0,0,0), "cm")))
  legend.diff<- get_legend(spring.res[[6]] + theme(legend.position = c(0.05, 0.13), plot.margin = unit(c(0,0,0,0), "cm")))
  legends<- plot_grid(legend.norm, legend.diff, nrow = 1, ncol = 2, align = "h")
  out.s2<- plot_grid(out.s, legends, nrow = 1, rel_widths = c(7, 0.9), rel_heights = c(7, 0.9), ncol = 2)
  ggsave(filename = paste("~/GitHub/SDMSeasonality/Results/Maps/", spp.use, "spring.jpg", sep = ""), out.s2, width = 14, height = 10)
  
}

# Save results
write.csv(mod.results, "~/GitHub/SDMSeasonality/Results/Tables/mod.results.table.csv")

#####
## Bhattacharyya's Index
#####
## Species data -- filtered
dat.path<- "~/GitHub/SDMSeasonality/Data/model.dat.rds"

# Read it in, filter to species with at least 250 observations and observed in at least 10 years, do some quick formatting to fit the GAM
dat<- readRDS(dat.path)

dat.filter.a<- dat %>%
  group_by(., SVSPP) %>%
  summarize_at(vars(PRESENCE), sum, na.rm = T) %>%
  filter(., PRESENCE >= 250)

dat.filter.b<- dat %>%
  filter(., SVSPP %in% dat.filter.a$SVSPP) %>%
  group_by(., SEASON, SVSPP) %>% 
  filter(., PRESENCE > 0) %>%
  summarize_at(vars(EST_YEAR), n_distinct, na.rm = T) %>%
  filter(., EST_YEAR >= 10) 
dat.filter.b$Count<- ifelse(dat.filter.b$EST_YEAR >= 10, 1, 0)
dat.filter.b<- dat.filter.b %>%
  group_by(., SVSPP) %>%
  summarize_at(vars(Count), sum, na.rm = T) %>%
  filter(., Count == 2)

dat<- dat %>%
  filter(., SVSPP %in% dat.filter.b$SVSPP)

# Add species common name
species.names<- read_csv("~/GitHub/SDMSeasonality/Data/svspp_spp names.csv") %>%
  dplyr::select(., -X4) 
species.names$SVSPP<- as.numeric(species.names$SVSPP)
dat<- dat %>%
  left_join(., species.names) 
dat<- dat[-which(grepl("UNCL", dat$COMNAME)),]
dat<- dat[-which(grepl("UNKNOWN 01", dat$COMNAME)),]
dat<- dat[-which(is.na(dat$COMNAME)),]

## KernelUD by season for each species 
# Nest the data
dat<- dat %>%
  group_by(., SEASON, COMNAME) %>%
  nest(., .key = "Data") %>%
  arrange(., COMNAME, SEASON)

# Home range estimates from KernelUD
# Need a spatial grid 
nelme<- st_read("~/GitHub/SDMSeasonality/Data/NELME_clipped.shp")
st_crs(nelme)<- "+init=epsg:4326"
grid.use<- st_make_grid(nelme, cellsize = 0.25, what = "centers")
grid.use<- as(st_sf(grid.use), "Spatial")
grid.use<- SpatialPixels(grid.use)

# Could we make a border?
kernel.border<- st_read("~/GitHub/SDMSeasonality/Data/NELME_KernelUD_Border.shp")
st_crs(kernel.border)<- "+init=epsg:4326"
kernel.border<- as(st_sf(kernel.border), "Spatial")

kernelUD_function<- function(df){
  if(FALSE){
    df<- dat$Data[[1]]
    h.use<- "LSCV"
    grid.use<- grid.use
    boundary.use<- kernel.border
  }
  
  # Spatial points
  dat.sp<- df %>%
    filter(., PRESENCE > 0) %>%
    dplyr::select(., DECDEG_BEGLON, DECDEG_BEGLAT) 
  coordinates(dat.sp)<- ~DECDEG_BEGLON+DECDEG_BEGLAT
  proj4string(dat.sp)<-  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  kern.temp<- kernelUD(xy = dat.sp, h = "LSCV", grid = grid.use)
  return(kern.temp)
}

# Run it
dat<- dat %>%
  mutate(., "KernelUD" = map(Data, possibly(kernelUD_function, NA)))

# Check em
if(FALSE){
  for (i in 1:nrow(dat)){
    plot(dat$KernelUD[[i]], main = paste(tolower(dat$SEASON[[i]]), "_", dat$COMNAME[[i]], sep = ""))
    grDevices::devAskNewPage(ask = TRUE)
  }
  grDevices::devAskNewPage(ask = FALSE)
}

# Calculating BA from seasonal home range estimates
dat.ba<- dat %>%
  dplyr::select(., SEASON, COMNAME, KernelUD) %>%
  spread(., SEASON, KernelUD)

ba_index_function<- function(fall.ud, spring.ud){
  if(FALSE){
    fall.ud<- dat.ba$FALL[[1]]
    spring.ud<- dat.ba$SPRING[[1]]
  }
  
  all.ud<- list(fall.ud, spring.ud)
  names(all.ud)<- c("Fall", "Spring")
  
  # Rasters 
  ud.rast<- stack(lapply(all.ud, raster))
  
  # Convert to estUDm class and then calculate kerneloverlap
  all.estUD<-lapply(lapply(unstack(ud.rast),as,"SpatialPixelsDataFrame"), new, Class='estUD', vol=F)
  class(all.estUD)<- 'estUDm'
  ba.index<- kerneloverlaphr(all.estUD, meth="BA", percent = 95, conditional=TRUE)
  return(mean(diag(apply(ba.index,2,rev))))
}

dat.ba<- dat.ba %>%
  mutate(., "BA.Index" = map2(FALL, SPRING, possibly(ba_index_function, NA)))
hist(unlist(dat.ba$BA.Index))
summary(unlist(dat.ba$BA.Index))

#####
## Results data summary
#####
# Read in results
res<- read_csv("~/GitHub/SDMSeasonality/Results/Tables/mod.results.table.csv")
res<- res[complete.cases(res),] %>%
  dplyr::select(., -X1)

# How many times is each model the best vs. the worst? For each species and season, assign rank?
# Don't need deviance explained...
devexp.index<- colnames(res)[which(grepl("DevExp", colnames(res)))]
res.table<- res %>%
  dplyr::select(., -one_of(devexp.index))

res.table.l<- res.table %>%
  gather(., "Model.Statistic", "Value", -Species)
res.table.l$Season<- str_extract(res.table.l$Model.Statistic, "[^.]*$")
res.table.l$Season<- factor(res.table.l$Season, levels = c("Fall", "Spring"))
res.table.l$Model<- ifelse(grepl("SeasonFactor", res.table.l$Model.Statistic), "SeasonFactor",
                        ifelse(grepl("SeasonInteraction", res.table.l$Model.Statistic), "SeasonInteraction",
                               ifelse(grepl("Independent.Fall", res.table.l$Model.Statistic), "SeasonIndependent.Fall",
                                      ifelse(grepl("Independent.Spring", res.table.l$Model.Statistic), "SeasonIndependent.Spring", "SeasonSST"))))
res.table.l$Model<- factor(res.table.l$Model, levels = c("SeasonSST", "SeasonFactor", "SeasonInteraction", "SeasonIndependent.Fall", "SeasonIndependent.Spring"))
res.table.l$Statistic<- ifelse(grepl("AUC", res.table.l$Model.Statistic), "AUC",
                               ifelse(grepl("Calib", res.table.l$Model.Statistic), "Calib",
                                      ifelse(grepl("RMSE.B", res.table.l$Model.Statistic), "RMSE.B", "Drop")))

# Only going to rank AUC and RMSE.B
res.table.l<- res.table.l %>%
  filter(., Statistic != "Drop")

#
if(FALSE){
  diffs<- res.table.l %>%
    filter(., Season == "Spring") %>%
    filter(., Model == "SeasonIndependent.Spring" | Model == "SeasonIndependent.Fall") %>%
    group_by(., Species, Statistic)
  
  diffs<- diffs %>%
    ungroup() %>%
    dplyr::select(., -Model, -Season, -Statistic) %>%
    spread(., Model.Statistic, Value)
  
  auc.diff<- diffs$SeasonalSST.Independent.Spring.AUC.Spring - diffs$SeasonalSST.Independent.Fall.AUC.Spring
  rmse.diff<- diffs$SeasonalSST.Independent.Spring.RMSE.B.Spring - diffs$SeasonalSST.Independent.Fall.RMSE.B.Spring
  diffs<- data.frame("Species" = diffs$Species, "AUC" = auc.diff, "RMSE" = rmse.diff)
  inds.better<- diffs[which(diffs$AUC > 0 & diffs$RMSE > 0),]
  indf.better<- diffs[which(diffs$AUC < 0 & diffs$RMSE < 0),]
  t<- rmse.diff[which(rmse.diff < 0)]
}

# Now rank
res.table.ranks.auc<- res.table.l %>%
  filter(., Statistic == "AUC") %>%
  group_by(., Species, Season) %>%
  mutate(., "Rank.AUC" = order(order(round(Value, 2), decreasing = T))) %>%
  arrange(., Species, Season, Rank.AUC) %>%
  dplyr::select(., -Value, -Model.Statistic, -Statistic)

res.table.ranks.calib<- res.table.l %>%
  filter(., Statistic == "Calib") %>%
  group_by(., Species, Season) %>%
  mutate(., "Rank.Calib" = order(order(round(Value, 2), decreasing = T))) %>%
  arrange(., Species, Season, Rank.Calib) %>%
  dplyr::select(., -Value, -Model.Statistic, -Statistic)

if(FALSE){
  temp<- res.table.ranks.calib[res.table.ranks.calib$Season == "Spring" & res.table.ranks.calib$Model == "SeasonIndependent.Spring" & res.table.ranks.calib$Rank.Calib == 5,]
}

res.table.ranks.rmse<- res.table.l %>%
  filter(., Statistic == "RMSE.B") %>%
  group_by(., Species, Season) %>%
  mutate(., "Rank.RMSE" = rank(round(Value, 2))) %>%
  arrange(., Species, Season, Statistic, Rank.RMSE) %>%
  dplyr::select(., -Value, -Model.Statistic, -Statistic)

# Lets join these...
res.table.both<- res.table.ranks.auc %>%
  left_join(., res.table.ranks.rmse) %>%
  left_join(., res.table.ranks.calib) %>%
  mutate(., "Sum_Ranks" = as.numeric(Rank.AUC) + as.numeric(Rank.Calib) + as.numeric(Rank.RMSE)) %>%
  dplyr::select(., Species, Season, Model, Sum_Ranks) %>%
  arrange(., Species, Season, Model, Sum_Ranks)

# Rank the sum of ranks?
res.table.summ<- res.table.both %>%
  mutate(., "Sum_Ranks_Rank" = rank(Sum_Ranks)) %>%
  mutate(., "Sum_Ranks_Rounded" = ifelse(Sum_Ranks_Rank == 1.50, 1, ifelse(Sum_Ranks_Rank == 4.50, 5, Sum_Ranks_Rank))) %>%
  group_by(., Season, Model, Sum_Ranks_Rounded) %>%
  summarize_at(., "Species", n_distinct) 

# Filter (only care about number of times it is the best or the worst)...
res.table.summ<- res.table.summ %>%
  filter(., Sum_Ranks_Rounded == 1.00 | Sum_Ranks_Rounded == 5.00) %>%
  arrange(., Season, Sum_Ranks_Rounded, Species) %>%
  mutate(., "Model.Rank" = paste(Model, Sum_Ranks_Rounded, sep = "_"))
res.table.summ$Model.Rank<- factor(res.table.summ$Model.Rank, levels = c("SeasonSST_1", "SeasonSST_5", "SeasonFactor_1", "SeasonFactor_5", "SeasonInteraction_1", "SeasonInteraction_5", "SeasonIndependent.Fall_1", "SeasonIndependent.Fall_5", "SeasonIndependent.Spring_1", "SeasonIndependent.Spring_5"))
res.table.summ<- res.table.summ %>%
  ungroup() %>%
  dplyr::select(., Season, Model.Rank, Species) %>%
  spread(., Model.Rank, Species)

# Easier to plot?
res.table.plot<- res.table.summ %>%
  gather(., Model.Rank, Species, -Season) %>%
  separate(., Model.Rank, c("Model", "Rank"), "_")

res.table.plot$Model.Form<- ifelse(grepl("SeasonFactor", res.table.plot$Model), "Both.SeasonFactor",
                             ifelse(grepl("SeasonInteraction", res.table.plot$Model), "Both.SeasonInteraction",
                                    ifelse(grepl("Independent.Fall", res.table.plot$Model), "FallOnly",
                                           ifelse(grepl("Independent.Spring", res.table.plot$Model), "SpringOnly", "Both.SeasonIgnored"))))
res.table.plot$Model.Form<- factor(res.table.plot$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))

res.rank.plot<- ggplot() +
  geom_col(data = res.table.plot, aes(x = Model.Form, y = Species, fill = Rank), position = "dodge") +
  scale_fill_manual(name = "Model Rank", values = c('#4daf4a', '#377eb8'), labels = c("Best", "Worst")) +
  facet_wrap(~Season) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())


# Panel for each stat?
res.table.both<- res.table.ranks.auc %>%
  left_join(., res.table.ranks.rmse) %>%
  left_join(., res.table.ranks.calib) %>%
  gather(., "Stat", "Value", -Species, -Season, -Model) %>%
  group_by(., Season, Model, Stat, Value) %>%
  summarize_at(., "Species", n_distinct) %>%
  dplyr::filter(., Value == 1 | Value == 5)
res.table.both$Value<- factor(res.table.both$Value, levels = c("1", "5"))

ggplot() +
  geom_col(data = res.table.both, aes(x = Model, y = Species, fill = Value), position = "dodge") +
  scale_fill_manual(name = "Model Rank", values = c('#4daf4a', '#377eb8'), labels = c("Best", "Worst")) +
  facet_wrap(~Season+Stat) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())


# Merge with BA index data
dat.ba.merge<- dat.ba %>%
  dplyr::select(., COMNAME, BA.Index) %>%
  as.data.frame()
colnames(dat.ba.merge)[1]<- "Species"
dat.ba.merge$BA.Index<- as.numeric(unlist(dat.ba.merge$BA.Index))

res.full<- res %>%
  left_join(., dat.ba.merge, by = "Species") %>%
  drop_na(., BA.Index)

# Code for relative migration
splits<- quantile(res.full$BA.Index, probs = c(0, 0.33, 0.67, 1))
res.full$MigrationIndexCode<- cut(res.full$BA.Index, splits, include.lowest = T)
levels(res.full$MigrationIndexCode)<- c("SM", "SS", "YRR")

#### Box and whisker plots: AUC, calibration RMSE
{
# PA AUC
yind<- "AUC"
yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
keep<- c("Species", "BA.Index", "MigrationIndexCode", yind.keep)
plot.dat<- res.full %>%
  dplyr::select(., keep) %>%
  gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode)
plot.dat$Season<- str_extract(plot.dat$Scenario, "[^.]*$")
plot.dat$Season<- factor(plot.dat$Season, levels = c("Fall", "Spring"))

plot.dat$Model.Form<- ifelse(grepl("SeasonFactor", plot.dat$Scenario), "Both.SeasonFactor",
                             ifelse(grepl("SeasonInteraction", plot.dat$Scenario), "Both.SeasonInteraction",
                                    ifelse(grepl("Independent.Fall", plot.dat$Scenario), "FallOnly",
                                           ifelse(grepl("Independent.Spring", plot.dat$Scenario), "SpringOnly", "Both.SeasonIgnored"))))
plot.dat$Model.Form<- factor(plot.dat$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))

auc.boxplot<- ggplot(plot.dat, aes(x = Model.Form, y = Value, fill = Model.Form)) + 
  geom_boxplot() +
  scale_fill_manual(name = "Model Form", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
  scale_x_discrete(levels(plot.dat$Model.Form), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
  scale_y_continuous(limits = c(0.5, 1)) +
  ylab("AUC") +
  facet_wrap(~Season) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())
ggsave("~/GitHub/SDMSeasonality/Results/Figures/AUCBoxPlot.jpg", plot = auc.boxplot, width = 11, height = 8)

# PA Calibration
yind<- "Calib"
yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
keep<- c("Species", "BA.Index", "MigrationIndexCode", yind.keep)
plot.dat<- res.full %>%
  dplyr::select(., keep) %>%
  gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode)
plot.dat$Season<- str_extract(plot.dat$Scenario, "[^.]*$")
plot.dat$Season<- factor(plot.dat$Season, levels = c("Fall", "Spring"))

plot.dat$Model.Form<- ifelse(grepl("SeasonFactor", plot.dat$Scenario), "Both.SeasonFactor",
                             ifelse(grepl("SeasonInteraction", plot.dat$Scenario), "Both.SeasonInteraction",
                                    ifelse(grepl("Independent.Fall", plot.dat$Scenario), "FallOnly",
                                           ifelse(grepl("Independent.Spring", plot.dat$Scenario), "SpringOnly", "Both.SeasonIgnored"))))
plot.dat$Model.Form<- factor(plot.dat$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))

calib.boxplot<- ggplot(plot.dat, aes(x = Model.Form, y = Value, fill = Model.Form)) + 
  geom_boxplot() +
  scale_fill_manual(name = "Model Form", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
  scale_x_discrete(levels(plot.dat$Model.Form), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
  ylab("Calibration statistic") +
  facet_wrap(~Season) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())
ggsave("~/GitHub/SDMSeasonality/Results/Figures/CalibBoxPlot.jpg", plot = calib.boxplot, width = 11, height = 8)

# Biomass
# RMSE
yind<- "RMSE.B"
yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
keep<- c("Species", "BA.Index", "MigrationIndexCode", yind.keep)
plot.dat<- res.full %>%
  dplyr::select(., keep) %>%
  gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode)
plot.dat$Log.Value<- log(plot.dat$Value)
plot.dat$Season<- str_extract(plot.dat$Scenario, "[^.]*$")
plot.dat$Season<- factor(plot.dat$Season, levels = c("Fall", "Spring"))

plot.dat$Model.Form<- ifelse(grepl("SeasonFactor", plot.dat$Scenario), "Both.SeasonFactor",
                             ifelse(grepl("SeasonInteraction", plot.dat$Scenario), "Both.SeasonInteraction",
                                    ifelse(grepl("Independent.Fall", plot.dat$Scenario), "FallOnly",
                                           ifelse(grepl("Independent.Spring", plot.dat$Scenario), "SpringOnly", "Both.SeasonIgnored"))))
plot.dat$Model.Form<- factor(plot.dat$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))

rmse.b.boxplot<- ggplot(plot.dat, aes(x = Model.Form, y = Log.Value, fill = Model.Form)) + 
  geom_boxplot() +
  scale_fill_manual(name = "Model Form", values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
  scale_x_discrete(levels(plot.dat$Model.Form), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
  ylab("Log(Biomass root mean squared error)") +
  facet_wrap(~Season) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())
ggsave("~/GitHub/SDMSeasonality/Results/Figures/BiomassRMSEBoxPlot.jpg", plot = rmse.b.boxplot, width = 11, height = 8)
}

#### Box and whisker plots: AUC, Calibration, RMSE vs. migration index 
{
  # PA AUC
  yind<- "AUC"
  yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
  keep<- c("Species", "BA.Index", "MigrationIndexCode", yind.keep)
  plot.dat<- res.full %>%
    dplyr::select(., keep) %>%
    gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode)
  plot.dat$Season<- str_extract(plot.dat$Scenario, "[^.]*$")
  plot.dat$Season<- factor(plot.dat$Season, levels = c("Fall", "Spring"))
  
  plot.dat$Model.Form<- ifelse(grepl("SeasonFactor", plot.dat$Scenario), "Both.SeasonFactor",
                               ifelse(grepl("SeasonInteraction", plot.dat$Scenario), "Both.SeasonInteraction",
                                      ifelse(grepl("Independent.Fall", plot.dat$Scenario), "FallOnly",
                                             ifelse(grepl("Independent.Spring", plot.dat$Scenario), "SpringOnly", "Both.SeasonIgnored"))))
  plot.dat$Model.Form<- factor(plot.dat$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))
  
  auc.migind.boxplot<- ggplot(plot.dat, aes(x = Model.Form, y = Value, color = MigrationIndexCode)) + 
    geom_boxplot(position = "dodge") +
    scale_color_manual(name = "Migration Index", values = c("firebrick3", "darkolivegreen4", "dodgerblue4"), labels = c("Seasonally Migrating", "Seasonally Shifting", "Year-round Resident")) +
    scale_x_discrete(levels(plot.dat$Model.Form), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
    scale_y_continuous(limits = c(0.5, 1)) +
    ylab("AUC") +
    facet_wrap(~Season) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank())
  ggsave("~/GitHub/SDMSeasonality/Results/Figures/AUCMigIndBoxPlot.jpg", plot = auc.migind.boxplot, width = 11, height = 8)  
  # PA Calibration
  yind<- "Calib"
  yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
  keep<- c("Species", "BA.Index", "MigrationIndexCode", yind.keep)
  plot.dat<- res.full %>%
    dplyr::select(., keep) %>%
    gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode)
  plot.dat$Season<- str_extract(plot.dat$Scenario, "[^.]*$")
  plot.dat$Season<- factor(plot.dat$Season, levels = c("Fall", "Spring"))
  
  plot.dat$Model.Form<- ifelse(grepl("SeasonFactor", plot.dat$Scenario), "Both.SeasonFactor",
                               ifelse(grepl("SeasonInteraction", plot.dat$Scenario), "Both.SeasonInteraction",
                                      ifelse(grepl("Independent.Fall", plot.dat$Scenario), "FallOnly",
                                             ifelse(grepl("Independent.Spring", plot.dat$Scenario), "SpringOnly", "Both.SeasonIgnored"))))
  plot.dat$Model.Form<- factor(plot.dat$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))
  
  calib.migind.boxplot<- ggplot(plot.dat, aes(x = Model.Form, y = Value, color = MigrationIndexCode)) + 
    geom_boxplot(position = "dodge") +
    scale_color_manual(name = "Migration Index", values = c("firebrick3", "darkolivegreen4", "dodgerblue4"), labels = c("Seasonally Migrating", "Seasonally Shifting", "Year-round Resident")) +
    scale_x_discrete(levels(plot.dat$Model.Form), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
    ylab("Calibration statistic") +
    facet_wrap(~Season) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank())
  ggsave("~/GitHub/SDMSeasonality/Results/Figures/CalibMigIndBoxPlot.jpg", plot = calib.migind.boxplot, width = 11, height = 8)
  
  # RMSE.B
  yind<- "RMSE.B"
  yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
  keep<- c("Species", "BA.Index", "MigrationIndexCode", yind.keep)
  plot.dat<- res.full %>%
    dplyr::select(., keep) %>%
    gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode)
  plot.dat$Season<- str_extract(plot.dat$Scenario, "[^.]*$")
  plot.dat$Season<- factor(plot.dat$Season, levels = c("Fall", "Spring"))
  plot.dat$Log.Value<- log(plot.dat$Value)
  
  plot.dat$Model.Form<- ifelse(grepl("SeasonFactor", plot.dat$Scenario), "Both.SeasonFactor",
                               ifelse(grepl("SeasonInteraction", plot.dat$Scenario), "Both.SeasonInteraction",
                                      ifelse(grepl("Independent.Fall", plot.dat$Scenario), "FallOnly",
                                             ifelse(grepl("Independent.Spring", plot.dat$Scenario), "SpringOnly", "Both.SeasonIgnored"))))
  plot.dat$Model.Form<- factor(plot.dat$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))
  
  rmse.b.migind.boxplot<- ggplot(plot.dat, aes(x = Model.Form, y = Log.Value, color = MigrationIndexCode)) + 
    geom_boxplot(position = "dodge") +
    scale_color_manual(name = "Migration Index", values = c("firebrick3", "darkolivegreen4", "dodgerblue4"), labels = c("Seasonally Migrating", "Seasonally Shifting", "Year-round Resident")) +
    scale_x_discrete(levels(plot.dat$Model.Form), labels = c("Fall and Spring\n Season Ignored", "Fall and Spring\n Season Factor", "Fall and Spring\n Season Interaction", "Fall Only", "Spring Only")) +
    ylab("Log(Root mean squared error)") +
    facet_wrap(~Season) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank())
  ggsave("~/GitHub/SDMSeasonality/Results/Figures/RMSEBMigIndBoxPlot.jpg", plot = rmse.b.migind.boxplot, width = 11, height = 8)
}

#### Center of gravity
yind<- "COG"
yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
keep<- c("Species", "BA.Index", "MigrationIndexCode", yind.keep)

res.cog<- res.full %>%
  dplyr::select(., keep)

# Now, we want to calculate the distance based of the baseline approach (ignoring season)...
res.cog.db<- data.frame("Species" = unique(res.cog$Species), "SeasonalSST.SeasonFactor.Dist.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.SeasonInteraction.Dist.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.Independent.Fall.Dist.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.Independent.Spring.Dist.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.SeasonFactor.Bear.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.SeasonInteraction.Bear.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.Independent.Fall.Bear.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.Independent.Spring.Bear.Fall" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.SeasonFactor.Dist.Spring" = rep(NA, length(unique(res.cog$Species))),  "SeasonalSST.SeasonInteraction.Dist.Spring" = rep(NA, length(unique(res.cog$Species))),  "SeasonalSST.Independent.Fall.Dist.Spring" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.Independent.Spring.Dist.Spring" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.SeasonFactor.Bear.Spring" = rep(NA, length(unique(res.cog$Species))),  "SeasonalSST.SeasonInteraction.Bear.Spring" = rep(NA, length(unique(res.cog$Species))),  "SeasonalSST.Independent.Fall.Bear.Spring" = rep(NA, length(unique(res.cog$Species))), "SeasonalSST.Independent.Spring.Bear.Spring" = rep(NA, length(unique(res.cog$Species))))

for(i in seq_along(res.cog$Species)){
  spp.dat<- res.cog[res.cog$Species == res.cog$Species[i], ]
  fall.pts<- data.frame("x" = c(spp.dat$SeasonalSST.COG.Long.Fall, spp.dat$SeasonalSST.SeasonFactor.COG.Long.Fall, spp.dat$SeasonalSST.SeasonInteraction.COG.Long.Fall, spp.dat$SeasonalSST.Independent.Fall.COG.Long.Fall, spp.dat$SeasonalSST.Independent.Spring.COG.Long.Fall),  "y" = c(spp.dat$SeasonalSST.COG.Lat.Fall, spp.dat$SeasonalSST.SeasonFactor.COG.Lat.Fall, spp.dat$SeasonalSST.SeasonInteraction.COG.Lat.Fall, spp.dat$SeasonalSST.Independent.Fall.COG.Lat.Fall, spp.dat$SeasonalSST.Independent.Spring.COG.Lat.Fall))
  fall.dist<- (round(distm(fall.pts), 3)/1000)[-1,1]
  fall.bear<- bearing(fall.pts[1,], fall.pts[-1,])
  
  spring.pts<- data.frame("x" = c(spp.dat$SeasonalSST.COG.Long.Spring, spp.dat$SeasonalSST.SeasonFactor.COG.Long.Spring, spp.dat$SeasonalSST.SeasonInteraction.COG.Long.Spring, spp.dat$SeasonalSST.Independent.Fall.COG.Long.Spring, spp.dat$SeasonalSST.Independent.Spring.COG.Long.Spring),  "y" = c(spp.dat$SeasonalSST.COG.Lat.Spring, spp.dat$SeasonalSST.SeasonFactor.COG.Lat.Spring, spp.dat$SeasonalSST.SeasonInteraction.COG.Lat.Spring, spp.dat$SeasonalSST.Independent.Fall.COG.Lat.Spring, spp.dat$SeasonalSST.Independent.Spring.COG.Lat.Spring))
  spring.dist<- (round(distm(spring.pts), 3)/1000)[-1,1]
  spring.bear<- bearing(spring.pts[1,], spring.pts[-1,])
  
  res.all<- c(fall.dist, fall.bear, spring.dist, spring.bear)
  res.cog.db[i,2:ncol(res.cog.db)]<- res.all
}
  
# How can we visualize this??
col.inds<- c("Species", colnames(res.cog.db)[which(grepl("Dist", colnames(res.cog.db)))])
res.cog.dist.plot<- res.cog.db %>%
  dplyr::select(one_of(col.inds)) %>%
  gather(., "Statistic", "Distance", -Species)

res.cog.dist.plot$Season<- str_extract(res.cog.dist.plot$Statistic, "[^.]*$")
res.cog.dist.plot$Season<- factor(res.cog.dist.plot$Season, levels = c("Fall", "Spring"))

res.cog.dist.plot$Model.Form<- ifelse(grepl("SeasonFactor", res.cog.dist.plot$Statistic), "Both.SeasonFactor",
                             ifelse(grepl("SeasonInteraction", res.cog.dist.plot$Statistic), "Both.SeasonInteraction",
                                    ifelse(grepl("Independent.Fall", res.cog.dist.plot$Statistic), "FallOnly",
                                           ifelse(grepl("Independent.Spring", res.cog.dist.plot$Statistic), "SpringOnly", "Both.SeasonIgnored"))))
res.cog.dist.plot$Model.Form<- factor(res.cog.dist.plot$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))
res.cog.dist.plot<- res.cog.dist.plot %>%
  dplyr::select(., -Statistic)

col.inds<- c("Species", colnames(res.cog.db)[which(grepl("Bear", colnames(res.cog.db)))])
res.cog.bear.plot<- res.cog.db %>%
  dplyr::select(one_of(col.inds)) %>%
  gather(., "Statistic", "Bearing", -Species)

res.cog.bear.plot$Season<- str_extract(res.cog.bear.plot$Statistic, "[^.]*$")
res.cog.bear.plot$Season<- factor(res.cog.bear.plot$Season, levels = c("Fall", "Spring"))

res.cog.bear.plot$Model.Form<- ifelse(grepl("SeasonFactor", res.cog.bear.plot$Statistic), "Both.SeasonFactor",
                                      ifelse(grepl("SeasonInteraction", res.cog.bear.plot$Statistic), "Both.SeasonInteraction",
                                             ifelse(grepl("Independent.Fall", res.cog.bear.plot$Statistic), "FallOnly",
                                                    ifelse(grepl("Independent.Spring", res.cog.bear.plot$Statistic), "SpringOnly", "Both.SeasonIgnored"))))
res.cog.bear.plot$Model.Form<- factor(res.cog.bear.plot$Model.Form, levels = c("Both.SeasonIgnored", "Both.SeasonFactor", "Both.SeasonInteraction", "FallOnly", "SpringOnly"))
res.cog.bear.plot<- res.cog.bear.plot %>%
  dplyr::select(., -Statistic)

res.cog.plot<- res.cog.dist.plot %>%
  left_join(., res.cog.bear.plot)

res.cog.index<- res.cog %>%
  dplyr::select(., Species, MigrationIndexCode)

res.cog.plot<- res.cog.plot %>%
  left_join(., res.cog.index)

res.cog.plot.out<- ggplot(res.cog.plot, aes(x = Bearing, y = Distance, color = MigrationIndexCode)) +
  geom_boxplot() +
  scale_color_manual(name = "Migration Index", values = c("firebrick3", "darkolivegreen4", "dodgerblue4"), labels = c("Seasonally Migrating", "Seasonally Shifting", "Year-round Resident")) +
  theme_bw() +
  coord_polar(start = pi, clip = "off") +
  scale_x_continuous(limits = c(-180, 180),
                     breaks =  seq(-180, 180, 45)) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500)) +
  ylab("Distance difference from SDM with seasonal SST only (km)") +
  xlab("Directional bearing from SDM with seasonal SST only") +
  facet_wrap(~Season+Model.Form, nrow = 2) +
  theme(strip.background = element_blank())
ggsave("~/GitHub/SDMSeasonality/Results/Figures/ResCOGPlot.jpg", plot = res.cog.plot.out, width = 11, height = 8)


##### Generalized Linear Mixed Effect Model of Metric ~ Model Form + Migration Index + Species Random Effect
{
# Presence Absence - AUC
res.full$trial.size<- rep(1, nrow(res.full))
res.full$BA.Index.CentScale<- scale(res.full$BA.Index)

yind<- "AUC"
yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
keep<- c("Species", "BA.Index", "MigrationIndexCode", "BA.Index.CentScale", yind.keep)
mod.dat<- res.full %>%
  dplyr::select(., keep) %>%
  gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode, -BA.Index.CentScale)
mod.dat$Season<- str_extract(mod.dat$Scenario, "[^.]*$")
mod.dat$Season<- factor(mod.dat$Season, levels = c("Fall", "Spring"))
mod.dat$Model.Form<- ifelse(grepl("SeasonFactor", mod.dat$Scenario), "SeasonFactor",
                             ifelse(grepl("SeasonInteraction", mod.dat$Scenario), "SeasonInteraction",
                                    ifelse(grepl("Independent.Fall", mod.dat$Scenario), "SeasonIndependent.Fall",
                                           ifelse(grepl("Independent.Spring", mod.dat$Scenario), "SeasonIndependent.Spring", "SeasonSST"))))
mod.dat$Model<- factor(mod.dat$Model.Form, levels = c("SeasonSST", "SeasonFactor", "SeasonInteraction", "SeasonIndependent.Fall", "SeasonIndependent.Spring"))
mod.dat$AUC.Mod<- asin(mod.dat$Value)

## Fall and spring prediction - AUC
mod.dat.f<- mod.dat[mod.dat$Season == "Fall",]
mod.par2.all.lm<- lm(AUC.Mod ~ Model + BA.Index.CentScale, data = mod.dat.f)
mod.par2.all.lmer<- lmer(AUC.Mod ~ Model + BA.Index.CentScale + (1|Species), data = mod.dat.f)
anova(mod.par2.all.lmer, mod.par2.all.lm) # Random effect significant

# Test significance of fixed effects with wald test
vc <- vcov(mod.par2.all.lmer, useScale = FALSE)
b <- fixef(mod.par2.all.lmer)
se <- as.numeric(sqrt(diag(vc)))
z <- as.numeric(as.numeric(b)) / sqrt(diag(vc))
P <- as.numeric(2 * (1 - pnorm(abs(z))))
fall.auc.table<- data.frame("Model" = names(b), "b" = as.numeric(b), cbind(se, z, P))
fall.auc.table$Stat.Season<- rep("Fall.AUC", nrow(fall.auc.table))

mod.dat.s<- mod.dat[mod.dat$Season == "Spring",]
mod.par2.all.lm<- lm(AUC.Mod ~ Model + BA.Index.CentScale, data = mod.dat.s)
mod.par2.all.lmer<- lmer(AUC.Mod ~ Model + BA.Index.CentScale + (1|Species), data = mod.dat.s)
anova(mod.par2.all.lmer, mod.par2.all.lm) # Random effect significant

# Test significance of fixed effects with wald test
vc <- vcov(mod.par2.all.lmer, useScale = FALSE)
b <- fixef(mod.par2.all.lmer)
se <- as.numeric(sqrt(diag(vc)))
z <- as.numeric(as.numeric(b) / sqrt(diag(vc)))
P <- as.numeric(2 * (1 - pnorm(abs(z))))
spring.auc.table<- data.frame("Model" = names(b), "b" = as.numeric(b), cbind(se, z, P))
spring.auc.table$Stat.Season<- rep("Spring.AUC", nrow(spring.auc.table))

auc.table.out<- bind_rows(fall.auc.table, spring.auc.table)
write.csv(auc.table.out, "~/GitHub/SDMSeasonality/Results/Tables/AUClmer.csv")

## Fall and Spring: Calibration
yind<- "Calib"
yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
keep<- c("Species", "BA.Index", "MigrationIndexCode", "BA.Index.CentScale", yind.keep)
mod.dat<- res.full %>%
  dplyr::select(., keep) %>%
  gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode, -BA.Index.CentScale)
mod.dat$Season<- str_extract(mod.dat$Scenario, "[^.]*$")
mod.dat$Season<- factor(mod.dat$Season, levels = c("Fall", "Spring"))
mod.dat$Model.Form<- ifelse(grepl("SeasonFactor", mod.dat$Scenario), "SeasonFactor",
                            ifelse(grepl("SeasonInteraction", mod.dat$Scenario), "SeasonInteraction",
                                   ifelse(grepl("Independent.Fall", mod.dat$Scenario), "SeasonIndependent.Fall",
                                          ifelse(grepl("Independent.Spring", mod.dat$Scenario), "SeasonIndependent.Spring", "SeasonSST"))))
mod.dat$Model<- factor(mod.dat$Model.Form, levels = c("SeasonSST", "SeasonFactor", "SeasonInteraction", "SeasonIndependent.Fall", "SeasonIndependent.Spring"))

# Some weirdness going on here as we have a few species with negative calibration values.
hist(mod.dat$Value)
temp<- mod.dat[mod.dat$Value <0,]
summary(temp)
unique(temp$Species)

best.norm<- bestNormalize(mod.dat$Value) # Picked the ordered, not a big fan of that one. What about second best, the Yeo Johnson one?
yeojohn<- yeojohnson(mod.dat$Value)
hist(yeojohn$x.t)
mod.dat$Calib.Mod<- yeojohnson(mod.dat$Value)$x.t

mod.dat.f<- mod.dat[mod.dat$Season == "Fall",]
mod.par2.all.lm<- lm(Calib.Mod ~ Model + BA.Index.CentScale, data = mod.dat.f)
mod.par2.all.lmer<- lmer(Calib.Mod ~ Model + BA.Index.CentScale + (1|Species), data = mod.dat.f)
qqnorm(resid(mod.par2.all.lmer))
qqline(resid(mod.par2.all.lmer))
anova(mod.par2.all.lmer, mod.par2.all.lm) # Random effect significant

# Test significance of fixed effects with wald test
vc <- vcov(mod.par2.all.lmer, useScale = FALSE)
b <- fixef(mod.par2.all.lmer)
se <- as.numeric(sqrt(diag(vc)))
z <- as.numeric(as.numeric(b) / sqrt(diag(vc)))
P <- as.numeric(2 * (1 - pnorm(abs(z))))
fall.calib.table<- data.frame("Model" = names(b), "b" = as.numeric(b), cbind(se, z, P))
fall.calib.table$Stat.Season<- rep("Fall.Calib", nrow(fall.calib.table))

mod.dat.s<- mod.dat[mod.dat$Season == "Spring",]
mod.par2.all.lm<- lm(Calib.Mod ~ Model + BA.Index.CentScale, data = mod.dat.s)
mod.par2.all.lmer<- lmer(Calib.Mod ~ Model + BA.Index.CentScale + (1|Species), data = mod.dat.s)
anova(mod.par2.all.lmer, mod.par2.all.lm) # Random effect significant

# Test significance of fixed effects with wald test
vc <- vcov(mod.par2.all.lmer, useScale = FALSE)
b <- fixef(mod.par2.all.lmer)
se <- as.numeric(sqrt(diag(vc)))
z <- as.numeric(as.numeric(b) / sqrt(diag(vc)))
P <- as.numeric(2 * (1 - pnorm(abs(z))))
spring.calib.table<- data.frame("Model" = names(b), "b" = as.numeric(b), cbind(se, z, P))
spring.calib.table$Stat.Season<- rep("Spring.Calib", nrow(spring.calib.table))

calib.table.out<- bind_rows(fall.calib.table, spring.calib.table)
write.csv(calib.table.out, "~/GitHub/SDMSeasonality/Results/Tables/Caliblmer.csv")

## Fall and Spring: RMSE Biomass
yind<- "RMSE.B" 
yind.keep<- colnames(res.full)[which(grepl(yind, names(res.full)))]
keep<- c("Species", "BA.Index", "MigrationIndexCode", "BA.Index.CentScale", yind.keep)
mod.dat<- res.full %>%
  dplyr::select(., keep) %>%
  gather(., "Scenario", "Value", -Species, -BA.Index, -MigrationIndexCode, -BA.Index.CentScale)
mod.dat$Season<- str_extract(mod.dat$Scenario, "[^.]*$")
mod.dat$Season<- factor(mod.dat$Season, levels = c("Fall", "Spring"))
mod.dat$Model.Form<- ifelse(grepl("SeasonFactor", mod.dat$Scenario), "SeasonFactor",
                            ifelse(grepl("SeasonInteraction", mod.dat$Scenario), "SeasonInteraction",
                                   ifelse(grepl("Independent.Fall", mod.dat$Scenario), "SeasonIndependent.Fall",
                                          ifelse(grepl("Independent.Spring", mod.dat$Scenario), "SeasonIndependent.Spring", "SeasonSST"))))
mod.dat$Model<- factor(mod.dat$Model.Form, levels = c("SeasonSST", "SeasonFactor", "SeasonInteraction", "SeasonIndependent.Fall", "SeasonIndependent.Spring"))
mod.dat$RMSE.Mod<- log(mod.dat$Value)

mod.dat.f<- mod.dat[mod.dat$Season == "Fall",]
mod.par2.all.lm<- lm(RMSE.Mod ~ Model + BA.Index.CentScale, data = mod.dat.f)
mod.par2.all.lmer<- lmer(RMSE.Mod ~ Model + BA.Index.CentScale + (1|Species), data = mod.dat.f)
anova(mod.par2.all.lmer, mod.par2.all.lm) # Random effect significant

# Test significance of fixed effects with wald test
vc <- vcov(mod.par2.all.lmer, useScale = FALSE)
b <- fixef(mod.par2.all.lmer)
se <- as.numeric(sqrt(diag(vc)))
z <- as.numeric(as.numeric(b) / sqrt(diag(vc)))
P <- as.numeric(2 * (1 - pnorm(abs(z))))
fall.rmse.table<- data.frame("Model" = names(b), "b" = as.numeric(b), cbind(se, z, P))
fall.rmse.table$Stat.Season<- rep("Fall.RMSE", nrow(fall.rmse.table))

mod.dat.s<- mod.dat[mod.dat$Season == "Spring",]
mod.par2.all.lm<- lm(RMSE.Mod ~ Model + BA.Index.CentScale, data = mod.dat.s)
mod.par2.all.lmer<- lmer(RMSE.Mod ~ Model + BA.Index.CentScale + (1|Species), data = mod.dat.s)
anova(mod.par2.all.lmer, mod.par2.all.lm) # Random effect significant

# Test significance of fixed effects with wald test
vc <- vcov(mod.par2.all.lmer, useScale = FALSE)
b <- fixef(mod.par2.all.lmer)
se <- as.numeric(sqrt(diag(vc)))
z <- as.numeric(as.numeric(b) / sqrt(diag(vc)))
P <- as.numeric(2 * (1 - pnorm(abs(z))))
spring.rmse.table<- data.frame("Model" = names(b), "b" = as.numeric(b), cbind(se, z, P))
spring.rmse.table$Stat.Season<- rep("Spring.RMSE", nrow(spring.rmse.table))

rmse.table.out<- bind_rows(fall.rmse.table, spring.rmse.table)
write.csv(rmse.table.out, "~/GitHub/SDMSeasonality/Results/Tables/RMSElmer.csv")
}


#####
## Study area map
#####
nelme<- st_read("~/GitHub/SDMSeasonality/Data/NELME_clipped.shp")
st_crs(nelme)<- "+init=epsg:4326"

nelme.sp<- as(st_sf(nelme), "Spatial")

bstrat<- st_read("~/GitHub/SDMSeasonality/Data/BottomTrawlStrata/BTS_Strata.shp") %>%
  st_transform(., "+init=epsg:4326")
strata.ca<- c(1351, 1310, 1320, 1410, 1420, 1490, 1990, 1410, 1420, 1490, 5440, 5480, 5430) # Canada

bstrat<- bstrat %>%
  filter(., !STRATA %in% strata.ca & STRATA <= 3990) %>%
  filter(., !as.numeric(STRATUMA) %in% strata.ca & as.numeric(STRATUMA) <= 3990)

# GoM
gom<- st_read("~/GitHub/SDMSeasonality/Data/GoMPhysioRegions/PhysioRegions_WGS84.shp")
st_crs(gom)<- "+init=epsg:4326"
gom<- gom[!gom$Region == "Seamount",] %>%
  st_union() 
gom.sp<- gom %>%
  st_sf()
gom.sp<- as(st_zm(gom), "Spatial")
gom.sp<- spTransform(gom.sp, proj.utm)
gom.sp.wgs<-  spTransform(gom.sp, st_crs(nelme)$proj4string)

# Buffer it a bit
gom.buff<- gBuffer(gom.sp, width = 7500)
gom.buff<- spTransform(gom.buff, st_crs(nelme)$proj4string)

# Southern regions
south<- erase(nelme.sp, gom.buff)

# Still a bit remaining...custom box to get rid of the rest of it
# Coordinates
ow<- data.frame("x" = c(-71, -71, -67, -67), "y" = c(42, 46, 46, 42))

# Convert coordinates to Spatial Polygons
ow.p<- Polygon(ow)
ow.ps<- Polygons(list(ow.p), 1)
ow.sp<- SpatialPolygons(list(ow.ps))
proj4string(ow.sp)<- proj4string(nelme.sp)
south2<- erase(south, ow.sp)
proj4string(south2)<- proj4string(nelme.sp)
south<- st_as_sf(south2)

# Spatial stuff -- gets us the states and shoreline
# Spatial projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

#Bounds
xlim.use<- c(-77, -65)
ylim.use<- c(35.05, 45.2)

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us<- raster::getData("GADM",country="USA",level=1)
us.states<- us[us$NAME_1 %in% states,]
us.states<- gSimplify(us.states, tol=0.01, topologyPreserve=TRUE)
us.states<- st_as_sf(us.states)
canada<- raster::getData("GADM",country="CAN",level=1)
ca.provinces<- canada[canada$NAME_1 %in% provinces,]
ca.provinces<- gSimplify(ca.provinces, tol=0.01, topologyPreserve=TRUE)
ca.provinces<- st_as_sf(ca.provinces)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

# Alright, plot time
plot.out<- ggplot() + 
  geom_sf(data = us.states, fill = "white", lwd = 0.4, show.legend = FALSE) +
  geom_sf(data = ca.provinces.f, fill = "white", lwd = 0.4, show.legend = FALSE) +
  geom_sf(data = bstrat, fill = "light gray", color = "black", show.legend = FALSE) +
  #geom_sf(data = gom, fill = NA, aes(color = "#377eb8"), lwd = 1.5, show.legend = TRUE) +
  #geom_sf(data = south, fill = NA, aes(color = "#ff7f00"), lwd = 1.5, show.legend = TRUE) +
  #scale_color_manual(name = "Region", values = c("#377eb8", "#ff7f00"), labels = c("Gulf of Maine", "Southern NELME")) +
  xlim(xlim.use) +
  ylim(ylim.use) +
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"))
ggsave("~/GitHub/SDMSeasonality/Results/Figures/StudyArea.jpg", plot.out, width = 8, height = 11)

us.large<- gSimplify(us, tol=0.1, topologyPreserve=TRUE)
us.large<- st_as_sf(us.large)

ca.large<- gSimplify(canada, tol=0.1, topologyPreserve=TRUE)
ca.large<- st_as_sf(ca.large)

plot.large<- ggplot() + 
  geom_sf(data = us.large, fill = "white", lwd = 0.7) +
  geom_sf(data = ca.large, fill = "white", lwd = 0.7) +
  xlim(c(-100, 0)) +
  ylim(c(25, 55)) +
  coord_sf(datum = NA) 
ggsave("~/GitHub/SDMSeasonality/Results/Figures/StudyAreaOverview.jpg", plot.large, width = 11, height = 8)
