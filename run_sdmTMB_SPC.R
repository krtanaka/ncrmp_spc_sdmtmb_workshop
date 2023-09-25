# Load necessary libraries
library(sdmTMB)      # For sdmTMB functions
library(dplyr)       # For data manipulation
library(ggplot2)     # For creating plots
library(rgdal)       # For spatial data operations
library(colorRamps)  # For color ramps
library(raster)      # For raster data manipulation
library(sf)          # For working with Simple Features
library(readr)       # For reading data
library(ggthemes)    # For additional ggplot themes
library(tidyr)       # For data tidying

# Clear the workspace
rm(list = ls())

# Select specific functions from the dplyr package
select = dplyr::select

# Load clean NCRMP survey data
load("data/clean_df.RData")

# Calculate the UTM zone based on the longitude of the first data point
# Project latitude and longitude coordinates to UTM
# Rename the columns for the UTM coordinates
# Add the UTM coordinates to the original data frame
zone <- (floor((df$lon[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("lon", "lat")]), paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("X", "Y")
df = cbind(df, xy_utm)
plot(xy_utm, pch = ".", bty = 'n')

# Read in Island Boundaries, also re-project latitude and longitude coordinates to UTM 
load('data/MHI_islands_shp.RData')
crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ISL_this = ISL_bounds[which(ISL_bounds$ISLAND %in% toupper(unique(df$island))),]
ISL_this_utm = spTransform(ISL_this,CRS(paste0("+proj=utm +units=km +zone=", zone)))
ISL_this_sf = st_transform(st_as_sf(ISL_this), crs = paste0("+proj=utm +units=km +zone=", zone))

rea_spde <- make_mesh(df, c("X", "Y"), cutoff  = 15, type = "cutoff")

# add physical barriers (i.e. coastlines and land mass) to mesh
rea_spde_coast = add_barrier_mesh(rea_spde , ISL_this_sf)
plot(rea_spde_coast$mesh, asp = 1, main = ""); axis(1); axis(2)
plot(ISL_this_utm, add = TRUE)
points(rea_spde_coast$loc_xy,col = "green", pch = ".", cex = 5)
bar_i = rea_spde_coast$barrier_triangles
norm_i = rea_spde_coast$normal_triangles
points(rea_spde_coast$spde$mesh$loc[,1], rea_spde_coast$spde$mesh$loc[,2], pch = ".", col = "black")
points(rea_spde_coast$mesh_sf$V1[bar_i], rea_spde_coast$mesh_sf$V2[bar_i], col = "red", pch = 20, cex = 0.5)
points(rea_spde_coast$mesh_sf$V1[norm_i], rea_spde_coast$mesh_sf$V2[norm_i], col = "blue", pch = 20, cex = 0.5)

# fit a spatiotemporal generalzied mixed additive model
fit <- sdmTMB(
  
  data = df, 
  formula = response ~ as.factor(year) + s(depth, k = 3),
  silent = F, 
  time = "year",
  mesh = rea_spde_coast,
  family = tweedie(link = "log")
  
); beepr::beep(2)

# check model fit & some diagnostics
fit
max(fit$gradients)
AIC(fit)
tidy(fit, conf.int = T)
tidy(fit, effects = "ran_pars", conf.int = TRUE)
sanity(fit)

visreg::visreg(fit, xvar = "depth", xlim = c(0, 30))
visreg::visreg(fit, xvar = "depth", scale = "response", xlim = c(0, 30), nn = 200)

# Make spatial prediction using a spatially resolved gridded data
load("data/mhi_grid.rdata")
p <- predict(fit, newdata = grid)

p %>% 
  group_by(X, Y, year) %>% 
  summarise(est = mean(est, na.rm = T)) %>% 
  ggplot(aes(X, Y, color = exp(est))) + 
  geom_point(shape = 22, alpha = 0.5, size = 0.5) + 
  coord_fixed() + 
  scale_color_viridis_c(trans = "sqrt") + 
  facet_wrap(~year)

# Get abundance/biomass index
p <- predict(fit, newdata = grid, return_tmb_object = T)
index <- get_index(p, area = rep(0.0081, nrow(grid)))

ggplot(index, aes(as.factor(year), est)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
  scale_x_discrete(limits = as.factor(2010:2019)) +
  labs(x = "Year", y = "Biomass (g)")

# Get Center of Gravity
cog <- get_cog(p, format = "wide")
ggplot(cog, aes(est_x, est_y, colour = year)) +
  # geom_point(data = df, aes(X, Y), color = "gray20") + 
  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y))
