# EXPLORE SATELLITE DATA FROM LANDSAT 7
# Proyecto de examen de medio termino de los chicos de RMS

# Install and load packages
install.packages("terra")
install.packages("raster")
install.packages("sf")
install.packages("readr")
install.packages("remotes")
remotes::install_github("bleutner/RStoolbox", force = TRUE)
library(terra)
library(raster)  
library(sf)             
library(RStoolbox)
library(readr)

#sumar lineas de dialogo 

rm(list = ls())

##########
# 1. Import meta-data Landsat 7 and bands based on MTL file for extracting information
setwd("D:/Process_landsat/010060_2001")

# 1. Leer el MTL como texto y extraer la zona UTM 
mtl_file <- "L71010060_06020011103_MTL.txt" 
mtl_text <- readLines(mtl_file, warn = FALSE)
print(mtl_text)

# Buscar la zona UTM (ej: "UTM_ZONE = 19")
utm_zone <- grep("UTM_ZONE", mtl_text, value = TRUE)
utm_zone <- as.numeric(gsub("[^0-9]", "", utm_zone))  # Extraer solo el número

# 2. Definir el CRS 
if (length(utm_zone) > 0) {
  crs_utm <- paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84 +units=m +no_defs")
} else {
  crs_utm <- "+proj=longlat +datum=WGS84"  # Fallback a WGS84 (lat/lon)
}

# 3. Cargar las bandas del Landsat 
band_paths <- list.files(pattern = "*_B[1-8]0\\.TIF$", full.names = TRUE)
# Reproyectar todas las bandas al CRS y resolución de la referencia
band_list_aligned <- lapply(band_paths, function(x) {
  project(rast(x), crs_utm)
})

stack_landsat <- rast(band_list_aligned[1:6])
crs(stack_landsat) <- crs_utm
# Definir una banda como referencia (ej: B20)
ref_band <- rast("L71010060_06020011103_B20.TIF")

# Verificar
print(crs(img))
plot(stack_landsat, 
     main = paste("Bandas Landsat -", names(stack_landsat)),
     nc = 3, nr = 3)

# 1. Create SpatRaster object from the bands we already loaded
# (Using stack_landsat from your original code instead of stackMeta)
ls <- stack_landsat
print(ls)

# Plot all bands
terra::plot(ls, main = "Landsat 7 Bands")

# 2. Spatial crop using the same extent from your code
print(terra::ext(ls))
e <- terra::ext(740000, 810000, -40000, 30000)
ls_crop <- terra::crop(ls, e)
terra::plot(ls_crop)

# 3. Saving results to disk
if(!dir.exists("output")) dir.create("output")
terra::writeRaster(ls_crop, filename="output/landsat7_crop.tif", overwrite=TRUE)
terra::writeRaster(ls, filename="output/landsat7.tif", overwrite=TRUE)

# 4. Perform DN to TOA reflectance (using metaData from text parsing)
# First we need to create a metadata list similar to what stackMeta would produce
meta_list <- list(
  RADIANCE_MAXIMUM_BAND = as.numeric(gsub(".*= ([0-9.]+).*", "\\1", 
                                          grep("LMAX_BAND", mtl_text, value = TRUE))),
  RADIANCE_MINIMUM_BAND = as.numeric(gsub(".*= ([0-9.]+).*", "\\1", 
                                          grep("LMIN_BAND", mtl_text, value = TRUE))),
  DATE_ACQUIRED = grep("ACQUISITION_DATE", mtl_text, value = TRUE) %>% 
    gsub(".*= (.*)", "\\1", .),
  SUN_ELEVATION = as.numeric(gsub(".*= ([0-9.]+).*", "\\1", 
                                  grep("SUN_ELEVATION", mtl_text, value = TRUE)))
)

ls_crop_ref <- ls_crop # In terra, radCor equivalent would be manual calculation
# Note: For precise atmospheric correction, consider using the landsat package or 
# manual calculation based on the metadata values

# Rename bands for Landsat 7
names(ls_crop_ref) <- c('blue', 'green', 'red', 'NIR', 'SWIR1', 'TIR1', 'SWIR2')
print(names(ls_crop_ref))

# Export raster
terra::writeRaster(ls_crop_ref, filename = "output/landsat7_crop_ref.tif", overwrite=TRUE)

# 5. Relation between bands
terra::pairs(ls_crop[[1:7]], main = "Bands relationship") 
terra::pairs(ls_crop[[c(3,4)]], main = "Red vs NIR relationship")

# 6. Stretching reflectance values
terra::hist(ls_crop_ref$red,
            main = "Distribution of reflectance values",
            xlab = "Radiance", ylab= "Frequency",
            col = "red", xlim = c(0, 1), breaks = 100, xaxt = 'n')
axis(side=1, at = seq(0,1, 0.1), labels = seq(0,1, 0.1))

red_strch <- terra::stretch(ls_crop_ref$red, minq=0.02, maxq=0.9)
green_strch <- terra::stretch(ls_crop_ref$green, minq=0.02, maxq=0.9)
blue_strch <- terra::stretch(ls_crop_ref$blue, minq=0.02, maxq=0.9)

# 7. False Color Composite plots
par(mfrow = c(1,2))
landsatRGB1 <- c(ls_crop_ref$red, ls_crop_ref$green, ls_crop_ref$blue)
landsatRGB2 <- c(red_strch, green_strch, blue_strch)
terra::plotRGB(landsatRGB1, axes = TRUE, main = "True Color - No Stretch") 
terra::plotRGB(landsatRGB2, axes = TRUE, main = "True Color - Stretched") 

terra::writeRaster(landsatRGB1, filename = "output/landsat7_crop_ref_RGB_nstch.tif", overwrite=TRUE)
terra::writeRaster(landsatRGB2, filename = "output/landsat7_crop_ref_RGB_stch.tif", overwrite=TRUE)

# 8. NDVI calculation
ndvi <- function(img, nir_band, red_band) {
  (img[[nir_band]] - img[[red_band]]) / (img[[nir_band]] + img[[red_band]])
}
ls_crop_ref_ndvi <- ndvi(ls_crop_ref, 4, 3) # NIR=4, Red=3 in our renamed bands
terra::plot(ls_crop_ref_ndvi, col = rev(terrain.colors(10)), main = "NDVI")
terra::writeRaster(ls_crop_ref_ndvi, filename = "output/ndvi_l7.tif", overwrite=TRUE)

# 9. NDSI calculation
ndsi <- function(img, green_band, swir_band) {
  (img[[green_band]] - img[[swir_band]]) / (img[[green_band]] + img[[swir_band]])
}
ls_crop_ref_ndsi <- ndsi(ls_crop_ref, 2, 5) # Green=2, SWIR=5 in our renamed bands
terra::plot(ls_crop_ref_ndsi, col = rev(terrain.colors(10)), main = "NDSI")
terra::writeRaster(ls_crop_ref_ndsi, filename = "output/ndsi_l7.tif", overwrite=TRUE)

# 10. Glacier classification
ls_crop_ref_ndsi_class <- terra::classify(ls_crop_ref_ndsi, 
                                          matrix(c(-Inf, 0.4, 0,
                                                   0.4, Inf, 1), 
                                                 ncol=3))
par(mfrow = c(1,2))
terra::plot(ls_crop_ref_ndsi, col = rev(terrain.colors(10)), main = 'NDSI')
terra::plot(ls_crop_ref_ndsi_class, col = c('grey90', 'blue'), main = 'Glacier Classification')
terra::writeRaster(ls_crop_ref_ndsi_class, filename = "output/landsat7_ndsi_class.tif", overwrite=TRUE)

# Reset graphical parameters
par(mfrow = c(1,1))