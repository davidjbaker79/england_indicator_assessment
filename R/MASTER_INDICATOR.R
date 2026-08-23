#=============================================================================#
#
# Title: Alternative interpretations of a biodiversity target reshape 
# conservation priorities
# 
# Analysis developed by David J Baker (2025-2026)
# Co-authors: Ilya MD Maclean, Kevin J Gaston, Marcus Rhodes
#
# Funding: This work was supported by the Natural Environment Research Council 
# funded “Renewing biodiversity through a people-in- nature approach (RENEW)” 
# project (no. NE/W004941/1). 
#
# Copyright: For the purpose of open access, the authors have applied a 
# Creative Commons Attribution (CC BY) license to outputs from this work.
#
#==============================================================================#

#================================== SUMMARY ===================================#
#
# The aim is to evaluate the trade-off between the spatial footprint of the
# conservation planning area required to meet England's statutory species 
# abundance targets and the constraints placed on the number of species meeting
# individual abundance targets.
#
# The analysis is theoretical, but grounded in the empirical system of England's
# All Species Abundance Indicator. The analysis uses a heuristic approach to 
# finding near-optimal minimum spatial footprint for conservation planning for 
# each constraint value (K = number of species also reaching +10% abundance 
# target).
#
# Thus, two simplifying assumptions are made:
#
# (1) species’ abundance is distributed across the occupied cells in proportion 
#     to the amount of suitable habitat currently present in each cell; and 
# (2) benefits of conservation action for a species only occur within occupied 
#     cells (no range expansion beyond the 10 km x 10 km scale of the analysis). 
# 
# We also assumed a linear and equal per-cell effect of treatment on abundance, 
# such that the post-treatment abundance multiplier for species is equal to 
# 1 + beta*P[cs], where Pcs = the proportion of a species range in cell c and 
# beta scales the response.
#
#==============================================================================# 

#==============================================================================#
# Setup
#==============================================================================#

# === Set library paths ===
.libPaths()

# === Define directories ===
sysPath <- Sys.getenv("OneDrive")
shStdDat <- file.path(sysPath, "RENEW/AA_RENEW_STANDARD_PRODUCTS/outputs")
shareDir <- file.path(sysPath,"RENEW/renew_data/")

# === Load packages ===
library(renewTools)
library(dplyr)
library(sf)
library(tidyr)
library(terra)
library(data.table)
library(exactextractr)
library(ggplot2)
library(readxl)
library(openxlsx)
library(stringr)
library(future)
library(future.apply)
library(lwgeom)
library(gtools)
library(scam)
library(scico)
library(ggspatial)
library(patchwork)
library(gratia)
library(cowplot)

# === Source optimisation functions ===
source("R/optimisation_functions.R")

# === Source plotting functions ===
source("R/build_pareto_front_plot.R")
source("R/plotting_functions.R")

# === GB boundary ===
gb <- st_read(file.path(shStdDat, "gb_mlwm_boundary.gpkg"))

# === Create GB OS Grid ===

# = Just make once
# GB10km <- make_os_grid(gb, 10000)
# GB10km <- st_intersection(GB10km, gb)
# st_write(GB10km, "outputs/GB10km.gpkg")

# = Load and process
GB10km <- st_read("outputs/GB10km.gpkg")
eng10km2 <- GB10km |>
  filter(CTRY22NM == "England") |>
  st_make_valid() |>
  st_collection_extract("POLYGON") |>
  st_cast("MULTIPOLYGON") |>
  group_by(os_gridref) |>
  summarise(.groups = "drop") |>
  select(os_gridref) |>
  mutate(area = as.numeric(st_area(geom)))

# === Total area km2 ===
total_area <- sum(eng10km2$area) / 1e+06

# === Cell by area ===
cellXarea <- eng10km2 |>
  st_drop_geometry() |> 
  mutate(area = area / 1e+06)

# === 10KM base raster ===
base_raster_10km <- rast(file.path(shStdDat,"BNG_10KM_GB_TEMPLATE.tif"))

# === Cell to os_gridref ===
cent <- terra::vect(sf::st_centroid(st_geometry(GB10km)))
GB10km$cell <- terra::cellFromXY(base_raster_10km, terra::crds(cent))
cell_to_os_gridref_lkup <- distinct(st_drop_geometry(GB10km[,c("cell","os_gridref")]))

# === Habitat data

# = Habitat type conversion lookup 
habitat_types_lookup <- tribble(
  ~ebrel_hab_code, ~habitat, ~narrow_classes,
  1, "deciduous_woodland", "BLWDY|BLWDM|SCRB",
  2, "wood_pasture", "WDPS",
  3, "wet_woodland", "WTWD",
  4, "neutral_grassland", "NUGR",
  5, "calcareous_grassland", "CAGR",
  6, "acid_grassland", "ACGR",
  7, "heathland", "HEATH|INRCK",
  8, "purple_moor_rush_pasture", "WTGR",
  9, "peatland", "BOG",
  10, "wetland", "FENS|FWLAK|FWPND|FWRIV",
  11, "sand_dunes", "CSTDUN|CSTSHG|CSTCLF",
  12, "intertidal", "CSTRCK|INTIDL|SLTW|SLTWLAG|SLTMA",
  #13, "arable_nature", "ARABNAT", # This is the old code for AES
  14, "arable", "ARAB",
  15, "improved_grassland", "IMGR",
  16, "coniferous_woodland", "COWDY|COWDM",
  17, "urban", "URBBLT|URBGSP"
)

# = Augmented CEH land cover for England
# Ebrel version matching habitat_types_lookup
habExist100m <- rast(file.path(shStdDat, "ebrel_landcover_100m_england_wales.tif"))


#==============================================================================#
# Compile species information to drive analysis
#==============================================================================#

# This is underpinned by RENEW trait database used to drive Ebrel

# === Species trait data for indicator species ===
# Used version current to 02_02_2026
species_info <- 
  fread("data/EBREL_MASTER_SPECIES_TRAITS_DATABASE_02_02_2026.csv")

# === Subset to ASAI ===
asai_info <- species_info[grep(species_info$LIST_INCL, pattern = "ASAI")]

# = Lesser spotted woodpecker name alignment is wrong across lists
idx <- agrep(asai_info$SCIENTIFIC_NAME, pattern = "Dendrocopos minor subsp. comminutus")
asai_info$SCIENTIFIC_NAME[idx] <- "Dryobates minor"

# === Name check
# Original BDI data used in creating species_info, but need to recheck for exact
# name matches and completeness to make sure names fully align
ASAI <- read_excel("data/Eng_BDI_Species_Master.xlsx", sheet = "ASAI") 
asai_info$BDI_NAME_CHECK <- NA

for (i in 1:nrow(asai_info)) {
  
  # Scientific name
  sci_nam <- asai_info$SCIENTIFIC_NAME[i]
  sci_nam
  
  # Check identical matching names
  exact_match <- sci_nam %in% ASAI$SCIENTIFIC_NAME
  
  # Either assign or check for approx match
  if (exact_match) {
    asai_info$BDI_NAME_CHECK[i] <- sci_nam
  } else {
    # Approximate match
    approx_match <- ASAI$SCIENTIFIC_NAME[agrep(ASAI$SCIENTIFIC_NAME, pattern = sci_nam)]
    
    # Assign if length 1
    if (length(approx_match) == 1) {
      asai_info$BDI_NAME_CHECK[i] <- approx_match
    } else {
      print(sci_nam)
      # Check for subspecies names and try sub-setting to just Genus species 
      gs_nam <- paste0(strsplit(sci_nam, split = " ")[[1]][1:2], collapse = " ")
      gs_match <- ASAI$SCIENTIFIC_NAME[agrep(ASAI$SCIENTIFIC_NAME, pattern = gs_nam)]
      print(gs_match)
      
      if (length(gs_match) == 1) {
        asai_info$BDI_NAME_CHECK[i] <-  gs_match
      } else {
       asai_info$BDI_NAME_CHECK[i] <- "CHECK"
      }
    }
  }
}

# = Identify name misalignments with trait database (i.e. CHECK)
which(asai_info$BDI_NAME_CHECK == "CHECK")
# [1] 389 422 888

# = Address name misalignments to adjust
# "Phengaris arion" = "Maculinea arion"
which(asai_info$SCIENTIFIC_NAME == "Maculinea arion") # Check exists
asai_info$BDI_NAME_CHECK[which(asai_info$SCIENTIFIC_NAME == "Phengaris arion")] <- "Maculinea arion"

# = "Anchoscelis litura" == "Agrochola litura" 
which(asai_info$SCIENTIFIC_NAME == "Agrochola litura") # Check exists
# This does, so just remove row of non-matching name
asai_info <- asai_info[!which(asai_info$SCIENTIFIC_NAME == "Anchoscelis litura"),] 

# = "Myotis mystacinus/brandtii" == "Myotis mystacinus" & "Myotis brandtii"
which(asai_info$SCIENTIFIC_NAME == "Myotis mystacinus/brandtii") # Check exists
# This doesn't exist as an aggregate but I think it should be kept separate because
# the mapping and traits are separate and distinct.
asai_info$BDI_NAME_CHECK[which(asai_info$SCIENTIFIC_NAME == "Myotis mystacinus")] <- "Myotis mystacinus"
asai_info$BDI_NAME_CHECK[which(asai_info$SCIENTIFIC_NAME == "Myotis brandtii")] <- "Myotis brandtii"

# = Miss-spelling
asai_info$SCIENTIFIC_NAME[which(asai_info$SCIENTIFIC_NAME == "Geum urbnum")] <- "Geum urbanum"

# = Recheck
which(asai_info$BDI_NAME_CHECK == "CHECK") # NONE

# = Check for missing TAXON_KEY
which(asai_info$TAXON_KEY == "") # NONE MISSING

# = CHECK incomplete habitat encoding
which(!asai_info$COMPLETED) # NONE MISSING

# ==== 
# === Now, need to check for species distribution mapping completeness
# ====

# 1) Check for mapping
# 2) Mark as completed
# 3) Move .tif to project data folder

# === tif check ===
asai_info$mapping_10km_complete <- FALSE
asai_info$mapping_type <- NA
asai_info$file_path <- NA

for (i in 1:nrow(asai_info)) { print(i)
  # Species
  SCIENTIFIC_NAME <- gsub(" |\\/", "_", gsub("\\.", "", asai_info$SCIENTIFIC_NAME[i]))
  SCIENTIFIC_NAME_ALT <- gsub(" |\\/", "_", gsub("\\.", "", asai_info$BDI_NAME_CHECK[i]))
  
  # Exact-match search string with .tif extension
  SCIENTIFIC_NAME_STR <- paste0("^", SCIENTIFIC_NAME, "\\.tif$")
  if (SCIENTIFIC_NAME != SCIENTIFIC_NAME_ALT) {
    SCIENTIFIC_NAME_STR <- paste0(
      "^(", SCIENTIFIC_NAME, "|", SCIENTIFIC_NAME_ALT, ")\\.tif$"
    )
  }
  
  # Check completed
  if (!asai_info$mapping_10km_complete[i]) {
    
    # First, check atlas data
    atlas_i <- character(0)
    file_path <- NULL
    if (asai_info$TAXON_GROUP_2[i] == "Bird") {
      file_path <- "BIOLOGICAL_RECORD_DATA/bird_atlas_2008_2011_10km_mapped"
      atlas_i <- list.files(file.path(shareDir, file_path), pattern = SCIENTIFIC_NAME_STR)
    } else if (asai_info$TAXON_GROUP_2[i] == "Butterfly") {
      file_path <- "BIOLOGICAL_RECORD_DATA/butterfly_atlas_2015_2019_10km_mapped"
      atlas_i <- list.files(file.path(shareDir, file_path), pattern = SCIENTIFIC_NAME_STR)
    } else if (asai_info$TAXON_GROUP_1[i] == "Vascular plant") {
      file_path <- "BIOLOGICAL_RECORD_DATA/plant_atlas_2020_10km_mapped"
      atlas_i <- list.files(file.path(shareDir, file_path), pattern = SCIENTIFIC_NAME_STR)
    }
    
    if (length(atlas_i) > 0) {
      # Check off tif
      asai_info$mapping_10km_complete[i] <- TRUE

      # Mark off mapping type
      asai_info$mapping_type[i] <- "Atlas"

      # File path 
      asai_info$file_path[i] <- file.path(shareDir, file_path, atlas_i)
      
      if (length(atlas_i) > 1) 
        break
    
      # Copy tif to project folder
      file.copy(from = file.path(shareDir, file_path, atlas_i), 
                to = file.path(getwd(), "data/asai_10km_mapping",  atlas_i))
      
    } else {
      # Check NBN 10KM mapping
      nbn_i <- character(0)
      file_path <- "BIOLOGICAL_RECORD_DATA/nbn_at_10km_mapped"
      pattern_i <- gsub("([][{}()+*^$.|\\\\?-])", "\\\\\\1", SCIENTIFIC_NAME)
      nbn_i <- list.files(file.path(shareDir, file_path), pattern = pattern_i)

      if (length(nbn_i) > 0) {
        # Check off tif
        asai_info$mapping_10km_complete[i] <- TRUE

        # Mark off mapping type
        asai_info$mapping_type[i] <- "NBN"

        # File path 
        asai_info$file_path[i] <- file.path(shareDir, file_path, nbn_i)
        
        # Copy tif to project folder
        file.copy(from = file.path(shareDir, file_path, nbn_i), 
                  to = file.path(getwd(), "data/asai_10km_mapping", nbn_i))
        
      } else {
        # Check off tif
        asai_info$mapping_10km_complete[i] <- FALSE
      }
    }
  }
}

missing_mapping <- asai_info[which( !asai_info$mapping_10km_complete),]
missing_mapping # All species complete

# === Remove duplicated rows (fix in original ebrel database script) ===
asai_info_compiled <- asai_info[!duplicated(asai_info$SCIENTIFIC_NAME),]
length(unique(asai_info_compiled$SCIENTIFIC_NAME))
nrow(asai_info_compiled)

# === Normalise names ===
normalize_taxon <- function(x) {
  x |>
    tolower() |>
    gsub("_", " ", x = _) |>
    gsub("[[:punct:]]+", " ", x = _) |>
    gsub("\\s+", " ", x = _) |>
    trimws()
}
asai_info_compiled$SCIENTIFIC_NAME_key <- normalize_taxon(asai_info_compiled$SCIENTIFIC_NAME)

# === Normalise BDI_NAME_CHECK ===
asai_info_compiled$BDI_NAME_CHECK <- normalize_taxon(asai_info_compiled$BDI_NAME_CHECK)

# === Save this subset ===
fwrite(asai_info_compiled, file = "outputs/bdi_asai_compiled_species_information.csv")

#==============================================================================#
# Create species x 10km OS grid cell matrix
#==============================================================================#

# === Create cell x species data.frame ===

# == Re-load compiled data
asai_info_compiled <- fread("outputs/bdi_asai_compiled_species_information.csv")

# == asai species list
f_asai <- list.files("data/asai_10km_mapping", full.names = FALSE)

# == Compile cell by species
asai_cxs <- lapply(1:nrow(asai_info_compiled), function(i){
  print(i)
  # Species name in the asai 
  spp_name_asai <- gsub(" ", "_", asai_info_compiled$BDI_NAME_CHECK[i])
  
  # Species name tif
  spp_name_file_path <- asai_info_compiled$file_path[i]
  spp_name_tif <- strsplit(spp_name_file_path, split = "/" )[[1]][6]
  
  # Load species tif
  tif <- rast(file.path("data/asai_10km_mapping", spp_name_tif))

  # Extract to dataframe with row number ids
  df <- terra::as.data.frame(tif, cell = TRUE, na.rm = TRUE)

  if (sum(df[,2]) == 0) {
    print(i)
    print("Zero")
  }
  
  # Name
  names(df)[2] <- spp_name_asai

  # Return
  return(df)

})
asai_cxs[sapply(asai_cxs, is.null)] <- NULL
asai_cxs <- Reduce(function(...) merge(..., by = "cell", all = TRUE), asai_cxs)

# === Remove Zero columns ===
zero_col <- colSums(asai_cxs, na.rm = TRUE) == 0
asai_cxs <- asai_cxs[,!zero_col]
dim(asai_cxs)
sum(asai_cxs$dryobates_minor, na.rm = TRUE)

# === Add in os_gridref ===
asai_cxs <- merge(asai_cxs, cell_to_os_gridref_lkup, by = c("cell"), all.x = TRUE)
asai_cxs <- asai_cxs[, c("cell", "os_gridref",sort(setdiff(names(asai_cxs), c("cell", "os_gridref"))))]

# === Replace NAs ===
asai_cxs[is.na(asai_cxs)] <- 0
asai_cxs[1:5,1:5]
dim(asai_cxs)

# === Save ===
fwrite(asai_cxs, file = "outputs/asai_cxs_gb.csv")

# ====
# === Process to England
# ====

# Just want England distributions for this analysis

# === Load cell x species ===
asai_cxs_gb <- fread("outputs/asai_cxs_gb.csv")
asai_cxs_gb[1:5,1:5]

# === Subset ASAI cells to England ===
asai_cxs_eng <- asai_cxs_gb[asai_cxs_gb$os_gridref %in% unique(eng10km$os_gridref), ]
asai_cxs_eng[, cell := NULL]
asai_cxs_eng[1:5,1:5]
dim(asai_cxs_eng)

# === Save ===
fwrite(asai_cxs_eng, file = "outputs/asai_cxs_eng.csv")


#==============================================================================#
# Extract LU/LC area by type by 10km
#==============================================================================#

# === Load species trait information ===
asai_cxs_eng <- fread("outputs/asai_cxs_eng.csv")
asai_info_compiled <- fread( "outputs/bdi_asai_compiled_species_information.csv")

# === Extract coverage
he10km_i <- exact_extract(habExist100m, eng10km2, coverage_area  = TRUE)
id_col <- "os_gridref"
hab_cov <- rbindlist(lapply(seq_along(he10km_i), function(i) {
  x <- as.data.table(he10km_i[[i]])
  x[, (id_col) := eng10km2[[id_col]][i]]
  x
}), use.names = TRUE, fill = TRUE)
hab_cov <- hab_cov[!is.na(value),
         .(area_m2 = sum(coverage_area, na.rm = TRUE)),
         by = c(id_col, "value")]

# === Calculate habitat proportions per grid cell ===

# = Cell areas in km2
cell_area <- setNames(as.numeric(st_area(eng10km2)), eng10km2[[id_col]])

# = Full cell area
A0 <- max(cell_area)   

# = Proportion relative to complete cell so not to overweight fractional cells
hab_cov[, prop := area_m2 / A0] 

# = Clamp here just in case value exceed 1 (they shouldn't)
hab_cov[prop > 1, prop := 1]   

# = Convert to wide format with raw proportions
LU_by_C <- dcast(hab_cov, as.formula(paste(id_col, "~ value")),
                 value.var = "prop", fill = 0)

# = Set habitat column names
hab_cols <- setdiff(names(LU_by_C), id_col)
setnames(LU_by_C, old = setdiff(names(LU_by_C), id_col),
         new = habitat_types_lookup$habitat)

#### Plot to check ######
LU_sf <- left_join(GB10km, LU_by_C)
ggplot() + 
  geom_sf(aes(fill = heathland, colour = heathland), data = LU_sf) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c()

# === Species x habitat (aligning habitat classes) ===

# = Split narrow classes into list
`%||%` <- function(x, y) if (is.null(x)) y else x
HT <- setDT(habitat_types_lookup)
HT[, narrow := strsplit(narrow_classes %||% "", "\\|")]  
HT[is.na(narrow_classes), narrow := list(list(character()))]

# = Ensure all referenced narrow columns exist in DT (create as 0 if missing)
all_narrow <- unique(unlist(HT$narrow))
all_narrow <- all_narrow[nzchar(all_narrow)]
missing <- setdiff(all_narrow, names(asai_info_compiled))
if (length(missing)) speciesHabitatTraits[, (missing) := 0L]

# = Build condensed columns (ANY of the narrow columns = 1)
out_names <- make.names(HT$habitat) 
for (i in seq_len(nrow(HT))) {
  sub <- HT$narrow[[i]]
  sub <- sub[nzchar(sub)]
  if (length(sub) == 0) {
    asai_info_compiled[, (out_names[i]) := 0L]
  } else {
    asai_info_compiled[, (out_names[i]) := as.integer(rowSums(.SD > 0) > 0), .SDcols = sub]
  }
}

# = Select only required columns
hab_cols <- c("deciduous_woodland","wood_pasture","wet_woodland",
              "neutral_grassland","calcareous_grassland","acid_grassland",
              "heathland","purple_moor_rush_pasture","peatland","wetland",
              "sand_dunes","intertidal","arable","improved_grassland",
              "coniferous_woodland","urban")
asai_info_aligned <- asai_info_compiled %>% 
  select("BDI_NAME_CHECK", hab_cols)

# ====
# === For each species, calculate proportion of cell with suitable habitat
# ====

# === Calculate proportions based on species habitat requirements ===
suit_sc_empty <- data.frame(os_gridref = LU_by_C$os_gridref, prop_suit = NA)
species_habitat_by_cell <- list()
species_asai <- names(asai_cxs_eng)[-c(1)]

for (i in seq_along(species_asai)) {
  # Species
  species_name <-  species_asai[i]
  species_name_str <- gsub("_", " ", species_name)

  # Species' habitat info
  hab_info_sp <- asai_info_aligned[asai_info_aligned$BDI_NAME_CHECK == species_name_str]
  hab_info_sp <- names(hab_info_sp)[which(hab_info_sp > 0)][-1]
  
  if (length(hab_info_sp) == 0) {
  
    # Get binomials for both asai_info_aligned$SCIENTIFIC_NAME and species_name
    db_names <-  sapply(asai_info_aligned$SCIENTIFIC_NAME_key, function(x){
      paste(strsplit(x, split = " ")[[1]][1:2], collapse = " ")
    })
 
     # Remove common infra-rank markers and everything after them
    species_name_binom <-  paste(strsplit(species_name, split = " ")[[1]][1:2], collapse = " ")
      
    hab_info_sp <- asai_info_aligned[db_names == species_name_binom]
    hab_info_sp <- names(hab_info_sp)[which(hab_info_sp > 0)]
    
  }
  
  if (length(hab_info_sp) == 0) {
    message("Species - ", species_name, ": i = ", i)
  }

  # Cell-wise proportion coverage 
  suit_sc <- suit_sc_empty
  for(c in 1:nrow(suit_sc)) {
    prop_suit <- sum(LU_by_C[c, names(LU_by_C) %in% hab_info_sp, with = FALSE])
    # Assuming that all cells will have some suitable habitat if the species is
    # thought to occur there and also unavailable habitat (e.g. roads etc) we 
    # rescale onto [0.05,0.95].
    suit_sc[c, "prop_suit"] <-  pmin(pmax(0.05, prop_suit), 0.95)
  }
  names(suit_sc)[2] <- species_name
  species_habitat_by_cell[[i]] <- suit_sc
  
}
prop_suit_cs <- Reduce(function(...) merge(..., by = "os_gridref", all = TRUE),
                       species_habitat_by_cell)

# === Save ===
fwrite(prop_suit_cs, file = "outputs/prop_suit_cs.csv")

#==============================================================================#
# Prepare model data
#==============================================================================#

# === Species by cell - occupancy ===
asai_cxs_eng <- fread("outputs/asai_cxs_eng.csv")
names(asai_cxs_eng) <- gsub(" ", "_", names(asai_cxs_eng))
asai_cxs_eng[1:5,1:5]
dim(asai_cxs_eng)

# === Remove No occurrences ===
cols0 <- names(asai_cxs_eng)[-1][colSums(asai_cxs_eng[, -1]) == 0]
asai_cxs_eng[, (cols0) := NULL]

# === Species by cell - suitable habitat area ===
prop_suit_cs <- fread("outputs/prop_suit_cs.csv")
names(prop_suit_cs) <- gsub(" ", "_", names(prop_suit_cs))
prop_suit_cs[1:5,1:5]
dim(prop_suit_cs)

# ====
# === Prepare optimisation data
# ====

# === Ensure cells included have habitat coverage ===
asai_cxs_eng <- asai_cxs_eng[asai_cxs_eng$os_gridref %in% prop_suit_cs$os_gridref, ]

# === Ensure the exact same species and order ===
asai_cxs_eng <- asai_cxs_eng[,names(asai_cxs_eng) %in% names(prop_suit_cs), with = FALSE]
prop_suit_cs <- prop_suit_cs[,names(prop_suit_cs) %in% names(asai_cxs_eng), with = FALSE]
table(names(asai_cxs_eng) == names(prop_suit_cs))

# === Match cell orders ===
asai_cxs_eng <- asai_cxs_eng[order(asai_cxs_eng$os_gridref),]
prop_suit_cs <- prop_suit_cs[order(prop_suit_cs$os_gridref),]

# = Check dims =
dim(asai_cxs_eng) == dim(prop_suit_cs)

# = Check data.tables =
asai_cxs_eng[1:5,1:5]
prop_suit_cs[1:5,1:5]
names(asai_cxs_eng)

# === Convert to matrix P[c, s] and H[c, s] ===

# = P[c, s]
P <- as.matrix(asai_cxs_eng[,2:ncol(asai_cxs_eng)])

# = H[c, s]
H <- as.matrix(prop_suit_cs[,2:ncol(prop_suit_cs)])

# = Check dims = 
dim(H) == dim(P)

# === Extract dimensions ===
nC <- nrow(P)  # number of 10 km cells
nS <- ncol(P)  # number of species


# === Calculate habitat-weighted presence (unnormalised) ===
# i.e., proportion of suitable habitat for each species where it occurs
PH <- P * H

# === Calculate the total area of suitable habitat across occupied cells ===
PD_dt <- cbind(asai_cxs_eng[,1], PH)
setDT(cellXarea)
species_cols <- setdiff(names(PD_dt), "os_gridref")
HArea <- cellXarea[PD_dt, on = "os_gridref"]
range_size <- data.frame(
  range_size = colSums(P * sweep(as.matrix(HArea[, ..species_cols]), 1, HArea$area, `*`))
)
fwrite(range_size, "outputs/range_sizes.csv")

cbind(species_cols,range_size)

summary(range_size$range_size)
7.255e+03

# === Produce species range histogram ===
species_habitat_range <- ggplot(data = range_size, aes(x = range_size)) +
  geom_histogram(
    bins = 100,
    fill = "grey80",
    colour = "black",
    linewidth = 0.3
  ) +
  labs(
    x = "Area of suitable habitat across occupied\nOS 10 km² grid cells (km²)",
    y = "Number of species"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    #axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    panel.grid = element_blank(),
    plot.margin = margin(10, 25, 10, 10)
  )
ggsave(filename = "outputs/species_range_habitat_weighted.png", 
       plot = species_habitat_range,
       width = 10,
       height = 6)



# === Denominator: total suitable habitat across occupied cells ===
R_s_weighted <- colSums(PH)

# = Check zeros
which(R_s_weighted == 0)

# = Normalised weights: each column (species) sums to 1
P_over_R <- sweep(PH, 2, R_s_weighted, "/")

# = Replace NA columns with zeros if needed (optional but often safer)
P_over_R[is.na(P_over_R)] <- 0

#### Plot to check ####
# sol_os <- cbind(asai_cxs_eng[,1], selected=P_over_R[,"passer_montanus"])
# sol_sf <- left_join(GB10km, sol_os)
# ggplot() + geom_sf(aes(fill = selected, colour = selected), data = sol_sf)
# sum(P_over_R[,"passer_montanus"])

# === Save all prepared data ===
optimisation_data <- list(
  nC = nC,
  nS = nS,
  asai_cxs_eng = asai_cxs_eng,
  prop_suit_cs = prop_suit_cs,
  P_over_R = P_over_R,
  PH = PH
)
saveRDS(optimisation_data, file = "outputs/optimisation_data.rds")

#==============================================================================#
# Analysis - optimisation
#==============================================================================#

# === Setup parellelisation ===
#plan(multisession, workers = 18)# max(1, parallel::detectCores() - 2))

# === Load optimisation data ===
optimisation_data <- readRDS("outputs/optimisation_data.rds")

# === os_gridref_inc ===
os_gridref_inc <- optimisation_data$asai_cxs_eng$os_gridref

# === GM targets ===
g    <- 0.10  # target +10% increase in the indicator
beta <- 0.5   # max proportional increase if full range treated 
GM_target <- log(1 + g)

# === K sequence search ===
K_seq <- seq(0, optimisation_data$nS, by = 1)

# === Costs with boundary adjustment ===
# Need boundary adjustments on costs as otherwise partial cells around the
# coast and borders and overly represented
cost_bound_adj <- make_boundary_adjusted_cost(
  cell_area = area_km2,
  full_cell_area = max(area_km2),
  eta = 1
)

# ====
# === Solve Pareto front for incrementally higher K = species +10%
# ====

pareto_front_on_K <- solve_pareto_front_incremental(
  P_over_R = optimisation_data$P_over_R,
  beta = beta,
  GM_target = GM_target,
  K_seq = K_seq,
  cost = cost_bound_adj,
  randomise_within_cost_prune = TRUE,
  top_species = 25,
  beam_width = 8,
  max_depth = 5,
  top_expand = 8,
  overshoot_lambda = 0.25,
  verbose = TRUE
)

# === Save ===
saveRDS(pareto_front_on_K, "outputs/pareto_front_on_K_0_to_1174_by_1_beta1_costEta1.rds")


#==============================================================================#
# Analysis - results
#==============================================================================#

# === Major results ===
#
# Scenarios:
# - 1) GM Constrained only
# - 2) GM Constrained + varying species % constraints
#
# 1) Pareto front Footprint vs. Species % +10%
#
#    Summary statistics:
#    * footprint at GM-constraint only, N85, N100
#
# 2) Spatial map showing selection across pareto front
#
#    Summary statistics:
#    * Spatial analysis of regional representation, habitat representation 
#
# 3) Figure showing range size of species vs proportion of range selected across pareto front
#
#    Summary statistics: 
#    * Summary statistics at GM-constraint only, N85, N100
#
# 4) Figure showing representation by taxonomic group across pareto front
#
#    Summary statistics: 
#    * Summary statistics at GM-constraint only, N85, N100

# ====
# === Optimisation data
# ====

# Load data
optimisation_data <-readRDS("outputs/optimisation_data.rds")

# Number of cells
nC <- optimisation_data$nC

# Number of species
nS <- optimisation_data$nS

# Species presence absence
asai_cxs_eng <- optimisation_data$asai_cxs_eng

# Proportion of cell suitable for each species
prop_suit_cs <- optimisation_data$prop_suit_cs

# Normalised weights
P_over_R <- optimisation_data$P_over_R

# Species names
species_names <- gsub("_", " ", names(asai_cxs_eng)[2:ncol(asai_cxs_eng)])

# Area
area_km2 <- cellXarea$area[which(cellXarea$os_gridref %in% asai_cxs_eng$os_gridref)]

# ====
# === Build pareto front (beta = 1 and eta = 1)
# ====

# Load optimisation results
pareto_front_on_K_b1_e1 <- readRDS("outputs/pareto_front_on_K_0_to_1174_by_1_beta1_costEta1.rds")

# Build plot
pareto_plot_b1_e1 <- build_pareto_front_plot(pareto_front_on_K_b1_e1)

# = Save
ggsave(filename = "outputs/pareto_front_panel_beta1.png", 
       plot = pareto_plot_b1_e1[[1]], width = 8, height = 6)

# === Extract summary statistics for N10s

# == Unconstrained 
min_footprint_beta1 <- sum(pareto_front_on_K_b1_e1[[1]]$best_state$x * area_km2)
min_footprint_beta1 
# 6697.649km2
(min_footprint_beta1  / total_area) * 100
# 5.0385%

# Percentage of species min footprint
pareto_front_on_K_b1_e1[[1]]$best_state$N10 / nS
#  0.4795571%

# Percentage of species +5% or more
length(which(pareto_front_on_K_b1_e1[[1]]$best_state$m_s >= 1.05)) / nS

# == N10 == 100%
full_K_solution <- pareto_front_on_K_b1_e1[[length(pareto_front_on_K_b1_e1)]]
max_N10_footprint_beta1 <- sum(full_K_solution$best_state$x * area_km2)
max_N10_footprint_beta1 
# 13828.09 km2
(max_N10_footprint_beta1  / total_area) * 100
# 10.40254%


# == Marginal gain as the additional species reaching +10% per extra area
# species_per_100km2_effective_footprint
deriv_df <- pareto_plot_b1_e1[[2]] |>
  mutate(
    species_per_km2 = abs(derivative) * nS / total_area,
    species_per_100km2 = species_per_km2 * 100,
    species_per_10km_cell = species_per_100km2
  )

# Across the entire pareto front
deriv_summary <- deriv_df |>
  summarise(
    mean_species_per_100km2 = mean(species_per_100km2, na.rm = TRUE),
    median_species_per_100km2 = median(species_per_100km2, na.rm = TRUE),
    min_species_per_100km2 = min(species_per_100km2, na.rm = TRUE),
    max_species_per_100km2 = max(species_per_100km2, na.rm = TRUE)
  )

# Up to the point where 85% of species meet targets.
deriv_summary_85 <- deriv_df |>
  filter(1 - y <= 0.85) |>
  summarise(
    mean_species_per_100km2 = mean(species_per_100km2, na.rm = TRUE),
    median_species_per_100km2 = median(species_per_100km2, na.rm = TRUE),
    min_species_per_100km2 = min(species_per_100km2, na.rm = TRUE),
    max_species_per_100km2 = max(species_per_100km2, na.rm = TRUE)
  )

# ====
# === 2) Plot minimal footprints at different N10
# ====

# We want to show how the spatial pattern of the footprint changes as you enforce
# N10 targets. To do this, plot the minimal spatial footprint (unconstrained) 
# and the footprint for N10 = 1174, and also a heat map of the frequency of selection.
# Also, add the 85%

# === For plotting, crop GB map
gb_eng_crop <- st_crop(gb, eng10km2)
eng_only <- gb[gb$CTRY22NM == "England",]

# === Biogeographic zones
bio_geo_zone_bd <- st_read("data/Biogeographic_Zones_Living_England_2021.gpkg")
bio_geo_zone <- st_intersection(eng10km2, bio_geo_zone_bd)
plot(bio_geo_zone[,3])

# === Load spatial protected area data
NP <- st_read("data/National_Parks_England.gpkg")
sum(NP$MEASURE)

# === Load spatial protected area data
SSSI <- st_read("data/SSSI_England.gpkg")

# === Load spatial protected area data
NL <- st_read("data/Areas_of_Outstanding_Natural_Beauty_England.gpkg")

# ====
# === Biogeographic representation across solution space
# ====

# === Zone area
zone_area_df <- bio_geo_zone %>%
  st_drop_geometry() %>%
  group_by(Zone) %>%
  summarise(
    zone_area_km2 = sum(area, na.rm = TRUE) / 1e6,
    .groups = "drop"
  )

# === Extract select by biogeographic region
footprint_per_zone <- lapply(c(seq(1, 1175, 10), 1175), function(i) {
  
  print(i)
  
  # Get temp x
  tmp <- pareto_front_on_K_b1_e1[[i]]$best_state$x
  
  # Map to actual area because best assumes equal cost for partial cells
  best_area_actual_size <- sum(tmp * area_km2, na.rm = TRUE)
  
  # Intersect with biogeographic zones
  tmp <- cbind(asai_cxs_eng[, 1], footprint = tmp)
  
  tmp_sf <- left_join(bio_geo_zone, tmp, by = join_by(os_gridref)) %>% 
    st_drop_geometry()
  
  # Summarise selected footprint in biogeographic zones
  bioGeoSum <- tmp_sf %>% 
    group_by(Zone) %>% 
    summarise(
      footprint_in_zone = sum(area[footprint == 1], na.rm = TRUE) / 1e6,
      .groups = "drop"
    ) %>%
    right_join(zone_area_df, by = "Zone") %>%
    mutate(
      footprint_in_zone = replace_na(footprint_in_zone, 0),
      
      k = pareto_front_on_K_b1_e1[[i]]$K,
      N10 = pareto_front_on_K_b1_e1[[i]]$best_state$N10,
      best_footprint_area_km = best_area_actual_size,
      
      total_zone_area_km2 = sum(zone_area_km2, na.rm = TRUE),
      
      expected_footprint_in_zone =
        best_footprint_area_km * zone_area_km2 / total_zone_area_km2,
      
      observed_prop =
        footprint_in_zone / sum(footprint_in_zone, na.rm = TRUE),
      
      expected_prop =
        zone_area_km2 / total_zone_area_km2,
      
      observed_expected_ratio =
        observed_prop / expected_prop,
      
      footprint_difference =
        footprint_in_zone - expected_footprint_in_zone,
      
      footprint_prop_difference =
        observed_prop - expected_prop
    )
  
  bioGeoSum
})
footprint_per_zone_df <- do.call(rbind, footprint_per_zone)

# === Plot of footprint per biogeographic zone by N10

# = Want the same colours on both panels
zone_levels_chr <- sort(unique(as.numeric(c(
  as.character(footprint_per_zone_df$Zone),
  as.character(bio_geo_zone_bd$Zone))))) %>% 
  as.character()
zone_cols <- setNames(
  hcl.colors(length(zone_levels_chr), palette = "Dark 3"),
  zone_levels_chr
)
footprint_plot_dat <- footprint_per_zone_df %>%
  mutate(
    Zone = factor(as.character(Zone), levels = zone_levels_chr),
    observed_expected_ratio = observed_expected_ratio
    )
map_dat <- bio_geo_zone_bd %>%
  mutate(
    Zone = factor(as.character(Zone), levels = zone_levels_chr)
  )

# = And, numeric region labels
zone_labels <- data.frame(
  Zone = zone_levels_chr,
  zone_lab = zone_levels_chr
)
# ... Label positions for inset map
zone_label_pts <- bio_geo_zone_bd %>%
  mutate(
    Zone = factor(as.character(Zone), levels = zone_levels_chr)
  ) %>%
  st_point_on_surface() %>%
  left_join(zone_labels, by = c("Zone" = "Zone"))

# = Footprint by N10
footprint_plot_dat <- footprint_per_zone_df %>%
  mutate(
    Zone = factor(
      as.character(Zone),
      levels = zone_levels_chr
    ),
    N10_pc = (N10 / nS) * 100
  )

# = Smoother
smooth_dat <- footprint_plot_dat %>%
  group_by(Zone) %>%
  group_modify(~ {
    
    mod <- gam(
      observed_expected_ratio ~ s(N10_pc),
      data = .x
    )
    
    x_new <- seq(
      min(.x$N10_pc, na.rm = TRUE),
      max(.x$N10_pc, na.rm = TRUE),
      length.out = 200
    )
    
    data.frame(
      N10_pc = x_new,
      observed_expected_ratio = predict(
        mod,
        newdata = data.frame(N10_pc = x_new)
      )
    )
  }) %>%
  ungroup()

# = End points for smoothed lines
zone_end_labels <- smooth_dat %>%
  group_by(Zone) %>%
  slice_max(N10_pc, n = 1, with_ties = FALSE) %>%
  ungroup()

# = Plot smooths
footprint_per_zone_by_N10 <- ggplot() +
  geom_point(
    data = footprint_plot_dat,
    aes(
      x = N10_pc,
      y = observed_expected_ratio,
      colour = Zone
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = smooth_dat,
    aes(
      x = N10_pc,
      y = observed_expected_ratio,
      colour = Zone,
      group = Zone
    ),
    linewidth = 1.5
  ) +
  ggrepel::geom_text_repel(
    data = zone_end_labels,
    aes(
      x = N10_pc,
      y = observed_expected_ratio,
      label = Zone,
      colour = Zone
    ),
    direction = "y",
    hjust = 0,
    nudge_x = 2,
    segment.colour = "grey60",
    box.padding = 0.2,
    point.padding = 0.2,
    min.segment.length = 0,
    seed = 1,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = zone_cols,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    plot.margin = margin(
      t = 10,
      r = 35,
      b = 10,
      l = 10
    )
  ) +
  labs(
    x = "Species reaching +10% abundance target (%)",
    y = "Observed / expected footprint\n(values > 1 = greater than expected)"
  )

# = Insert map with region number labels
bio_geo_zone_map_inset <- ggplot() +
  geom_sf(data = gb_eng_crop, fill = "grey95", colour = NA) +
  geom_sf(data = eng_only, fill = "grey80", colour = NA) +
  geom_sf(
    data = bio_geo_zone_bd %>%
      mutate(
        Zone = factor(
          as.character(Zone),
          levels = zone_levels_chr
        )
      ),
    aes(fill = Zone),
    colour = "grey30",
    linewidth = 0.15
  ) +
  geom_sf_text(
    data = zone_label_pts,
    aes(label = Zone),
    size = 3,
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_manual(values = zone_cols, drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(
      fill = "white",
      colour = "grey30",
      linewidth = 0.3
    ),
    plot.margin = margin(2, 2, 2, 2)
  )


# = Graph with map inset in top right
zones_panel <- footprint_per_zone_by_N10 +
  inset_element(
    bio_geo_zone_map_inset,
    left = 0.56,
    bottom = 0.54,
    right = 1.13,
    top = 0.985,
    align_to = "panel",
    clip = TRUE
  ) #+plot_layout(tag_level = "keep")


# = Save
ggsave(filename = "outputs/footprint_by_zone.png", 
       plot = zones_panel, width = 9, height = 6)

# ====
# === Selection frequency
# ====

# = Selection frequency
selection_freq <- lapply(1:length(pareto_front_on_K_b1_e1), function(x) {
  tmp <- pareto_front_on_K_b1_e1[[x]]$best_state$x
})
selection_freq <- rowSums(as.data.frame(do.call(cbind, selection_freq)))
selection_freq <- cbind(asai_cxs_eng[, 1], selection_freq = selection_freq)
selection_freq_sf <- left_join(eng10km2, selection_freq, by = join_by(os_gridref)) %>% 
  mutate(selection_freq = ifelse(selection_freq == 0, NA, (selection_freq / 1176) * 100))
# Selection frequency figure
cols <- c("#1A0500", "#5C0F00", "#8B1A1A", "#C5301F", "#E8502A", "#F59E00", "#FFE800")

# = Plot selection frequency
cell_selection_frequency <- ggplot(selection_freq_sf) +
  geom_sf(data = gb_eng_crop, fill = "grey95") +
  geom_sf(aes(fill = selection_freq, colour = selection_freq)) +
  scale_fill_gradientn(
    colours  = cols,
    values  = scales::rescale(c(0, 5, 20, 60, 100)),
    limits   = c(0, 100),
    na.value = "grey80",
    name     = "Selected (%)",
    guide    = guide_colourbar(
      barwidth     = 0.8,
      barheight    = 8,
      ticks.colour = "grey40"
    )
  ) +
  scale_colour_gradientn(
    colours  = cols,
    values  = scales::rescale(c(0, 5, 20, 60, 100)),
    limits   = c(0, 100),
    na.value = "grey80",
    name     = "Selected (%)",
    guide    = guide_colourbar(
      barwidth     = 0.8,
      barheight    = 8,
      ticks.colour = "grey40"
    )
  ) +
  theme_minimal(base_size = 16) +      
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.15, 0.55),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 10),
    plot.background   = element_rect(fill = NA, colour = NA),
    plot.margin       = margin(10, 10, 10, 10)
  ) +
  coord_sf(
    xlim = st_bbox(gb_eng_crop)[c("xmin", "xmax")],
    ylim = st_bbox(gb_eng_crop)[c("ymin", "ymax")],
    expand = FALSE
  )

# = Number selection 100%
select_freq <- selection_freq$selection_freq / nS
length(select_freq[select_freq == 1])
length(select_freq[select_freq > 0])

# = Combined panel

# Aspect ratio of the map extent
bb <- st_bbox(eng10km2)
map_aspect <- as.numeric(
  (bb["xmax"] - bb["xmin"]) /
    (bb["ymax"] - bb["ymin"])
)
zones_width <- 9
map_width   <- 6 * map_aspect
zones_panel <- zones_panel +
  labs(tag = "a)")

# Treat the complete plot + inset as ONE panel (necessary for annotations)
zones_panel_wrapped <- wrap_elements(
  full = zones_panel
)

# Create panel
Zonal_Select_Panel <-
  zones_panel_wrapped +
  cell_selection_frequency +
  plot_layout(
    widths = c(zones_width, map_width)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0, 1)
  )

# Save at exactly the dimensions implied by the panels
ggsave(
  "outputs/Zonal_Select_Panel.png",
  Zonal_Select_Panel,
  width = zones_width + map_width,
  height = 6
)

# === 
# === Uncontrained (K_target = 0)
# ===

# = Get the K_target = 0 solution
K0_solution <- pareto_front_on_K_b1_e1[[1]]$best_state$x

# = Convert to sf
unconst_footprint <- cbind(asai_cxs_eng[, 1], footprint = K0_solution)
unconst_footprint_sf <- left_join(eng10km2, unconst_footprint, by = join_by(os_gridref)) %>%
  mutate(footprint = factor(footprint)) %>%
  filter(!is.na(footprint))

# = Create map of selected cells
unconst_footprint_sol <- ggplot(unconst_footprint_sf[unconst_footprint_sf$footprint == 1, ]) +
  geom_sf(data = gb_eng_crop, fill = "grey95") +
  geom_sf(data = eng_only, fill = "grey80") +
  #geom_sf(data = bio_geo_zone, colour = "grey20", fill = NA) +
  geom_sf(data = NP, colour = NA, fill = "forestgreen", alpha = 0.2) +
  geom_sf(data = NL, colour = NA, fill = "blue", alpha = 0.2) +
  geom_sf(data = SSSI, colour = NA, fill = "orange", alpha = 0.5) +
  geom_sf(aes(fill = footprint, colour = footprint)) +
  scale_fill_manual(values = c("1" = "red", "0" = "grey80")) +
  scale_colour_manual(values = c("1" = "red", "0" = "grey80")) +
  annotation_scale(
    location   = "br",
    width_hint = 0.25,
    line_width = 2,
    text_cex = 1.2,
    style      = "ticks",
    line_col   = "grey10",
    text_col   = "grey10"
  ) +
  annotation_north_arrow(
    location = "tr",
    pad_y    = unit(0.5, "cm"),
    style    = north_arrow_orienteering(
      line_col = "grey10",
      text_col = "grey10",
      text_size = 12
    )
  ) +
  theme_minimal(base_size = 14) +      
  theme(legend.position = "none") +
  coord_sf(
    xlim = st_bbox(gb_eng_crop)[c("xmin", "xmax")],
    ylim = st_bbox(gb_eng_crop)[c("ymin", "ymax")],
    expand = FALSE
  )

# === 
# === N10 = 85% 
# ===

# === Which is the best 85% solution
which_best_85pc <- rbindlist(lapply(1:length(pareto_front_on_K_b1_e1), function(x){
  tmp <- pareto_front_on_K_b1_e1[[x]]
  data.frame(x = x, N10 = tmp$best_state$N10, best_area = tmp$best_area)
}))

# How many species is 85%
nSp85 <- ceiling(nS*0.85)

# Best index
best_85pc_ind <- which_best_85pc[which_best_85pc$N10 == nSp85,]

# === Get the N10 = 85% solution
N10_85_solution <- pareto_front_on_K_b1_e1[[best_85pc_ind$x]]$best_state$x

# = Convert to sf
pc85_footprint <- cbind(asai_cxs_eng[, 1], footprint = N10_85_solution)
pc85_footprint_sf <- left_join(eng10km2, pc85_footprint, by = join_by(os_gridref)) %>%
  mutate(footprint = factor(footprint)) %>%
  filter(!is.na(footprint))

# = Create map of selected cells
pc85_footprint_sol <- ggplot(pc85_footprint_sf[pc85_footprint_sf$footprint == 1, ]) +
  geom_sf(data = gb_eng_crop, fill = "grey95") +
  geom_sf(data = eng_only, fill = "grey80") +
  #geom_sf(data = bio_geo_zone, colour = "grey20", fill = NA) +
  geom_sf(data = NP, colour = NA, fill = "forestgreen", alpha = 0.2) +
  geom_sf(data = NL, colour = NA, fill = "blue", alpha = 0.2) +
  geom_sf(data = SSSI, colour = NA, fill = "orange", alpha = 0.5) +
  geom_sf(aes(fill = footprint, colour = footprint)) +
  scale_fill_manual(values = c("1" = "red", "0" = "grey80")) +
  scale_colour_manual(values = c("1" = "red", "0" = "grey80")) +
  theme_minimal(base_size = 14) +      
  theme(legend.position = "none") +
  coord_sf(
    xlim = st_bbox(gb_eng_crop)[c("xmin", "xmax")],
    ylim = st_bbox(gb_eng_crop)[c("ymin", "ymax")],
    expand = FALSE
  )

# === 
# === N10 = 1174
# ===

# = Get the K_target = 0 solution
N10_100_solution <- pareto_front_on_K_b1_e1[[length(pareto_front_on_K_b1_e1)]]$best_state$x

# = Convert to sf
pc100_footprint <- cbind(asai_cxs_eng[, 1], footprint = N10_100_solution)
pc100_footprint_sf <- left_join(eng10km2, pc100_footprint, by = join_by(os_gridref)) %>%
  mutate(footprint = factor(footprint)) %>%
  filter(!is.na(footprint))

# Map
pc100_footprint_sol <- ggplot(pc100_footprint_sf[pc100_footprint_sf$footprint == 1, ]) +
  geom_sf(data = gb_eng_crop, fill = "grey95") +
  geom_sf(data = eng_only, fill = "grey80") +
  #geom_sf(data = bio_geo_zone, colour = "grey20", fill = NA) +
  geom_sf(data = NP, colour = NA, fill = "forestgreen", alpha = 0.2) +
  geom_sf(data = NL, colour = NA, fill = "blue", alpha = 0.2) +
  geom_sf(data = SSSI, colour = NA, fill = "orange", alpha = 0.5) +
  geom_sf(aes(fill = footprint, colour = footprint)) +
  scale_fill_manual(values = c("1" = "red", "0" = "grey80")) +
  theme_minimal(base_size = 14) +      
  theme(legend.position = "none") +
  coord_sf(
    xlim = st_bbox(gb_eng_crop)[c("xmin", "xmax")],
    ylim = st_bbox(gb_eng_crop)[c("ymin", "ymax")],
    expand = FALSE
  )


# ====
# === Create solution mapping panel
# ====

# = Panel with min, max, and selection freq panels
footprint_panel <- 
  unconst_footprint_sol + pc85_footprint_sol +  pc100_footprint_sol  +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0, 1)
  )

# = Save
ggsave(filename = "outputs/footprint_panel_beta1_beta1_eta1_v3.png", 
       plot = footprint_panel, width = 15, height = 5)

  
# ====
# === Habitat representation
# ===

# == Ebrel version matching habitat_types_lookup
habExist100m <- rast(file.path(shStdDat, "ebrel_landcover_100m_england_wales.tif"))

# === Extract coverage
he10km_i <- exact_extract(habExist100m, eng10km2, coverage_area  = TRUE)
id_col <- "os_gridref"
hab_cov <- rbindlist(lapply(seq_along(he10km_i), function(i) {
  x <- as.data.table(he10km_i[[i]])
  x[, (id_col) := eng10km2[[id_col]][i]]
  x
}), use.names = TRUE, fill = TRUE)
hab_cov <- hab_cov[!is.na(value),
                   .(area_m2 = sum(coverage_area, na.rm = TRUE)),
                   by = c(id_col, "value")]

# = Summarise habitat by selected cells

# Convert to dt
habitat_types_lookup_dt <- as.data.table(habitat_types_lookup)

# Marcus prefers grass and rush pasture
habitat_types_lookup_dt$habitat[8] <- "Grass and rush pasture"

# Calculate expected habitat area across all available cells
hab_expected <- hab_cov[,.(expected_area = sum(area_m2, na.rm = TRUE)), by = value]

# Calculate observed habitat area for smallest / unconstrained footprint
selected_cells_s <- unconst_footprint[footprint == 1]$os_gridref
hab_by_select_s <- hab_cov[
  os_gridref %in% selected_cells_s,
  .(observed_area_s = sum(area_m2, na.rm = TRUE)),
  by = value
]

# Calculate observed habitat area for 85% footprint
selected_cells_85 <- pc85_footprint[footprint == 1]$os_gridref
hab_by_select_85 <- hab_cov[
  os_gridref %in% selected_cells_85,
  .(observed_area_85 = sum(area_m2, na.rm = TRUE)),
  by = value
]

# Calculate observed habitat area for largest / constrained footprint
selected_cells_100 <- pc100_footprint[footprint == 1]$os_gridref
hab_by_select_100 <- hab_cov[
  os_gridref %in% selected_cells_100,
  .(observed_area_100 = sum(area_m2, na.rm = TRUE)),
  by = value
]

# Merge expected and observed summaries
hab_compare <- Reduce(
  function(x, y) merge(x, y, by = "value", all = TRUE),
  list(hab_expected, hab_by_select_s, hab_by_select_85, hab_by_select_100)
)

# Replace missing area values with zero
area_cols <- c("expected_area", "observed_area_s", "observed_area_85", "observed_area_100")
hab_compare[,(area_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = area_cols]

# Calculate total areas
total_expected <- hab_compare[, sum(expected_area, na.rm = TRUE)]
total_s <- hab_compare[, sum(observed_area_s, na.rm = TRUE)]
total_85 <- hab_compare[, sum(observed_area_85, na.rm = TRUE)]
total_100 <- hab_compare[, sum(observed_area_100, na.rm = TRUE)]

# Calculate habitat proportions
hab_compare[, expected_prop := if (total_expected > 0) expected_area / total_expected else NA_real_]
hab_compare[, observed_prop_s := if (total_s > 0) observed_area_s / total_s else NA_real_]
hab_compare[, observed_prop_85 := if (total_85 > 0) observed_area_85 / total_85 else NA_real_]
hab_compare[, observed_prop_100 := if (total_100 > 0) observed_area_100 / total_100 else NA_real_]

# Compute observed / expected ratios
hab_compare[, ratio_s := fifelse(expected_prop > 0, observed_prop_s / expected_prop, NA_real_)]
hab_compare[, ratio_85 := fifelse(expected_prop > 0, observed_prop_85 / expected_prop, NA_real_)]
hab_compare[, ratio_100 := fifelse(expected_prop > 0, observed_prop_100 / expected_prop, NA_real_)]

# Log2 transform ratios
# Values of zero observed area become NA rather than -Inf
hab_compare[, log2_ratio_s := fifelse(ratio_s > 0, log2(ratio_s), NA_real_)]
hab_compare[, log2_ratio_85 := fifelse(ratio_85 > 0, log2(ratio_85), NA_real_)]
hab_compare[, log2_ratio_100 := fifelse(ratio_100 > 0, log2(ratio_100), NA_real_)]

# Mean enrichment across both scenarios
hab_compare[, mean_log2_ratio := rowMeans(.SD, na.rm = TRUE),.SDcols = c("log2_ratio_s", "log2_ratio_100")]

# Add habitat names
hab_compare <- merge(
  hab_compare,
  habitat_types_lookup_dt[, .(ebrel_hab_code, habitat)],
  by.x = "value",
  by.y = "ebrel_hab_code",
  all.x = TRUE
)

# Tidy habitat names before setting factor order
hab_compare[, habitat := gsub("_", " ", habitat)]
hab_compare[, habitat := tools::toTitleCase(habitat)]

# Order habitats so greatest positive differences appear at the top after coord_flip()
hab_order <- hab_compare[order(mean_log2_ratio), habitat]

# Reshape to long format for plotting
hab_long <- melt(
  hab_compare,
  id.vars = c("value", "habitat"),
  measure.vars = c("log2_ratio_s", "log2_ratio_85", "log2_ratio_100"),
  variable.name = "scenario",
  value.name = "log2_ratio"
)

# Rename scenarios
hab_long[
  ,
  scenario := fifelse(
    scenario == "log2_ratio_s",
    "GM target only\n(unconstrained)", fifelse(
      scenario == "log2_ratio_85",
      "GM target & 85% species reaching +10%\nabundance target",
      "GM target & 100% species reaching +10%\nabundance target"
    )
  )
]

# Apply habitat ordering
hab_long[, habitat := factor(habitat, levels = hab_order)]

# Plot observed vs expected habitat ratios
habitat_representation <- ggplot(hab_long, aes(x = habitat, y = log2_ratio, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = rev(
      colorspace::desaturate(
        hcl.colors(
          length(unique(range_plot_data$type)),
          palette = "Zissou 1"
        ),
        amount = 0.25
      )
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(log[2]("observed / expected habitat proportion")),
    fill = NULL
  )  +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.275, 0.8),
    legend.text = element_text(size = 14),
    legend.key.height = unit(40, "pt"),
    legend.spacing.y = unit(8, "pt"),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16)
  )

# = Save
ggsave("outputs/habitat_representation_figure.png", 
       plot = habitat_representation, width = 12, height = 8)


# ====
# === Species traits associated with selection
# ====

# Look at the traits associated with species selected early vs late

# === Species names
species_names <- gsub("_", " ", colnames(optimisation_data$P_over_R))

# === Re-load compiled data
asai_info_compiled <- fread("outputs/bdi_asai_compiled_species_information.csv")
asai_info_compiled_1174 <- asai_info_compiled[asai_info_compiled$BDI_NAME_CHECK %in% species_names, ]

# === Complete missing CLASS assignments
# Fish - Actinopterygii
asai_info_compiled_1174$CLASS[asai_info_compiled_1174$ORDER %in% c(
  "Anguilliformes", "Atheriniformes", "Clupeiformes", "Cypriniformes", "Esociformes",
  "Gasterosteiformes", "Mugiliformes", "Osmeriformes", "Perciformes",
  "Pleuronectiformes", "Salmoniformes", "Scorpaeniformes", "Syngnathiformes"
  )] <- "Actinopterygii"

asai_info_compiled_1174$CLASS[asai_info_compiled_1174$SCIENTIFIC_NAME %in% c(
  "Gymnocephalus cernuus", "Liza aurata", "Liza ramada"
  )] <- "Actinopterygii"

# Tubular worm
asai_info_compiled_1174$CLASS[asai_info_compiled_1174$ORDER %in% c("Tricladida")]  <-"Turbellaria"

# === Create plot of inclusion by class and K10 target

# = Find the N10 ranking for each list element
N10_list_index <- unlist(lapply(1:length(pareto_front_on_K_b1_e1), function(i) {
  pareto_front_on_K_b1_e1[[i]]$best_state$N10
}))

# = Add progressive targets to data.table
add_progressive_targets(
  asai_info_compiled = asai_info_compiled_1174,
  pareto_front_on_K = pareto_front_on_K_b1_e1,
  N10_list_index = N10_list_index,
  species_names = species_names
)

# = Create figure
prog_N10_class_plot <- plot_progressive_taxon_stack(
  asai_info_compiled_1174,
  taxon_col = "CLASS"
)
ggsave("outputs/progressive_N10_by_class.png", 
       plot = prog_N10_class_plot, 
       width = 8, height = 7)

# = Table for statistics
prog_N10_class_table <- plot_progressive_taxon_stack(
  asai_info_compiled_1174,
  taxon_col =   c("CLASS"),
  return_dt = TRUE
)

# = Target N30 taxon representation
prog_N10_class_table[,
  prop := N / sum(N), by = first_included][, 
  total_class_n := as.integer(sub(".*\\(([0-9]+)\\)$", "\\1", 
                                  as.character(taxon_label)))
  ][
    , prop_expect := total_class_n / 1174
  ][
    , delta_expect := (prop - prop_expect) 
  ][first_included == "N_30"]

# = Split taxon_label into CLASS, DOMAIN, and taxon_n
class_by_domain_by_N10 <- copy(prog_N10_class_table)
class_by_domain_by_N10[, c("CLASS", "domain_with_n") := tstrsplit(as.character(taxon_label)," \\| ")
   ][,DOMAIN := sub(" \\([0-9]+\\)$", "", domain_with_n)
     ][, taxon_n := as.integer(sub(".*\\(([0-9]+)\\)$", "\\1", domain_with_n))
       ][,domain_with_n := NULL]

# Count species per first_included class per DOMAIN
domain_first_included_summary <- class_by_domain_by_N10[,.(N = sum(N)), by = .(DOMAIN, first_included)]

# === Calculate proportion within each DOMAIN
domain_first_included_summary[,prop := N / sum(N), by = first_included]

# = N30 
n30_non_terrestrial <- domain_first_included_summary[first_included == "N_30" & 
                                !(DOMAIN %in% c("Terrestrial", "Mixed/Other"))]
sum(n30_non_terrestrial$N)/nS
# 47.96% of species select by n30 associated with Freshwater, marine, intertidal
# vs. 

# = Total domain counts / proportions
domain_proportions <- asai_info_compiled_1174 %>% 
  group_by(DOMAIN) %>% 
  summarise(N = n())

# Non-terrestrial
non_terr <- domain_proportions %>% 
  filter(!(DOMAIN %in% c("Terrestrial", "Mixed/Other")))
sum(non_terr$N) / nS
# 29.39% of the full species list

# ===
# == Range size by scenario
# ===

# === Load ranges
range_size <- fread("outputs/range_sizes.csv")

# == Unconstrained solution ranges
unconstrained_ranges <-
  range_size$range_size[which(pareto_front_on_K_b1_e1[[1]]$best_state$m_s >= 1.1)]
summary(unconstrained_ranges)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 4.578  3004.536  8438.707  7443.045 11299.746 19708.577 

# == 85% species representation
sol85_ranges <-
  range_size$range_size[which(pareto_front_on_K_b1_e1[[best_85pc_ind$x]]$best_state$m_s >= 1.1)]
summary(sol85_ranges)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 4.578  2827.642  6714.096  7558.253 10607.182 50484.072 

# == 100% species representation
summary(range_size$range_size)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.390e-01 3.187e+03 7.255e+03 1.471e+04 1.195e+04 1.203e+05 


# Number of species with large ranges than the largest in the unconstrainted solution
sum(range_size$range_size > 19708.577)

# == Test the expectation of inclusion for small range species 
# Use a permutation test to test the observed against expected inclusion of
# species with the smallest ranges (lowest 10%).

# Just for repeatability
set.seed(123)

# Set pools
selected_ranges <- unconstrained_ranges
pool_ranges     <- range_size$range_size
n_selected <- length(selected_ranges)

# Set smallest 10% of ranges
small_range_threshold <- quantile(
  pool_ranges,
  probs = 0.1,
  na.rm = TRUE
)

# Observed number with <= 10% range size
observed_n_small <- sum(selected_ranges <= small_range_threshold)

# Permutation selection of observed_n_small
n_iter = 100
null_n_small <- replicate(
  n_iter,
  {
    random_sample <- sample(
      pool_ranges,
      size = n_selected,
      replace = FALSE
    )
    
    sum(random_sample <= small_range_threshold)
  }
)

# P-value for this sample
p_value_small <- (
  sum(null_n_small >= observed_n_small) + 1
) / (n_iter + 1)

data.frame(
  threshold = small_range_threshold,
  observed_n_small = observed_n_small,
  expected_n_small = mean(null_n_small),
  observed_proportion = observed_n_small / n_selected,
  expected_proportion = mean(null_n_small) / n_selected,
  p_value = p_value_small
)

# === Create range size histogram
range_plot_data <- data.frame(
  range_size = c(
    unconstrained_ranges,
    sol85_ranges,
    range_size$range_size
  ),
  type = rep(
    c("GM target (unconstrained)", 
      "GM target & 85% species reaching +10% abundance target",
      "GM target & 100% species reaching +10% abundance target"),
    c(length(unconstrained_ranges), length(sol85_ranges), nrow(range_size))
  )
)

# Facet
range_size_by_solution <- ggplot(
  range_plot_data,
  aes(
    x = range_size,
    fill = type
  )
) +
  geom_histogram(
    bins = 100,
    alpha = 0.8,
    colour = "white",
    linewidth = 0.15
  ) +
  facet_wrap(
    ~ type,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = rev(
      colorspace::desaturate(
        hcl.colors(
          length(unique(range_plot_data$type)),
          palette = "Zissou 1"
        ),
        amount = 0.25
      )
    )
  ) +
  labs(
    x = expression("Range size (km"^2*")"),
    y = "Number of species"
  ) +
  guides(fill = "none") +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      hjust = 0
    ),
    plot.margin = margin(
      t = 5.5,
      r = 20,
      b = 5.5,
      l = 5.5,
      unit = "pt"
    )
  )

# = Save
ggsave(filename = "outputs/range_size_by_solution.png", 
       plot = range_size_by_solution, width = 8, height = 10)



#------------------------------------------------------------------------------#
#----------------------------------END-----------------------------------------#
#------------------------------------------------------------------------------#