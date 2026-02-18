library(move2)
library(dplyr)
library(purrr)
library(sf)

# -------------------------------------------------------------------
# get study info
# -------------------------------------------------------------------
studies <- movebank_download_study_info(i_have_download_access = TRUE) %>%
  filter(taxon_ids == "Limosa limosa") %>%
  filter(name != "Nora") %>%
  filter(license_type != "CC_0")

# -------------------------------------------------------------------
# download study data and attach deployment metadata
# -------------------------------------------------------------------
all_studies_data <- map(studies$id, movebank_download_study)
names(all_studies_data) <- as.character(studies$name)

for (i in seq_along(all_studies_data)) {
  deployment_info <- movebank_download_deployment(study_id = studies$id[i])
  all_studies_data[[i]] <- list(
    data = all_studies_data[[i]],
    deployment = deployment_info
  )
}

# -------------------------------------------------------------------
# flatten location data
# -------------------------------------------------------------------
events_all <- NULL
for (i in seq_along(all_studies_data)) {
  ev <- all_studies_data[[i]]$data
  # extract coordinates and bind to attribute table
  coords <- sf::st_coordinates(ev) %>%
    as.data.frame()
  attrs <- sf::st_drop_geometry(sf::st_as_sf(ev))
  ev_df <- cbind(coords, attrs)
  ev_df$study_name <- names(all_studies_data)[i]
  # stack across studies
  if (is.null(events_all)) {
    events_all <- ev_df
  } else {
    events_all <- bind_rows(events_all, ev_df)
  }
}

# -------------------------------------------------------------------
# flatten metadata
# -------------------------------------------------------------------
deployments_all <- NULL
for (i in seq_along(all_studies_data)) {
  dep <- all_studies_data[[i]]$deployment
  # add study name
  dep$study_name <- names(all_studies_data)[i]
  # stack across studies
  if (is.null(deployments_all)) {
    deployments_all <- dep
  } else {
    deployments_all <- bind_rows(deployments_all, dep)
  }
}

# -------------------------------------------------------------------
# join locations and metadata
# -------------------------------------------------------------------
events_joined <- events_all %>%
  left_join(
    deployments_all,
    by = c("study_name", "individual_local_identifier"),
    suffix = c("", "_dep")
  )

# -------------------------------------------------------------------
# clean up: remove duplicates/unwanted columns
# -------------------------------------------------------------------
events_joined <- events_joined %>%
  select(where(~ !all(is.na(.)))) %>%
  distinct(study_name, individual_local_identifier, timestamp, .keep_all = TRUE)

cols_to_remove <- c(
  "sensor_type_id","gps_satellite_count","light_level",
  "acceleration_raw_x","acceleration_raw_y","acceleration_raw_z",
  "gps_time_to_fix","magnetic_field_raw_x","magnetic_field_raw_y","magnetic_field_raw_z",
  "ornitela_transmission_protocol","argos_altitude","argos_best_level","argos_calcul_freq",
  "argos_iq","argos_nb_mes","argos_nb_mes_120","argos_nopc","argos_pass_duration",
  "argos_sensor_1","argos_sensor_2","argos_sensor_3","argos_sensor_4",
  "argos_valid_location_algorithm","argos_location_1","argos_location_2",
  "height_above_ellipsoid","argos_error_radius","argos_gdop","argos_orientation",
  "argos_semi_major","argos_semi_minor","gps_fix_type_raw","lotek_crc_status",
  "lotek_crc_status_text","manually_marked_valid","icarus_timestamp_accuracy",
  "icarus_timestamp_source","transmission_protocol",
  "location_accuracy_comments","deployment_id_dep","sensor_type_ids",
  "deploy_off_person","deploy_on_person","exact_date_of_birth",
  "group_id","mates","capture_method",
  "visible","barometric_height","capture_location",
  "deploy_on_location","deploy_off_location",
  "individual_number_of_deployments","mortality_location"
)

events_joined <- events_joined %>%
  select(-any_of(cols_to_remove))

events_joined <- events_joined %>%
  select(
    individual_local_identifier,
    X,
    Y,
    timestamp,
    study_name,
    everything()
  )

# -------------------------------------------------------------------
# remove bar-tailed/islandica
# -------------------------------------------------------------------

events_joined <- events_joined %>%
  filter(taxon_detail != "islandica" &
           taxon_detail != "Limosa lapponica" &
           taxon_detail != "Limosa limosa islandica")

# -------------------------------------------------------------------
# remove outliers
# -------------------------------------------------------------------

# calculate 5-point rolling averages and sd's for latitude and longitude
library(zoo)
events_joined <- events_joined %>%
  dplyr::group_by(individual_local_identifier) %>%
  dplyr::arrange(individual_local_identifier,timestamp) %>%
  dplyr::mutate(lat.mean.5d = zoo::rollapply(Y, FUN = mean, width=5, fill = NA),
                lat.sd.5d = zoo::rollapply(Y, FUN = sd, width=5, fill = NA),
                lon.mean.5d = zoo::rollapply(X, FUN = mean, width=5, fill = NA),
                lon.sd.5d = zoo::rollapply(X, FUN = sd, width=5, fill = NA)) %>%
  dplyr::ungroup()

# calculate DEV of each point to rolling average, square this TO SQDEV, and divide by SD
events_joined$lat.dev.to.roll <- abs(events_joined$Y - events_joined$lat.mean.5d)^2/events_joined$lat.sd.5d
events_joined$lon.dev.to.roll <- abs(events_joined$X - events_joined$lon.mean.5d)^2/events_joined$lon.sd.5d

# flag outliers if SQDEV of latitude or longitude > 10
events_joined$outlier <- ifelse(events_joined$lat.dev.to.roll > 10 | events_joined$lon.dev.to.roll > 10,'outlier','normal')
events_joined_outliers_removed <- subset(events_joined, outlier == "normal")

# export
write.csv(events_joined_outliers_removed,'all_locations.csv', row.names = FALSE)
