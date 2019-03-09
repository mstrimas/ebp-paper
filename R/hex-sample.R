hex_sample <- function(x, spacing = 5, 
                       regime = c("both", "positive", "negative")) {
  stopifnot(is.data.frame(x), 
            c("observation_date", "longitude", "latitude") %in% names(x))
  stopifnot(is.numeric(spacing), length(spacing) == 1, spacing > 0)
  regime <- match.arg(regime)
  
  # generate hexagonal grid
  dggs <- dggridR::dgconstruct(spacing = spacing)
  # get hexagonal cell id and week number for each checklist
  x_cell <- x %>% 
    mutate(cell = dggridR::dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
           week = lubridate::week(observation_date),
           obs = as.logical(species_observed))
  # sample one checklist per grid cell per week
  # sample detection/non-detection independently 
  if (regime == "both") {
    x_ss <- x_cell %>% 
      dplyr::group_by(obs, week, cell) %>% 
      dplyr::sample_n(size = 1) %>% 
      dplyr::ungroup()
  } else if (regime == "positive") {
    x_ss <- x_cell %>% 
      dplyr::filter(obs) %>% 
      dplyr::group_by(week, cell) %>% 
      dplyr::sample_n(size = 1) %>% 
      dplyr::ungroup()
    x_ss <- dplyr::bind_rows(x_ss, dplyr::filter(x_cell, !obs))
  } else if (regime == "negative") {
    x_ss <- x_cell %>% 
      dplyr::filter(!obs) %>% 
      dplyr::group_by(week, cell) %>% 
      dplyr::sample_n(size = 1) %>% 
      dplyr::ungroup()
    x_ss <- dplyr::bind_rows(x_ss, dplyr::filter(x_cell, obs))
  }
  dplyr::select(x_ss, -obs, -cell, -week)
}
