
# for example
#  data <- doses_per_combination(
#    treatment_1_IC50 = -7.25,   # 56 nM
#    treatment_1_IC50 = -7.25,   # 56 nM
#    ray_angle = pi/2,
#    n_doses_per_ray = 8,
#    dose_range = 3.5)
#  10^(data$treatment_1_doses_log_molar + 6) = c(1, 3, 10, 30, 100, 300, 1000, 3000) # nM
doses_per_combination <- function(
    treatment_1_IC50,
    treatment_2_IC50,
    ray_angle,
    n_doses_per_ray,
    dose_range = 3.5){

    # not confident in this math
    # the idea is given the angle, come up with the coordinates where it crosses the boarder
    # of the 1-1 square
    if(ray_angle > pi/2){
        ray_multiplier_1 <- 1/tan(pi/2 - ray_angle) # cotanget
        ray_multiplier_2 <- 1
    } else if(ray_angle == pi/2){
        ray_multiplier_1 <- 1
        ray_multiplier_2 <- 1
    } else {
        ray_multiplier_1 <- 1
        ray_multiplier_2 <- 1/tan(ray_angle)
    }
    data.frame(
        treatment_1_dose = seq(
            from = treatment_1_IC50 - dose_range / 2,
            to = (treatment_1_IC50 + dose_range / 2) * ray_multiplier_1,
            length.out = n_doses_per_ray),
        treatment_2_dose = seq(
            from = treatment_2_IC50 - dose_range / 2,
            to = (treatment_2_IC50 + dose_range / 2) * ray_multiplier_2,
            length.out = n_doses_per_ray))
}

# treatments is a data.frame with columns
#   [treatment_id, IC50, is_primary, is_secondary]
#   IC50 should be in log molar concentration e.g. 1 uM would be -6
generate_combination_plate_maps <- function(
    treatments,
    n_rays_per_combination,
    n_doses_per_ray,
    n_replicas_per_dose,
    dose_range = 3.5,
    n_rows_per_plate = 16,
    n_columns_per_plate = 24,
    n_pc_per_plate = 32,
    n_nc_per_plate = 32,
    verbose = FALSE) {

    # we already have the single agent rays
    ray_angles <- seq(0, 90, length.out = n_rays_per_combination + 2) / 360 * 2 * pi
    ray_angles <- ray_angles[-length(ray_angles)][-1]
    
    if (verbose) {
        cat(
            "Computing ", n_rays_per_combination, " rays per combination with angles [", paste0(ray_angles, collapse = ", "), "] radians.", sep = "")
    }
    
    n_primary <- treatments %>% dplyr::filter(is_primary) %>% nrow()
    n_secondary <- treatments %>% dplyr::filter(is_secondary) %>% nrow()
    if (verbose) {
        cat(
            "  N primary:   ", n_primary, "\n",
            "  N secondary: ", n_secondary, "\n",
            "  N rays per combination: ", n_rays_per_combination, "\n",
            "  N doses per ray: ", n_doses_per_ray, "\n", sep = "")
    }

    # compute the plate coordinates for each well treatment
    # as a shortcut
    #  assume the left two columns are PC
    #  assume the right two columns are NC
    assertthat::assert_that(n_pc_per_plate == 2 * n_rows_per_plate)
    assertthat::assert_that(n_nc_per_plate == 2 * n_rows_per_plate)
    n_pc_columns_per_plate <- n_pc_per_plate / n_rows_per_plate
    n_nc_columns_per_plate <- n_nc_per_plate / n_rows_per_plate
    n_treatment_columns_per_plate <-
        n_columns_per_plate - n_pc_columns_per_plate - n_nc_columns_per_plate


    # add control wells
    n_wells_per_plate <- n_rows_per_plate * n_columns_per_plate
    n_treatment_wells_per_plate <- n_wells_per_plate - n_pc_per_plate - n_nc_per_plate    
    
    # primary vs secondary excluding primary vs primary
    # times n rays per combination
    treatment_combination_rays <- expand.grid(
        treatment_1 = treatments %>%
            dplyr::filter(is_primary) %>%
            magrittr::extract2("treatment_id"),
        treatment_2 = treatments %>%
            dplyr::filter(is_secondary) %>%
            magrittr::extract2("treatment_id"),
        ray_angle = ray_angles,
        replica = 1:n_replicas_per_dose,
        stringsAsFactors = FALSE)


    
    # exclude combinations where treatment_1 == treatment_2
    treatment_combination_rays <- treatment_combination_rays %>%
        dplyr::filter(treatment_1 != treatment_2)


    # exclude combinations (treatment_2, treatment_1) where (treament_1, treatment_2) is tested
    # this happens exactly when the a compound is both a primary and secondary treatment and
    # treatment_1 > treatment_2
    primary_and_secondary_treatments <- treatments %>% dplyr::filter(is_primary, is_secondary)

    treatment_combination_rays <- treatment_combination_rays %>%    
        dplyr::anti_join(
            treatment_combination_rays %>%
            dplyr::semi_join(primary_and_secondary_treatments, by = c("treatment_1" = "treatment_id")) %>%
            dplyr::semi_join(primary_and_secondary_treatments, by = c("treatment_2" = "treatment_id")) %>%
            dplyr::filter(treatment_1 > treatment_2),
            by = c("treatment_1", "treatment_2"))

    # add treatment metadata
    treatment_combination_rays <- treatment_combination_rays %>%
        dplyr::left_join(
            treatments %>%
            dplyr::select(
                treatment_1 = treatment_id,
                treatment_1_IC50 = IC50),
            by = c("treatment_1")) %>%
        dplyr::left_join(
            treatments %>%
            dplyr::select(
                treatment_2 = treatment_id,
                treatment_2_IC50 = IC50),
            by = c("treatment_2"))

   # compute the doses for each ray
   well_treatments <- treatment_combination_rays %>%
       dplyr::rowwise() %>%
       dplyr::do({
           treatment_combination <- .
           doses_per_combination(
               treatment_combination$treatment_1_IC50,
               treatment_combination$treatment_2_IC50,
               treatment_combination$ray_angle,
               n_doses_per_ray,
               dose_range = 3.5) %>%
               dplyr::mutate(
                   treatment_1 = treatment_combination$treatment_1,
                   treatment_1_IC50 = treatment_combination$treatment_1_IC50,
                   treatment_2 = treatment_combination$treatment_2,
                   treatment_2_IC50 = treatment_combination$treatment_2_IC50,
                   ray_angle = treatment_combination$ray_angle,
                   replica = treatment_combination$replica)
       }) %>%
       dplyr::ungroup()

    # assign well coordinates to each well treatment
    well_treatments <- well_treatments %>%
        dplyr::arrange(
            desc(treatment_1 == "Lactoferrin"),
            desc(treatment_2 == "Lactoferrin"),
            treatment_1,
            treatment_2,
            replica) %>%
        dplyr::mutate(
            treatment_index = dplyr::row_number(),
            plate_index = floor((treatment_index - 1) / n_treatment_wells_per_plate) + 1) %>%
        dplyr::group_by(plate_index) %>%
        dplyr::mutate(
            per_plate_treatment_index = dplyr::row_number(),
            row_index =
                floor((per_plate_treatment_index - 1) / n_treatment_columns_per_plate) + 1,
            column_index =
                (per_plate_treatment_index - 1) %% n_treatment_columns_per_plate +
                1 + n_pc_columns_per_plate) %>%
        dplyr::mutate(
            well_id = paste0(LETTERS[row_index], column_index)) %>%
        dplyr::ungroup()

    well_pc <- expand.grid(
        plate_index = well_treatments %>%
            dplyr::distinct(plate_index) %>%
            magrittr::extract2("plate_index"),
        per_plate_pc_index = 1:n_pc_per_plate) %>%
        dplyr::mutate(
            row_index = floor((per_plate_pc_index - 1) / n_pc_columns_per_plate) + 1,
            column_index = (per_plate_pc_index - 1) %% n_pc_columns_per_plate + 1) %>%
        dplyr::mutate(
            well_id = paste0(LETTERS[row_index], column_index))
        

    well_nc <- expand.grid(
        plate_index = well_treatments %>%
            dplyr::distinct(plate_index) %>%
            magrittr::extract2("plate_index"),
        per_plate_nc_index = 1:n_nc_per_plate) %>%
        dplyr::mutate(
            row_index = floor((per_plate_nc_index - 1) / n_nc_columns_per_plate) + 1,
            column_index = (per_plate_nc_index - 1) %% n_nc_columns_per_plate + 1 +
                n_pc_columns_per_plate + n_treatment_columns_per_plate) %>%
        dplyr::mutate(
            well_id = paste0(LETTERS[row_index], column_index))

    wells <- dplyr::bind_rows(
        well_treatments %>%
        dplyr::mutate(
            condition = "Treatment") %>%
        dplyr::select(
            plate_index,
            row_index,
            column_index,
            well_id,
            condition,
            ray_angle,
            replica,
            treatment_1,
            treatment_1_IC50,
            treatment_2,
            treatment_2_IC50,
            treatment_1_dose,
            treatment_2_dose),
        well_pc %>%
        dplyr::mutate(
            condition = "PC",
            ray_angle = NA,
            replica = NA,
            treatment_1 = NA,
            treatment_1_IC50 = NA,
            treatment_2 = NA,
            treatment_2_IC50 = NA,
            treatment_1_dose = NA,
            treatment_2_dose = NA) %>%
        dplyr::select(
            plate_index,
            row_index,
            column_index,
            well_id,
            condition,
            treatment_1,
            treatment_1_IC50,
            treatment_2,
            treatment_2_IC50,
            treatment_1_dose,
            treatment_2_dose),            
        well_nc %>%
        dplyr::mutate(
            condition = "NC",
            ray_angle = NA,
            replica = NA,            
            treatment_1 = NA,
            treatment_1_IC50 = NA,
            treatment_2 = NA,
            treatment_2_IC50 = NA,
            treatment_1_dose = NA,
            treatment_2_dose = NA) %>%
        dplyr::select(
            plate_index,
            row_index,
            column_index,
            well_id,
            condition,
            treatment_1,
            treatment_1_IC50,
            treatment_2,
            treatment_2_IC50,
            treatment_1_dose,
            treatment_2_dose))
    wells %>%
        dplyr::arrange(
            plate_index,
            row_index,
            column_index) %>%
        dplyr::mutate(
            treatment_1_dose_nM = 10^(treatment_1_dose + 9),
            treatment_2_dose_nM = 10^(treatment_2_dose + 9))
}
