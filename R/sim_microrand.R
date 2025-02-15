
#' Title
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
skip.seed.streams <- function(n) {
  x <- .Random.seed
  for (i in seq_len(n))
    x <- parallel::nextRNGStream(x)
  assign('.Random.seed', x, pos=.GlobalEnv)
}

# start block
RNGkind("L'Ecuyer-CMRG")



#' Title
#'
#' @param N
#' @param capacity
#' @param return_dt_sim
#'
#' @return
#' @export
#'
#' @examples
sim_4t_microrand <- function(N, capacity = Inf, return_dt_sim = FALSE) {
  # N: scalar or vector of sample sizes to be simulated

  dt_aux0 <- lapply(N, \(x){
    data.frame(sample_size = x, pid = 1:x)
  }) %>% dplyr::bind_rows()

  # Generate baseline data
  dt_aux <- dt_aux0 %>%
    dplyr::group_split(sample_size, pid) %>%
    lapply(\(x) {
      dt <- data.frame(
        sample_size = x$sample_size[1],
        pid = x$pid[1],
        study_week = 13:52)

      dplyr::bind_rows(
        dt %>% dplyr::mutate(set = "A1", sd_g = .20, sd_e = .20/2, b3 = .005),
        dt %>% dplyr::mutate(set = "B1", sd_g = .16, sd_e = .16/2, b3 = .005),
        dt %>% dplyr::mutate(set = "A2", sd_g = .20, sd_e = .20/2, b3 = .0025),
        dt %>% dplyr::mutate(set = "B2", sd_g = .16, sd_e = .16/2, b3 = .0025),
        dt %>% dplyr::mutate(set = "A3", sd_g = .20, sd_e = .20/2, b3 = .001),
        dt %>% dplyr::mutate(set = "B3", sd_g = .16, sd_e = .16/2, b3 = .001),
        dt %>% dplyr::mutate(set = "A4", sd_g = .20, sd_e = .20/2, b3 = .0005),
        dt %>% dplyr::mutate(set = "B4", sd_g = .16, sd_e = .16/2, b3 = .0005),
        dt %>% dplyr::mutate(set = "A5", sd_g = .20, sd_e = .20/2, b3 = .0001),
        dt %>% dplyr::mutate(set = "B5", sd_g = .16, sd_e = .16/2, b3 = .0001),
        dt %>% dplyr::mutate(set = "C1", sd_g = .20, sd_e = .20/2, b3 = 0),
        dt %>% dplyr::mutate(set = "C2", sd_g = .16, sd_e = .16/2, b3 = 0)) %>%
        dplyr::group_by(set) %>%
        dplyr::mutate(
          time = study_week - 13,
          tide_week = (study_week - 12 + sample(1:4, 1) - 1) %% 4 + 1,

          gamma = rnorm(1, 0, sd_g^2),
          error = rnorm(length(study_week), 0, sd_e^2),

          treatment = 'Default',
          TIR = .75 + gamma - .005*time + error,
          iter = cumsum(as.numeric(tide_week == 1)),
          rand = 0) %>%
        dplyr::ungroup()
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(tide_week <= 2)

  # Iterate study week
  dt_sim <- dt_aux %>%
    dplyr::group_split(sample_size, set) %>%
    lapply(\(x) {
      aux_risk <- x %>%
        dplyr::filter(tide_week == 1, TIR < .65)

      if (nrow(aux_risk) == 0) return()

      min_k <- min(aux_risk$iter)
      dt_iter <- x
      for (k in min_k:10) {
        dt_iter <- dt_iter %>%
          dplyr::mutate(
            elig = ifelse(tide_week == 1, as.numeric(TIR < .65), NA),
            rand_trt = ifelse((iter == k)*elig == 1, sample(c(0, 1), nrow(.), replace = T), NA)) %>%
          dplyr::group_by(pid) %>%
          dplyr::mutate(
            rand_trt = ifelse(is.na(rand_trt), '', rand_trt),
            rand_trt = paste0(dplyr::lag(rand_trt, n = 1, default = ''),
                              dplyr::lag(rand_trt, n = 2, default = '')),
            rand_trt = ifelse(rand_trt == '', NA, as.numeric(rand_trt))) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            rand = rand + as.numeric(!is.na(rand_trt)),
            treatment = ifelse(!is.na(rand_trt) & rand_trt == 1, 'Addon', treatment),
            rand_TIR = .75 + gamma - .005*time + b3*time*rand_trt + error,
            TIR = ifelse(!is.na(rand_trt) & rand_trt == 1, rand_TIR, TIR))

        # print(dt_iter, n=100)
      }

      dt_iter %>%
        dplyr::mutate(
          sample_size = x$sample_size[1],
          set = x$set[1],
          rand_trt = NULL,
          rand_TIR = NULL) %>%
        dplyr::filter(iter > 1, tide_week == 1, rand == 1)
    }) %>%
    dplyr::bind_rows()

  # Impose capacity
  if (is.finite(capacity)) {
    dt_sim <- dt_sim %>%
      # dplyr::filter(sample_size==200, set=='B', time == 38) %>%
      dplyr::arrange(sample_size, set, time, treatment) %>%
      dplyr::group_by(sample_size, set, time) %>%
      dplyr::mutate(
        addon_cum = cumsum(as.numeric(treatment == 'Addon')),
        addon_cum = ifelse(treatment == 'Addon', addon_cum, NA)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        exceed_cap = ifelse(addon_cum > capacity, 1, 0),
        treatment = ifelse(!is.na(exceed_cap) & exceed_cap == 1, 'Default', treatment),
        TIR = ifelse(!is.na(exceed_cap) & exceed_cap == 1, .75 + gamma - .005*time + error, TIR)) %>%
      dplyr::arrange(sample_size, set, pid, time)
  }

  # Summary
  dt_summ <- dt_sim %>%
    dplyr::count(sample_size, set, time, treatment) %>%
    tidyr::pivot_wider(names_from = treatment, values_from = n) %>%
    dplyr::mutate(
      dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.na(.), 0, .)),
      atRisk = Addon + Default)

  if(return_dt_sim) return(dt_sim)

  # Fit LMM
  aux_time <- seq(8, 40, 8)
  safe_lmer <- purrr::possibly(purrr::quietly(lmer), otherwise = NULL)
  safe_broom_tidy <- purrr::possibly(\(m) {
    broom.mixed::tidy(m) %>%
      dplyr::filter(term == 'time:trt') %>%
      dplyr::select(-group)},
    otherwise = tidyr::tibble())

  dt_reg <- lapply(aux_time, \(t) {
    dt_mod = dt_sim %>%
      dplyr::filter(time <= t) %>%
      dplyr::mutate(trt = as.numeric(treatment == 'Addon')) %>%
      dplyr::nest_by(sample_size, set, b3) %>%
      dplyr::mutate(
        up_to_time = t,
        n_clust = length(unique(data$pid)),
        n_obs = nrow(data),
        mod = list(safe_lmer(TIR ~ time + time*trt + (1|pid), data = data)$result)) %>%
      dplyr::reframe(up_to_time, n_clust, n_obs, safe_broom_tidy(mod))
  }) %>% dplyr::bind_rows() %>% dplyr::filter(!is.na(term))

  return(list(lmer_reg = dt_reg, summary = dt_summ))
}


# sim_4t_microrand(N = 50)
# sim_4t_microrand(N = c(50, 100, 200))

# dt_sim = sim_4t_microrand(N = c(50, 100, 200), return_dt_sim = T)
#
# dt_summ <- dt_sim %>%
#   dplyr::count(sample_size, set, time, treatment) %>%
#   tidyr::pivot_wider(names_from = treatment, values_from = n) %>%
#   dplyr::mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)),
#          atRisk = Addon + Default)
#
# dt_sim %>%
#   dplyr::filter(sample_size==200, set=='B', time == 38) %>%
#   dplyr::arrange(sample_size, set, time, treatment) %>%
#   dplyr::group_by(sample_size, set, time) %>%
#   dplyr::mutate(
#     addon_cum = cumsum(as.numeric(treatment == 'Addon')),
#     addon_cum = ifelse(treatment == 'Addon', addon_cum, NA)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(
#     exceed_cap = ifelse(addon_cum > capacity, 1, 0),
#     treatment = ifelse(!is.na(exceed_cap) & exceed_cap == 1, 'Default', treatment),
#     TIR = ifelse(!is.na(exceed_cap) & exceed_cap == 1, .75 + gamma - .005*time + error, TIR)) %>%
#   dplyr::arrange(sample_size, set, pid, time)
#
#
# dt_summ %>%
#   dplyr::filter(sample_size==200, set=='B', time == 38)




