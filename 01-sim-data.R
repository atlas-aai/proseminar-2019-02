### Packages -------------------------------------------------------------------
needed_packages <- c("tidyverse", "glue", "mvtnorm", "bindata")
load_packages <- function(x) {
  if (!(x %in% installed.packages())) {
    install.packages(x, repos = "https://cran.rstudio.com/")
  }
  suppressPackageStartupMessages(require(x, character.only = TRUE))
}
vapply(needed_packages, load_packages, logical(1))


### Simulation parameters ------------------------------------------------------
nstu <- 1000
ndim <- 3
nitm <- ndim * 10

dim_cor <- matrix(data = c(1.0, 0.5, 0.2,
                           0.5, 1.0, 0.6,
                           0.2, 0.6, 1.0),
                  nrow = ndim, byrow = TRUE)

save(nstu, ndim, nitm, dim_cor, file = "output/data-sets/sim_params.rda")


### Functions ------------------------------------------------------------------
logit <- function(x) {
  log(x / (1 - x))
}
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}


### Make Q-matrix --------------------------------------------------------------
qmatrix <- map_dfc(seq_len(ndim), function(x, nitm, ndim) {
  att_name <- glue("dim_{x}")
  totalone <- nitm / ndim
  
  tibble(item_id = seq_len(nitm),
         section = rep(seq_len(ndim), each = totalone)) %>%
    mutate(!!att_name := case_when(section == x ~ 1L, TRUE ~ 0L)) %>%
    select(!!att_name)
},
                   nitm = nitm, ndim = ndim) %>%
  rowid_to_column(var = "item_id")


### Simulate MIRT data ---------------------------------------------------------
set.seed(9416)

# People
theta <- rmvnorm(n = nstu, mean = rep(0, ndim), sigma = dim_cor) %>%
  as_tibble(.name_repair = ~ glue("dim_{seq_len(ndim)}")) %>%
  rowid_to_column(var = "stu_id")

# Items
a_mean <- 1.25
a_sd <- 0.5
rl_loc <- log(a_mean^2 / sqrt(a_sd^2 + a_mean^2))
rl_shp <- sqrt(log(1 + (a_sd^2 / a_mean^2)))

mirt_params <- tibble(item_id = seq_len(nitm),
                      a = rlnorm(nitm, meanlog = rl_loc, sdlog = rl_shp),
                      b = rnorm(nitm, mean = 0, sd = 0.7))

# Response data
mirt_response <- theta %>%
  nest(-stu_id) %>%
  mutate(resp_data = map(data, function(x, mirt_params, qmatrix) {
    qmatrix %>%
      gather(key = "dim", value = "value", -item_id) %>%
      filter(value == 1L) %>%
      left_join(gather(x, key = "dim", value = "theta"), by = "dim") %>%
      left_join(mirt_params, by = "item_id") %>%
      mutate(prob = inv_logit(a * (theta - b)),
             rand = runif(n = nrow(.), min = 0, max = 1),
             score = case_when(rand <= prob ~ 1L, TRUE ~ 0L)) %>%
      select(item_id, score)
  },
                         mirt_params = mirt_params, qmatrix = qmatrix)) %>%
  select(-data) %>%
  unnest() %>%
  left_join(
    qmatrix %>%
      gather(key = "dim", value = "value", -item_id) %>%
      filter(value == 1L) %>%
      separate(dim, into = c("junk", "dim"), convert = TRUE) %>%
      select(item_id, dim),
    by = "item_id"
  )

save(theta, mirt_params, file = "output/data-sets/mirt_simvalues.rda")
write_rds(mirt_response, path = "output/data-sets/mirt_data.rds",
          compress = "gz")


### Simulate LCDM data ---------------------------------------------------------
set.seed(9416)

# People
profiles <- rmvbin(n = nstu, margprob = c(0.8, 0.5, 0.2), sigma = dim_cor) %>%
  as_tibble(.name_repair = ~ glue("dim_{seq_len(ndim)}")) %>%
  rowid_to_column(var = "stu_id")

# Items
mef_mean <- 3
mef_sd <- 1
rl_loc <- log(mef_mean^2 / sqrt(mef_sd^2 + mef_mean^2))
rl_shp <- sqrt(log(1 + (mef_sd^2 / mef_mean^2)))

lcdm_params <- tibble(item_id = seq_len(nitm),
                      int = rnorm(nitm, mean = -1.5, sd = 0.5),
                      mef = rlnorm(nitm, meanlog = rl_loc, sdlog = rl_shp))

# Response data
lcdm_response <- profiles %>%
  nest(-stu_id) %>%
  mutate(resp_data = map(data, function(x, lcdm_params, qmatrix) {
    qmatrix %>%
      gather(key = "dim", value = "value", -item_id) %>%
      filter(value == 1L) %>%
      left_join(gather(x, key = "dim", value = "mastery"), by = "dim") %>%
      left_join(lcdm_params, by = "item_id") %>%
      mutate(prob = inv_logit(int + (mef * mastery)),
             rand = runif(n = nrow(.), min = 0, max = 1),
             score = case_when(rand <= prob ~ 1L, TRUE ~ 0L)) %>%
      select(item_id, score)
  },
  lcdm_params = lcdm_params, qmatrix = qmatrix)) %>%
  select(-data) %>%
  unnest() %>%
  left_join(
    qmatrix %>%
      gather(key = "dim", value = "value", -item_id) %>%
      filter(value == 1L) %>%
      separate(dim, into = c("junk", "dim"), convert = TRUE) %>%
      select(item_id, dim),
    by = "item_id"
  ) %>%
  select(stu_id, item_id, dim, score) %>%
  arrange(stu_id, item_id)

save(profiles, lcdm_params, file = "output/data-sets/lcdm_simvalues.rda")
write_rds(lcdm_response, path = "output/data-sets/lcdm_data.rds",
          compress = "gz")
