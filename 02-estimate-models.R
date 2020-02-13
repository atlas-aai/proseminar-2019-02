### Packages -------------------------------------------------------------------
needed_packages <- c("tidyverse", "here", "glue", "rstan", "fs")
load_packages <- function(x) {
  if (!(x %in% installed.packages())) {
    install.packages(x, repos = "https://cran.rstudio.com/")
  }
  suppressPackageStartupMessages(require(x, character.only = TRUE))
}
vapply(needed_packages, load_packages, logical(1))

if (!dir_exists("output/estimated-models/")) {
  dir_create("output/estimated-models/")
}


### Get simulation parameters --------------------------------------------------
load(here("output/data-sets/sim_params.rda"))


### Estimate MIRT model --------------------------------------------------------
load(here("output/data-sets/mirt_simvalues.rda"))
mirt_response <- read_rds(here("output/data-sets/mirt_data.rds"))

# Stan data
stan_mirt <- list(
  I = nitm,
  J = nstu,
  D = ndim,
  N = nrow(mirt_response),
  ii = mirt_response$item_id,
  jj = mirt_response$stu_id,
  dd = mirt_response$dim,
  y = mirt_response$score
)

# Stan model
mirt_model <- stan(file = here("stan-models/mirt.stan"), data = stan_mirt,
                   chains = 4, iter = 2000, warmup = 1000, cores = 4,
                   control = list(adapt_delta = 0.95), seed = 1992)

write_rds(mirt_model, path = "output/estimated-models/mirt_stan.rds",
          compress = "gz")

# Check output
summary(mirt_model, pars = "theta")$summary %>%
  as_tibble(rownames = "param") %>%
  select(param, est = mean) %>%
  separate(param, into = c("param", "stu_id", "dim", "junk"),
           sep = "\\[|,|\\]", convert = TRUE) %>%
  select(stu_id, dim, est) %>%
  left_join(
    theta %>%
      gather(key = "dim", value = "true", -stu_id) %>%
      mutate(dim = as.integer(str_replace(dim, "dim_", ""))),
    by = c("stu_id", "dim")
  ) %>%
  ggplot(aes(x = true, y = est)) +
  facet_wrap(~ dim) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed()


### Estimate saturated LCDM ----------------------------------------------------
load(here("output/data-sets/lcdm_simvalues.rda"))
lcdm_response <- read_rds(here("output/data-sets/lcdm_data.rds"))
qmatrix <- read_rds(here("output/data-sets/qmatrix.rds"))

alpha_patt <- rep(list(c(0L, 1L)), ndim) %>%
  set_names(glue("dim_{seq_len(ndim)}")) %>%
  expand.grid() %>%
  as_tibble() %>%
  mutate(total = rowSums(.)) %>%
  select(total, everything()) %>%
  arrange_at(vars(total, desc(-one_of("total")))) %>%
  select(-total)

# Stan data
xi_sat <- matrix(0, nitm, nrow(alpha_patt))
for (i in 1:nitm) {
  for(c in 1:nrow(alpha_patt)) {
    xi_sat[i, c] <- prod(alpha_patt[c,]^select(qmatrix, -item_id)[i,])
  }
}

ragged_array <- lcdm_response %>%
  rowid_to_column() %>%
  group_by(stu_id) %>%
  summarize(start = min(rowid), num = n())

stan_lcdm_sat <- list(
  I = nitm,
  J = nstu,
  K = ndim,
  C = nrow(alpha_patt),
  N = nrow(lcdm_response),
  ii = lcdm_response$item_id,
  jj = lcdm_response$stu_id,
  y = lcdm_response$score,
  s = ragged_array$start,
  l = ragged_array$num,
  alpha = as.matrix(alpha_patt),
  xi = xi_sat
)

lcdm_inits <- map(seq_len(4), function(x, num) {
  list(
    intercept  = stats::runif(num, -2.25, -1.00),
    maineffect = stats::runif(num,  1.00,  4.50)
  )
},
                  num = nitm)

# Stan model
lcdm_model_sat <- stan(file = "stan-models/lcdm.stan", data = stan_lcdm_sat,
                       init = lcdm_inits, chains = 4, iter = 2000,
                       warmup = 1000, cores = 4, seed = 1992,
                       control = list(adapt_delta = 0.95))

write_rds(lcdm_model_sat, path = "output/estimated-models/lcdm_sat_stan.rds",
          compress = "gz")

# Check outputs
summary(lcdm_model_sat, pars = c("intercept", "maineffect"))$summary %>%
  as_tibble(rownames = "param") %>%
  select(param, est = mean) %>%
  separate(param, into = c("param", "item_id"), sep = "\\[") %>%
  mutate(param = case_when(param == "intercept" ~ "int", TRUE ~ "mef"),
         item_id = as.integer(str_replace(item_id, "\\]", ""))) %>%
  arrange(item_id, param) %>%
  left_join(
    lcdm_params %>%
      gather(key = "param", value = "true", -item_id),
    by = c("item_id", "param")
  ) %>%
  ggplot(aes(x = true, y = est)) +
  facet_wrap(~ param, scales = "free") +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

summary(lcdm_model_sat, pars = "prob_resp_attr")$summary %>%
  as_tibble(rownames = "param") %>%
  select(param, est = mean) %>%
  separate(param, into = c("param", "stu_id", "dim", "junk"),
           sep = "\\[|,|\\]", convert = TRUE) %>%
  select(stu_id, dim, est) %>%
  left_join(
    profiles %>%
      gather(key = "dim", value = "mastery", -stu_id) %>%
      mutate(dim = as.integer(str_replace(dim, "dim_", ""))),
    by = c("stu_id", "dim")
  )


### Estimate reduced LCDM ------------------------------------------------------
load(here("output/data-sets/lcdm_simvalues.rda"))
lcdm_response <- read_rds(here("output/data-sets/lcdm_data.rds"))
qmatrix <- read_rds(here("output/data-sets/qmatrix.rds"))

reduc_patt <- rep(list(c(0L, 1L)), ndim) %>%
  set_names(glue("dim_{seq_len(ndim)}")) %>%
  expand.grid() %>%
  mutate(total = rowSums(.)) %>%
  select(total, everything()) %>%
  arrange_at(vars(total, desc(-one_of("total")))) %>%
  rowid_to_column() %>%
  group_by(total) %>%
  top_n(-1, wt = rowid) %>%
  ungroup() %>%
  select(-rowid, -total)

# Stan data
xi_red <- matrix(0, nitm, nrow(reduc_patt))
for (i in 1:nitm) {
  for(c in 1:nrow(reduc_patt)) {
    xi_red[i, c] <- prod(reduc_patt[c,]^select(qmatrix, -item_id)[i,])
  }
}

ragged_array <- lcdm_response %>%
  rowid_to_column() %>%
  group_by(stu_id) %>%
  summarize(start = min(rowid), num = n())

stan_lcdm_red <- list(
  I = nitm,
  J = nstu,
  K = ndim,
  C = nrow(reduc_patt),
  N = nrow(lcdm_response),
  ii = lcdm_response$item_id,
  jj = lcdm_response$stu_id,
  y = lcdm_response$score,
  s = ragged_array$start,
  l = ragged_array$num,
  alpha = as.matrix(reduc_patt),
  xi = xi_red
)

lcdm_inits <- map(seq_len(4), function(x, num) {
  list(
    intercept  = stats::runif(num, -2.25, -1.00),
    maineffect = stats::runif(num,  1.00,  4.50)
  )
},
                  num = nitm)

# Stan model
lcdm_model_red <- stan(file = "stan-models/lcdm.stan", data = stan_lcdm_red,
                       init = lcdm_inits, chains = 4, iter = 2000,
                       warmup = 1000, cores = 4, seed = 1992,
                       control = list(adapt_delta = 0.95))

write_rds(lcdm_model_red, path = "output/estimated-models/lcdm_red_stan.rds",
          compress = "gz")


### Estimate separate latent class models --------------------------------------
lcdm_response <- read_rds(here("output/data-sets/lcdm_data.rds"))

lcdm_response %>%
  nest(-dim) %>%
  pwalk(.f = function(dim, data) {
    ragged_array <- data %>%
      rowid_to_column() %>%
      group_by(stu_id) %>%
      summarize(start = min(rowid), num = n())
    
    stan_lca <- list(
      I = length(unique(data$item_id)),
      J = length(unique(data$stu_id)),
      N = nrow(data),
      ii = data$item_id - ((dim - 1) * 10),
      jj = data$stu_id,
      y = data$score,
      s = ragged_array$start,
      l = ragged_array$num
    )
    
    lca_inits <- map(seq_len(4), function(x, num) {
      list(
        intercept  = stats::runif(num, -2.25, -1.00),
        maineffect = stats::runif(num,  1.00,  4.50)
      )
    },
                     num = stan_lca$I)
    
    lca_model <- stan(file = "stan-models/lca.stan", data = stan_lca,
                      init = lca_inits, chains = 4, iter = 2000,
                      warmup = 1000, cores = 4, seed = 1992,
                      control = list(adapt_delta = 0.95))
    
    write_rds(lca_model,
              path = here(glue("output/estimated-models/lca_dim_{dim}.rds")),
              compress = "gz")
  })
