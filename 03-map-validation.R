### Packages -------------------------------------------------------------------
needed_packages <- c("tidyverse", "here", "glue", "fs", "rstan", "hrbrthemes",
                     "colorblindr")
load_packages <- function(x) {
  if (!(x %in% installed.packages())) {
    install.packages(x, repos = "https://cran.rstudio.com/")
  }
  suppressPackageStartupMessages(require(x, character.only = TRUE))
}
vapply(needed_packages, load_packages, logical(1))


### Patterns of Mastery Profiles -----------------------------------------------
# Posterior predictive replications -----
lcdm_model_sat <- read_rds(here("output/estimated-models/lcdm_sat_stan.rds"))
lcdm_model_red <- read_rds(here("output/estimated-models/lcdm_red_stan.rds"))

lcdm_response <- read_rds(here("output/data-sets/lcdm_data.rds"))

y_rep_sat <- rstan::extract(lcdm_model_sat, pars = "y_rep",
                            permuted = TRUE)$y_rep %>%
  as_tibble(.name_repair = "universal") %>%
  rowid_to_column(var = "iteration") %>%
  group_by(iteration) %>%
  nest()

y_rep_red <- rstan::extract(lcdm_model_red, pars = "y_rep",
                            permuted = TRUE)$y_rep %>%
  as_tibble(.name_repair = "universal") %>%
  rowid_to_column(var = "iteration") %>%
  group_by(iteration) %>%
  nest()

all_rep_data <- map2_dfr(.x = y_rep_sat$data, .y = y_rep_red$data,
                         .f = function(x, y, orig_data) {
                           sat_scores <- x %>%
                             tidyr::gather() %>%
                             pull(value)
                           red_scores <- y %>%
                             tidyr::gather() %>%
                             pull(value)
                           
                           orig_data %>%
                             select(stu_id, item_id, dim) %>%
                             mutate(sat_score = sat_scores,
                                    red_score = red_scores)
                         },
                         orig_data = lcdm_response, .id = "iteration")

write_rds(all_rep_data, path = "output/model-fit/replicated-data.rds",
          compress = "gz")

# Posterior predictive model checks -----
all_rep_data <- read_rds(here("output/model-fit/replicated-data.rds"))
lcdm_response <- read_rds(here("output/data-sets/lcdm_data.rds"))

# Raw score distribution
stu_scores <- all_rep_data %>%
  rename(Saturated = sat_score, Reduced = red_score) %>%
  gather(key = "model", value = "score", Saturated:Reduced) %>%
  group_by(iteration, stu_id, model) %>%
  summarize(total_score = sum(score)) %>%
  ungroup()

score_counts <- stu_scores %>%
  count(iteration, model, total_score) %>%
  complete(iteration, model, total_score, fill = list(n = 0)) %>%
  mutate(model = factor(model, levels = c("Saturated", "Reduced")))

exp_counts <- score_counts %>%
  group_by(model, total_score) %>%
  summarize(exp = mean(n),
            lb = quantile(n, probs = 0.025),
            ub = quantile(n, probs = 0.975)) %>%
  ungroup()

obs_counts <- lcdm_response %>%
  group_by(stu_id) %>%
  summarize(total_score = sum(score)) %>%
  ungroup() %>%
  count(total_score) %>%
  complete(total_score = 0:30, fill = list(n = 0))

ggplot() +
  facet_wrap(~ model, ncol = 1) +
  geom_jitter(data = filter(score_counts, iteration %in% sample(1:4000, 1000)),
              aes(x = total_score, y = n), width = 0.3, height = 0,
              alpha = 0.2) +
  geom_line(data = exp_counts, aes(x = total_score, y = lb,
                                   color = "95% Credible Interval"),
            linetype = "dashed") +
  geom_line(data = exp_counts, aes(x = total_score, y = ub,
                                   color = "95% Credible Interval"),
            linetype = "dashed") +
  geom_line(data = obs_counts, aes(x = total_score, y = n,
                                   color = "Observed"),
            linetype = "solid") +
  geom_point(data = obs_counts, aes(x = total_score, y = n,
                                    color = "Observed"),
             size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("95% Credible Interval" = "#F2A900", "Observed" = "#99D2DE")) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = "Raw Score", y = "Number of Students", color = NULL) +
  theme_ipsum_ps() +
  theme(legend.position = "bottom") -> raw_dist

ggsave("raw_score_ppmc.png", plot = raw_dist, path = "figures/",
       width = 8, height = 8 * 0.618, units = "in", bg = "transparent",
       dpi = "retina")

# Chi-square
chi_square <-
  left_join(score_counts,
            select(exp_counts, model, total_score, exp),
            by = c("model", "total_score")) %>%
  mutate(piece = ((n - exp)^2) / exp) %>%
  group_by(model, iteration) %>%
  summarize(chisq = sum(piece))

obs_chi <-
  left_join(obs_counts,
            select(exp_counts, model, total_score, exp),
            by = c("total_score")) %>%
  mutate(piece = ((n - exp)^2) / exp) %>%
  group_by(model) %>%
  summarize(chisq = sum(piece))

text <- chi_square %>%
  left_join(rename(obs_chi, obschi = chisq), by = "model") %>%
  mutate(exceed = chisq > obschi) %>%
  summarize(ppp = mean(exceed)) %>%
  mutate(label = paste0("italic(ppp) == ", sprintf("%0.2f", ppp)))

ggplot() +
  facet_wrap(~ model, ncol = 1) +
  geom_histogram(data = chi_square, aes(x = chisq, y = ..density..,
                                        color = "Posterior Distribution",
                                        fill = "Posterior Distribution"),
                 alpha = 0.8, binwidth = 2, boundary = 0) +
  geom_histogram(data = chi_square, aes(x = chisq, y = ..density..,
                                        fill = "Observed"),
                 alpha = 0, binwidth = 2) +
  geom_density(data = chi_square, aes(x = chisq), color = "black",
               show.legend = FALSE) +
  geom_vline(data = obs_chi, aes(xintercept = chisq, color = "Observed"),
             linetype = "dashed", size = 1) +
  geom_text(data = text, aes(x = 85, y = 0, label = label),
            hjust = 0, vjust = 0, parse = TRUE, family = "IBMPlexSans") +
  expand_limits(x = 0) +
  scale_fill_manual(values = c("Posterior Distribution" = "#F2A900",
                               "Observed" = "#99D2DE"),
                    breaks = c("Posterior Distribution", "Observed")) +
  scale_color_manual(values = c("Posterior Distribution" = "#F2A900",
                                "Observed" = "#99D2DE")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  labs(x = expression(chi^2), y = "Posterior Density", fill = NULL) +
  theme_ipsum_ps() +
  theme(legend.position = "bottom") +
  guides(color = FALSE) -> chisq_dist

ggsave("chisq_ppmc.png", plot = chisq_dist, path = "figures/",
       width = 8, height = 8 * 0.618, units = "in", bg = "transparent",
       dpi = "retina")  

# Relative fit -----
lcdm_model_sat <- read_rds(here("output/estimated-models/lcdm_sat_stan.rds"))
lcdm_model_red <- read_rds(here("output/estimated-models/lcdm_red_stan.rds"))

loo_sat <- loo(lcdm_model_sat, cores = 7)
loo_red <- loo(lcdm_model_red, cores = 7)

compare(loo_sat, loo_red)

write_rds(loo_sat, "output/model-fit/loo_sat.rds", compress = "gz")
write_rds(loo_red, "output/model-fit/loo_red.rds", compress = "gz")


### Patterns of Attribute Mastery ----------------------------------------------
att_mastery <- dir_ls(here("output/estimated-models/"), regexp = "lca") %>%
  map(.f = function(x) {
    dim <- glue("dim_{str_sub(x, -5, -5)}")
    
    est_model <- read_rds(x)
    
    summary(est_model, pars = "prob_resp_class")$summary %>%
      as_tibble(rownames = "param") %>%
      select(param, !!dim := mean) %>%
      separate(param, into = c("param", "stu_id", "class", "junk"),
               sep = "\\[|,|\\]", convert = TRUE) %>%
      filter(class == 2) %>%
      select(stu_id, !!dim)
  }) %>%
  reduce(full_join, by = "stu_id") %>%
  mutate_if(is.double, list(~ case_when(. >= 0.8 ~ 1L, TRUE ~ 0L))) %>%
  mutate(rev1 = dim_2 > dim_1, rev2 = dim_3 > dim_2, bad = rev1 | rev2)

att_mastery %>%
  count(rev1, rev2) %>%
  filter(rev1 | rev2) %>%
  mutate(label = case_when(rev1 & !rev2 ~ "Initial/\nPrecursor",
                           !rev1 & rev2 ~ "Precursor/\nTarget")) %>%
  ggplot(aes(x = label, y = n)) +
  geom_col(fill = palette_OkabeIto[2]) +
  labs(x = "Reversal", y = "Number of Students") +
  theme_ipsum_ps() -> attribute_mastery

ggsave("attribute_mastery.png", plot = attribute_mastery, path = "figures/",
       width = 8 * 0.618, height = 8, units = "in", bg = "transparent",
       dpi = "retina")


### Patterns of Linkage Level Difficulty ---------------------------------------
load(here("output/data-sets/lcdm_simvalues.rda"))
lcdm_response <- read_rds(here("output/data-sets/lcdm_data.rds"))

profiles <- profiles %>%
  mutate(group = case_when(dim_1 == 0 & dim_2 == 0 & dim_3 == 0 ~ 0L,
                           dim_1 == 1 & dim_2 == 0 & dim_3 == 0 ~ 1L,
                           dim_1 == 0 & dim_2 == 1 & dim_3 == 0 ~ 1L,
                           dim_1 == 0 & dim_2 == 0 & dim_3 == 1 ~ 2L,
                           dim_1 == 1 & dim_2 == 1 & dim_3 == 0 ~ 2L,
                           dim_1 == 1 & dim_2 == 0 & dim_3 == 1 ~ 2L,
                           dim_1 == 0 & dim_2 == 1 & dim_3 == 1 ~ 3L,
                           dim_1 == 1 & dim_2 == 1 & dim_3 == 1 ~ 3L))

diff <- left_join(lcdm_response, select(profiles, stu_id, group), by = "stu_id") %>%
  group_by(group, dim, item_id) %>%
  summarize(pct_cor = mean(score),
            n = n(),
            se = sqrt((pct_cor * (1 - pct_cor)) / n),
            weight = (1 / (se^2))) %>%
  mutate(weight = weight / sum(weight)) %>%
  summarize(mean_pct_cor = sum(weight * pct_cor),
            se = sqrt(mean_pct_cor * (1 - mean_pct_cor) * sum(weight^2))) %>%
  nest(-group) %>%
  mutate(ordering = map(data, .f = function(x) {
    ret_frame <- tibble(compare = rep(NA_character_, nrow(x) - 1),
                        incorrect = rep(NA, nrow(x) - 1))
    for(i in 1:nrow(ret_frame)) {
      lower <- x$mean_pct_cor[i] + x$se[i]
      upper <- x$mean_pct_cor[i + 1] - x$se[i + 1]
      ret_frame$compare[i] <- glue("{x$dim[i]}_{x$dim[i + 1]}")
      ret_frame$incorrect[i] <- upper >= lower
    }
    
    flag <- any(ret_frame$incorrect)
    ret_frame %>%
      spread(key = compare, value = incorrect) %>%
      mutate(flag = flag)
  }))

diff %>%
  select(group, data) %>%
  unnest() %>%
  mutate(lb = mean_pct_cor - se,
         ub = mean_pct_cor + se,
         lb = case_when(lb < 0 ~ 0, TRUE ~ lb),
         ub = case_when(ub > 1 ~ 1, TRUE ~ ub),
         dim = factor(dim, levels = 1:3, labels = c("Initial", "Precursor",
                                                    "Target")),
         group = factor(group, levels = 0:3, labels = c("Group 1",
                                                        "Group 2", "Group 3",
                                                        "Group 4"))) %>%
  ggplot() +
  facet_wrap(~ group, ncol = 1) +
  geom_point(aes(x = mean_pct_cor, y = dim, color = dim)) +
  geom_errorbarh(aes(xmin = lb, xmax = ub, y = dim, color = dim), height = 0.4) +
  scale_color_manual(values = c("#F2A900", "#165B7E", "#8DD4BE")) +
  labs(x = expression(paste("Group ", italic("p"), "-value")), y = NULL) +
  theme_ipsum_ps() +
  guides(color = FALSE) -> pvalues

ggsave("pvalues.png", plot = pvalues, path = "figures/",
       width = 4, height = 6, units = "in", bg = "transparent",
       dpi = "retina")

diff %>%
  select(group, ordering) %>%
  unnest()
