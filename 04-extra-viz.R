### Packages -------------------------------------------------------------------
needed_packages <- c("tidyverse", "here", "glue", "rstan", "tidybayes",
                     "hrbrthemes", "colorblindr")
load_packages <- function(x) {
  if (!(x %in% installed.packages())) {
    install.packages(x, repos = "https://cran.rstudio.com/")
  }
  suppressPackageStartupMessages(require(x, character.only = TRUE))
}
vapply(needed_packages, load_packages, logical(1))


### MIRT visual ----------------------------------------------------------------
mirt_model <- read_rds(here("output/estimated-models/mirt_stan.rds"))

# Example graphic
pick <- summary(mirt_model, pars = "theta")$summary %>%
  as_tibble(rownames = "param") %>%
  select(param, est = mean) %>%
  separate(param, into = c("param", "stu_id", "dim", "junk"),
           sep = "\\[|,|\\]", convert = TRUE) %>%
  select(stu_id, dim, est) %>%
  mutate(dim = case_when(dim == 1 ~ "English Language Arts",
                         dim == 2 ~ "Mathematics",
                         dim == 3 ~ "Science")) %>%
  spread(key = dim, value = est) %>%
  filter(`English Language Arts` > Mathematics, Science > Mathematics) %>%
  mutate(diff1 = `English Language Arts` - Mathematics,
         diff2 = Science - Mathematics) %>%
  top_n(n = 20, wt = diff1) %>%
  top_n(n = 1, wt = diff2)

theta_draws <- spread_draws(mirt_model, theta[j,d])

theta_draws %>%
  ungroup() %>%
  filter(j == pick$stu_id) %>%
  mutate(d = case_when(d == 1 ~ "English Language Arts",
                       d == 2 ~ "Mathematics",
                       d == 3 ~ "Science"))  %>%
  ggplot(aes(x = theta)) +
  facet_wrap(~ d, ncol = 1) +
  geom_density(fill = palette_OkabeIto[2], color = NA, show.legend = FALSE) +
  expand_limits(x = c(-3, 3)) +
  scale_x_continuous(breaks = seq(-3, 3, 1)) +
  scale_fill_OkabeIto() +
  labs(x = expression(theta), y = NULL) +
  theme_ipsum_ps(axis_title_size = 15) +
  theme(axis.text.y = element_blank(),
        rect = element_rect(fill = "transparent")) -> theta_plot

ggsave("theta_plot.png", plot = theta_plot, path = "figures/",
       width = 8, height = 8 * 0.618, units = "in", bg = "transparent",
       dpi = "retina")


### LCDM example profile -------------------------------------------------------
lcdm_model_sat <- read_rds(here("output/estimated-models/lcdm_sat_stan.rds"))

profile_draws <- spread_draws(lcdm_model_sat, prob_resp_attr[j,d])

# Example graphic
pick <- summary(lcdm_model_sat, pars = "prob_resp_attr")$summary %>%
  as_tibble(rownames = "param") %>%
  select(param, est = mean) %>%
  separate(param, into = c("param", "stu_id", "dim", "junk"),
           sep = "\\[|,|\\]", convert = TRUE) %>%
  select(stu_id, dim, est) %>%
  mutate(dim = case_when(dim == 1 ~ "Initial",
                         dim == 2 ~ "Precursor",
                         dim == 3 ~ "Target")) %>%
  spread(key = dim, value = est) %>%
  filter(Initial > Precursor, Precursor > Target,
         Target > 0.2, Target < 0.5) %>%
  top_n(n = -1, wt = Precursor)

profile_draws %>%
  ungroup() %>%
  filter(j == pick$stu_id) %>%
  mutate(prob_resp_attr = round(prob_resp_attr, 5)) %>%
  count(d, prob_resp_attr) %>%
  complete(d, prob_resp_attr = seq(0, 1, 0.00001), fill = list(n = 0)) %>%
  arrange(d, prob_resp_attr) %>%
  group_by(d) %>%
  mutate(cumulative = cumsum(n),
         prop = cumulative / max(cumulative),
         trans = 1 - prop) %>%
  ungroup() %>%
  mutate(d = factor(d, levels = c(1, 2, 3), labels = c("Initial", "Precursor",
                                                       "Target"))) %>%
  ggplot(aes(x = as.factor(d), y = prob_resp_attr)) +
  geom_tile(aes(alpha = trans), width = 0.5, fill = palette_OkabeIto[2],
            show.legend = FALSE) +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_y_percent(breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Probability of Mastery") +
  theme_ipsum_ps() +
  theme(panel.grid.major.x = element_blank(),
        rect = element_rect(fill = "transparent")) -> profile_plot

ggsave("profile_plot.png", plot = profile_plot, path = "figures/",
       width = 5, height = 5 * 0.618, units = "in", bg = "transparent",
       dpi = "retina")


### LCDM profile assignment ----------------------------------------------------
load(here("output/data-sets/sim_params.rda"))
lcdm_model_sat <- read_rds(here("output/estimated-models/lcdm_sat_stan.rds"))

alpha_patt <- rep(list(c(0L, 1L)), ndim) %>%
  set_names(glue("dim_{seq_len(ndim)}")) %>%
  expand.grid() %>%
  as_tibble() %>%
  mutate(total = rowSums(.)) %>%
  select(total, everything()) %>%
  arrange_at(vars(total, desc(-one_of("total")))) %>%
  select(-total)

summary(lcdm_model_sat, pars = "prob_resp_class")$summary %>%
  as_tibble(rownames = "param") %>%
  select(param, est = mean) %>%
  separate(param, into = c("param", "stu_id", "class", "junk"),
           sep = "\\[|,|\\]", convert = TRUE) %>%
  select(stu_id, class, est) %>%
  group_by(stu_id) %>%
  top_n(1, wt = est) %>%
  ungroup() %>%
  count(class) %>%
  left_join(
    alpha_patt %>%
      rowid_to_column(var = "class") %>%
      gather(key = "dim", value = "master", -class) %>%
      group_by(class) %>%
      summarize(label = paste0("[", paste(master, collapse = ","), "]")),
    by = "class"
  ) %>%
  mutate(label = fct_reorder(label, class)) %>%
  ggplot(aes(x = label, y = n)) +
  geom_col(fill = "#9BD3DD") +
  labs(x = "Attribute Profile", y = "Students") +
  theme_ipsum_ps() -> profile_assignment

ggsave("profile_assgn.png", plot = profile_assignment, path = "figures/",
       width = 8, height = 8 * 0.618, units = "in", bg = "transparent",
       dpi = "retina")


### DLM aggregation ------------------------------------------------------------
skills <- tibble(
  xmin = 0, xmax = 2,
  ymin = c(0, 3, 6, 9, 12),
  ymax = c(2, 5, 8, 11, 14),
  level = "Skill"
)

stand <- tibble(
  xmin = 4, xmax = 6,
  ymin = 0, ymax = 14,
  level = "Standard"
)

stand_skills <- tibble(
  xmin = 4.25, xmax = 5.75,
  ymin = c(0.25, 3.25, 6.25, 9.25, 12.25),
  ymax = c(1.75, 4.75, 7.75, 10.75, 13.75),
  level = "Skill"
)

concep <- tibble(
  xmin = 8, xmax = 12,
  ymin = 0, ymax = 14,
  level = "Conceptual Area"
)

concep_stand <- tibble(
  xmin = c(8.25, 10.25), xmax = c(9.75, 11.75),
  ymin = 0.25, ymax = 13.75,
  level = "Standard"
)

concep_stand_skill <- bind_rows(
  tibble(
    xmin = 8.50, xmax = 9.50,
    ymin = c(0.50, 3.50, 6.50, 9.50, 12.50),
    ymax = c(1.50, 4.50, 7.50, 10.50, 13.50),
    level = "Skill"
  ),
  tibble(
    xmin = 10.50, xmax = 11.50,
    ymin = c(0.50, 3.50, 6.50, 9.50, 12.50),
    ymax = c(1.50, 4.50, 7.50, 10.50, 13.50),
    level = "Skill"
  )
)

content <- tibble(
  xmin = 14, xmax = 22,
  ymin = 0, ymax = 14,
  level = "Content"
)

content_concep <- tibble(
  xmin = c(14.25, 18.25), xmax = c(17.75, 21.75),
  ymin = 0.25, ymax = 13.75,
  level = "Conceptual Area"
)

content_concep_stand <- tibble(
  xmin = c(14.5, 16.25, 18.5, 20.25), xmax = c(15.75, 17.5, 19.75, 21.50),
  ymin = 0.5, ymax = 13.5,
  level = "Standard"
)

content_concep_stand_skill <- bind_rows(
  tibble(
    xmin = 14.75, xmax = 15.50,
    ymin = c(0.75, 3.69, 6.63, 9.57, 12.50),
    ymax = c(1.50, 4.44, 7.38, 10.32, 13.25),
    level = "Skill"
  ),
  tibble(
    xmin = 16.5, xmax = 17.25,
    ymin = c(0.75, 3.69, 6.63, 9.57, 12.50),
    ymax = c(1.50, 4.44, 7.38, 10.32, 13.25),
    level = "Skill"
  ),
  tibble(
    xmin = 18.75, xmax = 19.50,
    ymin = c(0.75, 3.69, 6.63, 9.57, 12.50),
    ymax = c(1.50, 4.44, 7.38, 10.32, 13.25),
    level = "Skill"
  ),
  tibble(
    xmin = 20.50, xmax = 21.25,
    ymin = c(0.75, 3.69, 6.63, 9.57, 12.50),
    ymax = c(1.50, 4.44, 7.38, 10.32, 13.25),
    level = "Skill"
  )
)

pld <- tibble(
  xmin = 24, xmax = 27,
  ymin = c(0, 3.66, 7.33, 11),
  ymax = c(3, 6.66, 10.33, 14),
  level = "Content"
)

all_data <- bind_rows(skills, stand, stand_skills, concep, concep_stand,
                      concep_stand_skill, content, content_concep, content_concep_stand,
                      content_concep_stand_skill, pld)

text_data = tibble(
  x = c(1, 5, 10, 18, 25.5),
  y = 14.50,
  label = c("Skills", "Standards", "Content Strand", "Subject", "Performance Level")
)


ggplot(all_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_rect(aes(fill = level), show.legend = FALSE) +
  geom_text(data = text_data, aes(x = x, y = y, label = label,
                                  family = "IBMPlexSans"),
            inherit.aes = FALSE, size = 6) +
  scale_fill_viridis_d(option = "D",
                       limits = rev(c("Skill", "Standard", "Conceptual Area",
                                      "Content"))) +
  expand_limits(x = c(0, 28)) +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank()) -> agg_levels

ggsave("agg_levels.png", plot = agg_levels, path = "figures/",
       width = 10, height = 8 * 0.618, units = "in", bg = "transparent",
       dpi = "retina")
