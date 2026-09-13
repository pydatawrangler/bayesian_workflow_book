library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
library(posterior)
library(brms)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 16))
library(marginaleffects)

print_stan_file <- function(file) {
  code <- readLines(file)
  if (isTRUE(getOption("knitr.in.progress")) &
        identical(knitr::opts_current$get("results"), "asis")) {
    # In render: emit as-is so Pandoc/Quarto does syntax highlighting
    block <- paste0("```stan", "\n", paste(code, collapse = "\n"), "\n", "```")
    knitr::asis_output(block)
  } else {
    writeLines(code)
  }
}

df_bioassay <- read.csv("data/bioassay.csv")

with(df_bioassay,
  plot(dose, deaths, xlab = "Dose log(g/ml)", ylab = "# of deaths",
  pch = 19, cex = 1.5, bty = "l"))

df_bioassay |> 
  ggplot(aes(x = dose, y = deaths)) +
  geom_point(size = 3) +
    labs(x = "Dose log(g/ml)", y = "# of deaths")

bioassay_data <- with(df_bioassay,
                      list(J = nrow(df_bioassay),
                           x = dose,
                           n = batch_size,
                           y = deaths))

bioassay_stan_file <- root("bioassay1.stan")
print_stan_file(bioassay_stan_file)

mod1 <- cmdstan_model(bioassay_stan_file, pedantic = TRUE)

fit1 <- mod1$sample(data = bioassay_data, refresh = 0)

fit1

draws1 <- fit1$draws(format="df")

par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
with(df_bioassay,
     plot(dose, deaths/batch_size, xlab = "Dose log(g/ml)", ylab = "Pr (death)",
          pch=19, cex=1.5, bty="l"))
invlogit <- plogis
for (s in sample(nrow(draws1), 20)) {
  curve(invlogit(draws1$a[s] + draws1$b[s] * x),
        col = "red", lwd = 0.5, add = TRUE)
}
curve(invlogit(mean(draws1$a) + mean(draws1$b) * x),
      col = "blue", lwd = 2, add = TRUE)

draws1 |>
  resample_draws(ndraws = 20) |>
  expand_grid(x=seq(-1, 1, length = 100)) |>
  mutate(y = plogis(a + b * x)) |>
  ggplot() +
  geom_point(data = df_bioassay, aes(x = dose, y = deaths / batch_size), size = 3) +
  geom_line(aes(x = x, y = y, group = .draw), alpha = .5, color = "red") +
  geom_function(fun = \(x) plogis(mean(draws1$a) + mean(draws1$b)*x),
                color = "blue",
                linewidth = 1) +
  labs(x = "Dose log(g/ml)", y = "Pr(death)")

draws1 <- draws1 |>
  mutate_variables(LD50_log_g_ml = -a / b,
                   LD50_mg_ml = 1000 * exp(LD50_log_g_ml))

draws1$LD50_log_g_ml <- -draws1$a /  draws1$b
draws1$LD50_mg_ml <- 1000 * exp(draws1$LD50_log_g_ml)

draws1 |>
  subset_draws(variable = "LD50_log_g_ml") |>
  summarize_draws()

ggplot_LD50_mg_ml <-
  ggplot(mapping = aes(x = LD50_mg_ml)) +
  scale_x_log10(limits = c(375,2700)) +
  labs(x = "LD50 mg/ml") +
  geom_vline(xintercept = c(500, 2000), alpha=0.1) +
  annotate(geom = "text", x = 500*0.96, y = .97, hjust = 1, label = "Category 3") +
  annotate(geom = "text", x = 1000, y = .97, label = "Category 4") +
  annotate(geom = "text", x = 2000*1.04, y = .97, hjust = 0, label = "Category 5") +
  theme_sub_axis_y(text = element_blank(),
                                        line = element_blank(),
                                        ticks = element_blank(),
                                        title = element_blank())

ggplot_LD50_mg_ml +
  stat_dots(data = draws1, quantiles = 100) +
  coord_cartesian(expand = c(bottom = FALSE))

ggplot_LD50_mg_ml +
  stat_slab(data = draws1, color = "gray", fill = NA) +
  coord_cartesian(expand = c(bottom = FALSE))

# brms model and inference

bfit1 <- brm(deaths | trials(batch_size) ~ dose,
             family = binomial(),
             prior = c(prior(normal(0, 5), class = Intercept),
                       prior(normal(0, 5), lb = 0, class = b)),
             data = df_bioassay,
             refresh = 0)

bfit1

bdraws1 <- as_draws_df(bfit1)

bdraws1 |>
  resample_draws(ndraws = 20) |>
  expand_grid(x=seq(-1, 1, length = 100)) |>
  mutate(y = plogis(b_Intercept + b_dose * x)) |>
  ggplot() +
  geom_line(aes(x = x, y = y, group = .draw), alpha = .5, color = "red") +
  geom_point(data = df_bioassay, aes(x = dose, y = deaths / batch_size), size = 3) +
  geom_function(fun = \(x) plogis(mean(bdraws1$b_Intercept) + mean(bdraws1$b_dose)*x),
                color = "blue",
                linewidth = 1) +
  labs(x = "Dose log(g/ml)", y = "Pr(death)")

p1 <- plot(conditional_effects(bfit1), plot=FALSE)[[1]] +
  labs(x = "Dose log(g/ml)", y = "Pr(death)")
p1

p1 +
  geom_point(data = df_bioassay,
  inherit.aes = FALSE,
aes(x = dose, y = deaths / batch_size),
size = 3)
