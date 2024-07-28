# Computes regressions evaluating:
#
# RQ3) the moderating effect OSA on cognitive performance in de novo PD
#
# The RQ is addressed via a series Bayesian test/re-test GLMMs with marginal means postprocessing of posteriors


rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

# load packages
library(here)
library(tidyverse)
library(brms)
library(priorsense)
library(bayestestR)
library(bayesplot)
library(ggh4x)
#library(patchwork)

theme_set( theme_bw() )
color_scheme_set("viridisA") # colour scheme

sapply( c("figures","tables"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders
psych <- read.csv( here("helpers","psychs.csv"), sep = ";") # read helper files with variable names

# in-house functions
source( here("scripts","utils.R") )


# PRE-PROCESSING  ----

# read the data
d0 <-
  
  # read and pre-process variables
  read.csv( here("_data","primary_dataset.csv"), sep = "," ) %>%
  mutate(
    id = Study.ID,
    sex = factor( case_when(GENDER == 1 ~ "male", GENDER == 0 ~ "female") ),
    group = factor( case_when(SUBJ == "PD" ~ "PD", SUBJ == "CON" ~ "HC") ),
    osa = factor( case_when(AHI.F == "H" ~ "+", AHI.F == "L" ~ "-") ),
    event = factor(
      case_when(event == "enrollment" ~ "enrollment", event == "followup_4" ~ "retest"),
      levels = c("enrollment", "retest"),
      ordered = T
    ),
    age = AGE,
    edu = EDU.Y
  ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(sex) <- contr.equalprior_pairs(2)
    contrasts(osa) <- contr.equalprior_pairs(2)
    contrasts(group) <- contr.equalprior_pairs(2)
    contrasts(event) <- contr.equalprior_pairs(2)
  } )

# extract response time vatiable names to be log-transformed
rt_vars <- subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable

# prepare data set for analyses
df <- d0 %>%
  
  select(id, event, sex, group, osa, event, age, edu, all_of(psych$variable) ) %>%
  filter( complete.cases(event) ) %>%
  mutate( across( all_of(rt_vars), ~ -log(.x) ) )

# amount of data points available
write.table(
  
  x = lapply(
    
    with( psych, setNames(variable,variable ) ),
    function(y)
      
      df[ complete.cases(df[, y]), c("event","group","osa") ] %>%
      table() %>%
      as.data.frame() %>%
      mutate(osa = paste0("OSA",osa) ) %>%
      mutate(
        domain = with( psych, domain[variable == y] ),
        y = with( psych, label[variable == y] ),
        .before = 1
      )
    
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    pivot_wider(names_from = c("event","group","osa"), names_sep = "_", values_from = Freq),
  
  file = here("tables","cognition_data_n.csv"),
  sep = ",",
  row.names = F,
  quote = F
    
)

# extract scaling values, i.e.,
# enrollment full sample means and SDs
scl <- sapply(
  
  c(psych$variable[-13], "age", "edu"),
  function(i)
    
    with(subset(df, event == "enrollment"), c(M = mean( get(i), na.rm = T ), SD = sd( get(i), na.rm = T ) ) )
  
) %>% t()

# scale the variables
df <- df %>%
  
  mutate(
    across(
      all_of( rownames(scl) ),
      ~ (.x - scl[cur_column(), "M"]) / scl[cur_column(), "SD"],
      .names = "{.col}_sc"
    )
  )


# MODEL FITTING ----

# set-up priors
prior <- c(
  
  prior("normal(0, 1)", class = "Intercept"),
  prior("normal(0, 1)", class = "b"),
  prior("exponential(1)", class = "sd", group = "id"),
  prior("exponential(1)", class = "sigma")
  
)

# fit the models
fit <- lapply(
  
  with( psych, setNames( paste0(psych$variable[-13],"_sc"), psych$variable[-13]) ),
  function(y) {
    
    data <- df %>% mutate( y = get(y) )
    
    brm(
      formula = bf( y ~ 1 + event * group * osa + (1 | id) ),
      data = data,
      prior = prior,
      family = gaussian(link = "identity"),
      seed = 87542
    )
    
  }
)


# MODEL CHECKING ----

# check prior sensitivity of parameters
lapply(
  
  names(fit),
  function(y) {
    
    print(y)
    print( powerscale_sensitivity(fit[[y]]) )
    Sys.sleep(3)
    
  }
)

# do posterior predictive checks of the mean for each event:group:osa combination
plot_ppc_stat(fit, df, stat = "mean") # looks very good

# do posterior predictive checks of the SD for each event*group*osa combination

# one-way comparisons
plot_ppc_stat(fit, df, stat = "sd", var = "event") # all ok
plot_ppc_stat(fit, df, stat = "sd", var = "group") # questionable cases: tmt_a, tmt_b, gpt_phk, and gpt_lhk
plot_ppc_stat(fit, df, stat = "sd", var = "osa") # all ok

# two-way interactions
plot_ppc_stat(fit, df, stat = "sd", var = "event_group") # tmt_b, and gpt_phk but most likely to be solved by group main effect
plot_ppc_stat(fit, df, stat = "sd", var = "event_osa") # all ok
plot_ppc_stat(fit, df, stat = "sd", var = "group_osa") # gpt_phk, and gpt_lhk but most likely to be solved by group main effect

# three-way interaction
plot_ppc_stat(fit, df, stat = "sd", var = "event_group_osa", sleep = 10) # all reasonably ok


# RESULTS EXTRACTION ----

# extract contrast codings
contr <- contr_extr(c("event","group","osa"), data = df)

# save posterior summaries
write.table(
  
  x = lapply(
    
    names(fit),
    function(y) {

      `2exp` <- ifelse(y %in% rt_vars, T, F) # back-transform reaction times via exponentiation
      return( contr_comp(fit[[y]], y = y, contr = contr, resc = scl, negexp = `2exp`, summarise = T) ) # compute the posterior summaries

    }
  ) %>% do.call( rbind.data.frame, . ),
  
  file = here("tables","cognition_posterior_summaries.csv"),
  sep = ";",
  row.names = F,
  quote = F
  
)


# VISUALISATION ----

## Interaction boxplots ----

d0 %>%
  
  filter( complete.cases(event) ) %>%
  select(event, group, osa, all_of(psych$variable) ) %>%
  pivot_longer(
    cols = all_of(psych$variable),
    values_to = "score",
    names_to = "test"
  ) %>%
  mutate(
    test = factor(
      unlist(
        sapply( 1:nrow(.), function(i) with( psych, label[variable == test[i]] ) ),
        use.names = F
      ),
      levels = psych$label,
      ordered = T
    ),
    Group = group,
    Event = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", event, perl = T),
    `OSA: ` = factor(
      if_else(osa == "+", "AHI ≥ 15", "AHI < 15"),
      levels = c("AHI < 15","AHI ≥ 15"),
      ordered = T
    )
  ) %>%
  
  ggplot() +
  aes(y = score, x = Group) +
  geom_boxplot( aes(fill = `OSA: `), width = .6, position = position_dodge(.7), linewidth = .75 ) +
  geom_dotplot( aes(fill = `OSA: `), binaxis = "y", stackdir = "center", position = position_dodge(.7), dotsize = 2.5 ) +
  labs( y = "Test score" ) +
  facet_grid2( test ~ Event, scales = "free_y", independent = F ) +
  scale_fill_manual( values = c("#64CDCC","#F9A729") ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","cognition_boxplots.jpg"),
  dpi = 300,
  width = 6,
  height = 12
)

