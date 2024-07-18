# Computes regressions evaluating:
#
# RQ4) the moderating effect OSA on longitudinal general cognitive performance in PD
# RQ5) the moderating effect OSA on re-test cognitive profile in PD
#
# Both RQs are addressed via Bayesian GLMMs


rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

# load packages
library(here)
library(tidyverse)
library(brms)
library(priorsense)
library(bayestestR)

sapply( c("figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders
psych <- read.csv( here("helpers","psychs.csv"), sep = ";") # read helper files with variable names


# UTILS ----


# PRE-PROCESSING  ----

# read the data
d0 <-
  read.csv( here("_data","primary_dataset.csv"), sep = "," ) %>%
  mutate(
    id = Study.ID,
    sex = factor( case_when(GENDER == 1 ~ "male", GENDER == 0 ~ "female") ),
    pd = factor( case_when(SUBJ == "PD" ~ 1, SUBJ == "CON" ~ 0) ),
    osa = factor( case_when(AHI.F == "H" ~ 1, AHI.F == "L" ~ 0) ),
    age = AGE,
    edu = EDU.Y
  )

# prepare MoCA data set for RQ4
d1 <- d0 %>%
  select(id, event, sex, pd, osa, age, edu, moca) %>%
  filter( complete.cases(moca) ) %>%
  mutate( event_num = as.numeric( as.factor(event) ) - 1 )

# prepare battery data set for RQ5
d2 <- d0 %>%
  select(id, event, sex, pd, osa, age, edu, all_of(psych$variable[-13]) ) %>%
  filter( event %in% c("enrollment","followup_4") ) %>%
  mutate(
    across(
      all_of( subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable ),
      ~ -log(.x)
    )
  )


# RQ4: DO OSA+/- PATIENTS SHOW DISTINCTIVE LONGITUDINAL COGNITIVE TRAJECTORY COMPARED TO HEALTHY CONTROLS? ----

# plot number of assessments per participant
table(d1$id) %>%
  
  as.data.frame() %>%
  mutate( group =  if_else( grepl("CON",Var1), "HC", "PD" ) ) %>%
  
  ggplot() +
  aes(y = reorder(Var1, Freq), x = Freq, fill = group) +
  geom_bar( stat = "identity" ) +
  labs(y = NULL, x = "Number of assessments") +
  scale_x_continuous( breaks = seq(0,9,1), labels = seq(0,9,1) ) +
  facet_wrap( ~group, scales = "free_y" ) +
  theme_bw() +
  theme(legend.position = "none")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figs","cognition_longitudinal_assessments.jpg"),
  dpi = 300,
  width = 8,
  height = 8
)

## ---- longitudinal GLM of MoCA ----

# set-up priors
prior1 <- c(
  
  prior("normal(0, 3)", class = "Intercept"),
  prior("normal(0, 3)", class = "b"),
  prior("normal(0, 3)", class = "sd", group = "id"),
  prior("lkj(2)", class = "cor")
  
)

# fit it
fit1 <- brm(
  
  moca | trials(30) ~ 1 + (1 + event_num | id) + event_num * pd * osa,
  data = d1,
  prior = prior1,
  family = binomial()
  
)

# check for prior sensitivity
powerscale_sensitivity(fit1)

# do some posterior checks
pp_check(fit1, type = "stat_2d", ndraws = NULL )
pp_check(fit1, type = "bars_grouped", ndraws = NULL, group = "event_num" )
pp_check(fit1, type = "bars_grouped", ndraws = NULL, group = "pd" )
pp_check(fit1, type = "bars_grouped", ndraws = NULL, group = "osa" )

# extract and save posterior summaries
describe_posterior(fit1, test = "p_direction") %>%
  as.data.frame() %>%
  write.csv(file = here("tabs","cognition_longitudinal.csv"), sep = ",", row.names = F, quote = F)


# RQ5: DO OSA+/- PATIENTS SHOW DISTINCTIVE RE-TEST COGNITIVE PROFILE COMPARED TO HEALTHY CONTROLS? ----

# amount of data points available
lapply(
  
  setNames( psych$variable[-13], psych$variable[-13] ),
  function(y)
    
    table( subset(d2, complete.cases( get(y) ) )$id) %>%
    as.data.frame() %>%
    mutate( Group =  if_else( grepl("CON",Var1), 0, 1 ) ) %>%
    mutate( OSA = unlist( sapply( 1:nrow(.), function(j) with(d1, unique( osa[id == Var1[j]] ) ) ) ) ) %>%
    select(OSA, Group, Freq) %>%
    table()
  
)

# extract scaling values, i.e.,
# enrollment full sample means and SDs
scl <- sapply(
  
  c(psych$variable[-13], "age", "edu"),
  function(i)
    
    with(subset(d2, event == "enrollment"), c(M = mean( get(i), na.rm = T ), SD = sd( get(i), na.rm = T ) ) )
  
) %>% t()

# scale the variables
d2 <- d2 %>%
  
  mutate(
    across(
      all_of( rownames(scl) ),
      ~ (.x - scl[cur_column(), "M"]) / scl[cur_column(), "SD"],
      .names = "{.col}_sc"
    )
  )

# set-up priors
prior2 <- c(
  
  prior("normal(0, 1)", class = "Intercept"),
  prior("normal(0, 1)", class = "b"),
  prior("exponential(1)", class = "sd", group = "id"),
  prior("exponential(1)", class = "sigma")
  
)

# fit the models
fit2 <- lapply(
  
  with( psych, setNames( paste0(psych$variable[-13],"_sc"), psych$variable[-13]) ),
  function(y)
    
    brm(
      formula = bf( paste0(y," ~ 1 + event * pd * osa + (1 | id)") ),
      data = d2,
      prior = prior2,
      family = gaussian(link = "identity"),
      seed = 87542
    )
  
)

do.call(
  rbind.data.frame,
  lapply(
    
    names(fit2),
    function(y)
      
      describe_posterior(fit2[[y]], test = "p_direction") %>%
      as.data.frame() %>%
      mutate(Outcome = y, .before = 1)
    
  )
)