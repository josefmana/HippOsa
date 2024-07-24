# Computes regressions evaluating:
#
# RQ6) the simple effect OSA on MDS-UPDRS III scores in de novo PD
#
# addressed via a multivariate regression of MDS-UPDRS III subscores on OSA status, age and sex
# after evaluating MDS-UPDRS III scales internal consistency,
# and a Bayesian IRT model of the full MDS.UPDRS III scale for item-level inference.


library(brms)
library(priorsense)



# read the data
d1 <-
  read.csv( here("_data","mds_updrs_iii.csv"), sep = "," ) %>%
  left_join(
    read.csv( here("_data","primary_dataset.csv"), sep = "," ) %>%
      filter(event == "enrollment")
  ) %>%
  mutate(
    AHI.F = as.factor(AHI.F),
    RESP = as.ordered(response_num),
    SEX = factor( case_when(GENDER == 1 ~ "male", GENDER == 0 ~ "female") ),
    ITEM = unlist(sapply( 1:nrow(.), function(i) paste( strsplit(item[i], ".", fixed = T)[[1]][1:2], collapse = "." ) ) ),
    NAME = unlist(sapply( 1:nrow(.), function(i) paste( strsplit(item[i], ".", fixed = T)[[1]][-c(1:2)], collapse = "." ) ) )
  ) %>%
  within( . , contrasts(AHI.F) <- -contr.sum(2)/2 )



# ITEM RESPONSE MODELLING ----

# prepare formulas
formula1 <- list(
  
  `1PL_fixed` = bf( RESP ~ 1 + AHI.F +  (1 + AHI.F | ITEM) + (1 | Study.ID) ),
  `1PL_flex` = bf( RESP | thres(gr = ITEM) ~ 1 + AHI.F + AHI.F:ITEM + (1 | Study.ID) )
  
)

# re-fit the models for TMT-B and GPT via brms with differing variance for the groups
prior1 <- list(
  
  `1PL_fixed` =
    prior("normal(0, 3)", class = "Intercept") +
    prior("normal(0, 3)", class = "b") +
    prior("normal(0, 3)", class = "sd", group = "ITEM") +
    prior("normal(0, 3)", class = "sd", group = "Study.ID") +
    prior("lkj(1)", class = "cor"),
  
  `1PL_flex` =
    prior("normal(0, 3)", class = "Intercept") +
    prior("normal(0, 3)", class = "b") +
    prior("normal(0, 3)", class = "sd", group = "Study.ID")
  
)

# fit the models
fit1 <- lapply(
  
  setNames( names(prior1), names(prior1) ),
  function(k)
    
    brm( formula = formula1[[k]],
         data = d1,
         family = cumulative(link = "logit"),
         prior = prior1[[k]],
         seed = 87542
    )
  
)

# compare models' expected predictive performance
with( fit1, loo(`1PL_fixed`, `1PL_flex`) ) # selecting the fixed model for further description of patterns in the data (based on the model)


#### ---- plots ----

# some posterior predictive checks
lapply(
  
  c("ITEM","Study.ID"),
  function(k) {
    
    pp_check(fit1$`1PL_fixed`, type = "bars_grouped", ndraws = NULL, group = k) +
      theme(legend.position = "bottom")
    
    ggsave(
      plot = last_plot(),
      filename = here( "figures", paste0("motor_ppc_",tolower(k),".jpg") ),
      dpi = 300,
      width = ifelse(k == "ITEM", 10, 12),
      height = ifelse(k == "ITEM", 10, 12)
    )
    
  }
)

# plot item parameters
as_draws_df(fit1$`1PL_fixed`) %>%
  select( starts_with("r_ITEM") ) %>%
  pivot_longer(everything(), names_to = "coefficient", values_to = "value") %>%
  mutate(
    coefficient = sub("r_ITEM[", "", sub("]", "", coefficient, fixed = T), fixed = T),
    item = sub(",.*", "", coefficient),
    parameter = sub(".*,", "", coefficient)
  ) %>%
  group_by(item, parameter) %>%
  summarise(
    Md = median(value),
    Q2.5 = quantile(value, probs = .025),
    Q97.5 = quantile(value, probs = .975)
  ) %>%
  ungroup() %>%
  mutate( `Parameter: ` = factor(parameter, levels = c("Intercept","AHI.F1"), ordered = T) ) %>%
  
  ggplot() +
  aes(x = Md, xmin = Q2.5, xmax = Q97.5, y = reorder(item, Md, decreasing = T), colour = `Parameter: `) +
  geom_point( size = 4, position = position_dodge(width = -.75) ) +
  geom_pointrange( linewidth = 1.25, position = position_dodge(width = -.75) ) +
  scale_colour_manual( values = c("navyblue","red") ) +
  labs(
    y = NULL,
    x = "Parameter value",
    title = "Inverse difficulty item parameters of MDS-UPDRS III",
    subtitle = "The parameters are based on Item Response Theory logistic ordered model\nwith fixed item thresholds and hierarchical priors"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5)
  )

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","motor_item_parameters.jpg"),
  dpi = 300,
  width = 6,
  height = 11
)


# Show descriptively the profile PD-OSA+ vs PD-OSA- in MDS-UPDRS III