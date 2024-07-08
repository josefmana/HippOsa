# Computes regressions evaluating the effect of hippocampal volume on cognitive performance conditional on AHI and diagnosis.

rm( list = ls() ) # clear environment

# load packages
library(here)
library(openxlsx)
library(tidyverse)
library(ggdag)
library(patchwork)
library(MatchIt)
library(car) # for easy to compute Type II and Type III sum of squares
library(performance)

# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )

# list all tests in the battery
# note: need to extract only those measures that are common across PD and CON, it ain't a level-II battery
psych <- read.csv("psychs.csv", sep = ";")


# UTILS ----

# print rounded number
rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) )

# delete leading zero
zerolead <- function(x, d = 3) ifelse( x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = T) )

# calculate and print mean and SD
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) )

# fit propensity scores weighted regression
fit_reg <- function(d, outcomes) lapply(
  
  setNames(outcomes, outcomes),
  function(y)
    
    lm(
      formula = as.formula( paste0(y," ~ SUBJ * AHI.F + SBTIV") ),
      data = d,
      weights = weights
    )
  
)

# extract linear regression model diagnostics
lm_dia <- function(fit, outcomes) sapply(
  
  outcomes,
  function(y)
    c( p_breusch_pagan = c( check_heteroscedasticity(fit[[y]]) ),
       n_cook = sum( check_outliers(fit[[y]]), na.rm = T ),
       p_shapiro_wilk = c( check_normality(fit[[y]]) ),
       p_durbin_watson = c( check_autocorrelation(fit[[y]]) )
    )
) %>%
  
  t() %>%
  as.data.frame() %>%
  mutate( heteroscedasticity = ifelse(p_breusch_pagan < .05, "!", ""), .after = p_breusch_pagan ) %>%
  mutate( nonnormality = ifelse(p_shapiro_wilk < .05, "!", ""), .after = p_shapiro_wilk ) %>%
  mutate( autocorrelation = ifelse(p_durbin_watson < .05, "!", ""), .after = p_durbin_watson ) %>%
  mutate( across( starts_with("p_"), zerolead ) ) %>%
  rownames_to_column("Outcome")

# extract ANOVAs and regression coefficients for interactions
lm_res <- function(fit, outcomes) sapply(
  
  outcomes,
  function(y)
    
    Anova(fit[[y]], type = 3)["SUBJ:AHI.F", ] %>%
    as.data.frame() %>%
    select( starts_with("F value"), Df, starts_with("Pr(>") ) %>%
    mutate( # add beta coefficient estimates
      beta = summary(fit[[y]])$coefficients[sub( "AHI.F", "AHI.F1", sub("SUBJ","SUBJ1","SUBJ:AHI.F") ), "Estimate"],
      b2.5 = confint(fit[[y]])[sub( "AHI.F", "AHI.F1", sub("SUBJ","SUBJ1","SUBJ:AHI.F") ), "2.5 %"],
      b97.5 = confint(fit[[y]])[sub( "AHI.F", "AHI.F1", sub("SUBJ","SUBJ1","SUBJ:AHI.F") ), "97.5 %"]
    )
  
) %>%
  
  t() %>%
  as.data.frame() %>%
  mutate_if(is.list, as.numeric) %>%
  rownames_to_column("structure")


# PRE-PROCESSING  ----

# read the data
d0 <- left_join(
  
  # read hippocampus data first
  read.xlsx( here("_raw","Tabhipposubfields_rhx.xlsx") ) %>% rename_with( ~ gsub("-", "_", .x) ),
  read.xlsx( here("_raw","Tabhipposubfields_lhx.xlsx") ) %>% rename_with( ~ gsub("-", "_", .x) ),
  by = "Study.ID",
  suffix = c("_rhx","_lhx") # add suffixes to differentiate betwen right and left sides

) %>%
  
  # add all other subscortical structures and cognitive data
  left_join(
    
    # start with the subcortical data as they include demography and such
    read.xlsx( here("_raw","Tabsubcortex1.xlsx") ) %>% rename_with( ~ gsub("-", "_", .x) ),
    . ,
    by = "Study.ID"

  ) %>%
  left_join(

    read.csv( here("_raw","20221120_redcap_export.csv"), sep = "," ) %>%
      mutate( Study.ID = sub("-","",study_id) ) %>%
      select( Study.ID, all_of(psych$variable) ),
    by = "Study.ID"

  ) %>%
  
  # correct formatting errors
  mutate( across( ends_with("_rhx"), as.numeric ) )

# format it for analyses
df <- d0 %>%
  
  # pre-process brain and demography variables
  mutate_if( is.character, as.factor ) %>%
  mutate(
    PD = if_else(SUBJ == "PD", 1, 0),
    GENDER = as.factor(GENDER),
    across( c("EDU.Y","AGE"), ~ as.numeric( scale(.x) ) )
  ) %>%
  mutate(
    SBTIV = as.numeric( scale(SBTIV) ),
    across(
      contains("L_") | contains("R_") |
        ends_with("_lhx") | ends_with("_rhx") |
        contains("Left_") | contains("Right_"),
      ~ as.numeric( scale(.x) )
    )
  ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(SUBJ) <- -contr.sum(2)/2 # CON = -0.5, PD = 0.5
    contrasts(AHI.F) <- -contr.sum(2)/2 # High = -0.5, Low = 0.5
    contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
  } ) %>%

  # scale neuropsychology variables (keeping MoCA non-transformed and will use binomial likelihood instead)
  mutate(
    across( all_of( subset(psych, domain == "Memory")$variable ), ~ as.numeric( scale(.x) ) ),
    across( all_of( subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable ), ~ -as.numeric( scale( log(.x) ) ) ),
    moca_max = 30
  )

# check for missing values
sapply( 1:nrow(df),  function(i) if ( sum( is.na(df[i, ]) ) > 0 ) df[i, ] ) %>%
  do.call( rbind.data.frame , . ) %>%
  mutate_if( is.numeric, round, 2 )


# CAUSAL ASSUMPTIONS ----

# set-up coordinates for nodes
coords <- data.frame(

  name = c("cog","AHI","PD","Hippo","Age","Edu","Sex","TIV"),
  x = c(2,0,1,1,0,2,1,2),
  y = c(0,0,1,2,-2,-2,-2,2)

)

# DAG with AHI ~ Hippocampus being confounded
dag1 <- dagify(
  
  cog ~ AHI + PD + Hippo + Age + Edu + Sex + TIV,
  AHI ~ PD + Age + Sex,
  PD ~ Age + Sex,
  Hippo ~ PD + Age + TIV,
  Edu ~ Sex,
  TIV ~ PD + Age + Sex,
  AHI ~~ Hippo,
  Sex ~~ Hippo,
  Edu ~~ Hippo,
  Edu ~~ TIV,
  coords = coords

)

# DAG with AHI causing Hippocampal volume
dag2 <- dagify(

  cog ~ AHI + PD + Hippo + Age + Edu + Sex + TIV,
  AHI ~ PD + Age + Sex,
  PD ~ Age + Sex,
  Hippo ~ AHI + PD + Age + TIV,
  Edu ~ Sex,
  TIV ~ PD + Age + Sex,
  Sex ~~ Hippo,
  Edu ~~ Hippo,
  Edu ~~ TIV,
  coords = coords

)

# set-up theme
theme_set( theme_dag() )

# set-up function for adjustment sets plotting
#dag_plot <- function(dag, iv, dv, eff = "direct", leg = "none") ggdag_adjustment_set(dag, exposure = iv, outcome = dv, effect = eff, type = "minimal", shadow = F) + theme(legend.position = leg)

# plot it
( ggdag(dag1) | ggdag(dag2)  ) +  plot_annotation(tag_levels = "A")
#( dag_plot(dag1, "PD", "cog", "total", "none") |  dag_plot(dag2, "PD", "cog", "total", "none") ) /
#( dag_plot(dag1, c("PD","AHI"), "cog", "total", "none") | dag_plot(dag2, c("PD","AHI"), "cog", "total", "none") ) /
#( dag_plot(dag1, c("PD","AHI","Hippo"), "cog", "direct", "none") | dag_plot(dag2, c("PD","AHI","Hippo"), "cog", "direct", "none") ) +

# save it
ggsave( here("figs","hippo_dags.jpg"), dpi = 300, width = 12, height = 7 )


## ---- Propensity score matching -----

# fit a logistic regression based propensity score model based on sex and age
# (add TIV as well?)
m.out <- matchit(
  formula = PD ~ AGE + GENDER,
  data = df,
  method = "full",
  distance = "glm"
)

# matching looks to have worked well with these data and variables
plot(m.out, type = "density", interactive = F)

# prepare a data set with propensity scores weights
df.matched <- match.data(m.out)


# RQ1: DOES AHI AFFECT SUBCORTICAL STRUCTURES VOLUME DIFFERENTLY IN PD PATIENTS COMPARED TO HEALTHY CONTROLS? ----

# extract subcortical areas' variable names
subcort <-
  read.xlsx( here("_raw","Tabsubcortex1.xlsx") ) %>%
  rename_with( ~ gsub("-", "_", .x) ) %>%
  select( contains("_"), contains("-") ) %>%
  select( !starts_with("c") ) %>%
  names()


## Interaction Boxplots ----

d0 %>%
  
  select(SUBJ, AHI.F, all_of( paste0("c",subcort) ) ) %>%
  pivot_longer(
    cols = all_of( paste0("c",subcort) ),
    values_to = "volume",
    names_to = "structure"
  ) %>%
  mutate(
    structure = factor(
      case_when(
        structure == "cL_amygdala" ~ "Left Amygdala", structure == "cR_amygdala" ~ "Right Amygdala",
        structure == "cL_hippocampus" ~ "Left Hippocampus", structure == "cR_hippocampus" ~ "Right Hippocampus",
        structure == "cLeft_Accumbens_area" ~ "Left Accumbens", structure == "cRight_Accumbens_area" ~ "Right Accumbens",
        structure == "cLeft_Caudate" ~ "Left Caudate", structure == "cRight_Caudate" ~ "Right Caudate",
        structure == "cLeft_Pallidum" ~ "Left Pallidum", structure == "cRight_Pallidum" ~ "Right Pallidum",
        structure == "cLeft_Putamen" ~ "Left Putamen", structure == "cRight_Putamen" ~ "Right Putamen",
        structure == "cLeft_Thalamus" ~ "Left Thalamus", structure == "cRight_Thalamus" ~ "Right Thalamus"
      ),
      levels = expand.grid(
        side = c("Left","Right"),
        struct = c("Amygdala","Hippocampus","Accumbens","Caudate","Pallidum","Putamen","Thalamus")
      ) %>%
        mutate( label = paste0(side," ",struct) ) %>%
        select(label) %>%
        unlist(use.names = F),
      ordered = T
    ),
    Diagnosis = if_else(SUBJ == "PD", "PD", "HC"),
    `AHI: ` = factor(
      if_else(AHI.F == "H", "Pathological", "Non-pathological"),
      levels = c("Pathological","Non-pathological"),
      ordered = T
    )
  ) %>%
  
  ggplot() +
  aes(y = volume, x = Diagnosis) +
  geom_boxplot( aes(fill = `AHI: `), width = .6, position = position_dodge(.7), linewidth = .75 ) +
  geom_dotplot( aes(fill = `AHI: `), binaxis = "y", stackdir = "center", position = position_dodge(.7), dotsize = 1.5 ) +
  labs( y = bquote("Standardized volume "(mm^3) ) ) +
  facet_wrap( ~ structure, ncol = 2, scales = "free_y" ) +
  scale_fill_manual( values = c("#64CDCC","#F9A729") ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figs","subcort_boxplots.jpg"),
  dpi = 300,
  width = 6,
  height = 10
)


## ---- Weighted linear regressions ----

# fit a series of propensity-score weighted regressions
fit_sc <- fit_reg(df.matched, subcort)
dia_sc <- lm_dia(fit = fit_sc, outcomes = subcort) # run diagnostics
res_sc <- lm_res(fit = fit_sc, outcomes = subcort) # extract ANOVAs and beta coefficients

# save the results
write.table(x = dia_sc, file = here("tabs","subcort_diagnostics.csv"), sep = ",", row.names = F, quote = F) # save diagnostics
write.table(x = res_sc, file = here("tabs","subcort_interactions.csv"), sep = ",", row.names = F, quote = F) # save ANOVAs


# RQ2: DOES AHI AFFECT HIPPOCAMPAL STRUCTURES VOLUME DIFFERENTLY IN PD PATIENTS COMPARED TO HEALTHY CONTROLS? ----

# extract hippocampal areas variable names
hippo <- names(d0)[ grepl( "rhx|lhx", names(d0) ) ]


## ---- Weighted linear regressions ----

fit_hippo <- fit_reg(df.matched, hippo)
dia_hippo <- lm_dia(fit = fit_hippo, outcomes = hippo) # run diagnostics
res_hippo <- lm_res(fit = fit_hippo, outcomes = hippo) # extract ANOVAs and beta coefficients

# save the results
write.table(x = dia_hippo, file = here("tabs","hippo_diagnostics.csv"), sep = ",", row.names = F, quote = F) # save diagnostics
write.table(x = res_hippo, file = here("tabs","hippo_interactions.csv"), sep = ",", row.names = F, quote = F) # save ANOVAs


# RQ3: DOES AHI/HIPPOCAMPAL STRUCTURES VOLUME DIFFERENTLY IN PD PATIENTS COMPARED TO HEALTHY CONTROLS? ----





# MODEL SUMMARIES ----

domcol <- c("#F9A729","#64CDCC","#A4A4D5","#5FBB68") # colours for cognitive domains
theme_set( theme_bw(base_size = 14) ) # set theme for plotting

# extract statistical tests results
tab2 <- lapply(
  
  # outer loop through inference types
  setNames( names(fit), names(fit) ),
  function(k)
    
    # middle loops through predictors
    lapply(
      
      with( preds, setNames(predictor,predictor) ),
      function(i)
        
        # inner loop for outcomes
        sapply(
          
          names(fit[[k]]),
          function(j) {
            
            print( paste0("outcome: ",j,", predictor: ",i,", inference type:",k) )
            
            # extract ANOVA results
            Anova(fit[[k]][[j]][[i]], type = 3)[i, ] %>%
              as.data.frame() %>%
              select( starts_with("F value"), Df, starts_with("Pr(>") ) %>%
              mutate( # add beta coefficient estimates
                beta = summary( fit[[k]][[j]][[i]] )$coefficients[ sub( "AHI.F", "AHI.F1", sub("SUBJ","SUBJ1",i) ), "Estimate" ],
                b2.5 = confint( fit[[k]][[j]][[i]] )[ sub( "AHI.F", "AHI.F1", sub("SUBJ","SUBJ1",i) ), "2.5 %" ],
                b97.5 = confint( fit[[k]][[j]][[i]] )[ sub( "AHI.F", "AHI.F1", sub("SUBJ","SUBJ1",i) ), "97.5 %" ]
              ) %>%
              return()
            
          }
        ) %>%
        
        t() %>%
        as.data.frame() %>%
        mutate( across( everything(), ~ unlist(.x, use.names = F) ) ) %>%
        rownames_to_column("Outcome") %>%
        mutate( Predictor = preds[ preds$predictor == i, "label2" ], .before = 1 ) %>%
        mutate( `Domain: ` = sapply( 1:nrow(.), function(j) psych[ psych$label == Outcome[j], "domain"] ), .after = 1 )
      
    ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    mutate( model_type = k ) %>%
    mutate( Predictor = factor(Predictor, levels = rev(preds$label2), ordered = T) ) %>%
    mutate( Outcome = factor(Outcome, levels = psych$label, ordered = T) ) %>%
    mutate( `Domain: ` = factor(`Domain: `, levels = unique(psych$domain), ordered = T) )
      
) %>%
  
  do.call( rbind.data.frame, . )

# write as csv
write.table( x = tab2, file = here("tabs","hippo_anovas.csv"), sep = ",", row.names = F, quote = F )

# plot beta weights
fig.beta <- lapply(
  
  with( tab2, setNames( unique(model_type), unique(model_type) ) ),
  function(i)
    
    subset(tab2, model_type == i) %>%
    ggplot() +
    aes(y = beta, ymin = b2.5, ymax = b97.5, x = Predictor, colour = `Domain: `) +
    geom_point(size = 5) +
    geom_linerange(linewidth = 1.4) +
    geom_hline( yintercept = 0, linetype = "solid", colour = "black", linewidth = .5 ) +
    scale_colour_manual( values = domcol ) +
    labs(y = "Standardized regression coefficent", x = NULL) + 
    coord_flip( ylim = c(-2,2) ) +
    facet_wrap( ~ Outcome ) +
    theme( legend.position = "bottom" )
  
)

# plot p-values
fig.pval <-lapply(
  
  with( tab2, setNames( unique(model_type), unique(model_type) ) ),
  function(i)
    
    subset(tab2, model_type == i) %>%
    ggplot() +
    aes(x = `Pr(>F)`, y = Predictor, fill = `Domain: `) +
    geom_bar(stat = "identity") +
    labs(x = "p-value", y = NULL) +
    geom_vline(xintercept = .05, linetype = "solid", linewidth = 1, colour = "black") +
    facet_wrap( ~ Outcome ) +
    scale_fill_manual( values = domcol ) +
    theme(legend.position = "bottom")
  
)

# plot interactions per outcome
fig.int <- lapply(
  
  with( psych, setNames(variable,label) ),
  function(i)
    
    d0 %>%
    select( SUBJ, AHI.F, all_of(i), all_of(sub("SUBJ:AHI.F:","c",preds$predictor)[-c(1:2)]) ) %>%
    pivot_longer( cols = starts_with("c"), names_to = "Structure", values_to = "Volume" ) %>%
    
    mutate( `Group: ` = case_when(SUBJ == "CON" ~ "HC", SUBJ == "PD" ~ "PD") ) %>%
    mutate( AHI = factor( case_when(AHI.F == "L" ~ "low AHI", AHI.F == "H" ~ "high AHI"), levels = paste0( c("low","high"), " AHI"), ordered = T ) ) %>%
    
    mutate( Structure = sub("c","SUBJ:AHI.F:",Structure) ) %>%
    mutate( Structure = sapply( 1:nrow(.), function(i) preds[ preds$predictor == Structure[i], "label2" ], USE.NAMES = F ) ) %>%
    mutate( Structure = factor( Structure, levels = preds$label2[-c(1:2)], ordered = T ) ) %>%
    
    ggplot() +
    aes(x = Volume, y = get(i), colour = `Group: `, fill = `Group: `) +
    geom_point(size = 3) +
    labs( x = bquote("Standardized volume"~("mm"^3) ), y = with( psych, label2[variable == i] ) ) +
    geom_smooth(method = "lm", linewidth = 1.5, alpha = .25) +
    ggh4x::facet_grid2( Structure ~ AHI, scales = "free_x", independent = "x" ) +
    theme_bw( base_size = 16 ) +
    theme(legend.position = "bottom")
  
)

# save them
for( i in names(fig.beta) ) ggsave( plot = fig.beta[[i]], file = here( "figs", paste0("hippo_betas_",i,".jpg") ), dpi = 300, width = 12.6, height = 13.3 )
for( i in names(fig.pval) ) ggsave( plot = fig.pval[[i]], file = here( "figs", paste0("hippo_pvalues_",i,".jpg") ), dpi = 300, width = 12.6, height = 13.3 )
for( i in names(fig.int) ) ggsave( plot = fig.int[[i]], file = here( "figs", paste0(i,"_intplot.jpg") ), dpi = 300, width = 9, height = 27 )


# BAYESIAN STATS ----

# Since all of TMT-A, TMT-B and GPT were potentially different conditionally on predictors, in this section I will model
# them better way in the Bayesian way.

# prepare data
d2 <- df %>%
  
