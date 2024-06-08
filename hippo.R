# Computes regressions evaluating the effect of hippocampal volume on cognitive performance conditional on AHI and diagnosis.

rm( list = ls() ) # clear environment

# load packages
library(here)
library(openxlsx)
library(tidyverse)
library(ggdag)
library(patchwork)
library(car) # for easy to compute Type II and Type III sum of squares
library(performance)

# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )

# list all tests in the battery
# note: need to extract only those measures that are common across PD and CON, it ain't a level-II battery
psych <- read.csv("psychs.csv", sep = ";")
preds <- read.csv("brains.csv", sep = ";")

# read the data set
d1 <- read.xlsx( here("_raw","TabHippAm1.xlsx") )
d2 <- read.csv( here("_raw","20221120_redcap_export.csv"), sep = "," )


# PRE-PROCESSING  ----

# prepare a data frame for analyses
d0 <- left_join( d1, d2 %>% mutate( Study.ID = sub("-","",study_id) ) %>% select( Study.ID, all_of(psych$variable) ), by = "Study.ID" )

# format it for analyses
df <- d0 %>%
  
  # pre-process predictor terms
  mutate_if( is.character, as.factor ) %>%
  mutate( TIV = as.numeric( scale(`SBTIVmm^3`) ), across( contains("L_") | contains("R_"), ~ as.numeric( scale(.x) ) ) ) %>%
  mutate( GENDER = as.factor(GENDER), across( c("EDU.Y","AGE"), ~ as.numeric( scale(.x) ) ) ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(SUBJ) <- -contr.sum(2)/2 # CON = -0.5, PD = 0.5
    contrasts(AHI.F) <- -contr.sum(2)/2 # High = -0.5, Low = 0.5
    contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
  } ) %>%

  # scale outcome variables
  mutate(
#    across( "avlt_1_5", ~ as.numeric( scale(.x) ) ),
    across( all_of( subset(psych, domain == "Memory")$variable ), ~ as.numeric( scale(.x) ) ),
    across( all_of( subset(psych, domain != "Memory")$variable ), ~ -as.numeric( scale( log(.x) ) ) )
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
  TIV ~ PD + Sex,
  AHI ~~ Hippo,
  coords = coords

)

# DAG with AHI causing Hippocampal volume
dag2 <- dagify(

  cog ~ AHI + PD + Hippo + Age + Edu + Sex + TIV,
  AHI ~ PD + Age + Sex,
  PD ~ Age + Sex,
  Hippo ~ AHI + PD + Age + TIV,
  Edu ~ Sex,
  TIV ~ PD + Sex,
  coords = coords

)

# set-up theme
theme_set( theme_dag() )

# set-up function for adjustment sets plotting
dag_plot <- function(dag, iv, dv, eff = "direct", leg = "none") ggdag_adjustment_set(dag, exposure = iv, outcome = dv, effect = eff, type = "minimal", shadow = F) + theme(legend.position = leg)

# plot it
( ggdag(dag1) | ggdag(dag2)  ) /
( dag_plot(dag1, "PD", "cog", "total", "none") |  dag_plot(dag2, "PD", "cog", "total", "none") ) /
( dag_plot(dag1, c("PD","AHI"), "cog", "total", "none") | dag_plot(dag2, c("PD","AHI"), "cog", "total", "none") ) /
( dag_plot(dag1, c("PD","AHI","Hippo"), "cog", "direct", "none") | dag_plot(dag2, c("PD","AHI","Hippo"), "cog", "direct", "none") ) +
  plot_annotation(tag_levels = "A")

# save it
ggsave( here("figs","hippo_dags.jpg"), dpi = 300, width = 10, height = 15 )
  

# STAT MODELS ----

# fit a series of univariate regressions
fit <- lapply(
  
  setNames( c("descriptive","causal"), c("descriptive","causal") ),
  function(type)
    
    lapply(
      
      setNames( 1:nrow(psych), psych$label ),
      function(i)
        
        lapply(
          
          setNames( 1:nrow(preds),preds$predictor ),
          function(j) {
            
            if(type == "causal") form <- paste0( psych$LHS[i], preds$RHS[j] )
            else form <- paste0( psych$LHS[i], gsub(":"," * ",preds$predictor[j] ) )

            return( lm( formula = as.formula(form), data = df ) )
            
          }
        )
    )
)


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

