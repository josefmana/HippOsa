# This scripts serves to explicitly state causal assumptions about relations between variables relevant to the project.

rm( list = ls() ) # clear environment

library(here)
library(tidyverse)
library(ggdag)
library(patchwork)

theme_set( theme_dag() ) # set-up theme for plotting
if( !dir.exists("figures") ) dir.create("figures") # prepare folder for figures

# set-up coordinates for nodes
coords <- data.frame(
  
  name = c("cog","AHI","PD","BS","Age","Edu","Sex","TIV"),
  x = c(2,0,1,1,0,2,1,2),
  y = c(0,0,1,2,-2,-2,-2,2)
  
)

# DAG with AHI ~ BS structure being confounded
dag1 <- dagify(
  
  cog ~ AHI + PD + BS + Age + Edu + Sex + TIV,
  AHI ~ PD + Age + Sex,
  PD ~ Age + Sex,
  BS ~ PD + Age + TIV,
  Edu ~ Sex,
  TIV ~ PD + Age + Sex,
  AHI ~~ BS,
  Sex ~~ BS,
  Edu ~~ BS,
  Edu ~~ TIV,
  coords = coords
  
)

# DAG with AHI causing BS structures volume
dag2 <- dagify(
  
  cog ~ AHI + PD + BS + Age + Edu + Sex + TIV,
  AHI ~ PD + Age + Sex,
  PD ~ Age + Sex,
  BS ~ AHI + PD + Age + TIV,
  Edu ~ Sex,
  TIV ~ PD + Age + Sex,
  Sex ~~ BS,
  Edu ~~ BS,
  Edu ~~ TIV,
  coords = coords
  
)

# plot it
( ggdag(dag1) | ggdag(dag2)  ) +  plot_annotation(tag_levels = "A")

# save it
ggsave( here("figures","dags.jpg"), dpi = 300, width = 12, height = 7 )
