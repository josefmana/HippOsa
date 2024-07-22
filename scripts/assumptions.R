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

# plot it
( ggdag(dag1) | ggdag(dag2)  ) +  plot_annotation(tag_levels = "A")

# save it
ggsave( here("figures","hippo_dags.jpg"), dpi = 300, width = 12, height = 7 )
