# Prepares interaction boxplots for brain or cognition variables

boxplots <- function(d0, df, help, scl, rt_vars, which = "brains") {
  
  # BRAINS ----
  if (which == "brains") {
    
    # prepare results for p-values reported in the boxplot
    p <-
      
      # fit model and extract contrasts
      fit_reg(d = df, outcomes = help$subco$scaled, X = "SUBJ * AHI.F", w = F) %>%
      mass() %>%
      
      # re-format for plotting
      mutate(
        p = paste0( zerolead(p.value), bh_adjust(p.value) ), # p-value labels
        y.position = unlist( # vertical displacement
          sapply(
            1:nrow(.),
            function(i)
              ifelse(
                term[i] == "PD - CON",
                max(df[ , y[i]]) + 3 * sd(df[ , y[i]]),
                max(df[ , y[i]]) + .5 * sd(df[ , y[i]])
              )
          )
        ),
        x = "Diagnosis",
        Diagnosis = if_else(SUBJ == "PD", "PD", "HC"),
        group1 = if_else(term == "PD - CON", "PD", Diagnosis),
        group2 = if_else(term == "PD - CON", "HC", Diagnosis),
        side = unlist(
          sapply(
            1:nrow(.),
            function(i)
              with( help$subco, side[scaled == y[i]] )
          ),
          use.names = F
        ),
        structure = factor(
          unlist(
            sapply( 1:nrow(.), function(i) with( help$subco, structure[scaled == y[i]] ) ),
            use.names = F
          ),
          levels = unique(help$subco$structure),
          ordered = T
        )
      )
    
    # plot it
    fig <- d0 %>%
      
      # prepare data
      select(SUBJ, AHI.F, all_of(help$subco$scaled) ) %>%
      pivot_longer( cols = all_of(help$subco$scaled), values_to = "volume", names_to = "struct" ) %>%
      mutate(
        side = unlist(
          sapply( 1:nrow(.), function(i) with( help$subco, side[scaled == struct[i]] ) ),
          use.names = F
        ),
        structure = factor(
          unlist(
            sapply( 1:nrow(.), function(i) with( help$subco, structure[scaled == struct[i]] ) ),
            use.names = F
          ),
          levels = unique(help$subco$structure),
          ordered = T
        ),
        Diagnosis = if_else(SUBJ == "PD", "PD", "HC"),
        `OSA: ` = factor(if_else(AHI.F == "H", "OSA+", "OSA-"), levels = c("OSA-","OSA+"), ordered = T)
      ) %>%
      
      # plotting proper
      ggplot() +
      aes(y = volume, x = Diagnosis) +
      geom_boxplot(
        aes(fill = `OSA: `),
        width = .6,
        position = position_dodge(.7),
        linewidth = .75
      ) +
      geom_dotplot(
        aes(fill = `OSA: `),
        binaxis = "y",
        stackdir = "center",
        position = position_dodge(.7),
        dotsize = 1.5
      ) +
      labs( y = bquote("Standardized volume "(mm^3) ) ) +
      facet_grid2( structure ~ side, scales = "free_y", independent = "y" ) +
      scale_fill_manual( values = c("#64CDCC","#F9A729") ) +
      stat_pvalue_manual( # p-values for simple main effects
        subset(p, term != "PD - CON"),
        label = "p",
        size = 3,
        position = position_dodge(.7),
        remove.bracket = T
      ) +
      stat_pvalue_manual( # p-values for interactions
        subset(p, term == "PD - CON"),
        label = "p",
        size = 3,
        tip.length = .1,
        vjust = -.33
      ) +
      scale_y_continuous( expand = expansion(mult = c(0, .19) ) ) +
      theme_bw(base_size = 10) +
      theme(
        legend.position = "bottom",
        panel.grid = element_blank()
      )
    
    # save it
    #ggsave(
    #  plot = last_plot(),
    #  filename = here("figures","subcortical_boxplots.jpg"),
    #  dpi = 300,
    #  width = 6,
    #  height = 10.5
    #)
    
  } else if (which == "cognition_1") {
    
    # COGNITION: PD/OSA interaction ----
    
    # prepare results for p-values reported in the boxplot
    p <-
      
      # fit model and extract contrasts
      fit_reg(d = df, outcomes = help$psych$variable, X = "SUBJ * AHI.F", w = F) %>%
      mass() %>%
      
      # re-format for plotting
      mutate(
        p = paste0( zerolead(p.value), bh_adjust(p.value) ), # p-value labels
        y.position = unlist(
          sapply(
            1:nrow(.),
            function(i) case_when(
              term[i] == "PD - CON" & !(y[i] %in% rt_vars) ~ max(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] + 2 * scl[y[i],"SD"],
              term[i] != "PD - CON" & !(y[i] %in% rt_vars) ~ max(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] + .3 * scl[y[i],"SD"],
              term[i] == "PD - CON" & y[i] %in% rt_vars ~ exp( -( min(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] - 1 * scl[y[i],"SD"]) ),
              term[i] != "PD - CON" & y[i] %in% rt_vars ~ exp( -( min(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] - .2 * scl[y[i],"SD"] ) )
            )
          )
        ), # vertical displacement
        x = "Diagnosis",
        Diagnosis = if_else(SUBJ == "PD", "PD", "HC"),
        group1 = if_else(term == "PD - CON", "PD", Diagnosis),
        group2 = if_else(term == "PD - CON", "HC", Diagnosis),
        test = factor(
          unlist(
            sapply( 1:nrow(.), function(i) with( help$psych, label[variable == y[i]] ) ),
            use.names = F
          ),
          levels = help$psych$label,
          ordered = T
        )
      )
    
    # plot it
    fig <- d0 %>%
      
      select(SUBJ, AHI.F, all_of(help$psych$variable) ) %>%
      pivot_longer(cols = all_of(help$psych$variable), values_to = "score", names_to = "test") %>%
      mutate(
        test = factor(
          unlist(
            sapply( 1:nrow(.), function(i) with( help$psych, label[variable == test[i]] ) ),
            use.names = F
          ),
          levels = help$psych$label,
          ordered = T
        ),
        Diagnosis = if_else(SUBJ == "PD", "PD", "HC"),
        `OSA: ` = factor(if_else(AHI.F == "H", "OSA+", "OSA-"), levels = c("OSA-","OSA+"), ordered = T)
      ) %>%
      
      ggplot() +
      aes(y = score, x = Diagnosis) +
      geom_boxplot( aes(fill = `OSA: `), width = .6, position = position_dodge(.7), linewidth = .75 ) +
      geom_dotplot(
        aes(fill = `OSA: `),
        binaxis = "y",
        stackdir = "center",
        position = position_dodge(.7),
        dotsize = 1.5
      ) +
      labs( y = "Test score" ) +
      facet_wrap( ~ test, scales = "free", nrow = 5 ) +
      scale_fill_manual( values = c("#64CDCC","#F9A729") ) +
      stat_pvalue_manual( # p-values for simple main effects
        subset(p, term != "PD - CON"),
        label = "p",
        size = 3,
        position = position_dodge(.7),
        remove.bracket = T
      ) +
      stat_pvalue_manual( # p-values for interactions
        subset(p, term == "PD - CON"),
        label = "p",
        size = 3,
        tip.length = .1,
        vjust = -.33
      ) +
      scale_y_continuous( expand = expansion(mult = c(0, .19) ) ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom",
        panel.grid = element_blank()
      )
    
    # save it
    #ggsave(
    #  plot = last_plot(),
    #  filename = here("figures","cognition_boxplots.jpg"),
    #  dpi = 300,
    #  width = 9,
    #  height = 9
    #)
    
  } else if (which == "cognition_2") {
    
    # COGNITION: PD main effects ----
    
    # prepare results for p-values reported in the boxplot
    p <-
      
      fit_reg(d = df, outcomes = help$psych$variable, X = "SUBJ * AHI.F", w = F) %>%
      mass(type = "fullest") %>%
      filter(term == "SUBJ") %>%
      
      # re-format for plotting
      mutate(
        p = paste0( zerolead(p.value), bh_adjust(p.value) ), # p-value labels
        y.position = unlist(
          sapply(
            1:nrow(.),
            function(i) case_when(
              is.na(AHI.F[i]) & !(y[i] %in% rt_vars) ~ max(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] + 1 * scl[y[i],"SD"],
              !is.na(AHI.F[i]) & !(y[i] %in% rt_vars) ~ max(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] + .2 * scl[y[i],"SD"],
              is.na(AHI.F[i]) & y[i] %in% rt_vars ~ exp( -( min(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] - .5 * scl[y[i],"SD"]) ),
              !is.na(AHI.F[i]) & y[i] %in% rt_vars ~ exp( -( min(df[ , y[i]], na.rm = T) * scl[y[i],"SD"] + scl[y[i],"M"] - .1 * scl[y[i],"SD"] ) )
            )
          )
        ), # vertical displacement
        x = "OSA",
        OSA = if_else(AHI.F == "H", "OSA+", "OSA-"),
        group1 = if_else(is.na(OSA), "OSA-", OSA),
        group2 = if_else(is.na(OSA), "OSA+", OSA),
        test = factor(
          unlist(
            sapply( 1:nrow(.), function(i) with( help$psych, label[variable == y[i]] ) ),
            use.names = F
          ),
          levels = help$psych$label,
          ordered = T
        )
      )
    
    # plot it
    fig <- d0 %>%
      
      select(SUBJ, AHI.F, all_of(help$psych$variable) ) %>%
      pivot_longer(cols = all_of(help$psych$variable), values_to = "score", names_to = "test") %>%
      mutate(
        test = factor(
          unlist(
            sapply( 1:nrow(.), function(i) with( help$psych, label[variable == test[i]] ) ),
            use.names = F
          ),
          levels = help$psych$label,
          ordered = T
        ),
        `Diagnosis: ` = if_else(SUBJ == "PD", "PD", "HC"),
        OSA = factor(if_else(AHI.F == "H", "OSA+", "OSA-"), levels = c("OSA-","OSA+"), ordered = T)
      ) %>%
      
      ggplot() +
      aes(y = score, x = OSA) +
      geom_boxplot( aes(fill = `Diagnosis: `), width = .6, position = position_dodge(.7), linewidth = .75 ) +
      geom_dotplot(
        aes(fill = `Diagnosis: `),
        binaxis = "y",
        stackdir = "center",
        position = position_dodge(.7),
        dotsize = 1.5
      ) +
      labs(y = "Test score", x = NULL) +
      facet_wrap( ~ test, scales = "free", nrow = 5 ) +
      scale_fill_manual( values = c("deepskyblue","red2") ) +
      stat_pvalue_manual( # p-values for simple main effects
        subset(p, !is.na(OSA) ),
        label = "p",
        size = 3,
        position = position_dodge(.7),
        remove.bracket = T
      ) +
      stat_pvalue_manual( # p-values for interactions
        subset(p, is.na(OSA) ),
        label = "p",
        size = 3,
        tip.length = 0,
        vjust = -.2
      ) +
      scale_y_continuous( expand = expansion(mult = c(0, .19) ) ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom",
        panel.grid = element_blank()
      )
    
    # save it
    #ggsave(
    #  plot = last_plot(),
    #  filename = here("figures","cognition_boxplots_group_contrast.jpg"),
    #  dpi = 300,
    #  width = 9,
    #  height = 9
    #)
    
  }
  
  return(fig)
  
}