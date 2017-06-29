#title: "Cardiac ACMG-ClinVar Penetrance Shiny App"
#author: "James Diao, under the supervision of Arjun Manrai"
#date: "June 29, 2017"

#rsconnect::deployApp('/Users/jamesdiao/Documents/Kohane_Lab/2017-ACMG-penetrance/Cardiac_Penetrance_App')

#Set working directory to file folder
#outdir <- getSrcDirectory(function(dummy) {dummy})
#setwd(outdir)

library(shiny)
library(rhandsontable)
library(plotly)
library(dplyr)
library(shinysky)
library(tidyr)
library(tibble)
# Packages to install
#pkg_list <- c("shiny","rhandsontable","plotly","ggplot2","dplyr")
#installed <- pkg_list %in% installed.packages()[,"Package"]
#if (!all(installed))
#  install.packages(pkg_list[!installed])
#sapply(pkg_list, require, character.only = T)

# load allele frequencies and other info
#setwd("/Users/jamesdiao/Documents/Kohane_Lab/2017-ACMG-penetrance/Penetrance_App")
super.levels <- c("AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH","SAS")
freq_gnomad.calc.gene <- readRDS(file = "freq_gnomad.calc.gene.RDS")

sample_size <- c(12020, 17210, 5076, 9435, 12897, 63369, 15391, 3234, 138632) %>% 
  setNames(c("AFR", "AMR", "ASJ","EAS", "FIN", "NFE", "SAS", "OTH", "GNOMAD"))

helpmsg <- 'Click "Make Plots" to generate/refresh figures (make take up to 5 seconds).'

# Read in the given .csv
DF <- read.csv(file = "Cardiac_Literature_Prevalence_Estimates.csv", 
               header = TRUE, stringsAsFactors = F, na.strings = "NA")
abbrev <- DF$Short_Name

DF <- DF %>% select(Evaluate, Abbreviation, Short_Name, 
                    Inverse_Prevalence, Inv_Prev_Lower, Inv_Prev_Upper, Prev_Confidence, 
                    Inverse_CAF, Inv_CAF_Lower, Inv_CAF_Upper, CAF_Confidence)


sample_beta_dist <- function(shapes) {
    rbeta(n = 10^4, shape = shapes[1], shape2 = shapes[2])
}

betaExpert <- function(best, lower, upper, p = 0.95) {
  if (missing(best)) 
    stop("'best' is missing")
  if (missing(lower) & missing(upper)) 
    stop("at least 'lower' or 'upper' must be specified")
  if (!missing(lower)) 
    if (lower > best) stop("'lower' cannot be greater than 'best'")
  if (!missing(upper)) 
    if (upper < best) stop("'upper' cannot be smaller than 'best'")
  if (!missing(lower) & !missing(upper)) 
    if (lower > upper) stop("'lower' cannot be greater than 'upper'")
  f_mode <- function(x, mode, p, target) {
    return(sum((qbeta(p = p, shape1 = x, shape2 = (x*(1-mode) + 2*mode - 1)/mode) - target)^2))
  }
  f_mode_zero <- function(x, p, target) {
    return((qbeta(p = p, shape1 = 1, shape2 = x) - target)^2)
  }
  f_mode_one <- function(x, p, target) {
    return((qbeta(p = p, shape1 = x, shape2 = 1) - target)^2)
  }
  f_mean <- function(x, mean, p, target) {
    return(sum((qbeta(p = p, shape1 = x, shape2 = (x*(1-mean))/mean) - target)^2))
  }
  if (!missing(lower) & missing(upper)) {
    target <- lower
    p <- 1-p
  }
  else if (!missing(upper) & missing(lower)) {
    target <- upper
  }
  else if (!missing(upper) & !missing(lower)) {
    target <- c(lower, upper)
    p <- c(0, p) + (1-p)/2
  }
  #method = mode
  if (best == 0) {
    a <- 1
    b <- optimize(f_mode_zero, c(0, 1000), p = p, target = target)$minimum
  }
  else if (best == 1) {
    a <- optimize(f_mode_one, c(0, 1000), p = p, target = target)$minimum
    b <- 1
  }
  else {
    a <- optimize(f_mode, c(0, 1000), mode = best, p = p, target = target)$minimum
    b <- (a * (1 - best) + 2 * best - 1)/best
  }
  out <- c(alpha = a, beta = b)
  return(out)
}

sample_beta_obs <- function(freq, n, trials) {
  if (missing(trials)) trials <- 10^4
  obs_success <- freq * n
  obs_failure <- n-obs_success
  out <- sapply(1:length(freq), function(i)
    rbeta(n = trials, shape = obs_success[i] + 1, shape2 = obs_failure[i] + 1) %>% 
      setNames(1:trials)
  )
  colnames(out) <- names(freq)
  rownames(out) <- NULL
  return(out)
}


### User Interface
ui <- shinyUI(fluidPage(
  titlePanel("Estimation of Penetrance from Population Parameters for Six Cardiac Phenotypes"),
  fluidPage(
    tabsetPanel(
      tabPanel("Table", h5(),
        helpText("Double-click on the table to edit values. Uncheck boxes to remove from analysis."), h5(),
        helpText('Move to the "Heatmap" or "Barplot" tabs to view figures.' ),
        rHandsontableOutput("hot"),
        h2(),
        actionButton("run", " Make Plots", icon = icon("bar-chart"), styleclass = "success"),
        h2(), h2(), 
        p('When editing values for prevalence or case allele frequency (CAF), you may change:'),
        tags$ol(
          tags$li('Point estimate'),
          tags$li('Lower bound'),
          tags$li('Upper bound'),
          tags$li('Confidence Level (that the true value lies between these two bounds)')
        ),
        p('A Beta distribution is fitted to these parameters, with the mode at the point estimate and 
          x% of the density between the lower and upper bounds, with x determined by the confidence level. '),
        h2(), em("James Diao, under the supervision of Dr. Arjun Manrai (29 June 2017)")
      ),
      tabPanel("Heatmap", h2(),
        helpText(helpmsg),
        actionButton("run1", " Make Plots", icon = icon("bar-chart"), styleclass = "success"),
        h4(),
        plotOutput("heatmap", height = "350px", width = "650px"),
        h2()
      ),
      tabPanel("Density plot", h2(),
        helpText(helpmsg), 
        actionButton("run2", " Make Plots", icon = icon("bar-chart"), styleclass = "success"),
        h4(),
        plotOutput("densityplot", height = "800px", width = "800px"),
        h2()
      ),
      tabPanel("Barplot", h2(),
        helpText(helpmsg), 
        actionButton("run3", " Make Plots", icon = icon("bar-chart"), styleclass = "success"),
        h4(),
        plotOutput("barplot", height = "600px", width = "800px"),
        h2()
      )
    )
  )
))

### Server
server <- shinyServer(function(input, output, session) {
  
  values <- reactiveValues()
  
  ## Handsontable
  observe({
    if (!is.null(input$hot)) {
      values[["previous"]] <- isolate(values[["DF"]])
      DF = hot_to_r(input$hot)
    } else {
      if (is.null(values[["DF"]]))
        DF <- DF
      else
        DF <- values[["DF"]]
    }
    values[["DF"]] <- DF
  })
  
  output$hot <- renderRHandsontable({
    DF <- values[["DF"]]
    if (!is.null(DF))
      rhandsontable(DF, useTypes = TRUE, 
                    readOnly = FALSE, stretchH = "all")
  })
  
  observeEvent(input$run | input$run1 | input$run2 | input$run3, {
    
    # Set Parameters
    finalDF <- isolate(values[["DF"]])
    keep <- finalDF$Evaluate %>% as.logical
    finalDF <- finalDF[keep,]
    dataset <- 'gnomAD'
    method <- 'Gene'
    abbrev <- finalDF$Short_Name
    named.freqs <- freq_gnomad.calc.gene[keep,]
    useDF <- finalDF %>% filter(Evaluate)
    
    sample.freqs <- lapply(1:length(sample_size), function(x) {
      subset <- names(sample_size)[x]
      col <- ifelse(subset=='GNOMAD', 'AF_GNOMAD', sprintf("AF_GNOMAD_%s",subset))
      out <- sample_beta_obs(freq = named.freqs[,col], n = 2*sample_size[x])
      colnames(out) <- abbrev
      rownames(out) <- NULL
      return(out)
    }) %>% setNames(names(sample_size))
    
    sample.prev <- sapply(1:nrow(useDF), function(i) {
      betaExpert(best = 1/useDF$Inverse_Prevalence[i], lower = 1/useDF$Inv_Prev_Lower[i], 
                 upper = 1/useDF$Inv_Prev_Upper[i], p = useDF$Prev_Confidence[i]) %>% 
        sample_beta_dist() 
    }) 
      
    sample.CAF <- sapply(1:nrow(useDF), function(i) {
      betaExpert(best = 1/useDF$Inverse_CAF[i], lower = 1/useDF$Inv_CAF_Lower[i], 
                 upper = 1/useDF$Inv_CAF_Upper[i], p = useDF$CAF_Confidence[i]) %>% 
        sample_beta_dist()
    }) 
    
    colnames(sample.prev) <- colnames(sample.CAF) <- abbrev
    
    sample.penetrance <- lapply(sample.freqs, function(freqs) {
      (sample.prev * sample.CAF / freqs) %>% pmin(1)
    })
    
    
    cred_intervals <- sapply(sample.penetrance, function(set){
      apply(set, 2, function(col) quantile(col, 0.95)) %>% pmin(1)
    })
    
    
    densityplot <- data.frame(Ancestry = rep(names(sample.penetrance), 
                                         each = nrow(sample.penetrance[[1]])),
               do.call('rbind', sample.penetrance)
    ) %>% 
      gather(key = Disease, value = Frequency, -Ancestry) %>%
      filter(Ancestry != 'GNOMAD') %>% 
      ggplot(aes(Frequency, color = Ancestry, fill = Ancestry)) + 
      geom_density(alpha = 0.3, adjust = 2) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
      facet_wrap(~Disease, nrow = 6, scales = 'free_y') + 
      ggtitle('Penetrance Posterior Distributions') + ylab('Density')
    
    col_label <- rep(1:length(abbrev), each = nrow(sample.prev))
    penetrance <- data.frame(Disease = rep(abbrev, each = 5),
                             sapply(sample.penetrance, c) %>% 
                               apply(2, function(col) {
                                 sapply(1:length(abbrev), function(x) {
                                   quantile(col[col_label == x], c(0.025, 0.025, 0.5, 0.975, 0.975))
                                 })
                               })
    ) %>% gather(key = Ancestry, value = Penetrance, -Disease) %>%
      mutate(Penetrance = pmin(1,Penetrance))
    
    #5%-95% confidence interval boxplots
    barplot <- penetrance %>% sample_n(min(nrow(penetrance),100000)) %>% filter(Ancestry != 'GNOMAD') %>% 
           ggplot(aes(x = Disease, y = Penetrance, fill = Ancestry)) +
           geom_boxplot(outlier.alpha = 0) + coord_flip() + 
           ggtitle("95% Credible Intervals for Penetrance Posterior Distributions") + ylim(0,1)
    
    
    penetrance_data <- data.frame(Disease = rep(rownames(cred_intervals), times = ncol(cred_intervals)),
                                  cred_intervals %>% as.data.frame %>% gather(key = Subset, value = Penetrance))
    penetrance_data$Subset <- penetrance_data$Subset %>% 
      replace(penetrance_data$Subset == 'GNOMAD','TOTAL')
    
    heatmap <- ggplot(aes(x=Disease, y = Subset), data = penetrance_data) + coord_flip() + 
           geom_tile(aes(fill = Penetrance), color = 'white') + xlab("Disease") + ylab("Ancestry") +
           scale_fill_gradient(low='white',high = 'darkblue', na.value = "grey50",
                               breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.50","0.75","1.00"), limits =c(0,1)) + 
           ggtitle("95% Upper Penetrance Bound by Ancestry") + 
           theme_minimal() + theme(axis.ticks = element_blank()) + 
           annotate("segment", y=c(0.5,8.5,9.5), yend=c(0.5,8.5,9.5), 
                    x=0.5, xend = length(abbrev)+0.5) +
           annotate("segment", y=0.5, yend=9.5, 
                    x=c(0.5,length(abbrev)+0.5), 
                    xend = c(0.5,length(abbrev)+0.5))
    
    # Barplot
    #penetrance_data <- gather(penetrance_data, Subset, Penetrance, -Disease)
    #barplot <- ggplot(aes(x=Disease, y=Penetrance), data = penetrance_data) + 
    #       geom_boxplot(position = 'identity', coef = 0, na.rm = T) + 
    #       facet_wrap(~Subset, ncol = 2) + coord_flip() + xlab(NULL) + 
    #       ggtitle(sprintf("Penetrance by Ancestry (%s)", dataset)) + 
    #       theme(axis.text.y=element_text(size=6), 
    #             axis.text.x = element_text(angle = -20, hjust = 0.4))
    #barplot <- ggplotly(barplot, height = 1200, width = 800)
    
    # Heatmap
    #m <- list(l = 170, r = 150, b = -50, t = 100, pad = 5)
    #vals <- unique(scales::rescale(c(penetrance_data$Penetrance))) %>% sort
    #setord <- order(vals, decreasing = FALSE)
    #cols <- scales::col_numeric("Blues", domain = NULL)(vals)
    #colz <- setNames(data.frame(vals[setord], cols[setord]), NULL)
    #heatmap <- plot_ly(
    #  x = factor(c(super.levels,"Total"), levels = c(super.levels,"Total")),
    #  y = factor(sapply(abbrev, function(x) rep(x,5)) %>% as.vector, levels = abbrev[ord]),
    #  z = penetrance_init[pos,][ord,] %>% as.matrix %>% signif(3), 
    #  type = "heatmap", height = 800, colorscale = colz
    #) %>% layout(autosize = T, margin = m, 
    #  title = sprintf("%s Penetrance by Ancestry (%s)", position, dataset)) 
    output$heatmap <- renderPlot({ heatmap })
    output$densityplot <- renderPlot({ densityplot })
    output$barplot <- renderPlot({ barplot })
  }) #Closes observeEvent "Make Plot"

})

# Run the application 
shinyApp(ui = ui, server = server)
