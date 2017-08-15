#title: "Cardiac ACMG-ClinVar Penetrance Shiny App"
#author: "James Diao, under the supervision of Arjun Manrai"
#date: "July 3, 2017"

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
super.levels <- c("AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH","SAS","GNOMAD")
freq_gnomad.calc.gene <- readRDS(file = "freq_gnomad.calc.gene.RDS")
sample_size <- c(12020, 17210, 5076, 9435, 12897, 63369, 15391, 3234, 138632) %>% 
  setNames(super.levels)
helpmsg <- 'Click "Make Plots" (in the Table tab) to generate/refresh figures (takes up to 5 seconds for 10,000 samples).'

# Read in the given .csv
DF <- read.csv(file = "Cardiac_Literature_Prevalence_Estimates.csv", 
               header = TRUE, stringsAsFactors = F, na.strings = "NA") %>% 
        select(Phenotype, Prevalence, Prev_Lower, Prev_Upper, #Prev_CL, 
               CAF, CAF_Lower, CAF_Upper) #, CAF_CL)
abbrev <- DF$Phenotype

sample_beta_dist <- function(shapes, trials) {
    rbeta(n = trials, shape = shapes[1], shape2 = shapes[2])
}

eval_frac <- function(frac) {
  sapply(frac, function(f){
    eval(parse(text=as.character(f)))
  }) %>% setNames(NULL)
}

betaExpert <- function(best, lower, upper, p = 0.95) {
  if (missing(best)) stop("The point estimate is missing")
  if (missing(lower) | missing(upper)) stop("Both lower and upper bounds must be specified")
  if (lower > best) stop("The lower bound cannot be greater than the point estimate")
  if (upper < best) stop("The upper bound cannot be smaller than the point estimate")
  if (lower > upper) stop("The lower bound cannot be greater than the upper bound")
  target <- c(lower, upper)
  p <- c(0, p) + (1-p)/2
  a <- b <- 1
  if (best == 0)
    b <- optimize(function(x, p, target) {
      return((qbeta(p = p, shape1 = 1, shape2 = x) - target)^2)
    }, c(0, 1000), p = p, target = target)$minimum
  else if (best == 1)
    a <- optimize(function(x, p, target) {
      return((qbeta(p = p, shape1 = x, shape2 = 1) - target)^2)
    }, c(0, 1000), p = p, target = target)$minimum
  else {
    a <- optimize(function(x, mode, p, target) {
      return(sum((qbeta(p = p, shape1 = x, shape2 = (x*(1-mode) + 2*mode - 1)/mode) - target)^2))
    }, c(0, 1000), mode = best, p = p, target = target)$minimum
    b <- (a * (1 - best) + 2 * best - 1)/best
  }
  return(c(alpha = a, beta = b))
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
        helpText("Double-click on the table to edit values.  
        Move to the 'Heatmap', 'Density Plot' or 'Barplot' tabs to view figures."),
        p('When editing values for prevalence or case allele frequency (CAF), you may change:'),
        tags$ol(
          tags$li('Point estimate (must be between both bounds)'),
          tags$li('Lower bound'),
          tags$li('Upper bound')
        ),
        rHandsontableOutput("hot"),
        h2(),
        actionButton("run", " Make Plots", icon = icon("bar-chart"), styleclass = "success"),
        h2(), h2(), 
        p('A Beta distribution is fitted to these parameters, with the mode at the point estimate and 
          x% of the density between the lower and upper bounds, with x determined by the confidence level. 
          The default value of x (95%) of may be modified in Options. Default table values are computed 
          using quantiles from posterior distributions inferred from epidemiological studies. '),
        conditionalPanel("input.plotparams",
          plotOutput("prev_prior", height = "200px", width = "700px"),
          plotOutput("CAF_prior", height = "200px", width = "700px")
        ),
        h2(), em("James Diao, under the supervision of Dr. Arjun Manrai (3 July 2017)")
      ),
      tabPanel("Heatmap", h2(),
        helpText(helpmsg),
        plotOutput("heatmap", height = "350px", width = "650px"),
        h2()
      ),
      tabPanel("Density Plot", h2(),
        helpText(helpmsg),
        plotOutput("densityplot", height = "800px", width = "800px"),
        h2()
      ),
      tabPanel("Barplot", h2(),
        helpText(helpmsg), 
        plotOutput("barplot", height = "600px", width = "800px"),
        h2()
      ),
      tabPanel("Options", h2(),
        selectizeInput("samples", label = "Number of samples drawn from posterior distribution", choices = 10^c(3:7), 
                       selected = 10^4, multiple = F),
        sliderInput("prev_cf", "Confidence Level in Prevalence Bounds", 
                    min = 0.05, max = 1, value = 0.95, step = 0.05),
        sliderInput("CAF_cf", "Confidence Level in CAF Bounds", 
                    min = 0.05, max = 1, value = 0.95, step = 0.05),
        checkboxInput("plotparams", label = "Plot parameter distributions with table", value = FALSE, width = NULL)
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
    
    if(input$plotparams==T) {
      lapply(1:nrow(DF), function(i) {
        betaExpert(best = eval_frac(DF$Prevalence[i]), 
                   lower = eval_frac(DF$Prev_Lower[i]), 
                   upper = eval_frac(DF$Prev_Upper[i]), 
                   p = input$prev_cf #p = useDF$Prev_CL[i]
        ) %>% as.list %>% setNames(c('shape1','shape2')) %>% return()
        #"stat_function(fun = dbeta, args = y[[%s]], aes(color = abbrev[%s]), size = 1)"
      }) -> y
      
      lapply(1:nrow(DF), function(i) {
        betaExpert(best = eval_frac(DF$CAF[i]), 
                   lower = eval_frac(DF$CAF_Lower[i]), 
                   upper = eval_frac(DF$CAF_Upper[i]), 
                   p = input$CAF_cf #p = useDF$CAF_CL[i]
        ) %>% as.list %>% setNames(c('shape1','shape2')) %>% return()
        #"stat_function(fun = dbeta, args = y[[%s]], aes(color = abbrev[%s]), size = 1)"
      }) -> z
      
      alpha <- 0.2; n = 600; geom = 'area'
      prev_prior <- ggplot(data.frame(x = c(0, min(1,1.5*sapply(y, function(i) (i$shape1-1)/(i$shape1+i$shape2-2)) %>% max))), aes(x = x)) +
        stat_function(fun = dbeta, args = y[[1]], aes(color = abbrev[1], fill = abbrev[1]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = y[[2]], aes(color = abbrev[2], fill = abbrev[2]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = y[[3]], aes(color = abbrev[3], fill = abbrev[3]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = y[[4]], aes(color = abbrev[4], fill = abbrev[4]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = y[[5]], aes(color = abbrev[5], fill = abbrev[5]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = y[[6]], aes(color = abbrev[6], fill = abbrev[6]), geom = geom, alpha = alpha, n = n) +
        labs(color = "Disease", fill = "Disease") + ggtitle('Prevalence Distribution') + 
        xlab('Prevalence') + ylab('Density') 
      
      CAF_prior <- ggplot(data.frame(x = c(0, 1)), aes(x = x)) +
        stat_function(fun = dbeta, args = z[[1]], aes(color = abbrev[1], fill = abbrev[1]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = z[[2]], aes(color = abbrev[2], fill = abbrev[2]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = z[[3]], aes(color = abbrev[3], fill = abbrev[3]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = z[[4]], aes(color = abbrev[4], fill = abbrev[4]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = z[[5]], aes(color = abbrev[5], fill = abbrev[5]), geom = geom, alpha = alpha, n = n) + 
        stat_function(fun = dbeta, args = z[[6]], aes(color = abbrev[6], fill = abbrev[6]), geom = geom, alpha = alpha, n = n) +
        labs(color = "Disease", fill = "Disease") + ggtitle('CAF Distribution') + 
        xlab('CAF') + ylab('Density')
      
      output$prev_prior <- renderPlot({ prev_prior })
      output$CAF_prior <- renderPlot({ CAF_prior })
    }
    
  })
  
  output$hot <- renderRHandsontable({
    DF <- values[["DF"]]
    if (!is.null(DF)) {
      rhandsontable(DF, useTypes = TRUE, 
         readOnly = FALSE, stretchH = "all") %>%
        hot_validate_numeric(col = 2:7, min = 0, max = 1) %>% 
        hot_cols(renderer = "function (instance, td, row, col, prop, value, cellProperties) {
                 Handsontable.renderers.TextRenderer.apply(this, arguments);
                 if (col > 1 && col < 4) {
                 td.style.background = 'lavenderblush';
                 } else if (col > 4) {
                 td.style.background = 'aliceblue';
                 } else if (col == 1) {
                 td.style.background = 'mistyrose';
                 } else if (col == 4) {
                 td.style.background = 'lightcyan';
                 }}")
        #hot_col(col = c(5,9), format = "0%") %>% 
        #hot_table(customBorders = list(
          #list(
          #  range = list(from = list(row = 0, col = 1),
          #               to = list(row = 5, col = 3)),
          #  top = list(width = 2, color = "black"),
          #  bottom = list(width = 2, color = "black")), 
          #list(
          #  range = list(from = list(row = 0, col = 4),
          #               to = list(row = 5, col = 6)),
          #  top = list(width = 2, color = "black"),
          #  left = list(width = 2, color = "black"),
          #  bottom = list(width = 2, color = "black"),
          #  right = list(width = 2, color = "black")),
          #list(
          #  range = list(from = list(row = 0, col = 0),
          #               to = list(row = 5, col = 0)),
          #  top = list(width = 2, color = "black"),
          #  left = list(width = 2, color = "black"),
          #  bottom = list(width = 2, color = "black"),
          #  right = list(width = 2, color = "black"))
          #))
    }
  })
  
  observeEvent(input$run, {
    
    # Set Parameters
    useDF <- isolate(values[["DF"]])
    keep <- rep(TRUE, nrow(useDF)) #finalDF$Evaluate %>% as.logical
    useDF <- useDF[keep,]
    abbrev <- useDF$Phenotype
    named.freqs <- freq_gnomad.calc.gene[keep,]
    trials <- input$samples
    
    sample.freqs <- lapply(1:length(sample_size), function(x) {
      subset <- names(sample_size)[x]
      col <- ifelse(subset=='GNOMAD', 'AF_GNOMAD', sprintf("AF_GNOMAD_%s",subset))
      out <- sample_beta_obs(freq = named.freqs[,col], n = 2*sample_size[x], trials = trials)
      colnames(out) <- abbrev
      rownames(out) <- NULL
      return(out)
    }) %>% setNames(names(sample_size))
    
    sample.prev <- sapply(1:nrow(useDF), function(i) {
      betaExpert(best = eval_frac(useDF$Prevalence[i]), 
                 lower = eval_frac(useDF$Prev_Lower[i]), 
                 upper = eval_frac(useDF$Prev_Upper[i]), 
                 p = input$prev_cf #p = useDF$Prev_CL[i]
                 ) %>% 
        sample_beta_dist(trials = trials) 
    }) 
      
    sample.CAF <- sapply(1:nrow(useDF), function(i) {
      betaExpert(best = eval_frac(useDF$CAF[i]), 
                 lower = eval_frac(useDF$CAF_Lower[i]), 
                 upper = eval_frac(useDF$CAF_Upper[i]), 
                 p = input$CAF_cf #useDF$CAF_CL[i]
                 ) %>% 
        sample_beta_dist(trials = trials)
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
    #       ggtitle("Penetrance by Ancestry (gnomAD)") + 
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
    #  title = sprintf("%s Penetrance by Ancestry (gnomAD)", position)) 
    output$heatmap <- renderPlot({ heatmap })
    output$densityplot <- renderPlot({ densityplot })
    output$barplot <- renderPlot({ barplot })
  }) #Closes observeEvent "Make Plot"

})

# Run the application 
shinyApp(ui = ui, server = server)