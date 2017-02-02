#title: "ACMG-ClinVar Penetrance Shiny App"
#author: "James Diao, under the supervision of Arjun Manrai"
#date: "January 26, 2017"

#rsconnect::deployApp('/Users/jamesdiao/Documents/Kohane_Lab/2017-ACMG-penetrance/Penetrance_App')

#Set working directory to file folder
#outdir <- getSrcDirectory(function(dummy) {dummy})
#setwd(outdir)

library(shiny)
library(rhandsontable)
library(plotly)
library(dplyr)
library(shinysky)
library(tidyr)
# Packages to install
#pkg_list <- c("shiny","rhandsontable","plotly","ggplot2","dplyr")
#installed <- pkg_list %in% installed.packages()[,"Package"]
#if (!all(installed))
#  install.packages(pkg_list[!installed])
#sapply(pkg_list, require, character.only = T)

# load allele frequencies and other info
#setwd("/Users/jamesdiao/Documents/Kohane_Lab/2017-ACMG-penetrance/Penetrance_App")
super.levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
load(file = "disease_level_AF.RData")
#save(freq_1000g.calc.gene, freq_1000g.count.gene, freq_exac.calc.gene, freq_gnomad.calc.gene, 
#     freq_1000g.calc.mim, freq_1000g.count.mim, freq_exac.calc.mim, freq_gnomad.calc.mim, 
#     freq_1000g.calc.medgen, freq_1000g.count.medgen, freq_exac.calc.medgen, freq_gnomad.calc.medgen, 
#     file = "Penetrance_App/disease_level_AF.RData")

# Read in the given .csv
DF <- read.csv(file = "Literature_Prevalence_Estimates.csv", 
               header = TRUE, stringsAsFactors = F, na.strings = "NA") %>%
  select(Evaluate, Abbreviation, Short_Name, Inverse_Prevalence, Case_Allele_Frequency)

### User Interface
ui <- shinyUI(fluidPage(
  titlePanel("Estimation of Penetrance from Population Parameters"),
  h4("James Diao - January 26, 2017"), h2(),
  sidebarLayout(
    sidebarPanel(
      h3("Control Panel"),
      wellPanel(
        radioButtons("dataset", "Select Dataset", c("gnomAD","ExAC", "1000G")),
        radioButtons("method", "Select Phenotype ID Method", c("Gene","MIM", "MedGen")),
        radioButtons("sort_by", "Select Row Order", c("Penetrance", "Disease Name")),
        #radioButtons("position", "Select Heatmap Values", c("Max", "Mean")),
        sliderInput("ah_range", "Case Allele Frequency Range", 
                    min = 0, max = 1, value = c(0.01,1), step = 0.05)
      ),
      wellPanel(
        h4("Generate Plots"),
        actionButton("run", " Make Plots", icon = icon("bar-chart"), styleclass = "success"),
        h4()
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Table", h5(),
          helpText("Double-click on the table to edit values. Uncheck boxes to remove from analysis."), h5(),
          helpText('Move to the "Heatmap" or "Barplot" tabs to view figures.' ),
          rHandsontableOutput("hot"),
          h2()
        ),
        tabPanel("Heatmap", h2(),
          helpText('Click "Make Plots" to generate/refresh figures.' ), h2(),
          plotlyOutput("heatmap", height = "800px"),
          h2()
        ),
        tabPanel("Barplot", h2(),
          helpText('Click "Make Plots" to generate/refresh figures.' ), h2(),
          plotlyOutput("barplot", height = "800px"),
          h2()
        )
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
      rhandsontable(DF, useTypes = TRUE, stretchH = "all")
  })
  
  observeEvent(input$run, {
    
    # Set Parameters
    finalDF <- isolate(values[["DF"]])
    keep <- finalDF$Evaluate %>% as.logical
    finalDF <- finalDF[keep,]
    ah_low <- input$ah_range[1]
    ah_high <- input$ah_range[2]
    dataset <- input$dataset
    method <- input$method
    range <- input$range
    #position <- input$position
    position <- "Max" #"Mean"
    pos <- replace(c(F,F,F,F,F), ifelse(position == "Max", 5, 3), T)
    abbrev <- finalDF$Short_Name
    acmg_ah <- finalDF$Case_Allele_Frequency %>% as.numeric %>% pmax(ah_low) %>% pmin(ah_high)
    
    sapply(c(super.levels,"Total"), function(superpop){
      # Map of disease name to disease tags
      find <- paste0("AF_", toupper(dataset))
      freq_name <- sprintf("freq_%s.calc.%s", tolower(dataset), tolower(method))
      named.freqs <- eval(parse(text=freq_name))[keep,]
      if (superpop != "Total") 
        find <- paste(find, superpop, sep = "_")
      named.freqs <- named.freqs[,find] %>% unlist %>% setNames(abbrev)
      allelic.het <- c(ah_low, ah_low, mean(c(ah_low, ah_high)), ah_high, ah_high) %>% 
        rep(nrow(finalDF)) %>% matrix(nrow = nrow(finalDF), byrow = T)
      allelic.het[,3] <- acmg_ah
      
      # Matrix of penetrance values for allelic het range, capped at 1
      set_to_na <- function(row) { replace(row, is.infinite(row), NA) %>% pmin(1)}
      apply(allelic.het / finalDF$Inverse_Prevalence / named.freqs, 1, set_to_na) %>% unlist
    }) %>% as.data.frame -> penetrance_init
    
    mat_pen <- matrix(penetrance_init$Total, ncol = 5, byrow = T)
    
    if (input$sort_by == "Penetrance") {
      if (position == "Max") {
        ord <- order(mat_pen[,5], mat_pen[,3], mat_pen[,1], decreasing = T)
      } else {
        ord <- order(mat_pen[,3], mat_pen[,5], mat_pen[,1], decreasing = T)
      }
    } else {
      ord <- length(abbrev):1
    }
    penetrance_data <- data.frame(penetrance_init, 
                                  "Disease" = factor(sapply(abbrev, 
                                                            function(x) rep(x,5)) %>% as.vector,
                                                     levels = abbrev[ord]) ) 
    
    # Barplot
    penetrance_data <- gather(penetrance_data, Subset, Penetrance, -Disease)
    barplot <- ggplot(aes(x=Disease, y=Penetrance), data = penetrance_data) + 
           geom_boxplot(position = 'identity', coef = 0, na.rm = T) + 
           facet_wrap(~Subset, ncol = 2) + coord_flip() + xlab(NULL) + 
           ggtitle(sprintf("Penetrance by Ancestry (%s)", dataset)) + 
           theme(axis.text.y=element_text(size=6), 
                 axis.text.x = element_text(angle = -20, hjust = 0.4))
    barplot <- ggplotly(barplot, height = 800)
    
    # Heatmap
    m <- list(l = 170, r = 150, b = -50, t = 100, pad = 5)
    vals <- unique(scales::rescale(c(penetrance_data$Penetrance))) %>% sort
    setord <- order(vals, decreasing = FALSE)
    cols <- scales::col_numeric("Blues", domain = NULL)(vals)
    colz <- setNames(data.frame(vals[setord], cols[setord]), NULL)
    heatmap <- plot_ly(
      x = factor(c(super.levels,"Total"), levels = c(super.levels,"Total")),
      y = factor(sapply(abbrev, function(x) rep(x,5)) %>% as.vector, levels = abbrev[ord]),
      z = penetrance_init[pos,][ord,] %>% as.matrix %>% signif(3), 
      type = "heatmap", height = 800, colorscale = colz
    ) %>% layout(autosize = T, margin = m, 
      title = sprintf("%s Penetrance by Ancestry (%s)", position, dataset)) 
    output$heatmap <- renderPlotly({ heatmap })
    output$barplot <- renderPlotly({ barplot })
  }) #Closes observeEvent "Make Plot"

})

# Run the application 
shinyApp(ui = ui, server = server)
