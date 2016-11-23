#title: "ACMG-ClinVar Penetrance Shiny App"
#author: "James Diao, under the supervision of Arjun Manrai"
#date: "November 23, 2016"

# Set working directory to file folder
outdir <- getSrcDirectory(function(dummy) {dummy})
setwd(outdir)

# Output name for saving modified csv
outfilename <- "Modified_Estimates"
# Packages to install
pkg_list <- c("knitr","shiny","shinysky","rhandsontable","plotly",
              "ggplot2","ggrepel","tibble","tidyr","dplyr")
installed <- pkg_list %in% installed.packages()[,"Package"]
if (!all(installed))
  install.packages(pkg_list[!installed])
sapply(pkg_list, require, character.only = T)

# load allele frequencies and other info
# setwd("/Users/jamesdiao/Documents/Kohane_Lab/2016-paper-ACMG-penetrance/Shiny_App")
load(file = "disease_level_AF.RData")
#save(freq_1000g.calc.gene, freq_1000g.count.gene, freq_exac.calc.gene, 
#     super.levels, file = "disease_level_AF.RData")

# Read in the given .csv
DF <- read.csv(file = "Literature_Prevalence_Estimates.csv", 
               header = TRUE, stringsAsFactors = F, na.strings = "NA") %>%
  select(Evaluate, Abbreviation, Short_Name, Inverse_Prevalence, Allelic_Heterogeneity)

# javascript function for log-scale ah_range
JS.logify <- "function logifySlider (sliderId) {
  $('#'+sliderId).data('ionRangeSlider').update({
  'prettify': function (num) { return (Math.pow(10, num)); }
  })
}"
JS.onload <- "$(document).ready(function() {
  setTimeout(function() {
  logifySlider('ah_range')
}, 5)}) "

### User Interface
ui <- shinyUI(fluidPage(
  titlePanel("Estimation of Penetrance from Population Parameters"),
  sidebarLayout(
    sidebarPanel(
      h3("Control Panel"),
      wellPanel(
        h4("Disease-Level Allele Frequencies"),
        radioButtons("dataset", "Select Dataset", 
                     c("ExAC", "1000 Genomes"))
      ),
      wellPanel(
        h4("Penetrance Estimate on Heatmap"),
        radioButtons("position", "Select Display", c("Max", "Mean"))
      ),
      wellPanel(
        h4("Allelic Heterogeneity Range"),
        tags$head(tags$script(HTML(JS.logify))),
        tags$head(tags$script(HTML(JS.onload))),
        sliderInput("ah_range", "Draws max and min plausible bounds", 
                    min = -4, max = 0, value = c(-3,0), step = 1)
      ),
      wellPanel(
        h4("Imputed Prevalence Range"),
        sliderInput("range", "Ex: With range 5, a point prevalence of 0.3 becomes the prevalence range: 0.1-0.5", 
                    min = 1, max = 30, value = 5, step = 1)
      ),
      wellPanel(
        h4("Generate Plots with Modified Table"),
        actionButton("run", " Make Plots", icon = icon("bar-chart"), styleclass = "success")
      ),
      wellPanel(
        h4("Save Modified Table"), 
        div(class='row', 
            div(class="col-sm-6", 
                actionButton("save", "Save Table", styleclass = "success"))
        ),
        shinyalert("saved", click.hide = TRUE, auto.close.after = 3)
      )
    ),
    
    mainPanel(
      wellPanel(
        helpText(sprintf("Working Directory: %s", outdir)),
        helpText("Double-click on the table to edit values. Uncheck boxes to remove from analysis.")
      ),
      rHandsontableOutput("hot"),
      br(),
      busyIndicator("In Progress: Please Wait", wait = 500),
      h2("Barplot"),
      plotOutput("barplot", height = "800px"),
      h2("Heatmap"),
      plotlyOutput("heatmap", height = "800px")
      
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
  
  ## Save 
  observeEvent(input$save, {
    finalDF <- isolate(values[["DF"]])
    write.csv(finalDF, file=file.path(outdir, sprintf("%s.csv", outfilename)), row.names = F, quote = F)
    showshinyalert(session, "saved", sprintf("File saved: %s.csv", outfilename), styleclass = "success")
  })
  
  observeEvent(input$run, {
    # Set Parameters
    finalDF <- isolate(values[["DF"]])
    keep <- finalDF$Evaluate %>% as.logical
    freq_1000g.calc <- freq_1000g.calc.gene[keep,]
    freq_1000g.count <- freq_1000g.count.gene[keep,]
    freq_exac.calc <- freq_exac.calc.gene[keep,]
    finalDF <- finalDF[keep,]
    ah_low = 10^input$ah_range[1]
    ah_high = 10^input$ah_range[2]
    dataset = input$dataset
    range = input$range
    position = input$position
    pos <- replace(c(F,F,F,F,F), ifelse(position == "Max", 5, 3), T)
    abbrev <- finalDF$Short_Name
    acmg_ah <- finalDF$Allelic_Heterogeneity %>% strsplit("-") %>% 
      sapply(as.numeric) %>% sapply(mean) %>% pmax(ah_low) %>% pmin(ah_high)
    
    sapply(c("Total",super.levels), function(superpop){
      # Map of disease name to disease tags
      if (dataset == "1000 Genomes") {
        find <- "AF_1000G"
        named.freqs <- freq_1000g.count
      }
      if (dataset == "ExAC") {
        find <- "AF_EXAC"
        named.freqs <- freq_exac.calc
      }
      if (superpop != "Total") 
        find <- paste(find, superpop, sep = "_")
      named.freqs <- named.freqs[,find] %>% unlist %>% setNames(abbrev)
      allelic.het <- c(ah_low, ah_low, mean(c(ah_low, ah_high)) %>% signif(3), ah_high, ah_high) %>% 
        rep(nrow(finalDF)) %>% matrix(nrow = nrow(finalDF), byrow = T)
      allelic.het[,3] <- acmg_ah
      sapply(finalDF$Inverse_Prevalence %>% strsplit("-"), function(prev){
        if (length(prev)>1) {
          temp <- 1/as.numeric(prev)
        } else {
          #take unique prev_1 values and compute a = 2k/(1+r) = lower value of a 5x range
          temp <- c(range,1)*2/(as.numeric(prev)*(1+range))
        }
        return(c(temp[2], temp[2], mean(c(temp[1], temp[2])), temp[1], temp[1]))
      }) %>% t -> prev_final
      # Matrix of penetrance values for allelic het range, capped at 1
      set_to_na <- function(row) { replace(row, is.infinite(row),NA) %>% pmin(1)}
      apply(prev_final*allelic.het / named.freqs, 1, set_to_na) %>% unlist
    }) %>% as.data.frame -> penetrance_init
    
    order(matrix(penetrance_init$Total, nrow = 5) %>% colSums, decreasing = T) -> ord
    penetrance_data <- data.frame(penetrance_init, 
                                  "Disease" = factor(sapply(abbrev, 
                                                            function(x) rep(x,5)) %>% as.vector,
                                                     levels = abbrev[ord]) ) 
    # Star/Radar Plot
    #temp <- penetrance_data[pos,] %>% select(-Disease, -Total)
    #rownames(temp) <- abbrev
    #order(rowSums(temp, na.rm = T), decreasing = T)[1:10] -> wanted
    #col <- 3
    #stars(temp[wanted,], scale = F, full = F, len = 1, nrow = ceiling(10/col), ncol = col, 
    #      flip.labels = F, key.loc = c(2*col,2), 
    #      main = sprintf("Radar Plot: %s Penetrance by Ancestry (%s)", position, dataset),
    #      draw.segments = T, col.segments = c('red','yellow','green','blue','purple'))
    #print("These are the top 10 diseases by summed allele frequencies. NULL values are not plotted.", quote = F)
    #print("Each radius is proportional to the penetrance of the disease in the given population.", quote = F)
    # Barplot
    penetrance_data <- gather(penetrance_data, Subset, Penetrance, -Disease)
    barplot <- ggplot(aes(x=Disease, y=Penetrance), data = penetrance_data) + 
           geom_boxplot(position = 'identity', coef = 0, na.rm = T) + coord_flip() + xlab(NULL) + 
           facet_wrap(~Subset, ncol = 2) + ggtitle(sprintf("Barplot: Penetrance by Ancestry (%s)", dataset)) + 
           theme(axis.text.y=element_text(size=6), axis.text.x = element_text(angle = -20, hjust = 0.4))
    #barplot <- ggplotly(barplot, height = 800)
    # Heatmap
    #heatmap <- ggplot(aes(x=Disease, y = Subset), data = penetrance_data[pos,]) + coord_flip() + 
    #       geom_tile(aes(fill = Penetrance), color = 'white') + xlab("Disease") + ylab("Ancestry") +
    #       scale_fill_gradient(low='white',high = 'darkblue', na.value = "grey50",
    #                           breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.50","0.75","1.00"), limits =c(0,1)) + 
    #       ggtitle(sprintf("Heatmap: %s Penetrance by Ancestry (%s)", position, dataset)) + 
    #       theme_minimal() + theme(axis.ticks = element_blank()) + 
    #       annotate("segment", y=c(0.5,5.5,6.5), yend=c(0.5,5.5,6.5), 
    #                x=0.5, xend = length(abbrev)+0.5) +
    #       annotate("segment", y=0.5, yend=6.5, 
    #                x=c(0.5,length(abbrev)+0.5), 
    #                xend = c(0.5,length(abbrev)+0.5))
    m <- list(l = 150, r = 150, b = -50, t = 100, pad = 5)
    heatmap <- plot_ly(
      x = factor(c(super.levels,"Total"), levels = c("Total",super.levels)),
      y = factor(sapply(abbrev, function(x) rep(x,5)) %>% as.vector, levels = abbrev[ord]),
      z = penetrance_init[pos,][ord,] %>% as.matrix, type = "heatmap", height = 800
    ) %>% layout(autosize = T, margin = m, 
      title = sprintf("%s Penetrance by Ancestry (%s)", position, dataset))
    output$barplot <- renderPlot({ barplot })
    output$heatmap <- renderPlotly({ heatmap })
    #cat("Dark gray boxes are NA: no associated variants discovered in that ancestral population.")
  }) #Closes observeEvent "Make Plot"
  
})

# Run the application 
shinyApp(ui = ui, server = server)

