#title: "ACMG-ClinVar Penetrance Shiny App"
#author: "James Diao, under the supervision of Arjun Manrai"
#date: "January 25, 2017"

# Set working directory to file folder
outdir <- getSrcDirectory(function(dummy) {dummy})
setwd(outdir)

# Output name for saving modified csv
outfilename <- "Modified_Estimates"
# Packages to install
pkg_list <- c("knitr","shiny","shinysky","rhandsontable","plotly",
              "ggplot2","tibble","tidyr","dplyr")
installed <- pkg_list %in% installed.packages()[,"Package"]
if (!all(installed))
  install.packages(pkg_list[!installed])
sapply(pkg_list, require, character.only = T)

# load allele frequencies and other info
# setwd("/Users/jamesdiao/Documents/Kohane_Lab/2016-paper-ACMG-penetrance/Shiny_App")
super.levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
load(file = "disease_level_AF.RData")
#save(freq_1000g.calc.gene, freq_1000g.count.gene, 
#     freq_exac.calc.gene, freq_gnomad.calc.gene, 
#     file = "disease_level_AF.RData")

# Read in the given .csv
DF <- read.csv(file = "../Supplementary_Files/Literature_Prevalence_Estimates.csv", 
               header = TRUE, stringsAsFactors = F, na.strings = "NA") %>%
  select(Evaluate, Abbreviation, Short_Name, Inverse_Prevalence, Case_Allele_Frequency)

### User Interface
ui <- shinyUI(fluidPage(
  titlePanel("Estimation of Penetrance from Population Parameters"),
  h4("James Diao | January 20, 2017"),
  sidebarLayout(
    sidebarPanel(
      h3("Control Panel"),
      wellPanel(
        radioButtons("dataset", "Select Dataset", c("gnomAD","ExAC", "1000G")),
        radioButtons("method", "Select Phenotype ID Method", c("Gene","MIM", "MedGen")),
        radioButtons("position", "Select Heatmap Values", c("Max", "Mean")),
        sliderInput("ah_range", "Case Allele Frequency Range", 
                    min = 0, max = 1, value = c(0.01,1), step = 0.05)
      ),
      wellPanel(
        h4("Generate Plots"),
        actionButton("run", " Make Plots", icon = icon("bar-chart"), styleclass = "success"),
        h4()
      ),
      wellPanel(
        h4("Save Modified Table"), 
        div(class='row', 
            div(class="col-sm-6", 
                actionButton("save", "Save Table", styleclass = "success"))
        )
      ),
      shinyalert("saved", click.hide = TRUE, auto.close.after = 3)
    ),
    
    mainPanel(
      wellPanel(
        helpText(sprintf("Working Directory: %s", outdir)),
        helpText("Double-click on the table to edit values. Uncheck boxes to remove from analysis.")
      ),
      rHandsontableOutput("hot"),
      br(),
      busyIndicator("In Progress: Please Wait", wait = 500),
      h2("Heatmap"),
      plotlyOutput("heatmap", height = "800px"),
      h2("Barplot"),
      plotlyOutput("barplot", height = "800px")
      
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
    #freq_1000g <- freq_1000g.calc.gene[keep,]
    #freq_exac <- freq_exac.calc.gene[keep,]
    #freq_gnomad <- freq_gnomad.calc.gene[keep,]
    finalDF <- finalDF[keep,]
    ah_low <- input$ah_range[1]
    ah_high <- input$ah_range[2]
    dataset <- input$dataset
    method <- input$method
    range <- input$range
    position <- input$position
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
      set_to_na <- function(row) { replace(row, is.infinite(row),NA) %>% pmin(1)}
      apply(allelic.het / finalDF$Inverse_Prevalence / named.freqs, 1, set_to_na) %>% unlist
    }) %>% as.data.frame -> penetrance_init
    
    mat_pen <- matrix(penetrance_init$Total, ncol = 5, byrow = T)
    if (position == "Max") {
      ord <- order(mat_pen[,5], mat_pen[,3], mat_pen[,1], decreasing = T)
    } else {
      ord <- order(mat_pen[,3], mat_pen[,5], mat_pen[,1], decreasing = T)
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
      #%>% layout(plot_bgcolor='rgb(40, 40, 80)') 
      #layout(plot_bgcolor='rgb(190, 190, 190)')
    output$heatmap <- renderPlotly({ heatmap })
    output$barplot <- renderPlotly({ barplot })
    #cat("Dark gray boxes are NA: no associated variants discovered in that ancestral population.")
  }) #Closes observeEvent "Make Plot"
  
})

# Run the application 
shinyApp(ui = ui, server = server)


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

#tags$head(tags$script(HTML(JS.logify))),
#tags$head(tags$script(HTML(JS.onload))),
# javascript function for log-scale ah_range
#JS.logify <- "function logifySlider (sliderId) {
#  $('#'+sliderId).data('ionRangeSlider').update({
# 'prettify': function (num) { return (Math.pow(10, num)); }
#  })
#}"
#JS.onload <- "$(document).ready(function() {
#  setTimeout(function() {
#  logifySlider('ah_range')
#}, 5)}) "