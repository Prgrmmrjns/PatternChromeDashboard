###########################
# Author: Jonas Wolber
# Description: This is the PatternChrome Dashboard. Type in a gene index number to receive the prediction
# and the prediction breakdown for the XGBoost model. The model is trained on cell line E003 using the training set.
###########################

#### Load libraries ####
library(xgboost)
library(data.table)
library(ggplot2)
library(shiny)
library(shinythemes)

#### Load data ####
#setwd("~/PatternChromeDashboard")
xgb_model <- xgb.load("xgb_model")
load("E003_gene_df.RData")
load("explainer.RData")
load("gene_symbols.RData")

#### Parameters ####
threshold = 0.05

#### XGB Explainer functions ####
showWaterfall = function(xgb.model, explainer, data.matrix, idx){
  valid = TRUE
  if (is.na(idx) | idx > nrow(data.matrix) | idx < 1){
    valid = FALSE
    idx <- 1
  }
  data <- slice(xgb.DMatrix(data = as.matrix(data.matrix[,-1])),as.integer(idx))
  nodes = predict(xgb.model,data,predleaf =TRUE)
  colnames = names(explainer)[1:(ncol(explainer)-2)]
  breakdown = data.table(matrix(0,nrow = nrow(nodes), ncol = length(colnames)))
  setnames(breakdown, colnames)
  num_trees = ncol(nodes)
  for (x in 1:num_trees){
    nodes_for_tree = nodes[,x]
    tree_breakdown = explainer[tree==x-1]
    
    breakdown_for_tree = tree_breakdown[match(nodes_for_tree, tree_breakdown$leaf),]
    breakdown = breakdown + breakdown_for_tree[,colnames,with=FALSE]
    
  }
  weight = rowSums(breakdown)
  pred = 1/(1+exp(-weight))
  breakdown_summary = as.matrix(breakdown)[1,]
  data_for_label = data.matrix[idx,]
  i = order(abs(breakdown_summary),decreasing=TRUE)
  breakdown_summary = breakdown_summary[i]
  data_for_label = data_for_label[i]
  intercept = breakdown_summary[names(breakdown_summary)=='intercept']
  data_for_label = data_for_label[names(breakdown_summary)!='intercept']
  breakdown_summary = breakdown_summary[names(breakdown_summary)!='intercept']
  
  i_other =which(abs(breakdown_summary)<threshold)
  other_impact = 0
  if (length(i_other > 0)){
    other_impact = sum(breakdown_summary[i_other])
    names(other_impact) = 'other'
    breakdown_summary = breakdown_summary[-i_other]
    data_for_label = data_for_label[-i_other]
  }
  
  if (abs(other_impact) > 0){
    breakdown_summary = c(intercept, breakdown_summary, other_impact)
    data_for_label = c("", data_for_label,"")
    labels = paste0(names(breakdown_summary)," = ", data_for_label)
    labels[1] = 'intercept'
    labels[length(labels)] = 'other'
  }else{
    breakdown_summary = c(intercept, breakdown_summary)
    data_for_label = c("", data_for_label)
    labels = paste0(names(breakdown_summary)," = ", data_for_label)
    labels[1] = 'intercept'
  }
  gene_expression_level <- data.matrix[idx, 1]
  annotation_height <- 0.05
  correct <- 1
  if (weight < 0){
    correct <- 0
  }
  if (weight > annotation_height){annotation_height <- weight + 0.2}
  else if (max(breakdown_summary) > annotation_height){annotation_height <- max(breakdown_summary) + 0.2}
  options(repr.plot.width = 5, repr.plot.height =3)
  if(valid){
    title <- paste("Gene ", rownames(data.matrix)[idx], " / ", gene_symbols[idx])
    subtitle <- paste("Actual gene expression level: ", gene_expression_level, "| PatternChrome prediction: ", correct)
    color <- "black"
  } else{
    title <- "Incorrect gene input"
    subtitle <- ""
    color <- "red"
  }
  waterfalls::waterfall(values = breakdown_summary,
                        rect_text_labels = round(breakdown_summary, 2),
                        labels = labels,
                        total_rect_text = round(weight, 2),
                        calc_total = TRUE,
                        total_axis_text = "Prediction")  +
    labs(title= title,
         subtitle = subtitle) +
    theme(plot.title = element_text(size=24, face="bold", colour = color),
          plot.subtitle = element_text(size=20, face="italic"))
}

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("superhero"),

    # Application title
    titlePanel("PatternChrome Prediction Breakdown"),
    tags$div("Find more information about PatternChrome on ", tags$a(href="https://github.com/Prgrmmrjns/PatternChromeDashboard", "Github") ), 
    fluidRow(
      # Sidebar with input
      column(3,
             wellPanel(
               selectInput("cell_line", "Cell line",
                           c("E003" = "e003")),
               textInput("gene", "Gene", value="SGIP1"),
               selectInput("method", "Input method", 
                   selected="gene_symbol", c("Gene symbol" = "gene_symbol", "GENCODE ID" = "gene_id", "Gene index in data matrix" = "index")),
               actionButton("submit", "Submit", class = "btn btn-primary")
             )       
      ),
      # Mainpanel with plot
      column(9,
             plotOutput("waterfall")
      )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  observeEvent(input$submit, {
    if(input$method == "index"){
      idx <- as.numeric(input$gene)
    }
    else{
      idx <- NA
      if(input$method == "gene_id"){
        if (input$gene %in% rownames(gene_df)){
          idx <- which(rownames(gene_df) == input$gene)
        }
      }
      else{
        if(input$gene %in% gene_symbols){
          idx <- which(gene_symbols == input$gene)
        }
      }
    }
    output$waterfall <- renderPlot({
      showWaterfall(xgb_model, explainer, gene_df,  idx)
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
