# app.R

# Load the packages from "install.r"
library(shiny)
library(Biostrings)
library(randomForest)
library(Peptides)
library(dplyr)
library(stringr)
library(ggplot2)
library(here)
library(DT)

# -------------------------
# Load trained model
rf_model <- readRDS("../model/rf_model.rds")

# -------------------------
# Feature computation function
compute_features <- function(seqs) {
  seqs <- seqs[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", seqs)] # Extract only sequences with standard amino acids
  feature_df <- data.frame(sequence = seqs, stringsAsFactors = FALSE) # Store as a data frame
  
  # Create feature data frame
  feature_df <- feature_df %>%
    mutate(
      length = lengthpep(sequence),
      pI = pI(sequence),
      hydrophobicity = hydrophobicity(sequence),
      charge = charge(sequence),
      amphipathicity = hmoment(sequence),
      aliphatic_index = aIndex(sequence),
      tiny_prop = str_count(sequence, "[ACGST]") / str_length(sequence),
      small_prop = str_count(sequence, "[ABCDGNPSTV]") / str_length(sequence),
      aromatic_prop = str_count(sequence, "[FHWY]") / str_length(sequence),
      pos_prop = str_count(sequence, "[HKR]") / str_length(sequence),
      neg_prop = str_count(sequence, "[DE]") / str_length(sequence),
      charged_prop = str_count(sequence, "[DEHKR]") / str_length(sequence),
      polar_prop = str_count(sequence, "[DEHKNQRSTZ]") / str_length(sequence)
    )
  
  # Add individual amino-acid proportions to the feature data frame 
  aa <- c("A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y")
  for (a in aa) {
    feature_df[[paste0(a, "_prop")]] <- str_count(feature_df$sequence, fixed(a)) / str_length(feature_df$sequence)
  }
  
  feature_matrix <- feature_df %>% select(-sequence)
  return(feature_matrix)
}

# -------------------------
# R SHINY
# -------------------------

# UI Elements
ui <- fluidPage(
  titlePanel("AMP_Predictor"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_file", "Upload FASTA File", accept = c(".fasta", ".fa")),
      actionButton("predict_btn", "Predict AMP"),
      br(), br(),
      downloadButton("download_results", "Download Results (CSV)")
    ),
    
    mainPanel(
      DT::DTOutput("results"),  # Table of sequence results on top
      br(),
      fluidRow(
        column(6, plotOutput("prediction_plot")),  # Histogram on the left
        column(6, plotOutput("pie_chart"))         # Pie chart on the right
      )
    )
  )
)


# -------------------------
# Server
server <- function(input, output) {
  
  # Read sequences and FASTA IDs
  fasta_data <- reactive({
    req(input$fasta_file)
    fasta <- readAAStringSet(input$fasta_file$datapath)
    seqs <- as.character(fasta)
    ids <- names(fasta)
    
    # Filter invalid sequences
    valid_idx <- !grepl("[^ACDEFGHIKLMNPQRSTVWY]", seqs)
    seqs <- seqs[valid_idx]
    ids <- ids[valid_idx]
    
    if(length(seqs) == 0){
      showNotification("No valid sequences found!", type = "error")
      return(NULL)
    }
    
    list(sequences = seqs, ids = ids)
  })
  
  # Compute features using the "compute_features" function
  features <- reactive({
    req(fasta_data())
    compute_features(fasta_data()$sequences)
  })
  
  # Make predictions
  predictions <- eventReactive(input$predict_btn, {
    req(features())
    prob_matrix <- predict(rf_model, features(), type = "prob")
    prob_amp <- prob_matrix[, "AMP"]
    class_label <- ifelse(prob_amp > 0.70, "AMP", "nonAMP") # If probability greater than 0.7 then classify as an AMP
    
    # Compute sequence lengths
    seq_lengths <- nchar(fasta_data()$sequences)
    
    # Create data frame of results
    data.frame(
      ID = fasta_data()$ids,
      Sequence = fasta_data()$sequences,
      Length = seq_lengths,
      AMP_Probability = round(prob_amp, 3),
      Prediction = class_label
    )
  })
  
  # Display results table on app
  output$results <- DT::renderDT({
    req(predictions())
    DT::datatable(
      predictions(),
      options = list(
        pageLength = 100,       # number of rows per page
        scrollY = "700px",      # vertical scroll
        scrollCollapse = TRUE,  # allow table to shrink if less rows
        autoWidth = TRUE        # adjust column widths
      ),
      rownames = FALSE          # hide row numbers
    )
  })
  
  # Display histogram of AMP probabilities
  output$prediction_plot <- renderPlot({
    req(predictions())
    df <- predictions()
    ggplot(df, aes(x = AMP_Probability)) +
      geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
      geom_vline(xintercept = 0.70, linetype = "dashed", color = "black", linewidth = 1) +
      annotate(
        "text",
        x = 0.70,
        y = Inf,
        label = "Predicted AMP (≥ 0.70)",
        hjust = -0.05, 
        vjust = 1.5,
        color = "black",
        size = 5
      ) +
      theme(axis.text=element_text(size=16),
                    axis.title=element_text(size=16)) +
      labs(title = "Distribution of AMP Prediction Scores", x = "AMP Probability", y = "Number of sequences")
  })
  
  # Display Pie chart of AMP status
  output$pie_chart <- renderPlot({
    req(predictions())
    
    # Count AMP vs nonAMP
    counts <- table(predictions()$Prediction)
    
    # Simple pie chart
    pie(counts,
        labels = paste0(names(counts), " (", counts, ")"),
        col = c("steelblue", "red"),
        main = "Proportion of AMP vs nonAMP")
  })

  
  # Download CSV button
  output$download_results <- downloadHandler(
    filename = function() { "AMP_predictions.csv" },
    content = function(file) {
      write.csv(predictions(), file, row.names = FALSE)
    }
  )
  
}

# -------------------------
# Run the app
shinyApp(ui = ui, server = server)
