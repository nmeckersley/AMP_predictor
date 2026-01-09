# app.R
library(shiny)
library(Biostrings)
library(randomForest)
library(Peptides)
library(dplyr)
library(stringr)
library(ggplot2)
library(here)  # robust paths

# -------------------------
# Load trained model
rf_model <- readRDS("model/rf_model.rds")

# -------------------------
# Feature computation function
compute_features <- function(seqs) {
  seqs <- seqs[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", seqs)]
  feature_df <- data.frame(sequence = seqs, stringsAsFactors = FALSE)
  
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
  
  aa <- c("A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y")
  for (a in aa) {
    feature_df[[paste0(a, "_prop")]] <- str_count(feature_df$sequence, fixed(a)) / str_length(feature_df$sequence)
  }
  
  feature_matrix <- feature_df %>% select(-sequence)
  return(feature_matrix)
}

# -------------------------
# UI
ui <- fluidPage(
  titlePanel("Antimicrobial Peptide Predictor"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_file", "Upload FASTA File", accept = c(".fasta", ".fa")),
      actionButton("predict_btn", "Predict AMP"),
      br(), br(),
      downloadButton("download_results", "Download Results (CSV)")
    ),
    
    mainPanel(
      tableOutput("results"),
      plotOutput("prediction_plot")
    )
  )
)

# -------------------------
# Server
server <- function(input, output) {
  
  # Reactive: Read sequences and FASTA IDs
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
  
  # Reactive: Compute features
  features <- reactive({
    req(fasta_data())
    compute_features(fasta_data()$sequences)
  })
  
  # Reactive: Make predictions
  predictions <- eventReactive(input$predict_btn, {
    req(features())
    prob_matrix <- predict(rf_model, features(), type = "prob")
    prob_amp <- prob_matrix[, "AMP"]
    class_label <- ifelse(prob_amp > 0.5, "AMP", "nonAMP")
    
    # Compute sequence lengths
    seq_lengths <- nchar(fasta_data()$sequences)
    
    data.frame(
      ID = fasta_data()$ids,
      Sequence = fasta_data()$sequences,
      Length = seq_lengths,
      AMP_Probability = round(prob_amp, 3),
      Prediction = class_label
    )
  })
  
  # Render results table
  output$results <- renderTable({
    req(predictions())
    predictions()
  })
  
  # Plot histogram of AMP probabilities
  output$prediction_plot <- renderPlot({
    req(predictions())
    df <- predictions()
    ggplot(df, aes(x = AMP_Probability)) +
      geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
      theme_minimal() +
      labs(title = "Distribution of AMP Prediction Scores", x = "Probability of being AMP", y = "Number of sequences")
  })
  
  # Download CSV
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
