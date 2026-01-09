#Load the required packages
library(shiny)
library(Biostrings) # For reading FASTA sequences
library(randomForest) # Using the model
library(ggplot2) # For plotting
library(dplyr)
library(stringr)
library(Peptides)  # For lengthpep(), pI(), hydrophobicity(), charge(), hmoment(), aIndex()



# Load the trained model
model_path <- readRDS("model/rf_model.rds")


###############################################################################

# Feature computation function
compute_features <- function(seqs) {
  
  # Remove sequences with non-standard amino acids
  seqs <- seqs[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", seqs)]
  
  # Create data frame
  feature_df <- data.frame(
    sequence = seqs,
    stringsAsFactors = FALSE
  )
  
  # Core peptide features
  feature_df <- feature_df %>%
    mutate(
      length = lengthpep(sequence),
      pI = pI(sequence),
      hydrophobicity = hydrophobicity(sequence),
      charge = charge(sequence),
      amphipathicity = hmoment(sequence),
      aliphatic_index = aIndex(sequence)
    )
  
  # Add a unique ID
  feature_df <- feature_df %>%
    mutate(id = sprintf("%05d", row_number())) %>%
    relocate(id, .before = 1)
  
  # Structural features
  feature_df <- feature_df %>%
    mutate(
      tiny_prop = str_count(sequence, "[ACGST]") / str_length(sequence),
      small_prop = str_count(sequence, "[ABCDGNPSTV]") / str_length(sequence),
      aromatic_prop = str_count(sequence, "[FHWY]") / str_length(sequence),
      pos_prop = str_count(sequence, "[HKR]") / str_length(sequence),
      neg_prop = str_count(sequence, "[DE]") / str_length(sequence),
      charged_prop = str_count(sequence, "[DEHKR]") / str_length(sequence),
      polar_prop = str_count(sequence, "[DEHKNQRSTZ]") / str_length(sequence)
    )
  
  # Amino acid proportions
  aa <- c("A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y")
  
  for (a in aa) {
    feature_df[[paste0(a, "_prop")]] <- str_count(feature_df$sequence, fixed(a)) / str_length(feature_df$sequence)
  }
  
  # Keep only numeric features for prediction
  feature_matrix <- feature_df %>% select(-sequence)
  
  return(feature_matrix)
}

#####################################################################################

# Define UI
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

# Define server
server <- function(input, output) {
  
  # Reactive: Read sequences from uploaded file
  sequences <- reactive({
    req(input$fasta_file)
    fasta <- readAAStringSet(input$fasta_file$datapath)
    as.character(fasta)
  })
  
  # Reactive: Compute features
  features <- reactive({
    req(sequences())
    compute_features(sequences())
  })
  
  # Reactive: Make predictions when button clicked
  predictions <- eventReactive(input$predict_btn, {
    req(features())
    predict(rf_model, features())
  })
  
  # Display results table
  output$results <- renderTable({
    req(predictions())
    data.frame(Sequence = sequences(), Prediction = predictions())
  })
  
  # Plot predictions
  output$prediction_plot <- renderPlot({
    req(predictions())
    df <- data.frame(Prediction = predictions())
    ggplot(df, aes(x = Prediction)) +
      geom_bar(fill = "steelblue") +
      theme_minimal() +
      labs(title = "AMP Predictions", y = "Count")
  })
  
  # Download results
  output$download_results <- downloadHandler(
    filename = function() { "AMP_predictions.csv" },
    content = function(file) {
      write.csv(data.frame(Sequence = sequences(), Prediction = predictions()), file, row.names = FALSE)
    }
  )
  
}

# Run the app
shinyApp(ui = ui, server = server)
