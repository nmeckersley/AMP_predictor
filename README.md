# AMP_Predictor

This Shiny app predicts **Antimicrobial Peptides (AMPs)** from protein sequence data (.fasta), using a trained Random Forest model. AMPs are short (10-50aa), naturally occurring proteins that play a key role in innate immunity by killing or inhibiting the growth of bacteria and other pathogens, often through disruption of microbial cell membranes. Identifying AMPs is important because they are promising candidates for new antimicrobial therapies, especially in the face of increasing antibiotic resistance.

## Note to the Marker

-   The process of training the model can be found in "notebooks/Coursework2.Rmd"

-   The trained model is saved in "model/rf_model.rds"

-   R shiny app is found in "R/app.r"

------------------------------------------------------------------------

## Features

-   Upload a **FASTA file** of protein sequences.
-   Computes sequence features automatically.
-   Predicts **AMP probability** (0–1) for each sequence.
-   Shows a **table** with:
    -   FASTA ID
    -   Sequence
    -   Sequence length
    -   AMP probability
    -   Prediction (AMP / nonAMP)\
-   Plots:
    -   **Histogram** of AMP probabilities
    -   **Pie chart** of AMP vs nonAMP counts\
-   Download results as **CSV**.

------------------------------------------------------------------------

## Installation

1.  Clone or download this repository.
2.  Open R or RStudio in the project folder.
3.  Run the install script to install dependencies:

source("install.R")

------------------------------------------------------------------------

## Usage

Open app.R in RStudio.

Click Run App.

If all packages were installed correctly, this will launch the Shiny app.

Click Browse to upload a FASTA file. There are three example files located in "data/":

1.  mock.fasta: 100 random short protein sequences from UniProt
2.  amps.fasta: 3,306 AMP sequences used for training
3.  non_amps.fasta: 11,000 sequences from UniProt filtered for non-AMP characteristics

Click "**Predict AMP**" to generate results.

The table of predictions can be downloaded as a .csv file for further analysis.

------------------------------------------------------------------------

## Notes

Sequences with non-standard amino acids will be ignored. FASTA IDs are preserved in the output table. Recommended for protein sequences containing only standard amino acid (ACDEFGHIKLMNPQRSTVWY).

------------------------------------------------------------------------

## Further Work

-   Seems to over predict sequences as AMP, will therefore require additional features in the training of the model.
