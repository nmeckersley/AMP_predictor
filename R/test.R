# Input AMP Data from a fasta file
upload_fasta <- read.fasta(file = "./Data/amps.fasta", seqtype = "AA")
seq_chars <- sapply(upload_fasta, function(x) paste(x, collapse = ""))

aa_df<- data.frame(
  sequence = seq_chars,
  stringsAsFactors = FALSE
)

# Remove non-standard amino acids
aa_df <- aa_df[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", aa_df$sequence), ]

feature_df <- aa_df %>%
  mutate(
    length = lengthpep(sequence),
    pI = pI(sequence),
    hydrophobicity = hydrophobicity(sequence),
    charge = charge(sequence),
    amphipathicity = hmoment(sequence),
    aliphatic_index = aIndex(sequence)

  )

# Creature a column for IDs
feature_df <- feature_df %>%
  mutate(id = sprintf("%05d", row_number())) %>%
  relocate(id, .before = 1)

# Additional structural features

# Tiny (A, C, G, S, T)
feature_df$tiny_prop <- str_count(feature_df$sequence, "[ACGST]") / str_length(feature_df$sequence)
# Small (A, B, C, D, G, N, P, S, T, V)
feature_df$small_prop <- str_count(feature_df$sequence, "[ABCDGNPSTV]")/ str_length(feature_df$sequence)
# Aromatic (F, H, W, Y)
feature_df$aromatic_prop <- str_count(feature_df$sequence, "[FHWY]")/ str_length(feature_df$sequence)
# Positive (H, K, R)
feature_df$pos_prop <- str_count(feature_df$sequence, "[HKR]")/ str_length(feature_df$sequence)
# Negative (D, E)
feature_df$neg_prop <- str_count(feature_df$sequence, "[DE]")/ str_length(feature_df$sequence)
# Charged (D, E, H, K, R)
feature_df$charged_prop <- str_count(feature_df$sequence, "[DEHKR]")/ str_length(feature_df$sequence)
# Polar (D, E, H, K, N, Q, R, S, T, Z)
feature_df$polar_prop <- str_count(feature_df$sequence, "[DEHKNQRSTZ]")/ str_length(feature_df$sequence)


# Add aa proportions
aa <- c("A","C","D","E","F","G","H","I","K","L",
        "M","N","P","Q","R","S","T","V","W","Y")

for (a in aa) {
  feature_df[[paste0(a, "_prop")]] <-
    str_count(feature_df$sequence, fixed(a)) /
    str_length(feature_df$sequence)
}

feature_df

