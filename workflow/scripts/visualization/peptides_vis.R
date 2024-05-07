library(Biostrings)
library(ggplot2)
library(cowplot)
library(tools)
library(patchwork)

#@ write outputs to file as log! 
log_file <- snakemake@log[['log_file']]
log <- file(log_file, open = "wt")
sink(log)
sink(log, type="message")

peptide_vis <- function(fasta_file, cut_off, sample_name){
  
  file_name <- file_path_sans_ext(basename(fasta_file))
  output_file <- paste0(file_name, "_short.fasta")
  writer <- file(output_file, "w")
  
  peptide_sequences <- readAAStringSet(fasta_file)
  
  under_cutoff <- 0
  upper_cutoff <- 0
  
  for (i in 1:length(peptide_sequences)) {
    sequence <- as.character(peptide_sequences[[i]])
    sequence_length <- nchar(sequence)
    sequence_name <- names(peptide_sequences)[i]
    complete_header <- paste0(">", sequence_name)
    
    if (sequence_length < cut_off) {

      under_cutoff <- under_cutoff + 1
      writeLines(complete_header, writer)
      writeLines(sequence, writer)
    } else {
      upper_cutoff <- upper_cutoff + 1
    }
  }

  close(writer)
  
  # Print the result
  cat("Total number of sequences       :", length(peptide_sequences), "\n")
  cat("Number of sequences upper cutoff:", upper_cutoff,"/ ",length(peptide_sequences),   "\n")
  cat("Number of sequences under cutoff:", under_cutoff,"/ ",length(peptide_sequences),  "\n")
  cat("Ratio: ", (under_cutoff / length(peptide_sequences)) * 100, "%\n")
  
  df <- data.frame(
    Category = c(sample_name),
    Value = c(under_cutoff, upper_cutoff)
  )
  
  p = ggplot(df, aes(x = Value, y =  Category, fill = ifelse(Value > under_cutoff, "Long", "Short"))) +
    geom_bar(stat = "identity", position = "stack", width = 0.3) +
    geom_text(aes(label = Value), position = position_stack(vjust = 0.6), color = "white") +
    scale_fill_manual(values = c("#89a8a2", "#EF6F61"), guide = "none") +
    labs(title = NULL,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank())
  
  return(p)
}

input_file <- snakemake@input[['input_peptides']]
k_mer <- snakemake@params[['k_mer']]
sample_name <- snakemake@params[['sample_name']]
output_file <- snakemake@output[['output_vis']]

p1 = peptide_vis(input_file, k_mer, paste0("Peptide Visualization of ", sample_name))
ggsave(output_file, plot = p1, width = 12, height = 7, units = "in")