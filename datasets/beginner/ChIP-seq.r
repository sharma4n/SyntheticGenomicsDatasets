# Load required libraries
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(Biostrings)

# Set seed for reproducibility
set.seed(123)

# Function to create synthetic ChIP-seq dataset
create_synthetic_chipseq <- function(num_peaks = 1000, 
                                     genome_size = 1e6,
                                     peak_width_mean = 200,
                                     peak_width_sd = 50,
                                     background_coverage = 5,
                                     peak_enrichment = 20) {
  
  # Create chromosome information
  chromosomes <- c("chr1", "chr2", "chr3")
  chr_sizes <- c(400000, 350000, 250000)  # sizes for our synthetic chromosomes
  
  # Initialize lists to store results
  all_peaks <- GRanges()
  coverage_data <- list()
  peak_sequences <- DNAStringSet()
  
  # Create peaks for each chromosome
  peak_counter <- 1
  for (chr_idx in 1:length(chromosomes)) {
    chr <- chromosomes[chr_idx]
    chr_size <- chr_sizes[chr_idx]
    
    # Determine number of peaks for this chromosome (proportional to size)
    num_chr_peaks <- round(num_peaks * (chr_size / sum(chr_sizes)))
    
    if (num_chr_peaks > 0) {
      # Generate random peak positions
      peak_starts <- sample(1:(chr_size - 500), num_chr_peaks, replace = FALSE)
      peak_widths <- round(rnorm(num_chr_peaks, peak_width_mean, peak_width_sd))
      peak_widths <- pmax(peak_widths, 50)  # ensure minimum width
      
      peak_ends <- peak_starts + peak_widths
      
      # Create GRanges object for peaks
      chr_peaks <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = peak_starts, end = peak_ends),
        strand = sample(c("+", "-"), num_chr_peaks, replace = TRUE),
        score = round(rnorm(num_chr_peaks, 100, 30)),  # peak scores
        pvalue = 10^(-runif(num_chr_peaks, 2, 10)),   # p-values
        qvalue = 10^(-runif(num_chr_peaks, 1, 8))     # q-values
      )
      
      all_peaks <- c(all_peaks, chr_peaks)
      
      # Generate coverage data (simulate read counts)
      coverage <- numeric(chr_size)
      
      # Add background coverage
      coverage[] <- rpois(chr_size, background_coverage)
      
      # Add enriched coverage at peaks
      for (i in 1:num_chr_peaks) {
        peak_range <- peak_starts[i]:peak_ends[i]
        # Ensure we don't go beyond chromosome boundaries
        peak_range <- peak_range[peak_range <= chr_size]
        coverage[peak_range] <- rpois(length(peak_range), peak_enrichment)
      }
      
      coverage_data[[chr]] <- coverage
    }
  }
  
  # Generate synthetic DNA sequences for peaks
  generate_peak_sequences <- function(peaks_granges) {
    sequences <- DNAStringSet()
    for (i in 1:length(peaks_granges)) {
      peak_width <- width(peaks_granges[i])
      # Generate random DNA sequence
      seq <- paste(sample(c("A", "C", "G", "T"), peak_width, replace = TRUE), 
                   collapse = "")
      sequences <- c(sequences, DNAString(seq))
    }
    names(sequences) <- paste0("peak_", 1:length(peaks_granges))
    return(sequences)
  }
  
  peak_sequences <- generate_peak_sequences(all_peaks)
  
  return(list(
    peaks = all_peaks,
    coverage = coverage_data,
    sequences = peak_sequences,
    chromosomes = data.frame(
      chr = chromosomes,
      size = chr_sizes
    )
  ))
}

# Function to simulate ChIP-seq experiment with input control
simulate_chipseq_experiment <- function(treatment_name = "H3K27ac", 
                                        num_replicates = 2) {
  
  # Simulate treatment data
  treatment_data <- create_synthetic_chipseq(
    num_peaks = 800,
    peak_enrichment = 25  # Higher enrichment for treatment
  )
  
  # Simulate input control (background)
  input_data <- create_synthetic_chipseq(
    num_peaks = 200,      # Fewer peaks in input
    peak_enrichment = 8   # Lower enrichment
  )
  
  # Create replicate data
  replicates <- list()
  for (rep in 1:num_replicates) {
    rep_data <- create_synthetic_chipseq(
      num_peaks = 700 + sample(-50:50, 1),  # Slight variation between replicates
      peak_enrichment = 20 + rnorm(1, 0, 2)
    )
    replicates[[paste0("rep", rep)]] <- rep_data
  }
  
  return(list(
    treatment = treatment_data,
    input = input_data,
    replicates = replicates,
    experiment_name = treatment_name
  ))
}

# Function to export data in common formats
export_chipseq_data <- function(chipseq_data, output_dir = "chipseq_synthetic_data") {
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # Export peaks as BED file
  peaks_bed <- data.frame(
    chr = as.character(seqnames(chipseq_data$peaks)),
    start = start(chipseq_data$peaks) - 1,  # BED is 0-based
    end = end(chipseq_data$peaks),
    name = paste0("peak_", 1:length(chipseq_data$peaks)),
    score = mcols(chipseq_data$peaks)$score,
    strand = as.character(strand(chipseq_data$peaks))
  )
  
  write.table(peaks_bed, file.path(output_dir, "peaks.bed"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Export peak sequences as FASTA
  writeXStringSet(chipseq_data$sequences, file.path(output_dir, "peak_sequences.fa"))
  
  # Export coverage data
  for (chr in names(chipseq_data$coverage)) {
    cov_df <- data.frame(
      position = 1:length(chipseq_data$coverage[[chr]]),
      coverage = chipseq_data$coverage[[chr]]
    )
    write.table(cov_df, file.path(output_dir, paste0("coverage_", chr, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  # Export metadata
  metadata <- data.frame(
    parameter = c("num_peaks", "genome_size", "experiment_type"),
    value = c(length(chipseq_data$peaks), 
              sum(sapply(chipseq_data$coverage, length)),
              "synthetic")
  )
  write.table(metadata, file.path(output_dir, "metadata.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Data exported to:", output_dir, "\n")
}

# Function to visualize the synthetic data
visualize_chipseq_data <- function(chipseq_data, chr_to_plot = "chr1") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)
  
  # Extract coverage for the specified chromosome
  coverage <- chipseq_data$coverage[[chr_to_plot]]
  positions <- 1:length(coverage)
  
  # Create coverage plot
  coverage_df <- data.frame(position = positions, coverage = coverage)
  
  p <- ggplot(coverage_df, aes(x = position, y = coverage)) +
    geom_line(alpha = 0.7, color = "blue") +
    labs(title = paste("Synthetic ChIP-seq Coverage -", chr_to_plot),
         x = "Genomic Position", y = "Coverage") +
    theme_minimal()
  
  print(p)
  
  # Create histogram of peak scores
  peak_scores <- mcols(chipseq_data$peaks)$score
  score_df <- data.frame(score = peak_scores)
  
  p2 <- ggplot(score_df, aes(x = score)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    labs(title = "Distribution of Peak Scores",
         x = "Peak Score", y = "Frequency") +
    theme_minimal()
  
  print(p2)
}

# Example usage and demonstration
if (TRUE) {  # Set to FALSE if you don't want to run the example
  # Create a synthetic ChIP-seq dataset
  cat("Creating synthetic ChIP-seq dataset...\n")
  chipseq_data <- create_synthetic_chipseq(num_peaks = 1000)
  
  # Display summary information
  cat("Dataset Summary:\n")
  cat("Number of peaks:", length(chipseq_data$peaks), "\n")
  cat("Chromosomes:", unique(as.character(seqnames(chipseq_data$peaks))), "\n")
  cat("Peak width statistics:\n")
  cat("  Mean:", mean(width(chipseq_data$peaks)), "\n")
  cat("  SD:", sd(width(chipseq_data$peaks)), "\n")
  cat("  Range:", range(width(chipseq_data$peaks)), "\n")
  
  # Show first few peaks
  cat("\nFirst 5 peaks:\n")
  print(head(chipseq_data$peaks, 5))
  
  # Simulate a complete experiment
  cat("\nSimulating complete ChIP-seq experiment...\n")
  experiment <- simulate_chipseq_experiment("H3K27ac", num_replicates = 2)
  
  # Export data
  cat("Exporting data...\n")
  export_chipseq_data(chipseq_data)
  
  # Visualize data
  cat("Creating visualizations...\n")
  visualize_chipseq_data(chipseq_data)
}

# Additional utility function to add motif occurrences
add_motifs_to_peaks <- function(chipseq_data, motif = "GATAA", probability = 0.3) {
  # This function adds specific motifs to peak sequences with given probability
  
  modified_sequences <- chipseq_data$sequences
  
  for (i in 1:length(modified_sequences)) {
    if (runif(1) < probability) {
      # Insert motif at random position within the peak
      seq_length <- length(modified_sequences[[i]])
      motif_length <- nchar(motif)
      
      if (seq_length > motif_length) {
        insert_pos <- sample(1:(seq_length - motif_length + 1), 1)
        
        # Create new sequence with inserted motif
        before_motif <- subseq(modified_sequences[[i]], 1, insert_pos - 1)
        after_motif <- subseq(modified_sequences[[i]], insert_pos + motif_length)
        new_seq <- c(before_motif, DNAString(motif), after_motif)
        
        modified_sequences[[i]] <- new_seq
      }
    }
  }
  
  chipseq_data$sequences <- modified_sequences
  return(chipseq_data)
}

# Example with motifs
cat("Adding transcription factor motifs...\n")
chipseq_with_motifs <- add_motifs_to_peaks(chipseq_data, motif = "GATAA", probability = 0.3)

cat("Synthetic ChIP-seq dataset creation complete!\n")
