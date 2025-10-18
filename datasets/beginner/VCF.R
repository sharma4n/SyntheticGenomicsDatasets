# Load required packages 
required_packages <- c("Biostrings", "dplyr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "Biostrings") {
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("Biostrings")
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Main function: Generate complete dataset 
generate_synthetic_genomics_dataset <- function(
    chromosome_length = 10000,
    num_variants = 200,
    output_prefix = "synthetic_dataset",
    seed = 1234) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  cat("=== Synthetic Genomics Dataset Generator ===\n")
  cat("Starting generation of", chromosome_length, "bp genome with", num_variants, "variants\n\n")
  
  # Step 1: Generate Reference Genome 
  cat("1. Generating reference genome...\n")
  bases <- c("A", "T", "C", "G")
  random_seq <- sample(bases, size = chromosome_length, replace = TRUE)
  synthetic_chr <- DNAString(paste(random_seq, collapse = ""))
  names(synthetic_chr) <- paste0("synthetic_chr_", chromosome_length, "bp")
  
  ref_file <- paste0(output_prefix, "_reference.fa")
  writeXStringSet(synthetic_chr, filepath = ref_file, format = "fasta")
  cat("   ✓ Reference genome saved:", ref_file, "\n")
  
  # Step 2: Generate Variants (VCF) 
  cat("2. Generating variants...\n")
  variant_positions <- sort(sample(2:(chromosome_length-1), num_variants))
  variant_type <- sample(c("SNP", "INDEL"), size = num_variants, replace = TRUE, prob = c(0.85, 0.15))
  
  vcf_df <- data.frame(
    CHROM = names(synthetic_chr),
    POS = variant_positions,
    ID = ".",
    REF = character(num_variants),
    ALT = character(num_variants),
    QUAL = 100,
    FILTER = "PASS",
    INFO = ".",
    stringsAsFactors = FALSE
  )
  
  # Create variants
  for (i in 1:num_variants) {
    pos <- variant_positions[i]
    ref_base <- as.character(subseq(synthetic_chr, pos, pos))
    
    if (variant_type[i] == "SNP") {
      possible_alt <- bases[bases != ref_base]
      vcf_df$REF[i] <- ref_base
      vcf_df$ALT[i] <- sample(possible_alt, 1)
    } else {
      if (sample(c(TRUE, FALSE), 1)) {
        # Insertion
        ins_length <- sample(1:3, 1)
        random_ins <- paste(sample(bases, ins_length, replace = TRUE), collapse = "")
        vcf_df$REF[i] <- ref_base
        vcf_df$ALT[i] <- paste0(ref_base, random_ins)
      } else {
        # Deletion
        del_length <- sample(1:2, 1)
        end_pos <- min(pos + del_length, chromosome_length)
        deleted_bases <- as.character(subseq(synthetic_chr, pos, end_pos))
        vcf_df$REF[i] <- deleted_bases
        vcf_df$ALT[i] <- ref_base
      }
    }
  }
  
  # Write VCF file
  vcf_header <- c(
    '##fileformat=VCFv4.2',
    '##fileDate=20240914',
    '##source=SyntheticGenomicsSuite_v1.0',
    paste0('##contig=<ID=', names(synthetic_chr), ',length=', chromosome_length, '>'),
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
  )
  
  vcf_file <- paste0(output_prefix, "_variants.vcf")
  writeLines(vcf_header, vcf_file)
  write.table(vcf_df, vcf_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE, append = TRUE)
  cat("   ✓ Variants saved:", vcf_file, "\n")
  
  # Step 3: Generate Summary Statistics 
  cat("3. Generating summary...\n")
  num_snps <- sum(variant_type == "SNP")
  num_indels <- sum(variant_type == "INDEL")
  
  cat("\n=== GENERATION COMPLETE ===\n")
  cat("Files created:\n")
  cat("  -", ref_file, "\n")
  cat("  -", vcf_file, "\n")
  cat("\nDataset Summary:\n")
  cat("  - Genome size:", chromosome_length, "bp\n")
  cat("  - Total variants:", num_variants, "\n")
  cat("  - SNPs:", num_snps, "\n")
  cat("  - INDELs:", num_indels, "\n")
  cat("  - Variant rate:", round(num_variants/chromosome_length * 1000, 2), "variants/kb\n")
  
  return(list(
    reference = synthetic_chr,
    variants = vcf_df,
    reference_file = ref_file,
    vcf_file = vcf_file
  ))
}

# Example usage function 
demo <- function() {
  cat("Running demo...\n")
  cat("This will create a 10,000 bp genome with 200 variants.\n")
  
  result <- generate_synthetic_genomics_dataset(
    chromosome_length = 10000,
    num_variants = 200,
    output_prefix = "demo_dataset",
    seed = 1234
  )
  
  cat("\nDemo complete! Check your files:\n")
  cat("- demo_dataset_reference.fa\n")
  cat("- demo_dataset_variants.vcf\n")
}

# FASTQ simulation function (for future use) 
simulate_fastq_reads <- function(dna_string, read_length = 150, coverage = 5, output_prefix = "sample") {
  cat("Simulating FASTQ reads...\n")
  num_reads <- ceiling((length(dna_string) * coverage) / read_length)
  start_pos <- sample(1:(length(dna_string) - read_length + 1), num_reads, replace = TRUE)
  
  reads <- subseq(dna_string, start = start_pos, width = read_length)
  
  # Generate quality scores
  generate_qual_string <- function(len) {
    qual_scores <- sample(30:40, len, replace = TRUE)
    qual_string <- rawToChar(as.raw(qual_scores + 33))
    return(qual_string)
  }
  
  qual_strings <- sapply(rep(read_length, num_reads), generate_qual_string)
  
  # Write FASTQ file
  fastq_file <- paste0(output_prefix, ".fastq")
  con <- file(fastq_file, "w")
  
  for (i in 1:num_reads) {
    writeLines(paste0("@read_", i), con)
    writeLines(as.character(reads[i]), con)
    writeLines("+", con)
    writeLines(qual_strings[i], con)
  }
  close(con)
  cat("✓ FASTQ file written:", fastq_file, "\n")
  return(fastq_file)
}

# Main execution block 
cat("Synthetic Genomics Dataset Generator loaded!\n")
cat("Available functions:\n")
cat("1. generate_synthetic_genomics_dataset() - Main function\n")
cat("2. demo() - Run a quick demo\n") 
cat("3. simulate_fastq_reads() - Generate FASTQ files\n")
cat("\nTo get started, run: demo()\n")
cat("Or create custom dataset: generate_synthetic_genomics_dataset(20000, 500, 'my_data')\n")
