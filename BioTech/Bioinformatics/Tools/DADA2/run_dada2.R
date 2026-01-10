#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- paste0("dada2_run_log_", timestamp, ".txt")

log_message <- function(msg) {
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# Clear previous log (if any)
write("", file = log_file)

if(length(args) < 1) {
  stop("Usage: Rscript run_dada2_script.R <base_path>")
}

base_path <- args[1]

taxonomy_ref <- "C:/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Tools/UNITE/QIIME/sh_qiime_release_19.02.2025/merged_sh_qiime_19.02.2025.fasta"
taxonomy_map <- "C:/Users/ketgl/OneDrive/Documents/GitHub/Python-Programs/BioTech/Bioinformatics/Tools/UNITE/QIIME/sh_qiime_release_19.02.2025/sh_taxonomy_qiime_ver10_dynamic_19.02.2025.txt"

log_message(paste0("Starting DADA2 run at ", Sys.time()))
log_message(paste0("Base path: ", base_path))
log_message(paste0("Taxonomy reference fasta: ", taxonomy_ref))
log_message(paste0("Taxonomy mapping file: ", taxonomy_map))

if(!dir.exists(base_path)) {
  stop(paste0("Base path does not exist: ", base_path))
}

if(!file.exists(taxonomy_ref)) {
  stop(paste0("Taxonomy reference fasta file does not exist: ", taxonomy_ref))
}

if(!file.exists(taxonomy_map)) {
  stop(paste0("Taxonomy mapping file does not exist: ", taxonomy_map))
}

if (!requireNamespace("dada2", quietly = TRUE)) {
  stop("dada2 package not installed. Please install it before running.")
}
library(dada2)

log_message("Loading taxonomy reference fasta once before barcode processing...")
tax_ref_loaded <- taxonomy_ref

# Get barcode directories
barcode_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
barcode_dirs <- barcode_dirs[grepl("barcode", basename(barcode_dirs), ignore.case = TRUE)]

# Filter only barcodes with numeric value > 19
barcode_dirs <- barcode_dirs[
  grepl("barcode\\d+", basename(barcode_dirs), ignore.case = TRUE) &
  as.integer(gsub("\\D+", "", basename(barcode_dirs))) > 19
]

if(length(barcode_dirs) == 0) {
  stop("No matching barcode directories found after filtering.")
}

log_message(paste("Filtered barcode directories:", paste(basename(barcode_dirs), collapse = ", ")))

for(barcode_dir in barcode_dirs) {
  log_message(paste0("Processing ", basename(barcode_dir)))

  # --------------------------
  # Recursively find _trimmed_cutadapt_PITS.fastq files
  # --------------------------
  fastq_files <- list.files(barcode_dir, pattern = "_trimmed_cutadapt_PITS\\.fastq$", 
                            full.names = TRUE, recursive = TRUE)
  
  if(length(fastq_files) == 0) {
    log_message(paste0("No _trimmed_cutadapt_PITS.fastq files found in ", basename(barcode_dir), ", skipping."))
    next
  }

  log_message(paste0("Found ", length(fastq_files), " FASTQ files for processing. First file: ", basename(fastq_files[1])))

  # --------------------------
  # Learning errors from FASTQ files (combined)
  # --------------------------
  log_message("Learning errors from FASTQ files...")

  err <- tryCatch({
    learnErrors(fastq_files, multithread = TRUE)
  }, error = function(e) {
    log_message(paste0("learnErrors error: ", e$message))
    NULL
  })

  if(is.null(err)) {
    next
  }

  sample_name <- tools::file_path_sans_ext(basename(fastq_files[1]))

  log_message("Running DADA2 core algorithm...")

  dada_out <- tryCatch({
    setDadaOpt(DETECT_SINGLETONS = TRUE)
    dada(fastq_files, err = err, multithread = TRUE)
  }, error = function(e) {
    log_message(paste0("dada error: ", e$message))
    NULL
  })

  if(is.null(dada_out)) {
    next
  }

  names(dada_out) <- sample_name

  seqtab <- makeSequenceTable(dada_out)
  log_message(paste0("Sequence table created with dimensions: ", paste(dim(seqtab), collapse = " x ")))

  log_message(paste0(Sys.time(), " - Assigning taxonomy..."))

  taxa <- tryCatch({
    assignTaxonomy(seqtab, tax_ref_loaded, multithread = TRUE, verbose = TRUE)
  }, error = function(e) {
    log_message(paste0("assignTaxonomy error: ", e$message))
    NULL
  })

  if(!is.null(taxa)) {
    log_message(paste0(Sys.time(), " - Taxonomy assignment completed successfully."))
  } else {
    next
  }

  taxa_df <- as.data.frame(taxa)
  output_csv <- file.path(barcode_dir, paste0(basename(barcode_dir), "_taxa.csv"))
  write.csv(taxa_df, output_csv, row.names = TRUE)
  log_message(paste0("Taxonomy results written to: ", output_csv))
}

log_message(paste0("DADA2 run finished at ", Sys.time()))
