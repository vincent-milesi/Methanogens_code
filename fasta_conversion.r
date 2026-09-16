# install.packages("CHNOSZ")
# install.packages("canprot")
# install.packages("xlsx")

library(CHNOSZ)
library(canprot)
library(xlsx)

# Read folder names from Excel file
folders <- xlsx::read.xlsx(r"(Q:\Documents\Projets\Methanogens\Metadata_methano_VM.xlsx)",
                            sheetName = "Feuil1", colIndex = 3)          
# Extract the column as a vector
folders <- folders[[1]]
# Store results for each organism
results <- list()
# Loop through each folder and process the protein.faa files
for (folder in folders[1:length(folders)]) {
    # Print a message indicating the folder being processed
    message("Processing: ", folder)
    # Path to NCBI dataset/data directory
    path <- file.path(r"(C:\Users\vmilesi\Documents\Omics_data\Methanogens)", folder, "ncbi_dataset", "ncbi_dataset", "data")
    # Find protein.faa files in the accession directories
    fasta_files <- list.files(path, pattern = "^protein\\.faa$", recursive = TRUE, full.names = TRUE)
    # Check if any protein.faa files were found
    if (length(fasta_files) == 0) {
        warning("No protein.faa found for: ", folder)
        next
    }

    # Read all FASTA files if there is more than one -> read both GCA and GCF: redundant? 
    # aa_list <- lapply(fasta_files, canprot::read_fasta)

    # Read only the first fasta file
    aa_list <- lapply(fasta_files[1], canprot::read_fasta)
    # Combine all proteins
    aa <- do.call(rbind, aa_list)
    # Average amino-acid composition per protein
    aa_mean <- canprot::sum_aa(aa, average = TRUE)
    # Modify metadata
    aa_mean$protein <- folder
    aa_mean$organism <- folder
    aa_mean$ref <- ''
    aa_mean$abbrv <- ''
    # Relative average amino-acid composition per protein (apply calculation only to numeric columns in the dataframe)
    aa_rel_mean <- aa_mean[sapply(aa_mean, is.numeric)]/CHNOSZ::protein.length(aa_mean)
    # Average protein length
    aa_mean$length <- CHNOSZ::protein.length(aa_mean)
    # Chemical formula of the average protein
    aa_mean$formula <- CHNOSZ::protein.formula(aa_mean)
    # Carbon oxidation state
    aa_mean$ZC <- CHNOSZ::ZC(aa_mean$formula)
    # Store result
    results[[folder]] <- cbind(aa_mean, rel = aa_rel_mean)
}
# Combine all organisms into one data frame
results <- do.call(rbind, results)
# Reset row names
rownames(results) <- NULL
# Save results to Excel file
xlsx::write.xlsx(results, file.path(r"(Q:\Documents\Projets\Methanogens\Methanogens_proteins.xlsx)"), row.names = FALSE)