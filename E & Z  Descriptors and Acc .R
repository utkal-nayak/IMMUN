###Load Required Package
library(seqinr)
###Load FASTA File
fasta_path <- file.choose()
fasta <- read.fasta(fasta_path, seqtype = "AA", as.string = FALSE, set.attributes = TRUE)

###Define E-Descriptor Table
e_table <- data.frame(
  aa = c("A","C","D","E","F","G","H","I","K","L",
         "M","N","P","Q","R","S","T","V","W","Y"),
  E1 = c(0.01, -0.13, 0.30, 0.22, -0.33, 0.22, 0.02, -0.35, 0.24, -0.27,
         -0.24, 0.26, 0.17, 0.15, 0.17, 0.20, 0.07, -0.27, -0.30, -0.14),
  E2 = c(0.13, 0.17, -0.06, -0.28, -0.02, 0.56, -0.18, 0.07, -0.34, 0.02,
         -0.14, 0.04, 0.29, -0.18, -0.36, 0.24, 0.15, 0.14, -0.19, -0.06),
  E3 = c(-0.48, 0.07, -0.01, -0.32, 0.07, -0.02, 0.04, -0.09, -0.04, -0.27,
         -0.16, 0.12, 0.41, -0.03, 0.11, -0.02, -0.02, -0.19, 0.39, 0.43),
  E4 = c(-0.04, 0.57, 0.23, 0.16, 0.00, 0.02, 0.28, -0.20, -0.33, -0.27,
         0.32, 0.12, -0.22, 0.04, -0.26, -0.07, -0.13, -0.20, 0.08, -0.10),
  E5 = c(0.18, -0.37, 0.16, 0.30, 0.21, 0.11, -0.02, -0.11, -0.03, 0.21,
         0.08, -0.06, 0.38, -0.11, -0.36, -0.20, -0.27, -0.30, 0.30, -0.09),
  stringsAsFactors = FALSE
)
rownames(e_table) <- e_table$aa
e_table$aa <- NULL

###Define Z-Descriptor Table
z_table <- data.frame(
  aa = c("A","C","D","E","F","G","H","I","K","L",
         "M","N","P","Q","R","S","T","V","W","Y"),
  Z1 = c(0.07, 0.71, 3.64, 3.08, -2.25, 0.00, 2.33, -1.56, 2.29, -1.81,
         -1.44, 2.05, 0.67, 1.65, 3.01, 0.64, 0.45, -1.52, -3.59, -2.96),
  Z2 = c(1.00, 0.25, 0.36, 0.39, -1.51, 1.00, 0.22, -0.97, 0.07, -1.15,
         -0.97, 0.62, 0.23, 0.98, 0.42, 0.04, 0.00, -0.63, -1.02, -0.77),
  Z3 = c(0.08, -1.18, -1.02, -1.50, 1.84, -2.33, -0.98, 1.52, 1.63, 1.14,
         0.95, 0.60, -1.52, 0.36, 0.00, -1.18, -0.52, 1.34, 1.67, 1.42),
  stringsAsFactors = FALSE
)
rownames(z_table) <- z_table$aa
z_table$aa <- NULL


###ACC Feature Extraction Function
acc_descriptor <- function(sequence, lag = 8, descriptor_table, prefix = "D") {
  aa_seq <- unlist(strsplit(sequence, split = ""))
  valid_aa <- aa_seq %in% rownames(descriptor_table)
  aa_seq <- aa_seq[valid_aa]
  
  L <- length(aa_seq)
  if (L < (lag + 1)) stop("Sequence too short for specified lag.")
  
  desc_mat <- t(sapply(aa_seq, function(aa) as.numeric(descriptor_table[aa, ])))
  means <- colMeans(desc_mat, na.rm = TRUE)
  
  d <- ncol(desc_mat)
  acc_feats <- numeric(d * d * lag)
  idx <- 1
  
  for (t in 1:lag) {
    for (i in 1:d) {
      for (j in 1:d) {
        acc_val <- sum((desc_mat[1:(L - t), i] - means[i]) * 
                         (desc_mat[(1 + t):L, j] - means[j])) / (L - t)
        acc_feats[idx] <- acc_val
        idx <- idx + 1
      }
    }
  }
  
  names(acc_feats) <- paste0(prefix, rep(1:d, each = d * lag), "_", 
                             rep(1:d, times = d * lag), "_lag", 
                             rep(1:lag, each = d * d))
  return(acc_feats)
}



####Extract features with accession numbers
all_features <- lapply(seq_along(fasta), function(i) {
  seq <- toupper(paste(fasta[[i]], collapse = ""))
  accession <- attr(fasta[[i]], "name") 
  tryCatch({
    e_feats <- acc_descriptor(seq, lag = 8, descriptor_table = e_table, prefix = "E")
    z_feats <- acc_descriptor(seq, lag = 8, descriptor_table = z_table, prefix = "Z")
    c(Accession = accession, e_feats, z_feats)
  }, error = function(e) {
    message(paste("Skipping sequence", i, "due to error:", e$message))
    return(NULL)
  })
})


###Combine into data frame
feature_df <- do.call(rbind, all_features[!sapply(all_features, is.null)])
feature_df <- as.data.frame(feature_df, stringsAsFactors = FALSE)
rownames(feature_df) <- NULL

###Move accession to first column and convert numeric columns
feature_df$Accession <- as.character(feature_df$Accession)
numeric_cols <- setdiff(names(feature_df), "Accession")
feature_df[numeric_cols] <- lapply(feature_df[numeric_cols], as.numeric)


###Save Data
write.csv(feature_df, "123&87 EZ_descriptor_ACC_features.csv", row.names = FALSE)
