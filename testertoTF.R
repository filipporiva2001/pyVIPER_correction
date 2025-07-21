# Set working directory
setwd("/Users/friva/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Desktop/SBP_CODE")
# Step 8: Compare Top 100 TFs between R and Python
# ------------------------------------------------

# Read both CSV files
r_df <- read.csv("R_top500_diffs.csv", stringsAsFactors = FALSE)
py_df <- read.csv("Python_top500_diffs.csv", stringsAsFactors = FALSE)

# Extract TF names
r_tfs <- r_df$TF
py_tfs <- py_df$TF

# Find overlap
overlap_tfs <- intersect(r_tfs, py_tfs)
unique_r <- setdiff(r_tfs, py_tfs)
unique_py <- setdiff(py_tfs, r_tfs)

# Print summary
cat("=== TF Overlap Summary ===\n")
cat("Number of TFs in R list: ", length(r_tfs), "\n")
cat("Number of TFs in Python list: ", length(py_tfs), "\n")
cat("Number of overlapping TFs: ", length(overlap_tfs), "\n\n")

# Print overlapping TFs
cat("Overlapping TFs:\n")
print(overlap_tfs)

# Optionally: Print TFs unique to each list
cat("\nTFs only in R list:\n")
print(unique_r)

cat("\nTFs only in Python list:\n")
print(unique_py)
