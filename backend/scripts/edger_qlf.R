# edgeR differential abundance test between two conditions; called by routers/differential_analysis.py
#
# Usage: Rscript edger_qlf.R <counts.bin> <n_rows> <n_cond1> <n_cond2> <lib_sizes> <out.bin>
#   counts.bin: int32 count matrix (n_rows x samples, column-major), condition 1 samples first
#   lib_sizes:  comma-separated total reads of each sample (whole file, before any filtering)
#   out.bin:    float64 matrix (n_rows x 4, column-major): kept (1/0), logFC, logCPM, FDR;
#               NA for sequences removed by filterByExpr, FDR NA when a condition has < 2 replicates

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
n_rows <- as.integer(args[2])
n_cond1 <- as.integer(args[3])
n_cond2 <- as.integer(args[4])
lib_sizes <- as.numeric(strsplit(args[5], ",")[[1]])
out_file <- args[6]

counts <- matrix(readBin(counts_file, "integer", n = n_rows * (n_cond1 + n_cond2), size = 4), nrow = n_rows)
group <- factor(rep(c("cond1", "cond2"), c(n_cond1, n_cond2)))
design <- model.matrix(~group)

# Standard edgeR workflow: filter low-abundance sequences, TMM normalization
y <- DGEList(counts = counts, lib.size = lib_sizes, group = group)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)

if (min(n_cond1, n_cond2) >= 2) {
    # Quasi-likelihood F-test (needs replicates to estimate dispersion)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    res <- glmQLFTest(fit, coef = 2)
    fdr <- p.adjust(res$table$PValue, method = "BH")
} else {
    # Without replicates no dispersion can be estimated: fold changes only, no test
    fit <- glmFit(y, design, dispersion = 0)
    res <- glmLRT(fit, coef = 2)
    fdr <- rep(NA_real_, nrow(y))
}

out <- matrix(NA_real_, n_rows, 4)
out[, 1] <- keep
out[keep, 2] <- res$table$logFC
out[keep, 3] <- res$table$logCPM
out[keep, 4] <- fdr
writeBin(as.vector(out), out_file, size = 8)
