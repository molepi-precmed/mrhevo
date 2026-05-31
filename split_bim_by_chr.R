## split_bim_by_chr.R
## One-time script: split refpop/kg.2020.hg38.eur.bim into per-chromosome RDS files.
## Output: refpop/bim_by_chr/chr_N.rds  (N = 1..22, X, Y)
## Each file is a data.table with cols: chr, rsid, pos, a1, a2

library(data.table)

bfile   <- "refpop/kg.2020.hg38.eur"
out_dir <- "refpop/bim_by_chr"
dir.create(out_dir, showWarnings = FALSE)

message("Loading full bim (may take ~60s)...")
t0 <- proc.time()["elapsed"]
bim <- fread(paste0(bfile, ".bim"), header = FALSE,
             select = c(1L, 2L, 4L, 5L, 6L),
             col.names = c("chr", "rsid", "pos", "a1", "a2"))
bim[, global_idx := .I]   # 1-based column index in BEDMatrix
message(sprintf("  Loaded %d rows in %.1fs", nrow(bim), proc.time()["elapsed"] - t0))

chrs <- unique(bim$chr)
message(sprintf("Writing %d per-chromosome RDS files...", length(chrs)))
for (ch in chrs) {
    out_file <- file.path(out_dir, sprintf("chr_%s.rds", ch))
    saveRDS(bim[chr == ch], out_file, compress = FALSE)
}
message(sprintf("Done. Files written to %s/", out_dir))
