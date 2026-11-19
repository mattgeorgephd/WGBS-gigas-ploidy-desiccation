S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
pnum <- function(v) as.numeric(sub("^< *", "", v))
for (s in c("ploidy", "heat", "shared")) for (o in c("BP", "CC", "MF")) {
  a <- readRDS(file.path(S, sprintf("gentable_%s_%s_legacy.rds", s, o))); b <- readRDS(file.path(S, sprintf("gentable_%s_%s_exact.rds", s, o)))
  pa <- pnum(a$classic); pb <- pnum(b$classic)
  sa <- a$GO.ID[pa < 0.01]; sb <- b$GO.ID[pb < 0.01]
  cat(sprintf("%-6s %s: whole table identical %s | terms with p < 0.01: legacy reader %d, complete read %d, in both %d | legacy-only %d, complete-only %d | top ten in common %d\n",
              s, o, identical(a, b), length(sa), length(sb), length(intersect(sa, sb)), length(setdiff(sa, sb)), length(setdiff(sb, sa)),
              length(intersect(a$GO.ID[1:10], b$GO.ID[1:10]))))
}
