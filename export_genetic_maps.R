args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(AlphaSimR))
files <- list.files(args[2],pattern="^founders[.]rds$",recursive=TRUE,full.names=TRUE)
for(file in files) {
  state <- readRDS(file);sp <- state$sp
  folder <- file.path(dirname(file),"initial")
  chip <- sp$snpChips[[1]]
  chr <- rep(seq_len(sp$nChr),chip@lociPerChr)
  position <- vapply(seq_along(chr),function(i) sp$genMap[[chr[i]]][chip@lociLoc[i]],numeric(1))
  markers <- data.frame(marker=paste0("SNP",seq_along(chr)),chromosome=chr,segregating_site=chip@lociLoc,
    position_morgan=position,coordinate_status="synthetic genetic map; not physical assembly bp")
  write.csv(markers,file.path(folder,"marker_map.csv"),row.names=FALSE)
  qtl <- do.call(rbind,lapply(seq_len(sp$nTraits),function(k) {
    trait <- sp$traits[[k]];chrom <- rep(seq_len(sp$nChr),trait@lociPerChr)
    data.frame(trait=unname(sp$traitNames[k]),chromosome=chrom,segregating_site=trait@lociLoc,
      position_morgan=vapply(seq_along(chrom),function(i) sp$genMap[[chrom[i]]][trait@lociLoc[i]],numeric(1)),
      additive_effect_internal=trait@addEff,dominance_effect_internal=trait@domEff,
      fixed_initial_cohort_sd=unname(state$reference_sd[k]),role="simulation truth only")
  }))
  write.csv(qtl,file.path(folder,"QTL_effects_TRUTH_ONLY.csv"),row.names=FALSE)
}
cat("Exported",length(files),"frozen genetic maps and truth-effect records\n")
