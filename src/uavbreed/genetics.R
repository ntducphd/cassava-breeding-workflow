args <- commandArgs(trailingOnly=TRUE)
action <- args[1]; root <- args[2]; config_path <- args[3]
crop <- args[4]; seed <- as.integer(args[5]); state_path <- args[6]
export_dir <- args[7]; rlib <- args[8]
.libPaths(c(rlib, .libPaths()))
suppressPackageStartupMessages(library(AlphaSimR))
suppressPackageStartupMessages(library(jsonlite))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
cfg <- fromJSON(config_path, simplifyVector=TRUE)
c <- cfg$crops[[crop]]
dir.create(export_dir,recursive=TRUE,showWarnings=FALSE)

founders <- function(n, inbred) {
  maps <- replicate(c$chromosomes,seq(0,c$map_morgan,length.out=cfg$genetics$sites_per_chr),simplify=FALSE)
  hap <- lapply(maps,function(map) {
    nr <- n * if(inbred) 1L else 2L
    x <- matrix(0L,nr,length(map)); x[,1] <- rbinom(nr,1,.5)
    for(j in 2:length(map)) {
      renew <- runif(nr) > exp(-diff(map)[j-1]/cfg$genetics$ld_block_morgan)
      x[,j] <- x[,j-1]; x[renew,j] <- rbinom(sum(renew),1,.5)
    }
    x
  })
  newMapPop(maps,hap,inbred=inbred,ploidy=2L)
}

hybrids <- function(parents, groups, sp) {
  a <- which(groups=="A"); b <- which(groups=="B")
  pairs <- as.matrix(expand.grid(a,b))
  p <- makeCross(parents,pairs,nProgeny=1L,simParam=sp)
  list(pop=p,family=paste0("A",pairs[,1]),mother=parents@id[pairs[,1]],father=parents@id[pairs[,2]],S=NULL)
}

outcross <- function(parents, sp, s, pairs=NULL) {
  n <- nInd(parents)
  if(is.null(pairs)) {
    pairs <- matrix(integer(),0,2)
    for(k in seq_len(cfg$genetics$families)) {
      options <- t(combn(seq_len(n),2))
      if(crop=="coffee") options <- options[vapply(seq_len(nrow(options)),function(i) any(!s[options[i,2],]%in%s[options[i,1],]),logical(1)),,drop=FALSE]
      if(!nrow(options)) stop("No compatible mating available")
      pairs <- rbind(pairs,options[sample.int(nrow(options),1),])
    }
  }
  p <- makeCross(parents,pairs,nProgeny=cfg$genetics$progeny_per_family,simParam=sp)
  repeat_pairs <- pairs[rep(seq_len(nrow(pairs)),each=cfg$genetics$progeny_per_family),,drop=FALSE]
  childS <- NULL
  if(crop=="coffee") childS <- t(vapply(seq_len(nrow(repeat_pairs)),function(i) {
    female <- repeat_pairs[i,1]; male <- repeat_pairs[i,2]
    accepted <- s[male,][!s[male,]%in%s[female,]]
    c(s[female,sample.int(2,1)],accepted[sample.int(length(accepted),1)])
  },integer(2)))
  list(pop=p,family=paste0("F",rep(seq_len(nrow(pairs)),each=cfg$genetics$progeny_per_family)),mother=parents@id[repeat_pairs[,1]],father=parents@id[repeat_pairs[,2]],S=childS)
}

if(action=="init") {
  inbred <- crop=="maize"
  f <- founders(if(inbred) 2L*cfg$genetics$maize_pool_size else cfg$genetics$founders,inbred)
  sp <- SimParam$new(f); sp$nThreads <- 1L
  sp$restrSegSites(minQtlPerChr=cfg$genetics$qtl_per_chr,minSnpPerChr=cfg$genetics$markers_per_chr,overlap=FALSE)
  sp$addTraitAD(cfg$genetics$qtl_per_chr,mean=rep(0,4),var=rep(1,4),
                meanDD=cfg$genetics$mean_dominance_degree,varDD=cfg$genetics$variance_dominance_degree,
                corA=diag(4),corDD=diag(4),useVarA=TRUE,name=c("yield","quality","water_response","resistance"))
  sp$addSnpChip(cfg$genetics$markers_per_chr)
  parents <- newPop(f,simParam=sp)
  if(inbred) {
    groups <- rep(c("A","B"),each=cfg$genetics$maize_pool_size)
    st <- hybrids(parents,groups,sp);st$parents <- parents;st$groups <- groups
  } else {
    s <- if(crop=="coffee") t(replicate(nInd(parents),sample.int(12,2))) else NULL
    st <- outcross(parents,sp,s)
  }
  st$sp <- sp;st$cycle <- 0L
  st$reference_mean <- colMeans(gv(st$pop))
  st$reference_sd <- apply(gv(st$pop),2,sd)
  st$base_marker_frequency <- colMeans(pullSnpGeno(st$pop,simParam=sp))/2
  st$seed <- seed
} else {
  st <- readRDS(state_path);sp <- st$sp
  choice <- fromJSON(file.path(export_dir,"parent_selection.json"))
  if(crop=="maize") {
    keep <- match(as.character(choice$parent_ids),st$parents@id)
    if(anyNA(keep)) stop("Unknown maize parent")
    old <- st$parents; pools <- st$groups
    nextp <- lapply(c("A","B"),function(group) {
      k <- keep[pools[keep]==group]
      if(length(k)<2L) stop("At least two parents per maize pool required")
      pairs <- matrix(sample(k,2L*cfg$genetics$maize_pool_size,replace=TRUE),ncol=2)
      for(i in seq_len(nrow(pairs))) while(pairs[i,1]==pairs[i,2]) pairs[i,2] <- sample(k,1)
      f1 <- makeCross(old,pairs,nProgeny=1L,simParam=sp)
      makeDH(f1,nDH=1L,simParam=sp)
    })
    parents <- c(nextp[[1]],nextp[[2]])
    groups <- rep(c("A","B"),each=cfg$genetics$maize_pool_size)
    newst <- hybrids(parents,groups,sp);newst$parents <- parents;newst$groups <- groups
  } else {
    k <- match(as.character(choice$parent_ids),st$pop@id)
    if(anyNA(k)) stop("Unknown selected parent")
    s <- if(crop=="coffee") st$S[k,,drop=FALSE] else NULL
    newst <- outcross(st$pop[k],sp,s)
  }
  newst$sp <- sp;newst$cycle <- st$cycle+1L
  newst$reference_mean <- st$reference_mean;newst$reference_sd <- st$reference_sd
  newst$base_marker_frequency <- st$base_marker_frequency;newst$seed <- st$seed
  st <- newst
}
saveRDS(st,state_path,version=3)
p <- st$pop
g <- gv(p); b <- bv(p,simParam=st$sp)
m <- pullSnpGeno(p,simParam=st$sp)
colnames(m) <- paste0("SNP",seq_len(ncol(m)))
traits <- data.frame(genotype_id=p@id,family=st$family,mother_id=st$mother,father_id=st$father)
for(i in seq_len(4)) {
  traits[[paste0("gv",i)]] <- (g[,i]-st$reference_mean[i])/st$reference_sd[i]
  traits[[paste0("bv",i)]] <- b[,i]/st$reference_sd[i]
}
if(!is.null(st$S)) {traits$S1<-st$S[,1];traits$S2<-st$S[,2]}
write.csv(traits,file.path(export_dir,"genetic_truth.csv"),row.names=FALSE)
write.csv(data.frame(genotype_id=p@id,m),file.path(export_dir,"markers.csv"),row.names=FALSE)
if(crop=="maize") {
  write.csv(data.frame(parent_id=st$parents@id,pool=st$groups),file.path(export_dir,"parents.csv"),row.names=FALSE)
  pm <- pullSnpGeno(st$parents,simParam=st$sp);colnames(pm)<-colnames(m)
  write.csv(data.frame(genotype_id=st$parents@id,pm),file.path(export_dir,"parent_markers.csv"),row.names=FALSE)
}
qtl <- st$sp$traits[[1]]; snp <- st$sp$snpChips[[1]]
qnames <- paste(rep(seq_len(st$sp$nChr),qtl@lociPerChr),qtl@lociLoc,sep="_")
snames <- paste(rep(seq_len(st$sp$nChr),snp@lociPerChr),snp@lociLoc,sep="_")
stopifnot(!length(intersect(qnames,snames)))
write_json(list(engine="AlphaSimR",version=as.character(packageVersion("AlphaSimR")),cycle=st$cycle,
  seed=seed,n_genotypes=nInd(p),chromosomes=st$sp$nChr,ploidy=2L,markers=ncol(m),qtl_per_trait=length(qnames),
  qtl_marker_overlap=length(intersect(qnames,snames)),founder_model="synthetic Markov-block haplotypes",
  reference_scaling="initial cohort fixed across cycles",genetic_map="synthetic Morgans; not physical assembly coordinates",
  R_version=R.version.string),file.path(export_dir,"genetics_audit.json"),auto_unbox=TRUE,pretty=TRUE)
