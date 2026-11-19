dir.create("analysis/env/Rlib",recursive=TRUE,showWarnings=FALSE)
lib <- normalizePath("analysis/env/Rlib",winslash="/",mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
deps <- c("Rcpp","Rdpack","R6")
missing <- deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos="https://cloud.r-project.org",lib=lib,type=if(.Platform$OS.type=="windows") "binary" else "source")
if(getRversion()<"4.4" || getRversion()>="4.5") stop("The pinned Windows binaries require R4.4.x; install matching source versions for other runtimes")
packages <- c(AlphaSimR="2.1.0",jsonlite="2.0.0")
for(p in names(packages)) {
  present <- requireNamespace(p,quietly=TRUE) && as.character(packageVersion(p))==packages[[p]]
  if(!present) {
    if(.Platform$OS.type!="windows") stop("Install the pinned source versions for this platform")
    url <- sprintf("https://cloud.r-project.org/bin/windows/contrib/4.4/%s_%s.zip",p,packages[[p]])
    install.packages(url,repos=NULL,type="win.binary",lib=lib)
  }
}
stopifnot(as.character(packageVersion("AlphaSimR"))=="2.1.0",as.character(packageVersion("jsonlite"))=="2.0.0")
writeLines(capture.output(sessionInfo()),"analysis/env/R_session.txt")
