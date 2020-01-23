getAndParseGSEMatrices=function (GEO, destdir, AnnotGPL, getGPL = TRUE, parseCharacteristics = TRUE) 
{
  GEO <- toupper(GEO)
  stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
  b = getDirListing(sprintf(gdsurl, stub, GEO))
  b=b[-1]
  message(sprintf("Found %d file(s)", length(b)))
  ret <- list()
  for (i in 1:length(b)) {
    message(b[i])
    destfile = file.path(destdir, b[i])
    if (file.exists(destfile)) {
      message(sprintf("Using locally cached version: %s", 
                      destfile))
    }
    else {
      download.file(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", 
                            stub, GEO, b[i]), destfile = destfile, mode = "wb", 
                    method = getOption("download.file.method.GEOquery"))
    }
    ret[[b[i]]] <- parseGSEMatrix(destfile, destdir = destdir, 
                                  AnnotGPL = AnnotGPL, getGPL = getGPL)$eset
  }
  return(ret)
}
environment(getAndParseGSEMatrices)<-asNamespace("GEOquery")
assignInNamespace("getAndParseGSEMatrices", getAndParseGSEMatrices, ns="GEOquery")
