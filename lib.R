

prepIntSites <- function(x){

  if(! 'timePointMonths' %in% names(x) | ! 'minTimePointDays' %in% names(x)) x <- bind_cols(x, expandTimePoints(x$timePoint))
  
  x <- subset(x, timePointDays >= minTimePointDays)
  if(nrow(x) == 0) stop('No sites remain after filtering on minTimePointDays')

  if(! 'patient' %in% names(x))   x$patient <- x$subject
  if(! 'trial' %in% names(x))     x$trial <- 'none'
  if(! 'seqnames' %in% names(x))  x$seqnames <- x$chromosome
  if(! 'start' %in% names(x))     x$start <- x$position
  if(! 'end' %in% names(x))       x$end <- x$position
  if(! 'posid' %in% names(x))     x$posid <- paste0(x$seqnames, x$strand, x$position)
  if(! 'siteLabel' %in% names(x)) x$siteLabel <- paste0(x$subject, '_', x$posid)
  if('chromosome' %in% names(x)) x$chromosome <- NULL
  
  x <- arrange(x, desc(estAbund))
  x <- x[! duplicated(x$posid),]
  
  x
}

expandTimePoints <- function (tps) {
  d <- tibble::tibble(tp = sub("_", ".", tps))
  d$n <- 1:nrow(d)
  d$timePointType <- stringr::str_match(base::toupper(d$tp),"[DMY]")
  d$timePointType[which(is.na(d$timePointType))] <- "X"
  
  d <- dplyr::bind_rows(lapply(split(d, d$timePointType), function(x) {
    n <- as.numeric(stringr::str_match(x$tp, "[\\d\\.]+")) * ifelse(grepl("\\-", x$tp), -1, 1)
    
    if (x$timePointType[1] == "D") {
      x$timePointMonths <- base::round(n/30.4167, digits = 0)
      x$timePointDays <- base::round(n, digits = 0)
    }
    else if (x$timePointType[1] == "M") {
      x$timePointMonths <- base::round(n, digits = 0)
      x$timePointDays <- base::round(n * 30.4167, digits = 0)
    }
    else if (x$timePointType[1] == "Y") {
      x$timePointMonths <- base::round(n * 12, digits = 0)
      x$timePointDays <- base::round(n * 365, digits = 0)
    }
    else {
      message("Warning - could not determine date unit for: ",
              paste0(unique(x$timePoint), collapse = ", "))
      x$timePointMonths <- n
      x$timePointDays <- n
    }
    x
  }))
  
  dplyr::arrange(d, n) %>% dplyr::select(timePointMonths, timePointDays)
}


createIntUCSCTrack <- function(d, abundCuts = c(5, 10, 20), 
                               posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
                               negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
                               title = 'intSites', outputFile = 'track.ucsc', visibility = 1, 
                               position = 'chr7:127471196-127495720', padSite = 0, siteLabel = NA){
  
  # Check function inputs.
  if(length(posColors) != length(negColors)) 
    stop('The pos and neg color vectors are not the same length.')
  
  if(length(abundCuts) != length(posColors) - 1) 
    stop('The number of aundance cut offs must be one less than the number of provided colors.')
  
  if(! all(c('start', 'end', 'strand', 'seqnames', 'estAbund') %in% names(d))) 
    stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.") 
  
  if(is.na(siteLabel) | ! siteLabel %in% names(d)) 
    stop('The siteLabel parameter is not defined or can not be found in your data.')
  
  # Cut the abundance data. Abundance bins will be used to look up color codes.
  cuts <- cut(d$estAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
  
  # Convert Hex color codes to RGB color codes. 
  posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
  negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
  
  # Create data fields needed for track table.
  d$score <- 0
  d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
  
  # Pad the site n NTs to increase visibility.
  if(padSite > 0){
    d$start <- floor(d$start - padSite/2)
    d$end   <- ceiling(d$end + padSite/2)
  }
  
  # Define track header.
  trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s\nbrowser position %s",
                       title, title, visibility, position)
  
  # Write out track table.
  write(trackHead, file = outputFile, append = FALSE)
  write.table(d[, c('seqnames', 'start', 'end', siteLabel, 'score', 'strand', 'start', 'end', 'color')], 
              sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
}
