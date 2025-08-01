library(dplyr)
library(data.table)
library(parallel)
library(GenomicRanges)
source('lib.R')

args <- commandArgs(trailingOnly=TRUE)

expSites <- args[1]
refSites <- args[2]

expSites <- 'data/lateSCID3.tsv.gz'
refSites <- 'data/safetyPaperIntSites.tsv.gz'

CPUs <- 60
minTimePointDays <- 180
minSitesPerPatientChromosome <- 100
window_width <- 10000
window_step <- 5000
n_draws <- 5000
min_eval <- 0.001
min_SitesPerCluster <- 10
chromosomeLengths <- 'data/hg38_chromosome_lengths.tsv'
outputFile <- 'output.tsv'

refSites <- readr::read_tsv(refSites, show_col_types = FALSE) %>% as.data.table()
expSites <- prepIntSites(readr::read_tsv(expSites, show_col_types = FALSE))
expSites <- mutate(expSites, k = paste0(trial, '~', patient, '~', seqnames, '~', strand)) %>% as.data.table()
chromosomeLengths <- readr::read_tsv(chromosomeLengths, show_col_types = FALSE) %>% as.data.table()

expSites <- rbindlist(lapply(split(expSites, expSites$k), function(x){
              if(n_distinct(x$posid) < minSitesPerPatientChromosome) x <- data.table()
              x
            }))

windows <- rbindlist(lapply(split(chromosomeLengths, 1:nrow(chromosomeLengths)), function(x){
             start <- seq(1, x$length, window_step)
             stop <- seq(window_step, x$length, window_step)
             data.table(chromosome = x$chromosome, start = start[1:(length(start)-1)], stop = stop)
            }))

refSites <- subset(refSites, seqnames %in% windows$chromosome)
expSites <- subset(expSites, seqnames %in% windows$chromosome)

cluster <- parallel::makeCluster(CPUs)

clusterExport(cluster, c('refSites', 'windows', 'minSitesPerPatientChromosome', 'n_draws', 'min_SitesPerCluster'))

worker <- function(x, useFocusedWindows = FALSE){
       library(dplyr)
       library(data.table)
       library(GenomicRanges)
    
       # Limit windows to those overlapping with query sites.
       q <- makeGRangesFromDataFrame(x)
       
       if(useFocusedWindows){
         if(x$k[1] %in% names(focusedWindows)){
           w <- focusedWindows[[x$k[1]]]
         } else {
           return(data.table())
         }
       } else {
         chromosome_windows <- windows[chromosome == x$seqnames[1]]
         w <- makeGRangesFromDataFrame(chromosome_windows)
         chromosome_windows <- chromosome_windows[unique(subjectHits(findOverlaps(q, w)))]
         w <- makeGRangesFromDataFrame(chromosome_windows)
       }
  
       # Create a set of comp sites for this chromosome, same orientation, excluding this patient.
       setkey(refSites, seqnames, strand, patient)
       compSites <- refSites[seqnames == x$seqnames[1] & strand == x$strand[1] & patient != x$patient[1]]

       o <- findOverlaps(q, w, ignore.strand = TRUE)
       m <- as.data.table(w[subjectHits(o)])[, c(-4,-5)]
       q_tab <- cbind(m[, .(count1 = .N), by = .(seqnames, start, end)])
       
       s_tab <- rbindlist(lapply(1:n_draws, function(draw){
                  set.seed(draw)
                  z <- makeGRangesFromDataFrame(compSites[sample(.N, n_distinct(x$posid))])
                  o <- findOverlaps(z, w, ignore.strand = TRUE)
                  m <- as.data.table(w[subjectHits(o)])[, c(-4,-5)]
                  cbind(draw  = draw, m[, .(count2 = .N), by = .(seqnames, start, end)])
               }))
       
       s_tab <- s_tab[complete.cases(s_tab)]
       
       setkey(s_tab, seqnames, start, end)
       
       j <- s_tab[q_tab, on = .(seqnames, start, end)]
       
       result <- j[, .(
         count = count1[1],
         sims = n_draws,
         numSimsMatching = sum(count2 >= count1[1]),
         pVal = (sum(count2 >= count1[1]) + 1) / (n_draws + 1)
       ), by = .(seqnames, start, end)]
       
       result <- result[count >= min_SitesPerCluster]
       
       if(nrow(result) > 0){
         return(cbind(k = x$k[1], result))
       } else {
         return(data.table())
       }
}


r <- rbindlist(parLapply(cluster, split(expSites, expSites$k), worker, useFocusedWindows = FALSE))
  

if(nrow(r) > 0){
    focusedWindows <- lapply(split(r, r$k), function(x){
      x <- subset(x, pVal <= min_eval & count >= min_SitesPerCluster)
  
      if(nrow(x) > 0){
        if(nrow(x) > 1){
          return(GenomicRanges::reduce(makeGRangesFromDataFrame(x)))
        } else {
          return(makeGRangesFromDataFrame(x))
        }
      } else {
        return(GRanges())
      }
    })

    clusterExport(cluster, c('focusedWindows'))

    expSites <- subset(expSites, k %in% names(focusedWindows))

    r2 <- rbindlist(parLapply(cluster, split(expSites, expSites$k), worker, useFocusedWindows = TRUE))
} else {
    r2 <- data.table()
}
  
if(nrow(r2) > 0){
    r2$range <- paste(round((r2$end - r2$start + 1) / 1000, 1), 'KB')
    r2 <- r2[complete.cases(r2)]
    r2 <- tidyr::separate(r2, k, c('trial', 'patient', 'chromosome', 'strand'), sep = '~')
}
  
parallel::stopCluster(cluster)

readr::write_tsv(r2, outputFile)

# Set hgt.customText to URL serving one of the bed files created below.
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.customText=http://bushmanlab.org/tracks/Boston_SCID2_p00009.bed

invisible(lapply(split(expSites, paste(expSites$trial, expSites$patient)), function(x){
  createIntUCSCTrack(as.data.frame(x), 
                     title = paste0(x$trial[1], '_', x$patient[1], '_SOCS'), 
                     siteLabel = 'siteLabel', 
                     outputFile = paste0(x$trial[1], '_', x$patient[1], '_SOCS.bed'))
  
  system(paste('gzip -9', paste0(x$trial[1], '_', x$patient[1], '_SOCS.bed')))
}))
