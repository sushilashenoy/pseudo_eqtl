
load('test_blat_matches.RData')
colnames(blat.matches)[2] <- 'chromosome'
load('entrez2genome.RData')

# Define functions

# Calculate distance (in Mb) to gene(s) as identified by row indices in e2g table
# Negative distances are 3', Positive distances are 5', Zero means gene overlap
find.distance.idx <- function(idx, chrom, min, max=min, verbose=FALSE) {
  # All input vectors should be the same length, or length 1
  stopifnot(length(chrom)==1 || length(idx)==length(chrom),
            length(min)==1 || length(idx)==length(min),
            length(max)==1 || length(idx)==length(max))
  
  same.chrom <- entrez2genome$chromosome[idx] == chrom
  dist <- with(entrez2genome[idx, ], cbind(min-start, max-start, min-end, max-end))
  
  if ( verbose ) 
    cat(dim(dist), '\n')
  
  direction <- sign(dist[, 1])
  strand <- entrez2genome$orientation[idx]
  orient <- c(-1, 1)[match(strand, c('-', '+'))]
  overlaps <- apply(dist, 1, function (d) any(sign(d[1]) != sign(d[2:4])))
  same.chrom <- entrez2genome$chromosome[idx] == chrom
  
  min.abs.dist <- apply(abs(dist), 1, min)
  
  distance <- min.abs.dist * direction * (!overlaps) * orient
  distance[!same.chrom] <- NA
  
  
  if ( verbose ) cat('Gene-genotype distances:', dist, '\tOverlap:', overlaps,
                     '\tDirection:', direction, '\tOrientation:', orient,
                     '\tChrom match:', same.chrom, '\n')
  
  return ( distance/1e6 )
}


# Calculate distance (in Mb) to alignments as identified by row indices in blat.matches table
find.match.distance.idx <- function(idx, chrom, min, max=min, verbose=FALSE) {
  # All input vectors should be the same length, or length 1
  stopifnot(length(chrom)==1 || length(idx)==length(chrom),
            length(min)==1 || length(idx)==length(min),
            length(max)==1 || length(idx)==length(max))
  
  same.chrom <- blat.matches$chromosome[idx] == chrom
  dist <- with(blat.matches[idx, ], cbind(min-start, max-start, min-end, max-end))
  
  if ( verbose ) 
    cat(dim(dist), '\n')
  
  direction <- sign(dist[, 1])
  overlaps <- apply(dist, 1, function (d) any(sign(d[1]) != sign(d[2:4])))
  same.chrom <- blat.matches$chromosome[idx] == chrom
  
  min.abs.dist <- apply(abs(dist), 1, min)
  
  distance <- min.abs.dist * direction * (!overlaps) 
  distance[!same.chrom] <- NA
  
  
  if ( verbose ) cat('Alignment-genotype distances:', dist, '\tOverlap:', overlaps,
                     '\tDirection:', direction, '\tOrientation:', orient,
                     '\tChrom match:', same.chrom, '\n')
  
  return ( distance/1e6 )
}


#' Identify trans eQTL that are cis to a region that aligns to the eQTL gene
#' Returns the index of the match, or a negative number if none is found.
#' -1 means the gene was not in the table
#' -2 means none of the alignments were cis to the eqtl gene
#' 
find.blat.pseudo.cis <- function(gene, chrom, min, max=min, cis.threshold=1, report=FALSE) {  
  
  # All input vectors should be the same length!
  stopifnot(length(gene)==length(chrom), length(gene)==length(min), length(gene)==length(max))
 
  # If we have multiple loci, deal with one at a time
  if ( length(gene) > 1 )
    return ( sapply(1:length(gene), function (i) find.blat.pseudo.cis(gene[i], chrom[i],
                                                                min[i], max[i],
                                                                cis.threshold=cis.threshold,
                                                                report=report)) )

  m.idx <- which(blat.matches$entrez==gene)
  
  # If there are no matches, abort
  if ( length(m.idx) < 1 )
    return ( -1 )
  
  
  # Find distances to all partner genes
  p.dist <- find.match.distance.idx(m.idx, chrom, min, max)
  
  if ( any(!is.na(p.dist) & abs(p.dist) < cis.threshold) ) {
    # Select the closest partner gene
    best.idx <- m.idx[which.min(abs(p.dist))]
    if ( report ) {
      cat('Found cis alignment for',
          entrez2genome$symbol[match(gene, entrez2genome$entrez)], 'trans hit (loci dist.', min(abs(p.dist)[!is.na(p.dist)]) , 'Mb).\n')
    }
    return ( best.idx )
  } 
  
  # If we get to this point, none of the related genes passed the cis threshold
  return ( -2 )
}



# Example run with test eqtl data

load('test_eqtl_data.RData')



e2g.idx <- match(eqtl.info$gene, entrez2genome$entrez)

eqtl.info$symbol <- entrez2genome$symbol[e2g.idx]

# Calculate distance for all eqtl loci
eqtl.info$distance <- with(eqtl.info, find.distance.idx(e2g.idx, chrom, min.pos, max.pos))

# Cis threshold is 1Mb 
eqtl.info$cis <- with(eqtl.info, !is.na(distance) & abs(distance) <= 1)

# Identify trans eQTL and look for pseudo cis genes
trans.idx <- which(!eqtl.info$cis)
pseudo.results <-with(eqtl.info[trans.idx, ],
                      find.blat.pseudo.cis(gene, chrom, min.pos, max.pos, report=TRUE))

# If pseudo.results is a postive number, we have a pseudo cis hit
eqtl.info$pseudo.cis <- FALSE
eqtl.info$pseudo.cis[trans.idx] <- pseudo.results > 0

# Select pseudo cis hits and calculate distances from loci to cis genes
pseudo.blat.eqtl <- subset(eqtl.info, pseudo.cis)
pseudo.blat.eqtl$m.idx <- pseudo.results[pseudo.results > 0]
pseudo.blat.eqtl$pseudo.dist <- with(pseudo.blat.eqtl, find.match.distance.idx(m.idx, chrom, min.pos, max.pos))


pseudo.blat.eqtl$pct.score <- blat.matches$pct.score[pseudo.blat.eqtl$m.idx]
pseudo.blat.eqtl$pct.query <- blat.matches$pct.query[pseudo.blat.eqtl$m.idx]
pseudo.blat.eqtl$pct.id <- blat.matches$pct.id[pseudo.blat.eqtl$m.idx]

