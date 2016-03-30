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

# Convenience function for when we haven't already found the gene in the table
find.distance <- function(gene, ...) {
  idx <- match(gene, entrez2genome$entrez)
  find.distance.idx(idx, ...)
}


#' Identify trans eQTL that are cis to a related gene
#' Returns the entrez ID of the cis related gene, or a negative number if none is found.
#' -1 means there were no related genes
#' -2 means there were related genes but they had no genomic coordinates
#' -3 means none of the related genes were cis to the eqtl gene

#' Side Note: For some reason there are a number of genes in the related genes
#' table that don't have any genomic coordinates. However, we can't just delete
#' them because there may be secondary linkages between genes that do have 
#' coordinates.
find.pseudo.cis <- function(idx, chrom, min, max=min, cis.threshold=1, report=FALSE) {  
  
  # All input vectors should be the same length!
  stopifnot(length(idx)==length(chrom), length(idx)==length(min), length(idx)==length(max))
 
  # If we have multiple loci, deal with one at a time
  if ( length(idx) > 1 )
    return ( sapply(1:length(idx), function (i) find.pseudo.cis(idx[i], chrom[i],
                                                                min[i], max[i],
                                                                cis.threshold=cis.threshold,
                                                                report=report)) )
  
  gene.id <- entrez2genome$entrez[idx]
  partners <- with(related.genes, partner[id==gene.id])
  
  # If there are no related genes, abort
  if ( length(partners) < 1 )
    return ( -2 )
  
  # Recursively collect partners of partners (excluding gene.id itself)
  plen <- 0
  while ( plen < (plen <- length(partners)) ) 
    partners <- setdiff(union(partners, with(related.genes, partner[id %in% partners])), gene.id)
  
  # Find partner genes in table
  p.idx <- match(partners, entrez2genome$entrez)
  p.idx <- p.idx[!is.na(p.idx)]
  
  # If no related genes have coordinates, abort
  if ( length(p.idx) < 1 )
    return ( -1 )
  
  # Find distances to all partner genes
  p.dist <- find.distance.idx(p.idx, chrom, min, max)
  
  if ( any(!is.na(p.dist) & abs(p.dist) < cis.threshold) ) {
    # Select the closest partner gene
    best.idx <- p.idx[which.min(abs(p.dist))]
    if ( report ) {
      cat('Found cis gene', entrez2genome$symbol[best.idx], 'for',
          entrez2genome$symbol[idx], 'trans hit (loci dist.', min(abs(p.dist)[!is.na(p.dist)]) , 'Mb).\n')
    }
    return ( entrez2genome$entrez[best.idx] )
  } 
  
  # If we get to this point, none of the related genes passed the cis threshold
  return ( -3 )
}



# Example run with test eqtl data

load('test_eqtl_data.RData')

load('related_genes.RData')
load('entrez2genome.RData')


e2g.idx <- match(eqtl.info$gene, entrez2genome$entrez)

# Calculate distance for all eqtl loci
eqtl.info$distance <- with(eqtl.info, find.distance.idx(e2g.idx, chrom, min.pos, max.pos))

# Cis threshold is 1Mb 
eqtl.info$cis <- with(eqtl.info, !is.na(distance) & abs(distance) <= 1)

# Identify trans eQTL and look for pseudo cis genes
trans.idx <- which(!eqtl.info$cis)
pseudo.results <-with(eqtl.info[trans.idx, ],
                      find.pseudo.cis(gene, chrom, min.pos, max.pos, report=TRUE))

# If pseudo.results is a postive number, we have a pseudo cis hit
eqtl.info$pseudo.cis <- FALSE
eqtl.info$pseudo.cis[trans.idx] <- pseudo.results > 0

# Select pseudo cis hits and calculate distances from loci to cis genes
pseudo.cis.eqtl <- subset(eqtl.info, pseudo.cis)
pseudo.cis.eqtl$cis.gene <- pseudo.results[pseudo.results > 0]
pseudo.cis.eqtl$pseudo.dist <- with(pseudo.cis.eqtl, find.distance(cis.gene, chrom, min.pos, max.pos))

print(pseudo.cis.eqtl)