



#' Previously I was treating the last two columns as "match length" and "match
#' score", but now after looking at the documentation for pslScore, these 
#' columns actually represent "score" followed by percentIdentity". Oops!
#' Anyway - for pct.score I need to divide by the length of the query which can
#' be guessed by looking at the max extent for a particular query sequence but
#' it is probably better to just get the lengths from the fasta file
#' 

# Read in a .score file
get.blat.scores <- function(scorefile, queryfile) {
  
  blat.cols <- c('chrom', 'start', 'end', 'query', 'score', 'pct.id')
  blat.score <- read.table(scorefile, sep='\t', col.names=blat.cols)
  query.cols <- c('id', 'length')
  query.info <- read.table(queryfile, sep='\t', col.names=query.cols)

  blat.query <- t(simplify2array(strsplit(blat.score$query, ':', fixed=TRUE)))
  query.pos <- t(simplify2array(strsplit(blat.query[, 2], '-', fixed=TRUE)))
  mode(query.pos) <- 'numeric'
  
  query.id <- t(simplify2array(strsplit(blat.query[, 1], '|', fixed=TRUE)))[, 4]
  
  q.size <- query.info$length[match(query.id, query.info$id)]
  
  blat.score$pct.query <- round(apply(query.pos, 1, diff)/q.size*100, 1)
  blat.score$pct.score <- round(blat.score$score/q.size*100, 1)
  blat.score$refseq <- query.id
  
  blat.score <- blat.score[with(blat.score, order(refseq, -pct.score, -pct.id)), ]
  
  # Some alignments are to patches, which I can't really use
  ucsc.chroms <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
  use.chroms <- c(1:22, 'X', 'Y', 'MT')
  
  blat.score <- subset(blat.score, chrom %in% ucsc.chroms)
  
  blat.score$chrom <- use.chroms[match(blat.score$chrom, ucsc.chroms)]
  
  rownames(blat.score) <- NULL
  return ( blat.score )
}

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
  
  return ( distance/1e6 )
}

# Convenience function for when we haven't already found the gene in the table
find.distance <- function(gene, ...) {
  idx <- match(gene, entrez2genome$entrez)
  find.distance.idx(idx, ...)
}


# Load scores from test.score
blat.scores <- get.blat.scores('test.score', 'all.human.rna.info')

# Remove matches where the score is less than 30% of the gene length
blat.matches <- subset(blat.scores, pct.score >= 30)

# Convert from refseq back to entrez genes
load('refseq2entrez.RData')
r2e.idx <- match(blat.matches$refseq, refseq2entrez$refseq)
blat.matches$entrez <- refseq2entrez$entrez[r2e.idx]
e2g.idx <- match(blat.matches$entrez, entrez2genome$entrez)

# Calculate distance
match.distance <- with(blat.matches, find.distance.idx(e2g.idx, chrom, start, end))
cis.match <- !is.na(match.distance) & abs(match.distance) == 0
blat.matches <- subset(blat.matches, !cis.match)

arrange.columns <- c('entrez', 'chrom', 'start', 'end', 'refseq', 'pct.id', 'pct.score', 'pct.query')

blat.matches <- blat.matches[, arrange.columns]
rownames(blat.matches) <- NULL

save(blat.matches, file='test_blat_matches.RData')