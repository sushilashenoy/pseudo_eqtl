

if ( !file.exists('all.human.rna.fasta') ) {
  #' There used to be a single file available for download with all (Human)
  #' sequences in refseq, but I can't find it at the moment - also it was
  #' missing a number of refseq sequences. There may be a way to do this that
  #' doesn't involve downloading 25 separate files but this seems to work so...
  require('RCurl')
  url <- 'ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/'
  
  message('Downloading file list..')
  filenamelist <- getURL(url, ftp.use.epsv=FALSE, ftplistonly=TRUE)
  filenames <- strsplit(filenamelist, '\n', fixed=TRUE)[[1]]
  
  # We only care about the RNA sequences
  sfx <- '.rna.fna.gz'
  choose.files <- substring(filenames, nchar(filenames)-nchar(sfx)+1)==sfx
  dl.urls <- paste0(url, filenames[choose.files])
  dl.files <- filenames[choose.files]
  
  message('Downloading ', sum(choose.files), ' files..')
  # Download the files
  mapply(curl_download, dl.urls, dl.files)
  
  #' This function removes all of the linebreaks in the unzipped file, then adds
  #' back only those linebreaks that separate sequences from descriptions. This
  #' makes it a lot easier to extract sequence IDs.
  magic <- function(infile, outfile) {
    if ( file.exists(outfile) ) return
    system(paste('gzip -dc', infile, '| sed -e "s/\\(^>.*$\\)/#\\1#/"',
                 '| tr -d "\\r" | tr -d "\\n" | sed -e "s/$/#/"',
                 '| tr "#" "\\n" | sed -e "/^$/d" > ', outfile))
  }
  fix.file <- function(prefix) magic(paste0(prefix, '.rna.fna.gz'), paste0(prefix, '.rna.fasta'))
  
  message('Reformatting files..')
  # Unzip and reformat the files slightly
  sapply(paste0('human.', 1:25), fix.file)
  
  # Combine all the files into a single file and then delete the individual & downloaded files
  message('Combining files..')
  system('cat human.*.rna.fasta > all.human.rna.fasta')
  
  message('Removing intermediate files..')
  system('rm human.*.rna.fasta')
  system('rm human.*.rna.fna.gz')
  
} else {
  message('Found previously downloaded all.human.rna.fasta, file age is ',
          format(Sys.time() - file.mtime('all.human.rna.fasta')))
}

if ( !file.exists('all.human.refseq.ids') ) {
  # Extract refseq IDs
  system('sed "n;d;" < all.human.rna.fasta | cut -d" " -f1 | cut -d"|" -f4 > all.human.refseq.ids')
}

if ( file.exists('refseq2entrez.RData') ) {
  load("refseq2entrez.RData")
  refseq2entrez <- subset(refseq2entrez, refseq != '-' & status != 'SUPPRESSED')
  all.refseq <- sort(unique(refseq2entrez$refseq))
  
  # Extract refseq IDs
  system('sed "n;d;" < all.human.rna.fasta | cut -d" " -f1 | cut -d"|" -f4 > all.human.refseq.ids')
  dl.refseq.ids <- scan('all.human.refseq.ids', character())
  
  message('Found ', sum(!dl.refseq.ids %in% all.refseq), ' sequences with no associated gene.')
  # Not worth the effort to remove these few sequences at the moment
}