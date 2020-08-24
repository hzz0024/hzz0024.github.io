# load the Fish's exact test outcome
ch_ref_name = 'CH_REF1.fet'  
ch_ref = read.delim(ch_ref_name, header = FALSE, sep='\t')
log_p <- ch_ref$V6
log_p <- as.numeric(gsub("1:2=","",log_p))

hist(log_p)

chr <- ch_ref$V1[log_p>2]
pos <- ch_ref$V2[log_p>2]
length(chr)
cat(paste0(chr, "\t", pos, "\n"))
