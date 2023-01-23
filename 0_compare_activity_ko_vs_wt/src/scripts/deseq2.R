library(DESeq2)

## good example ##
# https://lashlock.github.io/compbio/R_presentation.html

# get the counts file
args = commandArgs(trailingOnly=TRUE)
# check to see that two argument is given
if (length(args)!=4) {
  stop("Four arguments must be supplied (input file).n", call.=FALSE)
}

infilename = args[1]
designmatfilename = args[2]
outfilename = args[3]
contraststr = args[4]

# read count file
tabla <- read.table(infilename, sep=",", row.names=1, header=TRUE)

# design matrix or treatment data construction
design_table <- read.table(designmatfilename, sep=",", row.names=1, header=TRUE)
design_table[] <- lapply(design_table, factor)

# construct deseq2 object
dds <- DESeqDataSetFromMatrix(countData = tabla,
                              colData = design_table,
                              design = ~ condition)

# run DESeq pipeline, will normalize using median of ratios
dds <- DESeq(dds)

# get the results
contrast_vec = unlist(strsplit(contraststr, ","))
res <- results(dds, contrast=contrast_vec)

# save to file .. 
write.table(res, file=outfilename, sep=",", row.names=TRUE, col.names=TRUE)
