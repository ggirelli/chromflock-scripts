#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190708
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(require(argparser))
suppressMessages(require(data.table))
suppressMessages(require(pbapply))
suppressMessages(require(pdist))

setDTthreads(1)
pboptions(type = "timer")

# INPUT ========================================================================

script_name = 'merge_juicer_matrices.R'
parser = arg_parser('Merge observed contacts previously extracted by juicer into
a single RDS file.', name = script_name)

parser = add_argument(parser, 'intraDir',
	'Path to folder with intra-chromosomal contacts.')
parser = add_argument(parser, 'interDir',
	'Path to folder with inter-chromosomal contacts.')
parser = add_argument(parser, 'outDir',
	'Path to output folder.')
parser = add_argument(parser, 'beadSize', 'Bead size in nt.', type = class(0))
parser = add_argument(parser, 'normType',
	'HiC data normalization type. e.g., none, kr,...', type = class(""))

parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

p = parse_args(parser)
attach(p['' != names(p)])

cat(sprintf("
 # %s

  Intra folder : %s
  Inter folder : %s
 Output folder : %s
 Normalization : %s
         Beads : %e
       Threads : %d

\n",
	script_name, intraDir, interDir, outDir, normType, beadSize, threads))

# FUNCTIONS ====================================================================

read_intra_matrix = function(fpath) {
	meta = as.data.table(t(unlist(strsplit(fpath, "\\."))[c(1, 3, 4, 5, 6)]))
	setnames(meta, c("accession", "chrom", "dtype", "ntype", "bsize"))
	intraData = fread(file.path(intraDir, fpath),
		col.names = c("startA", "startB", "value"))
	intraData[, c("chromA", "chromB") := list(meta$chrom, meta$chrom)]
	intraData = cbind(intraData, meta[, .(accession, dtype, ntype, bsize)])
	intraData = intraData[, .(chromA, startA, chromB, startB, value,
		accession, dtype, ntype, bsize)]
	return(intraData)
}

read_inter_matrix = function(fpath) {
	meta = as.data.table(t(unlist(strsplit(fpath, "\\."))[c(1, 3, 4, 5, 6, 7)]))
	setnames(meta,
		c("accession", "chromA", "chromB", "dtype", "ntype", "bsize"))
	interData = fread(file.path(interDir, fpath),
		col.names = c("startA", "startB", "value"))
	interData = cbind(interData, meta)[, .(chromA, startA, chromB, startB,
		value, accession, dtype, ntype, bsize)]
	return(interData)
}

# RUN ==========================================================================

flist = list.files(intraDir,
	sprintf("observed\\.%s\\.%d\\.", normType, beadSize))
cat(sprintf(" Reading %d intra-contact files...\n", length(flist)))
intraData = rbindlist(pblapply(flist, read_intra_matrix, cl = threads))

flist = list.files(interDir,
	sprintf("observed\\.%s\\.%d\\.", normType, beadSize))
cat(sprintf(" Reading %d inter-contact files...\n", length(flist)))
interData = rbindlist(pblapply(flist, read_inter_matrix, cl = threads))

cat(" Merging contact files...\n")
hicData = rbind(intraData, interData)[order(chromA, startA, chromB, startB)]

cat(" Writing output...\n")
saveRDS(hicData, file.path(outDir, sprintf("hic.merged.%d.rds", beadSize)))

# END --------------------------------------------------------------------------

################################################################################
