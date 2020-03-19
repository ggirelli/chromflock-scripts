#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190705
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

# script_name = 'L2tsv.R'
# parser = arg_parser('Lorem ipsum dolor sit amet, consectetur adipisicing elit.
# Maiores id qui molestias quo ex ab sit fugiat, consequuntur doloremque magnam
# beatae fugit laborum ipsum quidem voluptatibus laudantium esse quae optio!',
# 	name = script_name)

# parser = add_argument(parser, 'labPath', 'Path to uint8 label file.')
# parser = add_argument(parser, 'outPath', 'Path to output tsv file.')

# # rootDir = "/mnt/data/chromflock/20190704-noy_new3Trans/10000G"
# # labPath = "/mnt/data/chromflock/20190704-noy_new3Trans/L_noy_new3Trans.uint8"
beadSize = 1e6
# # contactLab = "all"
# # description = NA
# # with_gpseq = T
# # threads = 20
# # sphere = 1
# # volume = .2
# # rthr = 4

# p = parse_args(parser)
# attach(p['' != names(p)])

# if ( is.na(description) ) description = basename(rootDir)

#labPath = "/mnt/data/chromflock/HAP1_100kb/1000G/8.000000_1000_MAX.L.uint8"
#outPath = "/mnt/data/chromflock/HAP1_100kb/beadLabels.tsv"

labPath = "/mnt/data/chromflock/20191213_100kb/1000/18.000000_1000_MAX.L.uint8"
outPath = "/mnt/data/chromflock/20191213_100kb/1000/beadLabels.tsv"

# FUNCTIONS ====================================================================

read_bead_labels = function(lpath) {
	LH = file(lpath, 'rb')
	dt = data.table(chromID = readBin(LH, integer(),
		n = 30430, size = 1, endian = "little"))
	close(LH)

	dt = rbindlist(by(dt, dt$chromID, FUN = function(ct) {
		ct$start = (1:nrow(ct) - 1) * beadSize
		ct$end = ct$start + beadSize
		return(ct)
	}))

	dt[,chrom := as.character(chromID)]
	dt[chromID == 23, chrom := "X"]
	dt[,chrom := paste0("chr", chrom)]
	dt[chrom == "chr9", chrom := "chr9:22"]
	dt[chrom == "chr22", chrom := "chr22:9"]

	return(dt[, .(chrom, start, end, chromID)])
}

# RUN ==========================================================================

cat("Retrieving bead labels...\n")
bLabs = read_bead_labels(labPath)
write.table(bLabs, outPath, quote=F, sep="\t", row.names=F, col.names=T)

# END --------------------------------------------------------------------------

################################################################################
