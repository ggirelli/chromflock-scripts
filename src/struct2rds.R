#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20190704
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(require(argparser))
suppressMessages(require(data.table))
suppressMessages(require(pbapply))

# INPUT ========================================================================

script_name = 'struct2rds.R'
parser = arg_parser('Read 3D structures and assemble them into a single RDS.
It expects as input the path to a folder with "cf_NNN" subfolders, each
containing a "coords.csv" file. Also, a uint8 label file is required, alongside
bead size in nt, and Hi-C contact type.', name = script_name)

parser = add_argument(parser, 'rootDir', 'Path to folder with structure data.')
parser = add_argument(parser, 'labPath', 'Path to uint8 label file.')
parser = add_argument(parser, 'beadSize', 'Bead size in nt.', type = class(0))
parser = add_argument(parser, 'contactLab',
	'Label for utilized contacts, e.g., all, intra, inter,...')

parser = add_argument(parser, '--description',
	'A short dataset label. Defaults to rootDir basename.', default = NA)
parser = add_argument(parser, "--with-gpseq",
	"GPSeq was included in the chromflock run.", flag = TRUE)
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

p = parse_args(parser)
attach(p['' != names(p)])

if ( is.na(description) ) description = basename(rootDir)

cat(sprintf("
 # %s

     Root : %s
    Label : %s
    Beads : %e
 Contacts : %s
    GPSeq : %s
  Threads : %d
   Descr. : %s
\n",
	script_name, rootDir, labPath, beadSize, contactLab,
	with_gpseq, threads, description))

# FUNCTIONS ====================================================================

read_structure = function(spath, lpath) {
	dt = fread(spath)
	colnames(dt) = c("x", "y", "z", "r", "i")[1:ncol(dt)]

	dt$chromID = readBin(file(lpath, 'rb'), integer(),
		n = 3043, size = 1, endian = "little")

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
	dt = dt[, .(chrom, start, end, x, y, z, r, chromID)]
	setkeyv(dt, c("chrom", "start", "end"))

	return(dt)
}

read_all_structures = function(dpath, lpath, nthreads = 1) {
	cat(sprintf(" Reading structures from '%s'\n", dpath))
	pboptions(type = "txt")
	out = rbindlist(pblapply(list.files(dpath, "cf_.*", full.names = T),
		function(dpath) {
			name = basename(dpath)
			spath = file.path(dpath, "coords.csv")
			s = read_structure(spath, lpath)
			s$structure = name
			return(s)
		}, cl = nthreads
	))
	setkeyv(out, c("chrom", "start", "end"))
	return(out)
}

# RUN ==========================================================================

structData = read_all_structures(rootDir, labPath, threads)[,
	.(rmean = mean(r, na.rm = T), rmedian = median(r, na.rm = T)),
	by = c("chrom", "start", "end", "chromID")]

structData[, n := beadSize]
structData[, gpseq := with_gpseq]
structData[, contacts := contactLab]
structData[, label := description]

cat(" Writing output...\n")
saveRDS(structData, file.path(dirname(rootDir),
	sprintf("%s.radial.rds", description)))

# END --------------------------------------------------------------------------

################################################################################
