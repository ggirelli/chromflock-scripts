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

script_name = 'struct2dist.R'
parser = arg_parser('Read 3D structures and extract mean pairwise bead distance.
It expects as input the path to a folder with "cf_NNN" subfolders, each
containing a "coords.csv" file. Also, a uint8 label file is required, alongside
bead size in nt', name = script_name)

parser = add_argument(parser, 'rootDir', 'Path to folder with structure data.')
parser = add_argument(parser, 'labPath', 'Path to uint8 label file.')
parser = add_argument(parser, 'beadSize', 'Bead size in nt.', type = class(0))
parser = add_argument(parser, 'contactLab',
	'Label for utilized contacts, e.g., all, intra, inter,...')
parser = add_argument(parser, 'nBeads', 'Number of beads', type=class(0))

parser = add_argument(parser, '--description',
	'A short dataset label. Defaults to rootDir basename.', default = NA)
parser = add_argument(parser, "--with-gpseq",
	"GPSeq was included in the chromflock run.", flag = TRUE)
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

p = parse_args(parser)
attach(p['' != names(p)])

# rootDir = "/mnt/data/chromflock/20190704-noy_new3Trans/10000G"
# labPath = "/mnt/data/chromflock/20190704-noy_new3Trans/L_noy_new3Trans.uint8"
# beadSize = 1e6
# contactLab = "all"
# description = NA
# with_gpseq = T
# threads = 20

outDir = NA
if ( is.na(description) ) description = basename(rootDir)

cat(sprintf("
 # %s

     Root : %s
    Label : %s
    Beads : %e (%d)
 Contacts : %s
    GPSeq : %s
  Threads : %d
   Descr. : %s
\n",
	script_name, rootDir, labPath, beadSize, nBeads, contactLab, with_gpseq,
	threads, description))

# FUNCTIONS ====================================================================

get_bead_distances = function(ssData, rootDir = NA) {
	pairIDs = data.table(expand.grid(1:nrow(ssData), 1:nrow(ssData)))[Var1 < Var2]
	pairIDs = pairIDs[order(Var2)][order(Var1)]
	dData = cbind(pairIDs, data.table(d3d = dist(ssData[, .(x, y, z)])))
	setnames(dData, c("Var1", "Var2"), c("A", "B"))
	remove("pairIDs")

	if ( !is.na(rootDir) )
		saveRDS(dData, file.path(rootDir, ssData[1, structure], "bead_dist.rds"))

	return(dData[, d3d])
}

read_structure = function(spath) {
	out = fread(spath)
	colnames(out) = c("x", "y", "z", "r", "i")[1:ncol(out)]
	return(out)
}

read_all_structures = function(dpath, nthreads = 1) {
	cat(sprintf(" Reading structures from '%s'\n", dpath))
	return(pblapply(list.files(dpath, "cf_.*", full.names = T),
		function(dpath) {
			name = basename(dpath)
			spath = file.path(dpath, "coords.csv")
			s = read_structure(spath)
			s$structure = name
			return(s)
		}, cl = nthreads
	))
}

read_bead_labels = function(lpath) {
	LH = file(lpath, 'rb')
	dt = data.table(chromID = readBin(LH, integer(),
		n = nBeads, size = 1, endian = "little"))
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

cat("Reading structures...\n")
sData = read_all_structures(rootDir, threads)

cat("Checking structures size...\n")

nBeads = unique(unlist(lapply(sData, nrow)))
stopifnot(length(nBeads) == 1)

cat("Calculating average pair-wise bead distance\n")
pb <- txtProgressBar(min = 1, max = length(sData), style = 3)
mDist = dist(sData[[1]][, .(x, y, z)])
setTxtProgressBar(pb, 1)
for (i in 2:length(sData)) {
	mDist = mDist + dist(sData[[i]][, .(x, y, z)])
	setTxtProgressBar(pb, i)
}
mDist = mDist / length(sData)
dData = data.table(d3dmean = mDist)
remove("pb", "mDist")

cat("\nRetrieving bead labels...\n")
bLabs = read_bead_labels(labPath)
pairIDs = data.table(expand.grid(1:nrow(bLabs), 1:nrow(bLabs))
	)[Var1 < Var2][order(Var1, Var2)]
dData = cbind(pairIDs, dData)
setnames(dData, c("Var1", "Var2"), c("A", "B"))
setkeyv(dData, c("A", "B"))

cat("Labeling Distances...\n")
distances = cbind(bLabs[pairIDs$Var1], bLabs[pairIDs$Var2], pairIDs)
setnames(distances,
	c(paste0(names(bLabs), "A"), paste0(names(bLabs), "B"), "A", "B"))
setkeyv(distances, c("A", "B"))
distances = distances[dData]

distances[, size := beadSize]
distances[, gpseq := with_gpseq]
distances[, contacts := contactLab]
distances[, label := description]
distances[, c("A", "B") := NULL]

cat("Writing output...\n")
saveRDS(distances, file.path(dirname(rootDir),
	sprintf("%s.dist.rds", description)))

# END --------------------------------------------------------------------------

################################################################################
