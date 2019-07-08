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

script_name = 'struct2rcontacts.R'
parser = arg_parser('Read 3D structures and extract pair-wise bead distances.
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
parser = add_argument(parser, "--save-single",
	"Save single structure distances.", flag = TRUE)
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

p = parse_args(parser)
attach(p['' != names(p)])

outDir = NA
if ( save_single ) outDir = rootDir
if ( is.na(description) ) description = basename(rootDir)

cat(sprintf("
 # %s

          Root : %s
         Label : %s
         Beads : %e
      Contacts : %s
         GPSeq : %s
       Threads : %d
   Save single : %s
        Descr. : %s
\n",
	script_name, rootDir, labPath, beadSize, contactLab,
	with_gpseq, threads, save_single, description))

# FUNCTIONS ====================================================================

get_bead_distances = function(ssData, rootDir = NA) {
	pairIDs = data.table(expand.grid(1:nrow(ssData), 1:nrow(ssData)))[Var1 < Var2]
	dData = cbind(pairIDs, data.table(d3d = dist(ssData[, .(x, y, z)])))
	setnames(dData, c("Var1", "Var2"), c("A", "B"))
	remove("pairIDs")

	if ( !is.na(rootDir) )
		saveRDS(dData, file.path(rootDir, ssData[1, structure], "bead_dist.rds"))

	return(dData[, d3d])
}

read_structure = function(spath) {
	return(fread(spath, col.names = c("x", "y", "z", "r", "i")))
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
		n = 3043, size = 1, endian = "little"))
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

cat("Calculating pair-wise distances...\n")
dData = do.call(cbind, pblapply(sData,
	get_bead_distances, outDir, cl = threads))

cat("Calculating mean/median pair-wise distances...\n")
dData = rbindlist(pbapply(dData, 1, function(x) data.table(
	d3dmean = mean(x, na.rm = T), d3dmedian = median(x, na.rm = T)),
	cl = threads))

cat("Retrieving bead labels...\n")
bLabs = read_bead_labels(labPath)
pairIDs = data.table(expand.grid(1:nrow(bLabs), 1:nrow(bLabs)))[Var1 < Var2]
dData = cbind(pairIDs, dData)
setnames(dData, c("Var1", "Var2"), c("A", "B"))
setkeyv(dData, c("A", "B"))

cat("Labeling Distances...\n")
distances = data.table(expand.grid(1:nrow(bLabs), 1:nrow(bLabs)))[Var1 < Var2]
distances = cbind(bLabs[distances$Var1], bLabs[distances$Var2], distances)
setnames(distances, c(names(bLabs), paste0(names(bLabs), 2), "A", "B"))
setkeyv(distances, c("A", "B"))
distances = distances[cData]

cat("Writing output...\n")
saveRDS(distances, file.path(dirname(rootDir),
	sprintf("%s.dist.rds", description)))

# END --------------------------------------------------------------------------

################################################################################
