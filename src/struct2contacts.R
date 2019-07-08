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
parser = arg_parser('Read 3D structures and extract their captured contacts.
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
parser = add_argument(parser, "--sphere",
	"Sphere radius. Default: 1", default = 1, type = class(0))
parser = add_argument(parser, "--volume", type = class(0),
	"Fraction of sphere volume occupied by beads. Default: 0.2", default = 0.2)
parser = add_argument(parser, "--rthr", type = class(0),
	"Radius threshold to consider two beads in contact. Default: 4", default = 4)
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
 Sphere radius : %.3f
 Vol. fraction : %.3f
   R threshold : %.3f
       Threads : %d
        Descr. : %s
\n",
	script_name, rootDir, labPath, beadSize, contactLab, with_gpseq,
	sphere, volume, rthr, threads, description))

# FUNCTIONS ====================================================================

get_captured_contacts = function(ssData, sphere, volume, rthr) {
	pairIDs = data.table(expand.grid(1:nrow(ssData), 1:nrow(ssData)))[Var1 < Var2]
	dData = data.table(d3d = dist(ssData[, .(x, y, z)]))

	beadRadius = ((sphere**3)*volume/nrow(ssData))**(1/3)
	rthr_effective = beadRadius * rthr

	contactIDs = dData$d3d <= rthr_effective
	cData = cbind(pairIDs[contactIDs], dData[contactIDs])
	setnames(cData, c("Var1", "Var2"), c("A", "B"))

	cData[, structure := ssData[1, structure]]
	return(cData)
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

cat("Identifying contacts...\n")
cData = rbindlist(pblapply(sData, get_captured_contacts,
	sphere, volume, rthr, cl = threads))

cat("Counting contacts...\n")
setDTthreads(threads)
cData = cData[, .N, by = c("A", "B")]
setkeyv(cData, c("A", "B"))
setDTthreads(1)

cat("Retrieving bead labels...\n")
bLabs = read_bead_labels(labPath)

cat("Labeling contacts...\n")
contacts = data.table(expand.grid(1:nrow(bLabs), 1:nrow(bLabs)))[Var1 < Var2]
contacts = cbind(bLabs[contacts$Var1], bLabs[contacts$Var2], contacts)
setnames(contacts, c(names(bLabs), paste0(names(bLabs), 2), "A", "B"))
setkeyv(contacts, c("A", "B"))
contacts = contacts[cData]

contacts[, size := beadSize]
contacts[, gpseq := with_gpseq]
contacts[, contacts := contactLab]
contacts[, label := description]
contacts[, sphere := sphere]
contacts[, volume := volume]
contacts[, rthr := rthr]

cat("Writing output...\n")
saveRDS(contacts, file.path(dirname(rootDir),
	sprintf("%s.contacts.rds", description)))

# END --------------------------------------------------------------------------

################################################################################
