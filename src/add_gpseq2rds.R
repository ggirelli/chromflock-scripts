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

script_name = 'add_gpseq2rds.R'
parser = arg_parser('Add GPSeq information to structure RDS',
	name = script_name)

parser = add_argument(parser, 'metaGPSeq', 'Path to GPSeq meta table.')
parser = add_argument(parser, 'rootGPSeq', 'Path to GPSeq score root folder.')
parser = add_argument(parser, 'outDir', 'Path to output folder.')

parser = add_argument(parser, '--rds',
	'REQUIRED. Path to structure RDS file(s).', nargs = Inf, type = class(""))
parser = add_argument(parser, '--step', 'GPSeq bin step. Default: 1e5',
	default = 1e5, type = class(0))
parser = add_argument(parser, "--interpolate",
	"Perform bead radial position interpolation based on GPSeq bins. (FISH)",
	flag = TRUE)
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

p = parse_args(parser)
attach(p['' != names(p)])

stopifnot(!is.na(rds))
stopifnot(0 != length(rds))
stopifnot(file.exists(metaGPSeq))
stopifnot(dir.exists(rootGPSeq))
stopifnot(dir.exists(outDir))

print(step)
cat(sprintf("
 # %s

         RDS : %s
  GPSeq root : %s
  GPSeq meta : %s
      Output : %s
        Step : %d
 Interpolate : %s
     Threads : %d
\n",
	script_name, paste(rds, collapse = "\n               "),
	rootGPSeq, metaGPSeq, outDir, step, interpolate, threads))

pboptions(type = "timer")

# FUNCTIONS ====================================================================

chrom2chromID = function(chrom, nchrom = 24, hetero = c("X", "Y")) {
	# Convert a chromosome name (e.g., "chr1", "chrX") to a numerical ID.
	if ( grepl(":", chrom) ) {
		return(floor(as.numeric(gsub(":", ".",
			substr(chrom, 4, nchar(chrom))))))
	} else {
		chromID = substr(chrom, 4, nchar(chrom))
		if ( chromID %in% hetero )
			chromID = nchrom - which(rev(hetero) == chromID) + 1
		chromID = as.numeric(chromID)
		stopifnot(!is.na(chromID))
		return(chromID)
	}
}

add_chrom_ID = function(data, key = "chrom") {
	# Add chromosome ID to a data.table.
	# key should be the name of the column with chromosome names.
	require(data.table)
	stopifnot(is.data.table(data))
	stopifnot("chrom" %in% colnames(data))

	cid_table = data.table(chrom = as.character(unique(data$chrom)))
	cid_table$chromID = unlist(lapply(cid_table$chrom, FUN = chrom2chromID))
	setkeyv(cid_table, "chrom")
	setkeyv(data, key)

	return(data[cid_table,, nomatch = 0])
}

average_and_rescale_gpseq = function(data) {
	# Averages and rescales GPSeq data
	suppressMessages(require(outliers))

	mData = data[, .(
			prob_g = mean(prob_g, na.rm = T),
			chromID = unique(chromID)
		), by = c("chrom", "start", "end")]

	mData[, prob_g := log2(prob_g)]
	mData[!is.na(prob_g), outlierStatus := scores(
		mData[!is.na(prob_g), prob_g], "iqr", lim=1.5)]

	lowest = mData[FALSE == outlierStatus, min(prob_g, na.rm = T)]
	mData[, prob_g := prob_g - lowest]
	highest = mData[FALSE == outlierStatus, max(prob_g, na.rm = T)]
	mData[, prob_g := prob_g / highest]

	mData[, prob_g := 2**prob_g]
	mData[, outlierStatus := NULL]

	return(mData)
}

read_gpseq_score_advanced = function(meta, path,
	binSize = NA, binStep = NA, groupSize = 0, csm = 3,
	custom = F, ftype = "rescaled", scores = c("prob_g")
) {
	# 
	# "meta" can be either a path to a metadata file or a data.table containing
	# the metadata. The following columns are required: dataset, flag, and
	# cFolder. "dataset" is the dataset identifier, it should be the name of a
	# subfolder in path. "cFolder" is the name of the subfolder in path/dataset
	# (e.g., "10-ON" or "all") which contains the GPSeq score tracks. "flag" is
	# used to distinguish different tracks generated with similar parameters.
	# 
	# "path" is the root GPSeq score folder path. In the scope of this repo, it
	# should be stored in config_paths$gpseq_centrality_path$local.
	# 
	# If only meta and path are provided, chromosome-wide scores are loaded.
	# Same behavior if either binSize or binStep are NA.
	# 
	# Specify both binSize and binStep to load that specific binning.
	# 
	# Set custom to TRUE to load a customBin track from the specified folder.
	# If multiple customBin were applied, they should be distinguished based on 
	# the flag column from meta.
	# 
	# "csm" stands for CutSite Mode and can be one of the following: 1 for
	# "universal", 2 for "union", 3 for "separate" (default), 4 for
	# "intersection". Check the help page of gpseqc_estimate for more details.
	# 
	# "ftype" can be one of the following: "combined" to load assembled bed
	# files with each region/cutsite repeated per condition, "estimated" to
	# load the GPSeq score track before rescaling, "rescaled" to load the
	# rescaled GPSeq score track.
	# 
	
	require(data.table)

	stopifnot(csm %in% 1:4)
	stopifnot(ftype %in% c("combined", "estimated", "rescaled"))

	if ( is.character(meta) & length(meta) == 1 ) {
		stopifnot(file.exists(meta))
		meta = fread(meta)
	}

	stopifnot(is.data.table(meta))

	required_columns = c("dataset", "flag", "cFolder")
	for (reqCol in required_columns) {
		if ( ! reqCol %in% colnames(meta) )
			stop(sprintf("Missing %s column from metadata.", reqCol))
	}

	binFlag = "bins"
	if ( custom ) {
		binFlag = "customBins"
	} else {
		if ( is.na(binSize + binStep) ) {
			binFlag = sprintf("%s.chrWide", binFlag)
		} else {
			binFlag = sprintf("%s.size%d.step%d", binFlag, binSize, binStep)
		}
	}

	groupFlag = ""
	if ( groupSize != 0 ) {
		groupFlag = sprintf(".group%d", groupSize)
	}

	meta$flag[is.na(meta$flag)] = ""

	gpsq = rbindlist(lapply(1:nrow(meta), FUN = function(i) {
		metai = meta[i, ]

		if ( 0 != nchar(metai$flag) )
			if ( !grepl("^\\.", metai$flag) )
				metai$flag = sprintf(".%s", metai$flag)

		cond_dir_path = file.path(path, metai$dataset, metai$cFolder, "/")
		if ( !dir.exists(cond_dir_path) ) {
			print(paste0("Skipped metadata row #", i,
				": condition folder not found ('", cond_dir_path, "')"))
			print(metai)
			return(NULL)
		}

		tablePath = file.path(path, metai$dataset, metai$cFolder, 
			sprintf("%s%s.%s.%s%s.csm%d.%s", metai$dataset, metai$flag, ftype,
				binFlag, groupFlag, csm, "rmOutliers_chi2.rmAllOutliers.tsv"))
		if ( !file.exists(tablePath) ) {
			print(paste0("Skipped metadata row #", i,
				": table not found ('", tablePath, "')"))
			print(metai)
			return(NULL)
		}

		data = fread(tablePath)
		if ( ! "combined" == ftype ) {
			selected_columns = c(1:3, which(colnames(data) %in% scores))
			data = data[, ..selected_columns]
		}

		data = add_chrom_ID(data)

		for ( col in colnames(metai) ) {
			eval(parse(text = sprintf("data$meta_%s = metai$%s", col, col)))
		}

		data = data.table(data)
		setkeyv(data, c("chrom", "start", "end"))

		data$csm = csm

		return(data)
	}))

	return(gpsq)
}

read_gpseq_score = function(meta, path,
	binSize = NA, binStep = NA, custom = F, ftype = "rescaled") {
	# Function to read GPSeq centrality scores.
	# By default, it assumes group size of 0 and separate cutsite mode.
	# To change group size or cutsite domain use read_gpseq_score_advanced.
	# 
	# See read_gpseq_score_advanced description for more details.
	return(read_gpseq_score_advanced(meta, path, binSize, binStep,
		custom = custom, ftype = ftype))
}

# RUN ==========================================================================

cat("Reading RDS...\n")
structData = rbindlist(pblapply(rds, readRDS, cl = threads))

req_cols = c("chrom", "start", "end", "chromID", "rmean", "rmedian",
	"n", "gpseq", "contacts", "label")
stopifnot(all(req_cols %in% names(structData)))

setkeyv(structData, c("chrom", "start", "end"))

if ( 1 != length(structData[, unique(n)]) )
	stop(sprintf("Multiple bead sizes detected: %s", structData[, unique(n)]))
beadSize = structData[, unique(n)]

cat("Reading GPSeq...\n")
	gMeta = fread(metaGPSeq, key = "condition")

if (interpolate) {
	gData = read_gpseq_score(gMeta, rootGPSeq, custom = T)
	gData = average_and_rescale_gpseq(gData)
	setkeyv(gData, c("chrom", "start", "end"))
	gData[, mid := (end+start)/2]

	cat("Combining (w/ interpolation)...\n")
	fData = rbindlist(pblapply(1:nrow(gData), function(rid) {
		row = gData[rid,]
		out = as.data.table(do.call(cbind, as.list(by(structData, structData$label,
			function(st) {
				ct = copy(st[chrom == row$chrom])
				ct[, mid := (start+end)/2]
				
				midDiff = ct$mid - row$mid
				borderIDs = c(last(which(midDiff < 0)), which(midDiff > 0)[1])

				deltaR = ct[borderIDs, diff(rmedian)]

				ct[borderIDs[1], rmedian]+deltaR/beadSize*midDiff[borderIDs[1]]
			}
		))))
		return(cbind(row, out))
	}, cl = threads))

	cat("Writing output...\n")
	saveRDS(structData, file.path(outDir, "radial.prob_g.interpol.rds"))
} else {
	gData = read_gpseq_score(gMeta, rootGPSeq, beadSize, step)
	gData = average_and_rescale_gpseq(gData)
	setkeyv(gData, c("chrom", "start", "end"))

	cat("Combining...\n")
	structData = gData[structData]

	cat("Writing output...\n")
	saveRDS(structData, file.path(outDir, "radial.prob_g.rds"))
}

# END --------------------------------------------------------------------------

################################################################################
