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
suppressMessages(require(cowplot))
suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(pbapply))
suppressMessages(require(viridis))

theme_set(theme_cowplot())
setDTthreads(1)
pboptions(type = "timer")

# INPUT ========================================================================

script_name = 'plot_struct.R'
parser = arg_parser('Generate plots comparing structures and GPSeq.',
	name = script_name)

parser = add_argument(parser, 'gRDS',
	'Path to input GPSeq rds generated with add_gpseq2rd.R')
parser = add_argument(parser, 'fRDS',
	'Path to input FISH rds generated with add_gpseq2rd.R')
parser = add_argument(parser, 'outDir', 'Path to output folder.')
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

p = parse_args(parser)
attach(p['' != names(p)])

cat(sprintf("
 # %s

 GPSeq RDS : %s
 iFISH RDS : %s
    Output : %s
   Threads : %d

", script_name, gRDS, fRDS, outDir, threads))

# FUNCTIONS ====================================================================

save_and_plot <- function(x, bname, width, height,
	dpi=300, use.cowplot=FALSE, ncol=1, nrow=1,
	base_height=4, base_width=NULL, base_aspect_ratio=1,
	plot=FALSE, formats = c("png", "pdf")){
  # Function to save the plot (and separately its data)
  # to file and show the plot in the notebook
  if( !use.cowplot ){
	if ( "png" %in% formats) {
		png(filename=file.path(paste0(bname, ".png")),
			units="in", height=height, width=width, res=dpi)
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "pdf" %in% formats) {
		cairo_pdf(filename=file.path(paste0(bname, "_cairo_pdf.pdf")),
			onefile = TRUE, height=height, width=width, family="Helvetica",
			pointsize=8, antialias="none")
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "eps" %in% formats) {
		cairo_ps(filename=file.path(paste0(bname, "_cairo_ps.eps")),
			onefile = TRUE, height=height, width=width, family="Helvetica",
			pointsize=8, antialias="none")
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
		postscript(file=file.path(paste0(bname, "_postscript.eps")),
			onefile = TRUE, paper="special", height=height, width=width,
			family="Helvetica", pointsize=8, horizontal=FALSE)
		print(x)
		while(length(dev.list())>0) invisible(dev.off())
	}
  }else{
	if ( "png" %in% formats) {
		save_plot(x, filename=file.path(paste0(bname, ".png")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			dpi=dpi)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "pdf" %in% formats) {
		save_plot(x, filename=file.path(paste0(bname, "_cairo_pdf.pdf")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			device=cairo_pdf)
		while(length(dev.list())>0) invisible(dev.off())
	}
	if ( "eps" %in% formats) {
		save_plot(x, filename=file.path(paste0(bname, "_cairo_ps.eps")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			device=cairo_ps)
		save_plot(x, filename=file.path(paste0(bname, "_postscript.eps")),
			ncol = ncol, nrow = nrow, base_height = base_height,
			base_width = base_width, base_aspect_ratio = base_aspect_ratio,
			device="ps")  
		while(length(dev.list())>0) invisible(dev.off())
	}
  }
  if( plot ) print(x)
}

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

plot_chrom_profile = function(cData, val.var = "rmedian") {
	corData = rbindlist(by(cData,
		cData$gpseqLab, function(gt) {
		rbindlist(by(gt, gt$contactLab, function(ct) {
			pcor = cor(ct$prob_g, 1-unlist(ct[, ..val.var]),
				use = "pairwise.complete.obs", method = "pearson")
			scor = cor(ct$prob_g, 1-unlist(ct[, ..val.var]),
				use = "pairwise.complete.obs", method = "spearman")
			data.table(pearson = pcor, spearman = scor,
				contactLab = ct[1, contactLab], gpseqLab = ct[1, gpseqLab])
		}))
	}))

	p = ggplot(cData, aes(x = (start+end)/2e6)
		) + geom_line(aes(y = prob_g, color = "GPSeq"), size = 1
		) + geom_line(aes(y = 1-rmedian, color = "structures"), size = 1
		) + geom_text(data = corData, aes(label = sprintf(
				"Pearson: %.3f\nSpearman: %.3f\nR2: %.3f", pearson, spearman, pearson**2
			)), x = 0, y = 1.25, hjust = 0, vjust = 1, size = 3, lineheight = .9
		) + facet_grid(~gpseqLab~contactLab
		) + scale_color_brewer(palette = "Paired"
		) + guides(color = guide_legend(title = "Distance from")
		) + theme(legend.position = "top",
			axis.text.x = element_text(angle = 45, hjust = 1)
		) + xlab("Bin midpoint genomic coordinate (Mb)"
		) + ylab("Distance from lamina (a.u.)"
		) + ggtitle(cData[1, chrom]
		) + ylim(0, 1.25) + xlim(0, cData[, max(end, na.rm = T)/1e6])
	#save_and_plot(p, file.path(outDir,
	#	sprintf("3Dstruct.GPSeq.profile.%s", cData[1, chrom])),
	#	format = "png", width = 10, height = 6)
	return(p)
}

# RUN ==========================================================================

# Prepare GPSeq data
	cat("Reading GPSeq data...\n")

	gsData = readRDS(gRDS)
	gsData[, prob_g := log2(prob_g)]
	gsData[prob_g < 0, prob_g := 0]
	gsData[prob_g > 1, prob_g := 1]

	gsData[, gpseqLab := "No GPSeq"]
	gsData[(gpseq), gpseqLab := "GPSeq included"]
	gsData[, gpseqLab := factor(gpseqLab,
		levels = c("No GPSeq", "GPSeq included"))]

	gsData[, contactLab := "Only intra contacts"]
	gsData["all" == contacts, contactLab := "All contacts"]
	gsData[, contactLab := factor(contactLab,
		levels = c("Only intra contacts", "All contacts"))]

	binSize = gsData[, unique(end-start)]
	stopifnot(1 == length(binSize))
	binStep = gsData[, .(step = unique(diff(start))),
		by = c("contactLab", "gpseqLab", "label", "chrom")] [, unique(step)]
	stopifnot(1 == length(binStep))

# Plot correlation scatter
	cat("GPSeq vs Structure correlation plot...\n")

	corData = rbindlist(by(gsData[!is.na(chromID)],
		gsData[!is.na(chromID)]$gpseqLab, function(gt) {
		rbindlist(by(gt, gt$contactLab, function(ct) {
			pcor = cor(ct$prob_g, 1-ct$rmedian,
				use = "pairwise.complete.obs", method = "pearson")
			scor = cor(ct$prob_g, 1-ct$rmedian,
				use = "pairwise.complete.obs", method = "spearman")
			data.table(pearson = pcor, spearman = scor,
				contactLab = ct[1, contactLab], gpseqLab = ct[1, gpseqLab])
		}))
	}))

	p = ggplot(gsData[!is.na(chromID)], aes(x = prob_g, y = 1-rmedian)
		) + geom_point(aes(color = reorder(chrom, chromID)), alpha = .5
		) + geom_smooth(method = "lm", linetype = "dashed", color = "red", fill = NA
		) + geom_text(data = corData, aes(label = sprintf(
				"Pearson: %.3f\nSpearman: %.3f\nR2: %.3f",
				pearson, spearman, pearson**2
			)), x = 0, y = 1, hjust = 0, vjust = 1
		) + facet_grid(~gpseqLab~contactLab
		) + scale_color_viridis(discrete = T
		) + guides(color = guide_legend("Chromosome", nrow = 3)
		) + theme(legend.position = "top",
			strip.background = element_rect(fill = NA),
			strip.text = element_text(face = "bold")
		) + xlab("log2(GPSeq score)") + ylab("Distance from lamina from structures"
		) + xlim(0, 1) + ylim(0, 1) + coord_fixed()
	save_and_plot(p, file.path(outDir, "3Dstruct.GPSeq.cor"),
		format = "png", height = 10, width = 10)

# Plot chromosome profiles
	cat("GPSeq vs Structure chromosome-profiles...\n")

	gsData2 = rbindlist(by(gsData, gsData$chrom, function(ct) {
		rbindlist(by(ct, paste0(ct$gpseq, "~", ct$label, "~", ct$contact),
			function(dt) {
				dt = dt[order(start), .(chrom, start, end, chromID,
						prob_g, rmean, rmedian,
						gpseqLab, contactLab, label)]
				starts = with(dt, seq(min(start, na.rm = T), max(start, na.rm = T),
					by = binStep))
				missing = starts[!starts %in% dt$start]
				if ( 0 == length(missing) ) {
					return(dt)
				}
				dt = rbind(dt, data.table(
					chrom = dt[1, chrom],
					start = missing, end = missing + binSize,
					chromID = dt[1, chromID],
					prob_g = NA, rmean = NA, rmedian = NA,
					gpseqLab = dt[1, gpseqLab],
					contactLab = dt[1, contactLab],
					label = dt[1, label]
				))
				dt[order(start)]
			}
		))
	}))
	gsData2$chromID = unlist(lapply(gsData2$chrom, chrom2chromID))

	pList = pblapply(split(gsData2, gsData2$chromID), plot_chrom_profile, cl = threads)
	pdf(file.path(outDir, "3Dstruct.GPSeq.profile.pdf"), width = 10, height = 10)
	l = lapply(pList, print)
	graphics.off()

# Prepare iFISH data
	cat("Reading iFISH data...\n")

	fsData = readRDS(fRDS)
	fsData[, prob_g := log2(prob_g)]
	fsData[prob_g < 0, prob_g := 0]
	fsData[prob_g > 1, prob_g := 1]

	fsData = melt(fsData,
		id.vars = c("chrom", "start", "end", "chromID", "mid", "prob_g"))

	fsData[, gpseqLab := "No GPSeq"]
	fsData[grepl("G", variable), gpseqLab := "GPSeq included"]
	fsData[, gpseqLab := factor(gpseqLab,
		levels = c("No GPSeq", "GPSeq included"))]

	fsData[, contactLab := "Only intra contacts"]
	fsData[!grepl("intra", variable), contactLab := "All contacts"]
	fsData[, contactLab := factor(contactLab,
		levels = c("Only intra contacts", "All contacts"))]

# Plot correlation scatter
	cat("iFISH vs Structure correlation plot...\n")

	corData = rbindlist(by(fsData, fsData$gpseqLab, function(gt) {
		rbindlist(by(gt, gt$contactLab, function(ct) {
			pcor = cor(ct$prob_g, 1-ct$value,
				use = "pairwise.complete.obs", method = "pearson")
			scor = cor(ct$prob_g, 1-ct$value,
				use = "pairwise.complete.obs", method = "spearman")
			data.table(pearson = pcor, spearman = scor,
				contactLab = ct[1, contactLab], gpseqLab = ct[1, gpseqLab])
		}))
	}))

	p = ggplot(fsData, aes(x = prob_g, y = 1-value)
		) + geom_point(aes(color = reorder(chrom, chromID)), alpha = .5
		) + geom_smooth(method = "lm", linetype = "dashed", color = "red", fill = NA
		) + geom_text(data = corData, aes(label = sprintf(
				"Pearson: %.3f\nSpearman: %.3f\nR2: %.3f", pearson, spearman, pearson**2
			)), x = 0, y = 1, hjust = 0, vjust = 1
		) + facet_grid(~gpseqLab~contactLab
		) + scale_color_viridis(discrete = T
		) + guides(color = guide_legend("Chromosome", nrow = 3)
		) + theme(legend.position = "top",
			strip.background = element_rect(fill = NA),
			strip.text = element_text(face = "bold")
		) + xlab("Median normalized distance from lamina (iFISH)"
		) + ylab("Distance from lamina from structures"
		) + xlim(0, 1) + ylim(0, 1) + coord_fixed()
	save_and_plot(p, file.path(outDir, "3Dstruct.iFISH.cor"),
		format = "png", height = 10, width = 10)

# END --------------------------------------------------------------------------

################################################################################
