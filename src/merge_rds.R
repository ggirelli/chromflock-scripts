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

script_name = 'merge_rds.R'
parser = arg_parser('Merge multiple RDS files', name = script_name)

parser = add_argument(parser, 'output', 'Path to output RDS file.')

parser = add_argument(parser, '--rds',
	'REQUIRED. Path to structure RDS file(s).', nargs = Inf, type = class(""))
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

p = parse_args(parser)
attach(p['' != names(p)])

stopifnot(!is.na(rds))

cat(sprintf("
 # %s

     RDS : %s
  Output : %s
 Threads : %d

\n",
	script_name, paste(rds, collapse = "\n          "), output, threads))

# RUN ==========================================================================

saveRDS(rbindlist(pblapply(rds, readRDS, cl = threads)), output)

# END --------------------------------------------------------------------------

################################################################################
