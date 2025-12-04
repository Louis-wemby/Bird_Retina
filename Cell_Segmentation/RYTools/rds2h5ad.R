library(argparse)
library(Seurat)
library(reticulate)
use_python("/usr/local/bin/python", required = TRUE)

parser <- ArgumentParser(description = "convert rds to h5ad")
parser$add_argument("--infile", type = "character", required = TRUE, help = "Path to the rds file")
parser$add_argument("--outfile", type = "character", required = TRUE, help = "Path to save the output file")
args <- parser$parse_args()

infile = args$infile
outfile = args$outfile

object <- readRDS(infile)
genes <- as.data.frame(rownames(object), row.names = rownames(object))
names(genes) <- "Gene"
cells <- as.data.frame(colnames(object), row.names = colnames(object))
names(cells) <- "CellID"
row <- object@images[["slice1"]]@coordinates$row
col <- object@images[["slice1"]]@coordinates$col
coordinates <- list(matrix(c(row, col), ncol = 2))
names(coordinates) <- "spatial"
anndata <- import("anndata", delay_load = TRUE)
ad <- anndata$AnnData(X = object@assays$Spatial@counts, obs = genes, var = cells, varm = coordinates)
ad <- ad$T
ad$write_h5ad(outfile)
