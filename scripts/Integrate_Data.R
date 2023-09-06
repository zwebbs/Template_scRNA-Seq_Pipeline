# File Name: Integrative_Analysis.R
# Created By: ZW
# Created On: 2023-07-04
# Purpose: TODO

# Read in Command-line arguments from optparse
# -----------------------------------------------------
library(optparse)
option_list <- list(
    make_option(c("-o", "--outfile"),
                help="output file path for RData object of integrated data"),
    make_option(c("-i", "--manifest"),
                help="input file manifest in the form 'sample,counts-dir'")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Library Imports
# -----------------------------------------------------
library(dplyr)
library(ggplot2)
library(Seurat)


# Read in Data Objects and hard filter
# -----------------------------------------------------
# specify file paths
fp.df <- read.csv(args$manifest, header=F)
colnames(fp.df) <- c("name", "path","scrublet")
scr <- as.vector(fp.df$scrublet)
fps <- as.vector(fp.df$path)
names(scr) <- as.vector(fp.df$name)
names(fps) <- as.vector(fp.df$name)
file.paths <- as.list(fps)
scrublet.paths <- as.list(scr)

# read in data objects, create seurat objects
cat("Reading filtered data from 10X result...\n")
dat.obj.raw <- lapply(names(file.paths), function(x) 
                       Read10X(data.dir=file.paths[[x]]) %>%
                           CreateSeuratObject(counts=.,project=x)
                )
names(dat.obj.raw) <- names(file.paths)

# read in scrublet data. append to seurat object, and subset
for(n in names(scrublet.paths)) {
    # read in scrublet datasets
    scrubs <- read.csv(scrublet.paths[[n]], header=FALSE)
    colnames(scrubs) <- c("barcode","doublet.score","is.doublet")
    rownames(scrubs) <- scrubs$barcode
    
    # append data to metadata for each dataset
    cells.to.keep <- scrubs$barcode[which(scrubs$is.doublet=="False")]
    dat.obj.raw[[n]] <- subset(dat.obj.raw[[n]], cells=cells.to.keep)
}


# create percent.mt feature
cat("Calculating Percent Mitochondrial Reads..\n")
dat.obj.prep <- list()
for (n in names(dat.obj.raw)) {
    s.obj <- dat.obj.raw[[n]]
    s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern="^MT-")
    dat.obj.prep[[n]] <- s.obj
}
rm(dat.obj.raw)

# filter on percent.mt < 50%, nCount_RNA > 1500
mt.filt <- 50
rna.filt <- 1500
cat(sprintf("Filtering cells based on percent MT: (< %d) and RNA reads (> %d)\n",mt.filt,rna.filt))
dat.obj.filt <- lapply(dat.obj.prep, function(x)
                          subset(x, subset = percent.mt < mt.filt & nCount_RNA > rna.filt)
                         )
names(dat.obj.filt) <- names(file.paths)
rm(dat.obj.prep)

# Normalize and Scale Data then perform basic clustering 
# ----------------------------------------------------
dat.obj.normd <- lapply(dat.obj.filt, function(x)
                            NormalizeData(x) %>%
                                FindVariableFeatures(
                                    ., selection.method="vst",
                                    nfeatures=2500
                                )
                       )
rm(dat.obj.filt)

# Integrate Datasets from the four experiments
# -----------------------------------------------------
features <- SelectIntegrationFeatures(object.list=dat.obj.normd)
dat.obj.preint <- lapply(dat.obj.normd, function(x)
                         ScaleData(x,features=features) %>%
                             RunPCA(.,features=features)
                        )
eht.anchors <- FindIntegrationAnchors(dat.obj.preint,
                                      anchor.features=features,
                                      reduction="rpca"
                                     )
dat.int <- IntegrateData(anchorset=eht.anchors)

# save integrated dataset 
save(dat.int,file=args$outfile)







# filter based on mitochontrial read percentage (45%)






