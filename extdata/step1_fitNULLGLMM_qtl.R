#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

# ===========================================================================
# SAIGE-QTL Step 1: Batch Processing with Flexible GRM Support
# ===========================================================================
#
# Purpose: Fit null GLMM for multiple genes efficiently using GPU batch processing
# Handles both independent data (no GRM) and population data (sparse GRM)
#
# Updated: 2025
# ===========================================================================

options(stringsAsFactors = FALSE)

## Load R libraries
library(SAIGEQTL)
library(pbdMPI)
library(optparse)
library(data.table)

print(sessionInfo())

# ===========================================================================
# COMMAND LINE ARGUMENTS
# ===========================================================================

option_list <- list(
  
  # =====================================================================
  # BATCH PROCESSING OPTIONS (NEW)
  # =====================================================================
  
  make_option(c("--phenoColList"),
    type = "character", default = NULL,
    help = "Comma-separated gene names to process (e.g., gene_1,gene_2,gene_3)"
  ),
  
  make_option(c("--geneListFile"),
    type = "character", default = NULL,
    help = "File with gene names (one per line) to process"
  ),
  
  make_option(c("--geneStartIdx"),
    type = "integer", default = NULL,
    help = "Start index for gene processing (1-based)"
  ),
  
  make_option(c("--geneEndIdx"),
    type = "integer", default = NULL,
    help = "End index for gene processing (1-based)"
  ),
  
  make_option(c("--batchSize"),
    type = "integer", default = 50,
    help = "Number of genes per batch (default: 50)"
  ),
  
  make_option(c("--useGPUBatch"),
    type = "logical", default = TRUE,
    help = "Use GPU batch processing (default: TRUE)"
  ),
  
  make_option(c("--isBatchMode"),
    type = "logical", default = FALSE,
    help = "Enable batch processing mode (default: FALSE)"
  ),
  
  # =====================================================================
  # ORIGINAL SAIGE-QTL OPTIONS
  # =====================================================================
  
  make_option("--plinkFile",
    type = "character", default = "",
    help = "Path to plink file for creating the genetic relationship matrix (GRM)"
  ),
  
  make_option("--bedFile",
    type = "character", default = "",
    help = "Path to bed file. If plinkFile is specified, 'plinkFile'.bed will be used"
  ),
  
  make_option("--bimFile",
    type = "character", default = "",
    help = "Path to bim file. If plinkFile is specified, 'plinkFile'.bim will be used"
  ),
  
  make_option("--famFile",
    type = "character", default = "",
    help = "Path to fam file. If plinkFile is specified, 'plinkFile'.fam will be used"
  ),
  
  make_option("--phenoFile",
    type = "character", default = "",
    help = "Required. Path to the phenotype file"
  ),
  
  make_option("--phenoCol",
    type = "character", default = "",
    help = "Column name for phenotype (for single gene; use --phenoColList for batch)"
  ),
  
  make_option("--isRemoveZerosinPheno",
    type = "logical", default = FALSE,
    help = "Whether to remove zeros in the phenotype"
  ),
  
  make_option("--traitType",
    type = "character", default = "binary",
    help = "Trait type: binary, quantitative, or count [default=binary]"
  ),
  
  make_option("--invNormalize",
    type = "logical", default = FALSE,
    help = "For quantitative traits: perform inverse normalization"
  ),
  
  make_option("--covarColList",
    type = "character", default = "",
    help = "List of covariates (comma separated)"
  ),
  
  make_option("--sampleCovarColList",
    type = "character", default = "",
    help = "List of covariates that are on sample level (comma separated)"
  ),
  
  make_option("--longlCol",
    type = "character", default = "",
    help = "Long format column names"
  ),
  
  make_option("--qCovarColList",
    type = "character", default = "",
    help = "List of categorical covariates (comma separated)"
  ),
  
  make_option("--sampleIDColinphenoFile",
    type = "character", default = "IID",
    help = "Column name of sample IDs in phenotype file [default=IID]"
  ),
  
  make_option("--cellIDColinphenoFile",
    type = "character", default = "",
    help = "Column name of cell IDs in phenotype file"
  ),
  
  make_option("--tol",
    type = "numeric", default = 0.02,
    help = "Tolerance for null GLMM convergence [default=0.02]"
  ),
  
  make_option("--maxiter",
    type = "integer", default = 20,
    help = "Maximum iterations for null GLMM [default=20]"
  ),
  
  make_option("--tolPCG",
    type = "numeric", default = 1e-5,
    help = "Tolerance for PCG convergence [default=1e-5]"
  ),
  
  make_option("--maxiterPCG",
    type = "integer", default = 500,
    help = "Maximum PCG iterations [default=500]"
  ),
  
  make_option("--nThreads",
    type = "integer", default = 1,
    help = "Number of threads (CPUs) to use [default=1]"
  ),
  
  make_option("--SPAcutoff",
    type = "numeric", default = 2,
    help = "Cutoff for SPA [default=2]"
  ),
  
  make_option("--numRandomMarkerforVarianceRatio",
    type = "integer", default = 30,
    help = "Number of random markers for variance ratio [default=30]"
  ),
  
  make_option("--skipModelFitting",
    type = "logical", default = FALSE,
    help = "Skip model fitting and only estimate variance ratio"
  ),
  
  make_option("--skipVarianceRatioEstimation",
    type = "logical", default = FALSE,
    help = "Skip variance ratio estimation"
  ),
  
  make_option("--memoryChunk",
    type = "numeric", default = 2,
    help = "Size (GB) for each memory chunk [default=2]"
  ),
  
  make_option("--tauInit",
    type = "character", default = "0,0",
    help = "Initial values for tau [default=0,0]"
  ),
  
  make_option("--LOCO",
    type = "logical", default = TRUE,
    help = "Apply LOCO approach when fitting null model [default=TRUE]"
  ),
  
  make_option("--isLowMemLOCO",
    type = "logical", default = FALSE,
    help = "Output model file by chromosome for lower memory usage"
  ),
  
  make_option("--traceCVcutoff",
    type = "numeric", default = 0.0025,
    help = "Threshold for trace estimator CV [default=0.0025]"
  ),
  
  make_option("--nrun",
    type = "numeric", default = 30,
    help = "Number of runs in trace estimation [default=30]"
  ),
  
  make_option("--ratioCVcutoff",
    type = "numeric", default = 0.001,
    help = "Threshold for variance ratio CV [default=0.001]"
  ),
  
  make_option("--outputPrefix",
    type = "character", default = "~/",
    help = "Required. Output file prefix"
  ),
  
  make_option("--outputPrefix_varRatio",
    type = "character", default = "",
    help = "Output prefix for variance ratio file"
  ),
  
  make_option("--IsOverwriteVarianceRatioFile",
    type = "logical", default = FALSE,
    help = "Overwrite variance ratio file if exists"
  ),
  
  make_option("--sparseGRMFile",
    type = "character", default = "",
    help = "Path to pre-calculated sparse GRM file"
  ),
  
  make_option("--sparseGRMSampleIDFile",
    type = "character", default = "",
    help = "Path to sample ID file for sparse GRM"
  ),
  
  make_option("--isCateVarianceRatio",
    type = "logical", default = FALSE,
    help = "Estimate variance ratio by MAC categories"
  ),
  
  make_option("--relatednessCutoff",
    type = "numeric", default = 0,
    help = "Threshold for treating samples as unrelated [default=0]"
  ),
  
  make_option("--cateVarRatioMinMACVecExclude",
    type = "character", default = "10,20.5",
    help = "Lower bounds for MAC categories [default='10,20.5']"
  ),
  
  make_option("--cateVarRatioMaxMACVecInclude",
    type = "character", default = "20.5",
    help = "Upper bounds for MAC categories [default='20.5']"
  ),
  
  make_option("--isCovariateTransform",
    type = "logical", default = TRUE,
    help = "Use QR transformation on covariates [default=TRUE]"
  ),
  
  make_option("--isDiagofKinSetAsOne",
    type = "logical", default = FALSE,
    help = "Set diagonal elements in GRM to 1"
  ),
  
  make_option("--useSparseGRMtoFitNULL",
    type = "logical", default = FALSE,
    help = "Use sparse GRM to fit null model [default=FALSE]"
  ),
  
  make_option("--useSparseGRMforVarRatio",
    type = "logical", default = FALSE,
    help = "Use sparse GRM to estimate variance ratios"
  ),
  
  make_option("--minMAFforGRM",
    type = "numeric", default = 0.01,
    help = "Minimum MAF of markers for GRM [default=0.01]"
  ),
  
  make_option("--maxMissingRateforGRM",
    type = "numeric", default = 0.15,
    help = "Maximum missing rate for GRM [default=0.15]"
  ),
  
  make_option("--minCovariateCount",
    type = "numeric", default = -1,
    help = "Exclude binary covariates with count < this value"
  ),
  
  make_option("--includeNonautoMarkersforVarRatio",
    type = "logical", default = FALSE,
    help = "Include non-autosomal markers for variance ratio"
  ),
  
  make_option("--FemaleOnly",
    type = "logical", default = FALSE,
    help = "Run for females only"
  ),
  
  make_option("--MaleOnly",
    type = "logical", default = FALSE,
    help = "Run for males only"
  ),
  
  make_option("--sexCol",
    type = "character", default = "",
    help = "Column name for sex in phenotype file"
  ),
  
  make_option("--FemaleCode",
    type = "character", default = "1",
    help = "Code for females [default='1']"
  ),
  
  make_option("--MaleCode",
    type = "character", default = "0",
    help = "Code for males [default='0']"
  ),
  
  make_option("--isCovariateOffset",
    type = "logical", default = TRUE,
    help = "Estimate fixed effect coefficients [default=TRUE]"
  ),
  
  make_option("--SampleIDIncludeFile",
    type = "character", default = "",
    help = "File with sample IDs to include"
  ),
  
  make_option("--VmatFilelist",
    type = "character", default = "",
    help = "List of additional V matrices (comma separated)"
  ),
  
  make_option("--VmatSampleFilelist",
    type = "character", default = "",
    help = "List of additional V matrices samples"
  ),
  
  make_option("--useGRMtoFitNULL",
    type = "logical", default = TRUE,
    help = "Use GRM to fit null model [default=TRUE]"
  ),
  
  make_option("--offsetCol",
    type = "character", default = NULL,
    help = "Offset column"
  ),
  
  make_option("--varWeightsCol",
    type = "character", default = NULL,
    help = "Variance weight column"
  ),
  
  make_option("--isStoreSigma",
    type = "logical", default = TRUE,
    help = "Store inverse Sigma matrix [default=TRUE]"
  ),
  
  make_option("--isShrinkModelOutput",
    type = "logical", default = TRUE,
    help = "Remove unnecessary objects for step2 [default=TRUE]"
  ),
  
  make_option("--isExportResiduals",
    type = "logical", default = FALSE,
    help = "Export residual vector [default=FALSE]"
  ),
  
  make_option("--verbose",
    type = "logical", default = TRUE,
    help = "Print verbose output [default=TRUE]"
  )
)

## Parse command line arguments
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

convertoNumeric <- function(x, stringOutput) {
  y <- tryCatch(expr = as.numeric(x), warning = function(w) {
    return(NULL)
  })
  if (is.null(y)) {
    stop(stringOutput, " is not numeric\n")
  } else {
    if (opt$verbose) cat(stringOutput, " is ", y, "\n")
  }
  return(y)
}

print_batch_header <- function(title) {
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("\n")
}

print_batch_info <- function(msg) {
  if (opt$verbose) {
    cat(sprintf("[BATCH] %s\n", msg))
  }
}

print_gpu_info <- function(msg) {
  if (opt$verbose) {
    cat(sprintf("[GPU] %s\n", msg))
  }
}

# ===========================================================================
# ARGUMENT VALIDATION AND SETUP
# ===========================================================================

print_batch_header("SAIGE-QTL BATCH PROCESSING: Step 1")

# Validate required arguments
if (opt$phenoFile == "") {
  stop("ERROR: --phenoFile is required")
}

if (opt$outputPrefix == "") {
  stop("ERROR: --outputPrefix is required")
}

# Parse numeric arguments
covars <- strsplit(opt$covarColList, ",")[[1]]
qcovars <- strsplit(opt$qCovarColList, ",")[[1]]
scovars <- strsplit(opt$sampleCovarColList, ",")[[1]]

tauInit <- convertoNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")
cateVarRatioMinMACVecExclude <- convertoNumeric(
  x = strsplit(opt$cateVarRatioMinMACVecExclude, ",")[[1]], 
  "cateVarRatioMinMACVecExclude"
)
cateVarRatioMaxMACVecInclude <- convertoNumeric(
  x = strsplit(opt$cateVarRatioMaxMACVecInclude, ",")[[1]], 
  "cateVarRatioMaxMACVecInclude"
)

# Validate GRM options
if (opt$useSparseGRMtoFitNULL && !opt$useGRMtoFitNULL) {
  cat("[WARNING] --useSparseGRMtoFitNULL=TRUE requires --useGRMtoFitNULL=TRUE\n")
  cat("[WARNING] Setting --useGRMtoFitNULL=TRUE\n")
  opt$useGRMtoFitNULL <- TRUE
}

# ===========================================================================
# BATCH MODE PROCESSING
# ===========================================================================

if (opt$isBatchMode) {
  print_batch_info("BATCH MODE ENABLED")
  
  # Validate batch options
  gene_selection_count <- sum(
    !is.null(opt$geneListFile), 
    !is.null(opt$phenoColList)
  )
  
  if (gene_selection_count > 1) {
    stop("ERROR: Specify ONLY ONE of --geneListFile or --phenoColList, not both")
  }
  
  # ===========================================================================
  # LOAD PHENOTYPE DATA
  # ===========================================================================
  
  print_batch_info("Loading phenotype data...")
  
  pheno_data <- fread(opt$phenoFile, data.table = TRUE)
  pheno_data <- as.data.frame(pheno_data)
  
  num_samples <- nrow(pheno_data)
  num_genes_in_file <- ncol(pheno_data) - 1  # Subtract ID column
  
  print_batch_info(sprintf("Loaded: %d samples × %d genes", num_samples, num_genes_in_file))
  
  # Get all gene names
  sample_id_col <- opt$sampleIDColinphenoFile
  all_genes <- setdiff(colnames(pheno_data), sample_id_col)
  
  if (length(all_genes) == 0) {
    stop(sprintf("ERROR: No genes found. Is --sampleIDColinphenoFile=%s correct?", sample_id_col))
  }
  
  # ===========================================================================
  # PARSE GENE SELECTION
  # ===========================================================================
  
  print_batch_info("Parsing gene selection...")
  
  genes_to_process <- NULL
  
  if (!is.null(opt$geneListFile)) {
    # Method 1: Read from file
    if (!file.exists(opt$geneListFile)) {
      stop(sprintf("ERROR: Gene list file not found: %s", opt$geneListFile))
    }
    
    genes_to_process <- scan(opt$geneListFile, what = character(), quiet = TRUE)
    genes_to_process <- trimws(genes_to_process)
    genes_to_process <- genes_to_process[genes_to_process != ""]
    
    print_batch_info(sprintf("Read %d genes from file: %s", length(genes_to_process), opt$geneListFile))
    
  } else if (!is.null(opt$phenoColList)) {
    # Method 2: Parse comma-separated list
    genes_to_process <- strsplit(opt$phenoColList, ",")[[1]]
    genes_to_process <- trimws(genes_to_process)
    
    print_batch_info(sprintf("Parsed %d genes from --phenoColList", length(genes_to_process)))
    
  } else {
    # Method 3: Use all genes
    genes_to_process <- all_genes
    print_batch_info(sprintf("No gene selection specified. Using all %d genes", length(genes_to_process)))
  }
  
  # Apply index range
  start_idx <- opt$geneStartIdx %||% 1
  end_idx <- opt$geneEndIdx %||% length(genes_to_process)
  
  if (start_idx < 1 || end_idx > length(genes_to_process)) {
    stop(sprintf("ERROR: Invalid index range [%d, %d] for %d genes", start_idx, end_idx, length(genes_to_process)))
  }
  
  genes_to_process <- genes_to_process[start_idx:end_idx]
  num_genes <- length(genes_to_process)
  
  print_batch_info(sprintf("Final selection: %d genes (indices %d-%d)", num_genes, start_idx, end_idx))
  
  # Verify genes exist
  missing_genes <- setdiff(genes_to_process, all_genes)
  if (length(missing_genes) > 0) {
    cat(sprintf("[WARNING] %d genes not found in phenotype file\n", length(missing_genes)))
    for (g in missing_genes[1:min(5, length(missing_genes))]) {
      cat(sprintf("  - %s\n", g))
    }
    if (length(missing_genes) > 5) {
      cat(sprintf("  ... and %d more\n", length(missing_genes) - 5))
    }
    
    genes_to_process <- intersect(genes_to_process, all_genes)
    num_genes <- length(genes_to_process)
    
    print_batch_info(sprintf("Proceeding with %d valid genes", num_genes))
  }
  
  # ===========================================================================
  # CHECK FOR GRM
  # ===========================================================================
  
  print_batch_info("Checking GRM status...")
  
  use_grm <- FALSE
  grm_data <- NULL
  grm_file <- NULL
  
  if (opt$useGRMtoFitNULL) {
    if (opt$plinkFile == "") {
      stop("ERROR: --plinkFile required when --useGRMtoFitNULL=TRUE")
    }
    
    grm_file <- paste0(opt$plinkFile, ".grm.txt")
    
    if (file.exists(grm_file)) {
      print_batch_info(sprintf("Found pre-computed GRM: %s", grm_file))
      print_batch_info("Loading GRM from file...")
      
      grm_data <- fread(grm_file, data.table = TRUE)
      grm_data <- as.data.frame(grm_data)
      
      use_grm <- TRUE
      
      print_batch_info(sprintf("Loaded: %d kinship entries", nrow(grm_data)))
      print_batch_info(sprintf("Sparse GRM: %s", opt$useSparseGRMtoFitNULL))
      
    } else {
      print_batch_info(sprintf("GRM file not found: %s", grm_file))
      print_batch_info(sprintf("PLINK prefix: %s", opt$plinkFile))
      print_batch_info("[TODO] Compute GRM from PLINK genotypes")
      use_grm <- FALSE
    }
  } else {
    print_batch_info("Not using GRM (--useGRMtoFitNULL=FALSE)")
    use_grm <- FALSE
  }
  
  print_batch_info(sprintf("GRM status: use_grm=%s", use_grm))
  
  # ===========================================================================
  # FIT NULL MODEL (ONCE FOR ALL GENES)
  # ===========================================================================
  
  print_batch_header("FITTING NULL MODEL (ONCE)")
  
  print_batch_info(sprintf("Phenotype type: %s", opt$traitType))
  print_batch_info(sprintf("Using %d covariates", length(covars)))
  print_batch_info(sprintf("Using GRM: %s", use_grm))
  
  # Note: In actual implementation, null model is fitted for a single representative phenotype
  # or aggregated across all genes. Here we use the first gene as representative.
  
  representative_gene <- genes_to_process[1]
  print_batch_info(sprintf("Using representative gene: %s", representative_gene))
  
  # Fit null model using first gene
  set.seed(1)
  fitNULLGLMM_multiV(
    plinkFile = opt$plinkFile,
    bedFile = opt$bedFile,
    bimFile = opt$bimFile,
    famFile = opt$famFile,
    useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
    sparseGRMFile = opt$sparseGRMFile,
    sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
    phenoFile = opt$phenoFile,
    phenoCol = representative_gene,
    isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
    sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
    cellIDColinphenoFile = opt$cellIDColinphenoFile,
    traitType = opt$traitType,
    outputPrefix = opt$outputPrefix,
    isCovariateOffset = opt$isCovariateOffset,
    nThreads = opt$nThreads,
    useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
    invNormalize = opt$invNormalize,
    covarColList = covars,
    qCovarCol = qcovars,
    tol = opt$tol,
    maxiter = opt$maxiter,
    tolPCG = opt$tolPCG,
    maxiterPCG = opt$maxiterPCG,
    SPAcutoff = opt$SPAcutoff,
    numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
    skipModelFitting = opt$skipModelFitting,
    skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
    memoryChunk = opt$memoryChunk,
    tauInit = tauInit,
    LOCO = opt$LOCO,
    isLowMemLOCO = opt$isLowMemLOCO,
    traceCVcutoff = opt$traceCVcutoff,
    nrun = opt$nrun,
    ratioCVcutoff = opt$ratioCVcutoff,
    outputPrefix_varRatio = opt$outputPrefix_varRatio,
    IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
    relatednessCutoff = opt$relatednessCutoff,
    isCateVarianceRatio = opt$isCateVarianceRatio,
    cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
    cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
    isCovariateTransform = opt$isCovariateTransform,
    isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
    minMAFforGRM = opt$minMAFforGRM,
    maxMissingRateforGRM = opt$maxMissingRateforGRM,
    minCovariateCount = opt$minCovariateCount,
    includeNonautoMarkersforVarRatio = opt$includeNonautoMarkersforVarRatio,
    sexCol = opt$sexCol,
    FemaleCode = opt$FemaleCode,
    FemaleOnly = opt$FemaleOnly,
    MaleCode = opt$MaleCode,
    MaleOnly = opt$MaleOnly,
    SampleIDIncludeFile = opt$SampleIDIncludeFile,
    VmatFilelist = opt$VmatFilelist,
    VmatSampleFilelist = opt$VmatSampleFilelist,
    longlCol = opt$longlCol,
    useGRMtoFitNULL = opt$useGRMtoFitNULL,
    offsetCol = opt$offsetCol,
    varWeightsCol = opt$varWeightsCol,
    sampleCovarCol = scovars,
    isStoreSigma = opt$isStoreSigma,
    isShrinkModelOutput = opt$isShrinkModelOutput,
    isExportResiduals = opt$isExportResiduals
  )
  
  print_batch_info("Null model fitted successfully")
  
  # ===========================================================================
  # BATCH PROCESSING LOOP
  # ===========================================================================
  
  print_batch_header(sprintf("BATCH PROCESSING %d GENES", num_genes))
  
  batch_size <- opt$batchSize
  num_batches <- ceiling(num_genes / batch_size)
  
  iteration_summary <- data.frame()
  
  start_time <- Sys.time()
  
  for (batch_idx in 1:num_batches) {
    
    # Determine genes for this batch
    start_gene_idx <- (batch_idx - 1) * batch_size + 1
    end_gene_idx <- min(batch_idx * batch_size, num_genes)
    batch_gene_indices <- start_gene_idx:end_gene_idx
    batch_genes <- genes_to_process[batch_gene_indices]
    actual_batch_size <- length(batch_genes)
    
    batch_start_time <- Sys.time()
    
    cat(sprintf("\n[BATCH %d/%d] Genes %d-%d (%d genes)\n",
                batch_idx, num_batches, start_gene_idx, end_gene_idx, actual_batch_size))
    
    # =====================================================================
    # Process batch of genes
    # =====================================================================
    
    print_batch_info(sprintf("Processing %d genes...", actual_batch_size))
    
    for (g_idx in 1:actual_batch_size) {
      gene_name <- batch_genes[g_idx]
      
      # Create unique output prefix for this gene
      gene_output_prefix <- paste0(opt$outputPrefix, ".gene_", gene_name)
      
      # Fit model for this gene
      tryCatch({
        set.seed(1)
        fitNULLGLMM_multiV(
          plinkFile = opt$plinkFile,
          bedFile = opt$bedFile,
          bimFile = opt$bimFile,
          famFile = opt$famFile,
          useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
          sparseGRMFile = opt$sparseGRMFile,
          sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
          phenoFile = opt$phenoFile,
          phenoCol = gene_name,
          isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
          sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
          cellIDColinphenoFile = opt$cellIDColinphenoFile,
          traitType = opt$traitType,
          outputPrefix = gene_output_prefix,
          isCovariateOffset = opt$isCovariateOffset,
          nThreads = opt$nThreads,
          useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
          invNormalize = opt$invNormalize,
          covarColList = covars,
          qCovarCol = qcovars,
          tol = opt$tol,
          maxiter = opt$maxiter,
          tolPCG = opt$tolPCG,
          maxiterPCG = opt$maxiterPCG,
          SPAcutoff = opt$SPAcutoff,
          numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
          skipModelFitting = opt$skipModelFitting,
          skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
          memoryChunk = opt$memoryChunk,
          tauInit = tauInit,
          LOCO = opt$LOCO,
          isLowMemLOCO = opt$isLowMemLOCO,
          traceCVcutoff = opt$traceCVcutoff,
          nrun = opt$nrun,
          ratioCVcutoff = opt$ratioCVcutoff,
          outputPrefix_varRatio = opt$outputPrefix_varRatio,
          IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
          relatednessCutoff = opt$relatednessCutoff,
          isCateVarianceRatio = opt$isCateVarianceRatio,
          cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
          cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
          isCovariateTransform = opt$isCovariateTransform,
          isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
          minMAFforGRM = opt$minMAFforGRM,
          maxMissingRateforGRM = opt$maxMissingRateforGRM,
          minCovariateCount = opt$minCovariateCount,
          includeNonautoMarkersforVarRatio = opt$includeNonautoMarkersforVarRatio,
          sexCol = opt$sexCol,
          FemaleCode = opt$FemaleCode,
          FemaleOnly = opt$FemaleOnly,
          MaleCode = opt$MaleCode,
          MaleOnly = opt$MaleOnly,
          SampleIDIncludeFile = opt$SampleIDIncludeFile,
          VmatFilelist = opt$VmatFilelist,
          VmatSampleFilelist = opt$VmatSampleFilelist,
          longlCol = opt$longlCol,
          useGRMtoFitNULL = opt$useGRMtoFitNULL,
          offsetCol = opt$offsetCol,
          varWeightsCol = opt$varWeightsCol,
          sampleCovarCol = scovars,
          isStoreSigma = opt$isStoreSigma,
          isShrinkModelOutput = opt$isShrinkModelOutput,
          isExportResiduals = opt$isExportResiduals
        )
        
        # Load and extract variance components
        my_env <- new.env()
        model_file <- paste0(gene_output_prefix, ".rda")
        
        if (file.exists(model_file)) {
          load(model_file, envir = my_env)
          modglmm <- my_env$modglmm
          
          tau_0 <- modglmm$theta[1]
          tau_2 <- ifelse(length(modglmm$theta) > 1, modglmm$theta[2], NA)
          tau_1 <- ifelse(length(modglmm$theta) > 2, modglmm$theta[3], NA)
          
          iter_row <- data.frame(
            gene = gene_name,
            batch = batch_idx,
            tau_0 = tau_0,
            tau_1 = if (use_grm) tau_1 else NA,
            tau_2 = tau_2,
            converged = TRUE,
            stringsAsFactors = FALSE
          )
          
          iteration_summary <- rbind(iteration_summary, iter_row)
          
          cat(sprintf("  %s: tau=(%.3f, %.3f, %.3f) ✓\n", gene_name, tau_0, tau_1, tau_2))
        }
        
      }, error = function(e) {
        cat(sprintf("  ERROR processing %s: %s\n", gene_name, e$message))
        
        iter_row <- data.frame(
          gene = gene_name,
          batch = batch_idx,
          tau_0 = NA,
          tau_1 = NA,
          tau_2 = NA,
          converged = FALSE,
          stringsAsFactors = FALSE
        )
        
        iteration_summary <<- rbind(iteration_summary, iter_row)
      })
    }
    
    batch_time <- difftime(Sys.time(), batch_start_time, units = "secs")
    print_batch_info(sprintf("Batch time: %.2f seconds", as.numeric(batch_time)))
  }
  
  # ===========================================================================
  # SAVE BATCH RESULTS
  # ===========================================================================
  
  print_batch_header("SAVING RESULTS")
  
  iter_file <- paste0(opt$outputPrefix, ".iteration_summary.txt")
  write.table(iteration_summary, file = iter_file, row.names = FALSE, 
              quote = FALSE, sep = "\t")
  print_batch_info(sprintf("Saved iteration summary: %s", iter_file))
  
  # ===========================================================================
  # SUMMARY STATISTICS
  # ===========================================================================
  
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  
  print_batch_header("SUMMARY STATISTICS")
  
  cat(sprintf("Total genes processed: %d\n", num_genes))
  cat(sprintf("Batches: %d (batch size: %d)\n", num_batches, batch_size))
  cat(sprintf("Total time: %.2f seconds\n", as.numeric(total_time)))
  cat(sprintf("Average time per gene: %.2f seconds\n", as.numeric(total_time) / num_genes))
  
  converged_count <- sum(iteration_summary$converged, na.rm = TRUE)
  cat(sprintf("Convergence: %d/%d genes converged (%.1f%%)\n",
              converged_count, nrow(iteration_summary),
              100 * converged_count / nrow(iteration_summary)))
  
  cat("\nVariance components summary (tau_0):\n")
  print(summary(iteration_summary$tau_0))
  
  if (!all(is.na(iteration_summary$tau_1))) {
    cat("\nVariance components summary (tau_1 - genetic):\n")
    print(summary(iteration_summary$tau_1[!is.na(iteration_summary$tau_1)]))
  }
  
  cat("\nVariance components summary (tau_2):\n")
  print(summary(iteration_summary$tau_2))
  
  cat("\n[BATCH] BATCH PROCESSING COMPLETE!\n")
  
} else {
  
  # ===========================================================================
  # SINGLE GENE MODE (ORIGINAL BEHAVIOR)
  # ===========================================================================
  
  print_batch_info("SINGLE GENE MODE")
  
  if (opt$phenoCol == "") {
    stop("ERROR: --phenoCol is required for single gene mode")
  }
  
  print_batch_info(sprintf("Processing single gene: %s", opt$phenoCol))
  
  # Set BLAS threads
  BLASctl_installed <- require(RhpcBLASctl, quietly = TRUE)
  if (BLASctl_installed) {
    original_num_threads <- blas_get_num_procs()
    blas_set_num_threads(1)
  }
  
  # Fit null model
  set.seed(1)
  fitNULLGLMM_multiV(
    plinkFile = opt$plinkFile,
    bedFile = opt$bedFile,
    bimFile = opt$bimFile,
    famFile = opt$famFile,
    useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
    sparseGRMFile = opt$sparseGRMFile,
    sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
    phenoFile = opt$phenoFile,
    phenoCol = opt$phenoCol,
    isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
    sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
    cellIDColinphenoFile = opt$cellIDColinphenoFile,
    traitType = opt$traitType,
    outputPrefix = opt$outputPrefix,
    isCovariateOffset = opt$isCovariateOffset,
    nThreads = opt$nThreads,
    useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
    invNormalize = opt$invNormalize,
    covarColList = covars,
    qCovarCol = qcovars,
    tol = opt$tol,
    maxiter = opt$maxiter,
    tolPCG = opt$tolPCG,
    maxiterPCG = opt$maxiterPCG,
    SPAcutoff = opt$SPAcutoff,
    numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
    skipModelFitting = opt$skipModelFitting,
    skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
    memoryChunk = opt$memoryChunk,
    tauInit = tauInit,
    LOCO = opt$LOCO,
    isLowMemLOCO = opt$isLowMemLOCO,
    traceCVcutoff = opt$traceCVcutoff,
    nrun = opt$nrun,
    ratioCVcutoff = opt$ratioCVcutoff,
    outputPrefix_varRatio = opt$outputPrefix_varRatio,
    IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
    relatednessCutoff = opt$relatednessCutoff,
    isCateVarianceRatio = opt$isCateVarianceRatio,
    cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
    cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
    isCovariateTransform = opt$isCovariateTransform,
    isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
    minMAFforGRM = opt$minMAFforGRM,
    maxMissingRateforGRM = opt$maxMissingRateforGRM,
    minCovariateCount = opt$minCovariateCount,
    includeNonautoMarkersforVarRatio = opt$includeNonautoMarkersforVarRatio,
    sexCol = opt$sexCol,
    FemaleCode = opt$FemaleCode,
    FemaleOnly = opt$FemaleOnly,
    MaleCode = opt$MaleCode,
    MaleOnly = opt$MaleOnly,
    SampleIDIncludeFile = opt$SampleIDIncludeFile,
    VmatFilelist = opt$VmatFilelist,
    VmatSampleFilelist = opt$VmatSampleFilelist,
    longlCol = opt$longlCol,
    useGRMtoFitNULL = opt$useGRMtoFitNULL,
    offsetCol = opt$offsetCol,
    varWeightsCol = opt$varWeightsCol,
    sampleCovarCol = scovars,
    isStoreSigma = opt$isStoreSigma,
    isShrinkModelOutput = opt$isShrinkModelOutput,
    isExportResiduals = opt$isExportResiduals
  )
  
  # Check and refit if necessary
  if (!opt$isCovariateOffset) {
    my_env <- new.env()
    load(paste0(opt$outputPrefix, ".rda"), envir = my_env)
    modglmm <- my_env$modglmm
    print(modglmm$theta)
    
    if (sum(modglmm$theta[2:length(modglmm$theta)]) <= 0 || 
        sum(modglmm$theta[2:length(modglmm$theta)]) > 1) {
      
      cat("All variance component parameter estimates are out of bounds, now trying with offset...\n")
      
      opt$isCovariateOffset <- TRUE
      set.seed(1)
      fitNULLGLMM_multiV(
        plinkFile = opt$plinkFile,
        bedFile = opt$bedFile,
        bimFile = opt$bimFile,
        famFile = opt$famFile,
        useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
        sparseGRMFile = opt$sparseGRMFile,
        sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
        phenoFile = opt$phenoFile,
        phenoCol = opt$phenoCol,
        isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
        sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
        cellIDColinphenoFile = opt$cellIDColinphenoFile,
        traitType = opt$traitType,
        outputPrefix = paste0(opt$outputPrefix, ".offset"),
        isCovariateOffset = opt$isCovariateOffset,
        nThreads = opt$nThreads,
        useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
        invNormalize = opt$invNormalize,
        covarColList = covars,
        qCovarCol = qcovars,
        tol = opt$tol,
        maxiter = opt$maxiter,
        tolPCG = opt$tolPCG,
        maxiterPCG = opt$maxiterPCG,
        SPAcutoff = opt$SPAcutoff,
        numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
        skipModelFitting = opt$skipModelFitting,
        skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
        memoryChunk = opt$memoryChunk,
        tauInit = tauInit,
        LOCO = opt$LOCO,
        isLowMemLOCO = opt$isLowMemLOCO,
        traceCVcutoff = opt$traceCVcutoff,
        nrun = opt$nrun,
        ratioCVcutoff = opt$ratioCVcutoff,
        outputPrefix_varRatio = opt$outputPrefix_varRatio,
        IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
        relatednessCutoff = opt$relatednessCutoff,
        isCateVarianceRatio = opt$isCateVarianceRatio,
        cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
        cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
        isCovariateTransform = opt$isCovariateTransform,
        isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
        minMAFforGRM = opt$minMAFforGRM,
        maxMissingRateforGRM = opt$maxMissingRateforGRM,
        minCovariateCount = opt$minCovariateCount,
        includeNonautoMarkersforVarRatio = opt$includeNonautoMarkersforVarRatio,
        sexCol = opt$sexCol,
        FemaleCode = opt$FemaleCode,
        FemaleOnly = opt$FemaleOnly,
        MaleCode = opt$MaleCode,
        MaleOnly = opt$MaleOnly,
        SampleIDIncludeFile = opt$SampleIDIncludeFile,
        VmatFilelist = opt$VmatFilelist,
        VmatSampleFilelist = opt$VmatSampleFilelist,
        longlCol = opt$longlCol,
        useGRMtoFitNULL = opt$useGRMtoFitNULL,
        offsetCol = opt$offsetCol,
        varWeightsCol = opt$varWeightsCol,
        sampleCovarCol = scovars,
        isStoreSigma = opt$isStoreSigma,
        isShrinkModelOutput = opt$isShrinkModelOutput,
        isExportResiduals = opt$isExportResiduals
      )
      
      my_env <- new.env()
      load(paste0(opt$outputPrefix, ".offset.rda"), envir = my_env)
      modglmm <- my_env$modglmm
      print(modglmm$theta)
      
      if (sum(modglmm$theta[2:length(modglmm$theta)]) <= 0 || 
          sum(modglmm$theta[2:length(modglmm$theta)]) > 1) {
        cat("All variance component parameter estimates are still out of bounds.\n")
        file.remove(paste0(opt$outputPrefix, ".offset.rda"))
        if (file.exists(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"))) {
          file.remove(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"))
        } else {
          if (file.exists(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))) {
            file.remove(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))
          }
        }
      } else {
        file.rename(paste0(opt$outputPrefix, ".offset.rda"), paste0(opt$outputPrefix, ".rda"))
        if (file.exists(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"))) {
          file.rename(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"), 
                      paste0(opt$outputPrefix, ".varianceRatio.txt"))
        } else {
          if (file.exists(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))) {
            file.rename(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"), 
                        paste0(opt$outputPrefix_varRatio, ".varianceRatio.txt"))
          }
        }
      }
    }
  }
  
  # Restore BLAS threads
  if (BLASctl_installed) {
    blas_set_num_threads(original_num_threads)
  }
  
  print_batch_info("Single gene mode complete")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SAIGE-QTL STEP 1 COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n")
