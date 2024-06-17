
require(dplyr)
require(rlang)

vand <- function(...) Reduce(`&`, list(...))
vany <- function(...) Reduce(`|`, list(...))
excluded_table <- function(geno_tbl) {
     select(geno_tbl$excluded, scan_id, .reason, .exclude_call) %>% rename(reason = .reason)
}

sum_excluded <- function(geno_tbl) {
    total_count <- nrow(geno_tbl$included) + nrow(geno_tbl$excluded)
    data <- select(geno_tbl$excluded, scan_id, .reason, .order_excluded, .exclude_call) %>%
        rename(reason = .reason, call = .exclude_call)

    top_row <- tibble::tibble(.order_excluded = 0,
                              reason = 'Not Excluded',
                              call = NA,
                              number_excluded = 0,
                              cumulative_excluded = 0
                              )

    table <- data %>%
        dplyr::group_by(.order_excluded) %>%
        dplyr::summarize(reason = unique(reason), call = unique(call), number_excluded = n()) %>%
        dplyr::bind_rows(top_row, .) %>%
        arrange(.order_excluded) %>%
        dplyr::mutate(
                   cumulative_excluded = cumsum(number_excluded),
                   number_included = total_count - cumulative_excluded,
               ) %>%
        select(-.order_excluded)
    table
}

geno_tibble <- function(included, excluded) {
    if (!('tbl' %in% class(included)) | !('tbl' %in% class(excluded))) {
        stop('geno_tibble constructor type error')
    }
    out <- list(included = included, excluded = excluded)
    class(out) <- 'geno_tbl'
    out
}

scan_annotation_tbl <- function(geno_tibble) {
    class(geno_tibble) <- c('scan_annotation_tbl', 'geno_tbl')
    geno_tibble
}

snp_annotation_tbl <- function(geno_tibble) {
    class(geno_tibble) <- c('snp_annotation_tbl', 'geno_tbl')
    geno_tibble
}

as_geno_tibble <- function(.data) {
    UseMethod("as_geno_tibble", .data)
}

as_scan_annotation <- function(.data) {
    UseMethod("as_scan_annotation", .data)
}
as_snp_annotation <- function(.data) {
    UseMethod("as_snp_annotation", .data)
}

all_scans <- function(.data) {
    dplyr:::bind_rows(.data$included, .data$excluded) %>%
        arrange(.row_number) %>%
        select(-.row_number)
}

as_geno_tibble.tbl <- function(.data) {
    included <- .data %>%
        dplyr:::mutate(
                    .row_number = dplyr::row_number()
                )
    excluded <- dplyr:::slice.tbl_df(.data, 0) %>%
        dplyr:::mutate.tbl_df( .reason = character(), .order_excluded = numeric(), .row_number = numeric(), .exclude_call = character())
    geno_tibble(included, excluded)
}

as_scan_annotation.tbl <- function(.data) {
    if (!('scan_id' %in% colnames(.data))) {
        .data <- mutate(.data, scan_id = dplyr::row_number())
    }
    .data <- as_geno_tibble(.data)
    class(.data) <- c('snp_annotation_tbl', 'geno_tbl')
    .data
}

as_snp_annotation.tbl <- function(.data) {
    if (!('snp_id' %in% colnames(.data))) {
        .data <- mutate(.data, snp_id = dplyr::row_number())
    }
    .data <- as_geno_tibble(.data)
    class(.data) <- c('snp_annotation_tbl', 'geno_tbl')
    .data
}

select.geno_tbl <- function(.data, ...) {
    dots <- rlang::quos(...)
    include <- .data$include %>% select(!!!dots, .row_number)
    exclude <- .data$exclude %>% select(!!!dots, .row_number, .order_excluded, .exclude_call)
    geno_tibble(include, exclude)
}
select.scan_annotation_tbl <- function(.data, ...) {
    scan_annotation_tbl(select.geno_tbl(.data, !!!rlang::quos(...)))
}
select.snp_annotation_tbl <- function(.data, ...) {
    snp_annotation_tbl(select.geno_tbl(.data, !!!rlang::quos(...)))
}

arrange.geno_tbl <- function(.data, ...) {
    .data$included <- dplyr::arrange(.data$included, !!!rlang::quos(...))
    .data
}

mutate.geno_tbl <- function(.data, ..., on_excluded = F) {
    dots <- rlang::quos(...)
    .data$included <- dplyr::mutate(.data$included, !!!dots)
    if (!on_excluded) {
        dots <- purrr::map(dots, ~set_expr(., NA))
    }
    .data$excluded <- dplyr::mutate(.data$excluded, !!!dots)
    .data
}

mutate.snp_annotation_tbl <- function(.data, ..., on_excluded = F) {
    snp_annotation_tbl(mutate.geno_tbl(.data, !!!rlang::quos(...), on_excluded = on_excluded))
}
mutate.scan_annotation_tbl <- function(.data, ..., on_excluded = F) {
    scan_annotation_tbl(mutate.geno_tbl(.data, !!!rlang::quos(...),on_excluded = on_excluded))
}


exclude.geno_tbl <- function(.data, ..., reason) {
    dots <- rlang::quos(...)


    if (nrow(.data$excluded) == 0) {
        order_excluded <- 1
    } else {
        order_excluded <- max(.data$excluded$.order_excluded) + 1
    }
    exclude_call <- paste(purrr::map(rlang::quos(...), rlang::get_expr), collapse = ', ')
    to_exclude <- dplyr::filter(.data$included, !!!dots)
    included <- setdiff(.data$included, to_exclude)
    excluded <- to_exclude %>%
        dplyr:::mutate.tbl_df(.reason = reason, .order_excluded = order_excluded, .exclude_call = exclude_call ) %>%
        dplyr:::bind_rows(.data$excluded, .)
    geno_tibble(included, excluded)
}

exclude.scan_annotation_tbl <- function(.data, ..., reason) {
    scan_annotation_tbl(exclude.geno_tbl(.data, !!!rlang::quos(...), reason = reason))
}
exclude.snp_annotation_tbl <- function(.data, ..., reason) {
    snp_annotation_tbl(exclude.geno_tbl(.data, !!!rlang::quos(...), reason = reason))
}

print.geno_tbl <- function(.data, ...) {
    to_print <- dplyr:::select.data.frame(.data$included, -.row_number)
    tibble:::print.tbl_df(to_print, ...)
}

print.scan_annotation_tbl <- function(.data, ...) {
    cat(paste0("# Scan Annotation\n"))
    print.geno_tbl(.data, ...)
    cat(paste0("Excluded ", nrow(.data$excluded), " scans total.\n"))
}

print.snp_annotation_tbl <- function(.data, ...) {
    cat(paste0("# SNP Annotation\n"))
    print.geno_tbl(.data, ...)
    cat(paste0("Excluded ", nrow(.data$excluded), " snps total.\n"))
}




### R code from vignette source 'DataCleaning.Rnw'

###################################################
### code chunk number 1: DataCleaning.Rnw:124-153
###################################################
library(GWASTools)
library(GWASdata)
# Load the SNP annotation (simple data frame)
data(illumina_snp_annot)
# Create a SnpAnnotationDataFrame
snpAnnot <- SnpAnnotationDataFrame(illumina_snp_annot)
# names of columns
varLabels(snpAnnot)
# data
head(pData(snpAnnot))
# Add metadata to describe the columns
meta <- varMetadata(snpAnnot)
meta[c("snpID", "chromosome", "position", "rsID", "alleleA", "alleleB",
  "BeadSetID", "IntensityOnly", "tAA", "tAB", "tBB", "rAA", "rAB", "rBB"),
  "labelDescription"] <- c("unique integer ID for SNPs",
  paste("integer code for chromosome: 1:22=autosomes,",
   "23=X, 24=pseudoautosomal, 25=Y, 26=Mitochondrial, 27=Unknown"),
  "base pair position on chromosome (build 36)",
  "RS identifier",
  "alelleA", "alleleB",
  "BeadSet ID from Illumina",
  "1=no genotypes were attempted for this assay",
  "mean theta for AA cluster",
  "mean theta for AB cluster",
  "mean theta for BB cluster",
  "mean R for AA cluster",
  "mean R for AB cluster",
  "mean R for BB cluster")
varMetadata(snpAnnot) <- meta


###################################################
### code chunk number 2: DataCleaning.Rnw:159-168
###################################################
snpID <- snpAnnot$snpID
snpID <- getSnpID(snpAnnot)
chrom <- snpAnnot[["chromosome"]]
chrom <- getChromosome(snpAnnot)
table(chrom)
chrom <- getChromosome(snpAnnot, char=TRUE)
table(chrom)
position <- getPosition(snpAnnot)
rsID <- getVariable(snpAnnot, "rsID")


###################################################
### code chunk number 3: DataCleaning.Rnw:191-201
###################################################
tmp <- snpAnnot[,c("snpID", "chromosome", "position")]
snp <- getAnnotation(tmp)
snp$flag <- sample(c(TRUE, FALSE), nrow(snp), replace=TRUE)
pData(tmp) <- snp
meta <- getMetadata(tmp)
meta["flag", "labelDescription"] <- "flag"
varMetadata(tmp) <- meta
getVariableNames(tmp)
varLabels(tmp)[4] <- "FLAG"
rm(tmp)


###################################################
### code chunk number 4: DataCleaning.Rnw:229-256
###################################################
# Load the scan annotation (simple data frame)
data(illumina_scan_annot)
# Create a ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(illumina_scan_annot)
# names of columns
varLabels(scanAnnot)
# data
head(pData(scanAnnot))
# Add metadata to describe the columns
meta <- varMetadata(scanAnnot)
meta[c("scanID", "subjectID", "family", "father", "mother",
  "CoriellID", "race", "sex", "status", "genoRunID", "plate",
  "batch", "file"), "labelDescription"] <-
   c("unique ID for scans",
  "subject identifier (may have multiple scans)",
  "family identifier",
  "father identifier as subjectID",
  "mother identifier as subjectID",
  "Coriell subject identifier",
  "HapMap population group",
  "sex coded as M=male and F=female",
  "simulated case/control status" ,
  "genotyping instance identifier",
  "plate containing samples processed together for genotyping chemistry",
  "simulated genotyping batch",
  "raw data file")
varMetadata(scanAnnot) <- meta


###################################################
### code chunk number 5: DataCleaning.Rnw:263-268
###################################################
scanID <- scanAnnot$scanID
scanID <- getScanID(scanAnnot)
sex <- scanAnnot[["sex"]]
sex <- getSex(scanAnnot)
subjectID <- getVariable(scanAnnot, "subjectID")


###################################################
### code chunk number 6: DataCleaning.Rnw:336-384
###################################################
# Define a path to the raw data files
path <- system.file("extdata", "illumina_raw_data", package="GWASdata")

geno.file <- "tmp.geno.gds"

# first 3 samples only
scan.annotation <- illumina_scan_annot[1:3, c("scanID", "genoRunID", "file")]
names(scan.annotation)[2] <- "scanName"

snp.annotation <- illumina_snp_annot[,c("snpID", "rsID", "chromosome", "position")]
# indicate which column of SNP annotation is referenced in data files
names(snp.annotation)[2] <-  "snpName"

col.nums <- as.integer(c(1,2,12,13))
names(col.nums) <- c("snp", "sample", "a1", "a2")

diag.geno.file <- "diag.geno.RData"
diag.geno <- createDataFile(path = path,
  filename = geno.file,
  file.type = "gds",
  variables = "genotype",
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation,
  sep.type = ",",
  skip.num = 11,
  col.total = 21,
  col.nums = col.nums,
  scan.name.in.file = 1,
  diagnostics.filename = diag.geno.file,
  verbose = FALSE)
# Look at the values included in the "diag.geno" object which holds
#   all output from the function call
names(diag.geno)
# `read.file' is a vector indicating whether (1) or not (0) each file
#   specified in the `files' argument was read successfully
table(diag.geno$read.file)
# `row.num' is a vector of the number of rows read from each file
table(diag.geno$row.num)
# `sample.match' is a vector indicating whether (1) or not (0)
#   the sample name inside the raw text file matches that in the
#   sample annotation data.frame
table(diag.geno$sample.match)
# `snp.chk' is a vector indicating whether (1) or not (0)
#   the raw text file has the expected set of SNP names
table(diag.geno$snp.chk)
# `chk' is a vector indicating whether (1) or not (0) all previous
#   checks were successful and the data were written to the data file
table(diag.geno$chk)


###################################################
### code chunk number 7: DataCleaning.Rnw:390-411
###################################################
check.geno.file <- "check.geno.RData"
check.geno <- checkGenotypeFile(path = path,
  filename = geno.file,
  file.type = "gds",
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation,
  sep.type = ",",
  skip.num = 11,
  col.total = 21,
  col.nums = col.nums,
  scan.name.in.file = 1,
  check.scan.index = 1:3,
  n.scans.loaded = 3,
  diagnostics.filename = check.geno.file,
  verbose = FALSE)
# Look at the values included in the "check.geno" object which holds
#   all output from the function call
names(check.geno)
# 'geno.chk' is a vector indicating whether (1) or not (0) the genotypes
#   match the text file
table(check.geno$geno.chk)


###################################################
### code chunk number 8: DataCleaning.Rnw:421-431
###################################################
(gds <- GdsGenotypeReader(geno.file))
nscan(gds)
nsnp(gds)
head(getScanID(gds))
head(getSnpID(gds))
head(getChromosome(gds))
head(getPosition(gds))
# genotypes for the first 3 samples and  the first 5 SNPs
getGenotype(gds, snp=c(1,5), scan=c(1,3))
close(gds)


###################################################
### code chunk number 9: DataCleaning.Rnw:459-478
###################################################
qxy.file <- "tmp.qxy.gds"

col.nums <- as.integer(c(1,2,5,16,17))
names(col.nums) <- c("snp", "sample", "quality", "X", "Y")

diag.qxy.file <- "diag.qxy.RData"
diag.qxy <- createDataFile(path = path,
  filename = qxy.file,
  file.type = "gds",
  variables = c("quality","X","Y"),
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation,
  sep.type = ",",
  skip.num = 11,
  col.total = 21,
  col.nums = col.nums,
  scan.name.in.file = 1,
  diagnostics.filename = diag.qxy.file,
  verbose = FALSE)


###################################################
### code chunk number 10: DataCleaning.Rnw:484-499
###################################################
check.qxy.file <- "check.qxy.RData"
check.qxy <- checkIntensityFile(path = path,
  filename = qxy.file,
  file.type = "gds",
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation,
  sep.type = ",",
  skip.num = 11,
  col.total = 21,
  col.nums = col.nums,
  scan.name.in.file = 1,
  check.scan.index = 1:3,
  n.scans.loaded = 3,
  diagnostics.filename = check.qxy.file,
  verbose = FALSE)


###################################################
### code chunk number 11: DataCleaning.Rnw:508-522
###################################################
(gds <- GdsIntensityReader(qxy.file))
nscan(gds)
nsnp(gds)
head(getScanID(gds))
head(getSnpID(gds))
head(getChromosome(gds))
head(getPosition(gds))
# quality score for the first 3 samples and  the first 5 SNPs
getQuality(gds, snp=c(1,5), scan=c(1,3))
# X intensity for the first 3 samples and  the first 5 SNPs
getX(gds, snp=c(1,5), scan=c(1,3))
# Y intensity for the first 3 samples and  the first 5 SNPs
getY(gds, snp=c(1,5), scan=c(1,3))
close(gds)


###################################################
### code chunk number 12: DataCleaning.Rnw:595-614
###################################################
bl.file <- "tmp.bl.gds"

col.nums <- as.integer(c(1,2,20,21))
names(col.nums) <- c("snp", "sample", "BAlleleFreq", "LogRRatio")

diag.bl.file <- "diag.bl.RData"
diag.bl <- createDataFile(path = path,
  filename = bl.file,
  file.type = "gds",
  variables = c("BAlleleFreq","LogRRatio"),
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation,
  sep.type = ",",
  skip.num = 11,
  col.total = 21,
  col.nums = col.nums,
  scan.name.in.file = 1,
  diagnostics.filename = diag.bl.file,
  verbose = FALSE)


###################################################
### code chunk number 13: DataCleaning.Rnw:621-625
###################################################
(gds <- GdsIntensityReader(bl.file))
getBAlleleFreq(gds, snp=c(1,5), scan=c(1,3))
getLogRRatio(gds, snp=c(1,5), scan=c(1,3))
close(gds)


###################################################
### code chunk number 14: DataCleaning.Rnw:629-632
###################################################
file.remove(geno.file, qxy.file, bl.file)
file.remove(diag.geno.file, diag.qxy.file, diag.bl.file)
file.remove(check.geno.file, check.qxy.file)


###################################################
### code chunk number 15: DataCleaning.Rnw:648-657
###################################################
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(genofile)
# only GDS file
genoData <- GenotypeData(gds)
# with scan annotation
genoData <- GenotypeData(gds, scanAnnot=scanAnnot)
# with scan and SNP annotation
genoData <- GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
genoData


###################################################
### code chunk number 16: DataCleaning.Rnw:664-683
###################################################
nsnp(genoData)
nscan(genoData)
# scan annotation
range(getScanID(genoData))
hasSex(genoData)
table(getSex(genoData))
hasScanVariable(genoData, "subjectID")
head(getScanVariable(genoData, "subjectID"))
getScanVariableNames(genoData)
# snp annotation
range(getSnpID(genoData))
table(getChromosome(genoData, char=TRUE))
head(getPosition(genoData))
hasSnpVariable(genoData, "rsID")
head(getSnpVariable(genoData, "rsID"))
getSnpVariableNames(genoData)
# genotypes
getGenotype(genoData, snp=c(1,5), scan=c(1,5))
close(genoData)


###################################################
### code chunk number 17: DataCleaning.Rnw:690-707
###################################################
# quality score, X and X intensity
qxyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
gds <- GdsIntensityReader(qxyfile)
qxyData <- IntensityData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
qxyData
getQuality(qxyData, snp=c(1,5), scan=c(1,5))
getX(qxyData, snp=c(1,5), scan=c(1,5))
getY(qxyData, snp=c(1,5), scan=c(1,5))
close(qxyData)
# BAF/LRR
blfile <- system.file("extdata", "illumina_bl.gds", package="GWASdata")
gds <- GdsIntensityReader(blfile)
blData <- IntensityData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
blData
getBAlleleFreq(blData, snp=c(1,5), scan=c(1,5))
getLogRRatio(blData, snp=c(1,5), scan=c(1,5))
close(blData)


###################################################
### code chunk number 18: DataCleaning.Rnw:754-767
###################################################
# open the GDS file and create a GenotypeData object
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(gdsfile)
# sex is required for this function, so we need the scan annotation
genoData <-  GenotypeData(gds, scanAnnot=scanAnnot)
# Calculate the number of missing calls for each snp over all samples
#     for each sex separately
miss <- missingGenotypeBySnpSex(genoData)
# Examine the results
names(miss)
head(miss$missing.counts)
miss$scans.per.sex
head(miss$missing.fraction)


###################################################
### code chunk number 19: DataCleaning.Rnw:776-783
###################################################
# Make sure ordering matches snp annotation
allequal(snpAnnot$snpID, as.numeric(names(miss$missing.fraction)))
snpAnnot$missing.n1 <- miss$missing.fraction
varMetadata(snpAnnot)["missing.n1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples",
  "except that females are excluded for Y chr SNPs")
summary(snpAnnot$missing.n1)


###################################################
### code chunk number 20: DataCleaning.Rnw:790-793
###################################################
hist(snpAnnot$missing.n1, ylim=c(0,100),
     xlab="SNP missing call rate",
     main="Missing Call Rate for All Probes")


###################################################
### code chunk number 21: DataCleaning.Rnw:796-807
###################################################
# Find the number of SNPs with every call missing
length(snpAnnot$missing.n1[snpAnnot$missing.n1 == 1])
# Fraction of autosomal SNPs with missing call rate < 0.05
x <- snpAnnot$missing.n1[snpAnnot$chromosome < 23]
length(x[x < 0.05]) / length(x)
# Fraction of X chromosome SNPs with missing call rate < 0.05
x <- snpAnnot$missing.n1[snpAnnot$chromosome == 23]
length(x[x < 0.05]) / length(x)
# Fraction of Y chromosome SNPs with missing call rate < 0.05
x <- snpAnnot$missing.n1[snpAnnot$chromosome == 25]
length(x[x < 0.05]) / length(x)


###################################################
### code chunk number 22: DataCleaning.Rnw:825-839
###################################################
# Want to exclude all SNP probes with 100% missing call rate
# Check on how many SNPs to exclude
sum(snpAnnot$missing.n1 == 1)
# Create a variable that contains the IDs of these SNPs to exclude
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n1 == 1]
length(snpexcl)
# Calculate the missing call rate per sample
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
names(miss)
head(miss$missing.counts)
head(miss$snps.per.chr)
# Check to make sure that the correct number of SNPs were excluded
sum(miss$snps.per.chr)
nrow(snpAnnot) - sum(miss$snps.per.chr)


###################################################
### code chunk number 23: DataCleaning.Rnw:846-855
###################################################
head(miss$missing.fraction)
# Check the ordering matches the sample annotation file
allequal(names(miss$missing.fraction), scanAnnot$scanID)
# Add the missing call rates vector to the sample annotation file
scanAnnot$missing.e1 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n1<1",
  "except that Y chr SNPs are excluded for females")
summary(scanAnnot$missing.e1)


###################################################
### code chunk number 24: DataCleaning.Rnw:867-870
###################################################
hist(scanAnnot$missing.e1,
     xlab="Fraction of missing calls over all probes",
     main="Histogram of Sample Missing Call Rate for all Samples")


###################################################
### code chunk number 25: DataCleaning.Rnw:873-879
###################################################
# Look at missing.e1 for males
summary(scanAnnot$missing.e1[scanAnnot$sex == "M"])
# Look at missing.e1 for females
summary(scanAnnot$missing.e1[scanAnnot$sex == "F"])
# Number of samples with missing call rate > 5%
sum(scanAnnot$missing.e1 > 0.05)


###################################################
### code chunk number 26: DataCleaning.Rnw:889-911
###################################################
auto <- colnames(miss$missing.counts) %in% 1:22
missa <- rowSums(miss$missing.counts[,auto]) / sum(miss$snps.per.chr[auto])
summary(missa)
missx <- miss$missing.counts[,"X"] / miss$snps.per.chr["X"]
summary(missx)
# check they match sample annotation file
allequal(names(missa), scanAnnot$scanID)
allequal(names(missx), scanAnnot$scanID)
# Add these separate sample missing call rates to the sample
# annotation
scanAnnot$miss.e1.auto <- missa
scanAnnot$miss.e1.xchr <- missx
# Order scanAnnot by missing.e1 so duplicate subjectIDs
# with a higher missing rate are marked as duplicates
scanAnnot <- scanAnnot[order(scanAnnot$subjectID, scanAnnot$missing.e1),]
scanAnnot$duplicated <- duplicated(scanAnnot$subjectID)
table(scanAnnot$duplicated, useNA="ifany")
# Put scanAnnot back in scanID order; this is very important!!
scanAnnot <- scanAnnot[order(scanAnnot$scanID),]
allequal(scanAnnot$scanID, sort(scanAnnot$scanID))
varMetadata(scanAnnot)["duplicated", "labelDescription"] <-
  "TRUE for duplicate scan with higher missing.e1"


###################################################
### code chunk number 27: DataCleaning.Rnw:928-938
###################################################
# Find the samples with missing.e1 > .05 and make a vector of
# scanID to exclude from the calculation
scan.exclude <- scanAnnot$scanID[scanAnnot$missing.e1 > 0.05]
# Call missingGenotypeBySnpSex and save the output
miss <- missingGenotypeBySnpSex(genoData, scan.exclude=scan.exclude)
snpAnnot$missing.n2 <- miss$missing.fraction
varMetadata(snpAnnot)["missing.n2", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples with missing.e1<0.05",
  "except that females are excluded for Y chr SNPs")
summary(snpAnnot$missing.n2)


###################################################
### code chunk number 28: DataCleaning.Rnw:948-958
###################################################
# Create a vector of the SNPs to exclude.
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n2 >= 0.05]
length(snpexcl)
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
# Add the missing call rates vector to the sample annotation file
scanAnnot$missing.e2 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e2", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n2<0.05",
  "except that Y chr SNPs are excluded for females")
summary(scanAnnot$missing.e2)


###################################################
### code chunk number 29: DataCleaning.Rnw:964-967
###################################################
hist(scanAnnot$missing.e2, xlab="Fraction of missing calls over all probes
     with missing call rate < 0.05",
     main="Histogram of Sample Missing Call Rate for all Samples")


###################################################
### code chunk number 30: DataCleaning.Rnw:980-984
###################################################
varLabels(scanAnnot)
# Check how many batches exist and how many samples are in each batch
length(unique(scanAnnot$batch))
table(scanAnnot$batch, useNA="ifany")


###################################################
### code chunk number 31: DataCleaning.Rnw:987-991
###################################################
# Plot the distribution of the number of samples per batch.
barplot(table(scanAnnot$batch),
        ylab="Number of Samples", xlab="Batch",
        main="Distribution of Samples per Batch")


###################################################
### code chunk number 32: DataCleaning.Rnw:994-1003
###################################################
# Examine the mean missing call rate per batch for all SNPs
batches <- unique(scanAnnot$batch)
bmiss <- rep(NA,length(batches)); names(bmiss) <- batches
bn <- rep(NA,length(batches)); names(bn) <- batches
for(i in 1:length(batches)) {
  x <- scanAnnot$missing.e1[is.element(scanAnnot$batch, batches[i])]
  bmiss[i] <- mean(x)
  bn[i] <- length(x)
}


###################################################
### code chunk number 33: DataCleaning.Rnw:1014-1016
###################################################
y <- lm(bmiss ~ bn)
anova(y)


###################################################
### code chunk number 34: DataCleaning.Rnw:1019-1023
###################################################
plot(bn, bmiss,
 xlab="Number of samples per batch", ylab="Mean missing call rate",
 main="Mean Missing Call Rate vs\nSamples per Batch")
abline(y$coefficients)


###################################################
### code chunk number 35: DataCleaning.Rnw:1044-1052
###################################################
res <- batchChisqTest(genoData, batchVar="batch", return.by.snp=TRUE)
close(genoData)
# chi-square values for each SNP
dim(res$chisq)
# genomic inflation factor
res$lambda
# average chi-square test statistic for each of the batches
res$mean.chisq


###################################################
### code chunk number 36: DataCleaning.Rnw:1057-1075
###################################################
x <- table(scanAnnot$race, useNA="ifany")
x
x[1] / sum(x)
x[2] / sum(x)
x <- table(scanAnnot$race, scanAnnot$batch)
x
# Run an approximate chi-square test to see if there are ethnic effects
chisq <- chisq.test(x)
chisq$p.value
# Calculate the fraction of samples in each batch that are CEU
batches <- unique(scanAnnot$batch)
eth <- rep(NA,length(batches)); names(eth) <- sort(batches)
for(i in 1:length(batches)){
  x <- scanAnnot$race[is.element(scanAnnot$batch, batches[i])]
  xl <- length(x[x == "CEU"])
  eth[i] <- xl / length(x)
}
allequal(names(eth), names(res$mean.chisq))


###################################################
### code chunk number 37: DataCleaning.Rnw:1078-1085
###################################################
# Plot the average Chi-Square test statistic against the
#     fraction of samples that are CEU
plot(eth, res$mean.chisq, xlab="Fraction of CEU Samples per Batch",
  ylab="Average Chi-square Test Statistic",
  main="Fraction of CEU Samples per Batch
  vs Average Chi-square Test Statistic")
abline(v=mean(eth), lty=2, col="red")


###################################################
### code chunk number 38: DataCleaning.Rnw:1116-1124
###################################################
qxyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
qualGDS <- GdsIntensityReader(qxyfile)
qualData <- IntensityData(qualGDS, scanAnnot=scanAnnot)
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genoGDS <- GdsGenotypeReader(genofile)
genoData <- GenotypeData(genoGDS, scanAnnot=scanAnnot)
qual.results <- qualityScoreByScan(qualData, genoData)
close(qualData)


###################################################
### code chunk number 39: DataCleaning.Rnw:1131-1133
###################################################
hist(qual.results[,"median.quality"], main="Median Genotype Quality Scores
  of Samples", xlab="Median Quality")


###################################################
### code chunk number 40: DataCleaning.Rnw:1164-1171
###################################################
blfile <- system.file("extdata", "illumina_bl.gds", package="GWASdata")
blGDS <- GdsIntensityReader(blfile)
blData <- IntensityData(blGDS, scanAnnot=scanAnnot)
nbins <- rep(12, 3)
slidingBAF12 <- sdByScanChromWindow(blData, genoData, nbins=nbins)
names(slidingBAF12)
dim(slidingBAF12[["21"]])


###################################################
### code chunk number 41: DataCleaning.Rnw:1179-1182
###################################################
sds.chr <- meanSdByChromWindow(slidingBAF12, scanAnnot$sex)
sds.chr[["21"]]
sds.chr[["X"]]


###################################################
### code chunk number 42: DataCleaning.Rnw:1188-1192
###################################################
res12bin4sd <- findBAFvariance(sds.chr, slidingBAF12, scanAnnot$sex,
                               sd.threshold=4)
head(res12bin4sd)
table(res12bin4sd[, "chromosome"])


###################################################
### code chunk number 43: DataCleaning.Rnw:1199-1205
###################################################
scanID <- as.integer(res12bin4sd[, "scanID"])
chrom <- as.integer(res12bin4sd[, "chromosome"])
chrom[res12bin4sd[, "chromosome"] == "X"] <- 23
bincode <- paste("Bin", res12bin4sd[, "bin"], sep = " ")
chromIntensityPlot(blData, scanID, chrom, info=bincode, ideogram=FALSE)
close(blData)


###################################################
### code chunk number 44: DataCleaning.Rnw:1231-1235
###################################################
miss <- missingGenotypeByScanChrom(genoData)
miss.rate <- t(apply(miss$missing.counts, 1, function(x) {
  x / miss$snps.per.chr}))
miss.rate <- as.data.frame(miss.rate)


###################################################
### code chunk number 45: DataCleaning.Rnw:1239-1242
###################################################
cols <- names(miss.rate) %in% c(1:22, "X", "XY")
boxplot(miss.rate[,cols], main="Missingness by Chromosome",
  ylab="Proportion Missing", xlab="Chromosome")


###################################################
### code chunk number 46: DataCleaning.Rnw:1245-1248
###################################################
boxplot(miss.rate$X ~ scanAnnot$sex,
  main="X Chromosome Missingness by Sex",
  ylab="Proportion Missing")


###################################################
### code chunk number 47: DataCleaning.Rnw:1255-1267
###################################################
# Calculate heterozygosity by scan by chromosome
het.results <- hetByScanChrom(genoData)
close(genoData)
# Ensure heterozygosity results are ordered correctly
allequal(scanAnnot$scanID, rownames(het.results))
# Write autosomal and X chr heterozygosity to sample annot
scanAnnot$het.A <- het.results[,"A"]
scanAnnot$het.X <- het.results[,"X"]
varMetadata(scanAnnot)["het.A", "labelDescription"] <-
  "fraction of heterozygotes for autosomal SNPs"
varMetadata(scanAnnot)["het.X", "labelDescription"] <-
  "fraction of heterozygotes for X chromosome SNPs"


###################################################
### code chunk number 48: DataCleaning.Rnw:1275-1277
###################################################
boxplot(scanAnnot$het.A ~ scanAnnot$race,
  main="Autosomal Heterozygosity")


###################################################
### code chunk number 49: DataCleaning.Rnw:1280-1283
###################################################
female <- scanAnnot$sex == "F"
boxplot(scanAnnot$het.X[female] ~ scanAnnot$race[female],
  main="X Chromosome Heterozygosity in Females")


###################################################
### code chunk number 50: DataCleaning.Rnw:1320-1325
###################################################
qxyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
intenGDS <- GdsIntensityReader(qxyfile)
inten.by.chrom <- meanIntensityByScanChrom(intenGDS)
close(intenGDS)
names(inten.by.chrom)


###################################################
### code chunk number 51: DataCleaning.Rnw:1333-1343
###################################################
mninten <- inten.by.chrom[[1]]  # mean intensities
dim(mninten)
# Check to be sure sample ordering is consistent
allequal(scanAnnot$scanID, rownames(mninten))
# Assign each sex a color
xcol <- rep(NA, nrow(scanAnnot))
xcol[scanAnnot$sex == "M"] <- "blue"
xcol[scanAnnot$sex == "F"] <- "red"
nx <- sum(snpAnnot$chromosome == 23)
ny <- sum(snpAnnot$chromosome == 25)


###################################################
### code chunk number 52: DataCleaning.Rnw:1353-1385
###################################################
#All intensities
x1 <-mninten[,"X"]; y1 <- mninten[,"Y"]
main1 <- "Mean X vs \nMean Y Chromosome Intensity"
#Het on X vs X intensity
x2 <- mninten[,"X"]; y2 <- scanAnnot$het.X
main2 <- "Mean X Chromosome Intensity vs
  Mean X Chromosome Heterozygosity"
# Het on X vs Y intensity
y3 <- mninten[,"Y"]; x3 <- scanAnnot$het.X
main3 <- "Mean X Chromosome Heterozygosity vs
  Mean Y Chromosome Intensity"
# X vs A het
x4 <- scanAnnot$het.A[scanAnnot$sex == "F"]
y4 <- scanAnnot$het.X[scanAnnot$sex == "F"]
main4 <- "Mean Autosomal Heterozygosity vs
  Mean X Chromosome Heterozygosity"
cols <- c("blue","red")
mf <- c("male", "female")
xintenlab <- paste("X intensity (n=", nx, ")", sep="")
yintenlab <- paste("Y intensity (n=", ny, ")", sep="")
pdf("DataCleaning-sex.pdf")
par(mfrow=c(2,2))
plot(x1, y1, xlab=xintenlab, ylab=yintenlab,
  main=main1, col=xcol, cex.main=0.8)
legend("topright",mf,col=cols,pch=c(1,1))
plot(x2, y2, col=xcol, xlab=xintenlab,
  ylab="X heterozygosity", main=main2, cex.main=0.8)
plot(x3, y3, col=xcol, ylab=yintenlab,
  xlab="X heterozygosity", main=main3, cex.main=0.8)
plot(x4,y4, col="red", xlab="Autosomal heterozygosity",
  ylab="X heterozygosity", main=main4, cex.main=0.8)
dev.off()


###################################################
### code chunk number 53: DataCleaning.Rnw:1412-1420
###################################################
library(SNPRelate)
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gdsobj <- snpgdsOpen(gdsfile)
ibdobj <- snpgdsIBDKING(gdsobj)
snpgdsClose(gdsobj)
names(ibdobj)
dim(ibdobj$kinship)
ibdobj$kinship[1:5,1:5]


###################################################
### code chunk number 54: DataCleaning.Rnw:1431-1436
###################################################
ped <- pData(scanAnnot)[,c("family", "subjectID", "father", "mother", "sex")]
dim(ped)
names(ped) <- c("family", "individ", "father", "mother", "sex")
ped[1:5, ]
(chk <- pedigreeCheck(ped))


###################################################
### code chunk number 55: DataCleaning.Rnw:1442-1445
###################################################
dups <- chk$duplicates
uni.ped <- pedigreeDeleteDuplicates(ped, dups)
(chk <- pedigreeCheck(uni.ped))


###################################################
### code chunk number 56: DataCleaning.Rnw:1450-1456
###################################################
ni <- chk$parent.no.individ.entry
parent <- data.frame(family=ni$family, individ=ni$parentID,
                     father=0, mother=0, sex="F",
                     stringsAsFactors=FALSE)
ped.complete <- rbind(uni.ped, parent)
(chk <- pedigreeCheck(ped.complete))


###################################################
### code chunk number 57: DataCleaning.Rnw:1464-1472
###################################################
ped.complete <- ped.complete[ped.complete$family != 58,]
subf <- chk$subfamilies.ident
table(subf$family)
subf.ids <- subf$individ[subf$subfamily == 2]
newfam <- ped.complete$individ %in% subf.ids
ped.complete$family[newfam] <- paste0(ped.complete$family[newfam], "-2")
table(ped.complete$family)
pedigreeCheck(ped.complete)


###################################################
### code chunk number 58: DataCleaning.Rnw:1482-1487
###################################################
rels <- pedigreePairwiseRelatedness(ped.complete)
length(rels$inbred.fam)
relprs <- rels$relativeprs
relprs[1:5,]
table(relprs$relation)


###################################################
### code chunk number 59: DataCleaning.Rnw:1497-1515
###################################################
samp <- pData(scanAnnot)[,c("scanID", "subjectID")]
samp <- samp[match(ibdobj$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")
ibd <- snpgdsIBDSelection(ibdobj, kinship.cutoff=1/32)
ibd <- merge(ibd, samp, by.x="ID1", by.y="scanID")
ibd <- merge(ibd, samp, by.x="ID2", by.y="scanID", suffixes=c("1","2"))
ibd$ii <- pasteSorted(ibd$Individ1, ibd$Individ2)

relprs$ii <- pasteSorted(relprs$Individ1, relprs$Individ2)
ibd <- merge(ibd, relprs[,c("ii","relation")], all.x=TRUE)
names(ibd)[names(ibd) == "relation"] <- "exp.rel"
ibd$exp.rel[ibd$Individ1 == ibd$Individ2] <- "Dup"
ibd$exp.rel[is.na(ibd$exp.rel)] <- "U"
table(ibd$exp.rel, useNA="ifany")

# assign observed relationships
ibd$obs.rel <- ibdAssignRelatednessKing(ibd$IBS0, ibd$kinship)
table(ibd$obs.rel, useNA="ifany")


###################################################
### code chunk number 60: DataCleaning.Rnw:1521-1533
###################################################
## thresholds for assigning relationships using kinship coefficients
## in table 1 of Manichaikul (2010)
cut.dup <- 1/(2^(3/2))
cut.deg1 <- 1/(2^(5/2))
cut.deg2 <- 1/(2^(7/2))
cut.deg3 <- 1/(2^(9/2))
cols <- c(Dup="magenta", PO="cyan", U="black")

plot(ibd$IBS0, ibd$kinship, col=cols[ibd$exp.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=names(cols), col=cols, pch=1)


###################################################
### code chunk number 61: DataCleaning.Rnw:1562-1585
###################################################
filt <- get(data(pcaSnpFilters.hg18))
chrom <- getChromosome(snpAnnot)
pos <- getPosition(snpAnnot)
snpID <- getSnpID(snpAnnot)
snp.filt <- rep(TRUE, length(snpID))
for (f in 1:nrow(filt)) {
  snp.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos
           & pos < filt$end.base[f]] <- FALSE
}
snp.sel <- snpID[snp.filt]
length(snp.sel)

sample.sel <- scanAnnot$scanID[!scanAnnot$duplicated]
length(sample.sel)

gdsobj <- snpgdsOpen(gdsfile)
snpset <- snpgdsLDpruning(gdsobj, sample.id=sample.sel, snp.id=snp.sel,
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(0.1))

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned)


###################################################
### code chunk number 62: DataCleaning.Rnw:1591-1595
###################################################
pca <- snpgdsPCA(gdsobj, sample.id=sample.sel, snp.id=snp.pruned)
names(pca)
length(pca$eigenval)
dim(pca$eigenvect)


###################################################
### code chunk number 63: DataCleaning.Rnw:1601-1609
###################################################
# Calculate the percentage of variance explained
# by each principal component.
pc.frac <- pca$eigenval/sum(pca$eigenval)
lbls <- paste("EV", 1:4, "\n", format(pc.frac[1:4], digits=2), sep="")
samp <- pData(scanAnnot)[match(pca$sample.id, scanAnnot$scanID),]
cols <- rep(NA, nrow(samp))
cols[samp$race == "CEU"] <- "blue"
cols[samp$race == "YRI"] <- "red"


###################################################
### code chunk number 64: DataCleaning.Rnw:1612-1614
###################################################
pairs(pca$eigenvect[,1:4], col=cols, labels=lbls,
  main = "CEU: blue, YRI: red")


###################################################
### code chunk number 65: DataCleaning.Rnw:1625-1638
###################################################
par.coord <- pca$eigenvect
rangel <- apply(par.coord, 2, function(x) range(x)[1])
rangeh <- apply(par.coord, 2, function(x) range(x)[2])
std.coord <- par.coord
for (i in 1:14)
  std.coord[,i] <- (par.coord[,i] - rangel[i])/(rangeh[i]-rangel[i])
plot(c(0,15), c(0,1), type = 'n', axes = FALSE, ylab = "", xlab = "",
 main = "Parallel Coordinates Plot
 CEU: blue, YRI: red")
for (j in 1:13)
  for (i in sample(1:nrow(std.coord)) )
    lines(c(j,j+1), std.coord[i,c(j,j+1)], col=cols[i], lwd=0.25)
axis(1, at = 1:14, labels = paste("PC",1:14, sep = "."))


###################################################
### code chunk number 66: DataCleaning.Rnw:1648-1659
###################################################
corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:4)
snpgdsClose(gdsobj)
snp <- snpAnnot[match(corr$snp.id, snpID),]
chrom <- getChromosome(snp, char=TRUE)
pdf("DataCleaning-corr.pdf")
par(mfrow=c(4,1))
for (i in 1:4) {
  snpCorrelationPlot(abs(corr$snpcorr[i,]), chrom,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()


###################################################
### code chunk number 67: DataCleaning.Rnw:1680-1685
###################################################
princomp <- as.data.frame(pca$eigenvect)
samples.nodup <- pData(scanAnnot)[!scanAnnot$duplicated,]
princomp$scanID <- as.factor(samples.nodup$scanID)
princomp$case.ctrl.status <- as.factor(samples.nodup$status)
princomp$race <- as.factor(samples.nodup$race)


###################################################
### code chunk number 68: DataCleaning.Rnw:1691-1698
###################################################
pc.percent <- 100 * pca$eigenval[1:32]/sum(pca$eigenval)
pc.percent
lbls <- paste("EV", 1:3, "\n", format(pc.percent[1:3], digits=2), "%", sep="")
table(samples.nodup$status)
cols <- rep(NA, nrow(samples.nodup))
cols[samples.nodup$status == 1] <- "green"
cols[samples.nodup$status == 0] <- "magenta"


###################################################
### code chunk number 69: DataCleaning.Rnw:1708-1710
###################################################
pairs(pca$eigenvect[,1:3], col=cols, labels=lbls,
  main = "First Three EVs by Case-Control Status")


###################################################
### code chunk number 70: DataCleaning.Rnw:1713-1715
###################################################
boxplot(princomp[, 1] ~ princomp$case.ctrl.status,
  ylab = "PC1", main = "PC1 vs. Case-control Status")


###################################################
### code chunk number 71: DataCleaning.Rnw:1718-1720
###################################################
boxplot(princomp[, 2] ~ princomp$case.ctrl.status,
  ylab = "PC2", main = "PC2 vs. Case-control Status")


###################################################
### code chunk number 72: DataCleaning.Rnw:1723-1725
###################################################
boxplot(princomp[, 3] ~ princomp$case.ctrl.status,
  ylab = "PC3", main = "PC3 vs. Case-control Status")


###################################################
### code chunk number 73: DataCleaning.Rnw:1728-1737
###################################################
aov.p1 <- aov(princomp[,1] ~ princomp$race *
  princomp$case.ctrl.status, princomp)
summary(aov.p1)
aov.p2 <- aov(princomp[,2] ~ princomp$race *
  princomp$case.ctrl.status, princomp)
summary(aov.p2)
aov.p3 <- aov(princomp[,3] ~ princomp$race *
  princomp$case.ctrl.status, princomp)
summary(aov.p3)


###################################################
### code chunk number 74: DataCleaning.Rnw:1757-1759
###################################################
lm.all <- lm(scanAnnot$missing.e1 ~ scanAnnot$status)
summary(aov(lm.all))


###################################################
### code chunk number 75: DataCleaning.Rnw:1762-1764
###################################################
boxplot(scanAnnot$missing.e1 ~ scanAnnot$status, ylab =
  "Mean missing call rate", main="Mean missing call rate by case status")


###################################################
### code chunk number 76: DataCleaning.Rnw:1780-1787
###################################################
blfile <- system.file("extdata", "illumina_bl.gds", package="GWASdata")
blgds <- GdsIntensityReader(blfile)
blData <-  IntensityData(blgds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genogds <- GdsGenotypeReader(genofile)
genoData <-  GenotypeData(genogds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)


###################################################
### code chunk number 77: DataCleaning.Rnw:1792-1795
###################################################
baf.sd <- sdByScanChromWindow(blData, genoData, var="BAlleleFreq")
med.baf.sd <- medianSdOverAutosomes(baf.sd)
low.qual.ids <- med.baf.sd$scanID[med.baf.sd$med.sd > 0.05]


###################################################
### code chunk number 78: DataCleaning.Rnw:1800-1817
###################################################
chrom <- getChromosome(snpAnnot, char=TRUE)
pos <- getPosition(snpAnnot)
hla.df <- get(data(HLA.hg18))
hla <- chrom == "6" & pos >= hla.df$start.base & pos <= hla.df$end.base
xtr.df <- get(data(pseudoautosomal.hg18))
xtr <- chrom == "X" & pos >= xtr.df["X.XTR", "start.base"] &
  pos <= xtr.df["X.XTR", "end.base"]
centromeres <- get(data(centromeres.hg18))
gap <- rep(FALSE, nrow(snpAnnot))
for (i in 1:nrow(centromeres)) {
  ingap <- chrom == centromeres$chrom[i] & pos > centromeres$left.base[i] &
    pos < centromeres$right.base[i]
  gap <- gap | ingap
}
ignore <- snpAnnot$missing.n1 == 1 #ignore includes intensity-only and failed snps
snp.exclude <- ignore | hla | xtr | gap
snp.ok <- snpAnnot$snpID[!snp.exclude]


###################################################
### code chunk number 79: DataCleaning.Rnw:1822-1827
###################################################
scan.ids <- scanAnnot$scanID[1:10]
chrom.ids <- 21:23
baf.seg <- anomSegmentBAF(blData, genoData, scan.ids=scan.ids,
  chrom.ids=chrom.ids, snp.ids=snp.ok, verbose=FALSE)
head(baf.seg)


###################################################
### code chunk number 80: DataCleaning.Rnw:1832-1838
###################################################
baf.anom <- anomFilterBAF(blData, genoData, segments=baf.seg,
  snp.ids=snp.ok, centromere=centromeres, low.qual.ids=low.qual.ids,
  verbose=FALSE)
names(baf.anom)
baf.filt <- baf.anom$filtered
head(baf.filt)


###################################################
### code chunk number 81: DataCleaning.Rnw:1850-1856
###################################################
loh.anom <- anomDetectLOH(blData, genoData, scan.ids=scan.ids,
  chrom.ids=chrom.ids, snp.ids=snp.ok, known.anoms=baf.filt,
  verbose=FALSE)
names(loh.anom)
loh.filt <- loh.anom$filtered
head(loh.filt)


###################################################
### code chunk number 82: DataCleaning.Rnw:1865-1879
###################################################
# create required data frame
baf.filt$method <- "BAF"
if (!is.null(loh.filt)) {
  loh.filt$method <- "LOH"
  cols <- intersect(names(baf.filt), names(loh.filt))
  anoms <- rbind(baf.filt[,cols], loh.filt[,cols])
} else {
  anoms <- baf.filt
}
anoms$anom.id <- 1:nrow(anoms)

stats <- anomSegStats(blData, genoData, snp.ids=snp.ok, anom=anoms,
                      centromere=centromeres)
names(stats)


###################################################
### code chunk number 83: DataCleaning.Rnw:1885-1888
###################################################
snp.not.ok <- snpAnnot$snpID[snp.exclude]
anomStatsPlot(blData, genoData, anom.stats=stats[1,],
  snp.ineligible=snp.not.ok, centromere=centromeres, cex.leg=1)


###################################################
### code chunk number 84: DataCleaning.Rnw:1898-1900
###################################################
lrr.sd <- sdByScanChromWindow(blData, var="LogRRatio", incl.hom=TRUE)
med.lrr.sd <- medianSdOverAutosomes(lrr.sd)


###################################################
### code chunk number 85: DataCleaning.Rnw:1906-1908
###################################################
baf.seg.info <- baf.anom$seg.info
loh.seg.info <- loh.anom$base.info[,c("scanID", "chromosome", "num.segs")]


###################################################
### code chunk number 86: DataCleaning.Rnw:1918-1923
###################################################
snpAnnot$eligible <- !snp.exclude
baf.low.qual <- anomIdentifyLowQuality(snpAnnot, med.baf.sd, baf.seg.info,
  sd.thresh=0.1, sng.seg.thresh=0.0008, auto.seg.thresh=0.0001)
loh.low.qual <- anomIdentifyLowQuality(snpAnnot, med.lrr.sd, loh.seg.info,
  sd.thresh=0.25, sng.seg.thresh=0.0048, auto.seg.thresh=0.0006)


###################################################
### code chunk number 87: DataCleaning.Rnw:1928-1930
###################################################
close(blData)
close(genoData)


###################################################
### code chunk number 88: DataCleaning.Rnw:1945-1957
###################################################
# anomalies to filter
anom.filt <- stats[,c("scanID", "chromosome", "left.base", "right.base")]
# whole.chrom column is required and can be used for sex chromosome
#   anomalies such as XXX
anom.filt$whole.chrom <- FALSE
# select unique subjects
subj <- scanAnnot$scanID[!scanAnnot$duplicated]
subj.filt.file <- "subj_filt.gds"
setMissingGenotypes(genofile, subj.filt.file, anom.filt,
                    file.type="gds", sample.include=subj, verbose=FALSE)
(gds <- GdsGenotypeReader(subj.filt.file))
close(gds)


###################################################
### code chunk number 89: DataCleaning.Rnw:1982-1984
###################################################
scan.excl <- scanAnnot$scanID[scanAnnot$missing.e1 >= 0.05]
length(scan.excl)


###################################################
### code chunk number 90: DataCleaning.Rnw:1995-2026
###################################################
snp.excl <- snpAnnot$snpID[snpAnnot$missing.n1 == 1]
length(snp.excl)
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genoGDS <- GdsGenotypeReader(genofile)
genoData <- GenotypeData(genoGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
dupdisc <- duplicateDiscordance(genoData, subjName.col="subjectID",
  scan.exclude=scan.excl, snp.exclude=snp.excl)
names(dupdisc)
head(dupdisc$discordance.by.snp)
length(dupdisc$discordance.by.subject)
dupdisc$discordance.by.subject[[2]]
# each entry is a 2x2 matrix, but only one value of each
# is important since these are all pairs
npair <- length(dupdisc$discordance.by.subject)
disc.subj <- rep(NA, npair)
subjID <- rep(NA, npair)
race <- rep(NA, npair)
for (i in 1:npair) {
  disc.subj[i] <- dupdisc$discordance.by.subject[[i]][1,2]
  subjID[i] <- names(dupdisc$discordance.by.subject)[i]
  race[i] <- scanAnnot$race[scanAnnot$subjectID == subjID[i]][1]
}
dat <- data.frame(subjID=subjID, disc=disc.subj, pop=race,
                  stringsAsFactors=FALSE)
summary(dat$disc)
# Assign colors for the duplicate samples based on population group.
dat$col <- NA
dat$col[dat$pop == "CEU"] <- "blue"
dat$col[dat$pop == "YRI"] <- "red"
dat <- dat[order(dat$disc),]
dat$rank <- 1:npair


###################################################
### code chunk number 91: DataCleaning.Rnw:2029-2034
###################################################
# Plot the sample discordance rates color coded by race.
plot(dat$disc, dat$rank, col=dat$col, ylab="rank",
  xlab="Discordance rate between duplicate samples",
  main="Duplicate Sample Discordance by Continental Ancestry")
legend("bottomright", unique(dat$pop), pch=rep(1,2), col=unique(dat$col))


###################################################
### code chunk number 92: DataCleaning.Rnw:2058-2059
###################################################
duplicateDiscordanceProbability(npair)


###################################################
### code chunk number 93: DataCleaning.Rnw:2080-2092
###################################################
men.list <- with(pData(scanAnnot), mendelList(family, subjectID,
  father, mother, sex, scanID))
res <- mendelListAsDataFrame(men.list)
head(res)
dim(res)
# Only want to use SNPs with missing.n1 < 0.05
snp.excl <- snpAnnot$snpID[snpAnnot$missing.n1 >= 0.05]
length(snp.excl)
mend <- mendelErr(genoData, men.list, snp.exclude=snp.excl)
names(mend)
head(mend$trios)
names(mend$snp)


###################################################
### code chunk number 94: DataCleaning.Rnw:2102-2105
###################################################
# Calculate the error rate
err <- mend$snp$error.cnt / mend$snp$check.cnt
table(err == 0, useNA="ifany")


###################################################
### code chunk number 95: DataCleaning.Rnw:2108-2110
###################################################
plot(err, rank(err), xlab="Error Rate (fraction)",
  ylab="rank", main="Mendelian Error Rate per SNP, ranked")


###################################################
### code chunk number 96: DataCleaning.Rnw:2118-2130
###################################################
fam <- mend$snp$error.cnt
n <- mend$snp$check.cnt
summary(fam)
# SNPs with errors
length(fam[n > 0 & fam > 0])
# SNPs for which more than one family has an error
length(fam[n > 0 & fam > 1])
# Get the SNPs with valid trios for error detection
val <- length(fam[n > 0])
noerr <- length(fam[n > 0 & fam == 0])
# Divide to get fraction with no errors
noerr / val


###################################################
### code chunk number 97: DataCleaning.Rnw:2137-2150
###################################################
snp.sel <- match(names(mend$snp$error.cnt), snpAnnot$snpID)
snpAnnot$mendel.err.count[snp.sel] <- mend$snp$error.cnt
snpAnnot$mendel.err.sampsize[snp.sel] <- mend$snp$check.cnt
allequal(snpAnnot$snpID, sort(snpAnnot$snpID))
# The high number of NA values is due to the filtering out of SNPs
#    before the Mendelian error rate calculation
sum(is.na(snpAnnot$mendel.err.count))
sum(is.na(snpAnnot$mendel.err.sampsize))
varMetadata(snpAnnot)["mendel.err.count", "labelDescription"] <-
  paste("number of Mendelian errors detected in trios averaged over",
        "multiple combinations of replicate genotyping instances")
varMetadata(snpAnnot)["mendel.err.sampsize", "labelDescription"] <-
  "number of opportunities to detect Mendelian error in trios"


###################################################
### code chunk number 98: DataCleaning.Rnw:2162-2179
###################################################
# Get a vector of SNPs to check
snp <- pData(snpAnnot)
snp$err.rate <- snp$mendel.err.count /
  snp$mendel.err.sampsize
snp <- snp[order(snp$err.rate, decreasing=TRUE),]
snp <- snp[1:9,]
xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
xyGDS <- GdsIntensityReader(xyfile)
xyData <- IntensityData(xyGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
pdf(file="DataCleaning-mendel.pdf")
par(mfrow = c(3,3))
mtxt <- paste("SNP", snp$rsID, "\nMendelian Error Rate",
   format(snp$err.rate, digits=5))
genoClusterPlot(xyData, genoData, snpID=snp$snpID, main.txt=mtxt,
                cex.main=0.9)
dev.off()
close(xyData)


###################################################
### code chunk number 99: DataCleaning.Rnw:2192-2209
###################################################
# Calculate the fraction of SNPs with an error for each trio
trios <- mend$trios
trios$Mend.err <- trios$Men.err.cnt/trios$Men.cnt
summary(trios$Mend.err)
# Start by pulling out the vectors needed from `trios'
tmp <- trios[, c("fam.id", "Mend.err")]; dim(tmp)
# Change fam.id to match the sample annotation column name
names(tmp) <- c("family", "Mend.err.rate.fam")
# Merge the variables into the sample annotation file
scanAnnot$mend.err.rate.fam <- NA
for (i in 1:nrow(tmp)) {
  ind <- which(is.element(scanAnnot$family, tmp$family[i]))
  scanAnnot$mend.err.rate.fam[ind] <- tmp$Mend.err.rate.fam[i]
}
head(scanAnnot$mend.err.rate.fam)
varMetadata(scanAnnot)["mend.err.rate.fam", "labelDescription"] <-
  "Mendelian error rate per family"


###################################################
### code chunk number 100: DataCleaning.Rnw:2217-2228
###################################################
# Get the families that have non-NA values for the family
#     Mendelian error rate
fams <- pData(scanAnnot)[!is.na(scanAnnot$mend.err.rate.fam) &
  !duplicated(scanAnnot$family), c("family",
  "mend.err.rate.fam", "race")]
dim(fams)
table(fams$race, useNA="ifany")
# Assign colors for the different ethnicities in these families
pcol <- rep(NA, nrow(fams))
pcol[fams$race == "CEU"] <- "blue"
pcol[fams$race == "YRI"] <- "red"


###################################################
### code chunk number 101: DataCleaning.Rnw:2231-2236
###################################################
plot(fams$mend.err.rate.fam*100, rank(fams$mend.err.rate.fam),
  main="Mendelian Error rate per Family, ranked",
  xlab="Mendelian error rate per family (percent)",
  ylab="rank", col=pcol)
legend("bottomright", c("CEU", "YRI"), pch=c(1,1), col=c("blue", "red"))


###################################################
### code chunk number 102: DataCleaning.Rnw:2254-2269
###################################################
head(pData(scanAnnot)[,c("father", "mother")])
nonfounders <- scanAnnot$father != 0 &
               scanAnnot$mother != 0
table(nonfounders)
scan.excl <- scanAnnot$scanID[scanAnnot$race != "CEU" |
   nonfounders | scanAnnot$duplicated]
length(scan.excl)

chr <- getChromosome(genoData)
auto <- range(which(chr %in% 1:22))
X <- range(which(chr == 23))
hwe <- exactHWE(genoData, scan.exclude=scan.excl, snpStart=auto[1], snpEnd=auto[2])
hweX <- exactHWE(genoData, scan.exclude=scan.excl, snpStart=X[1], snpEnd=X[2])
hwe <- rbind(hwe, hweX)
close(genoData)


###################################################
### code chunk number 103: DataCleaning.Rnw:2275-2285
###################################################
names(hwe)
dim(hwe)
# Check on sample sizes for autosomes and X chromosome
hwe$N <- hwe$nAA + hwe$nAB + hwe$nBB
summary(hwe$N[is.element(hwe$chr,1:22)])
summary(hwe$N[is.element(hwe$chr,23)])
hwe$pval[1:10]
sum(is.na(hwe$pval[hwe$chr == 23])) # X
hwe$MAF[1:10]
hwe[1:3, c("nAA", "nAB", "nBB")]


###################################################
### code chunk number 104: DataCleaning.Rnw:2292-2293
###################################################
summary(hwe$f)


###################################################
### code chunk number 105: DataCleaning.Rnw:2296-2298
###################################################
hist(hwe$f, main="Histogram of the Inbreeding Coefficient
  For CEU Samples", xlab="Inbreeding Coefficient")


###################################################
### code chunk number 106: DataCleaning.Rnw:2301-2304
###################################################
# Check the MAF of those SNPs with f=1
chkf <- hwe[!is.na(hwe$f) & hwe$f==1,]; dim(chkf)
summary(chkf$MAF)


###################################################
### code chunk number 107: DataCleaning.Rnw:2311-2322
###################################################
hwe.0 <- hwe[hwe$MAF > 0,]; dim(hwe.0)
# Only keep the autosomal SNPs for first plot
pval <- hwe.0$pval[is.element(hwe.0$chr, 1:22)]
length(pval)
pval <- pval[!is.na(pval)]
length(pval)
# X chromosome SNPs for plot 2
pval.x <- hwe.0$pval[is.element(hwe.0$chr, 23)]
length(pval.x)
pval.x <- pval.x[!is.na(pval.x)]
length(pval.x)


###################################################
### code chunk number 108: DataCleaning.Rnw:2324-2331
###################################################
pdf(file = "DataCleaning-hwe.pdf")
par(mfrow=c(2,2))
qqPlot(pval=pval, truncate = FALSE, main="Autosomes, all")
qqPlot(pval=pval, truncate = TRUE, main="Autosomes, truncated")
qqPlot(pval=pval.x, truncate = FALSE, main="X chromosome, all")
qqPlot(pval=pval.x, truncate = TRUE, main="X chromosome, truncated")
dev.off()


###################################################
### code chunk number 109: DataCleaning.Rnw:2337-2340
###################################################
plot(hwe.0$MAF, -log10(hwe.0$pval),
  xlab="Minor Allele Frequency", ylab="-log(p-value)",
  main="Minor Allele Frequency vs\nP-value")


###################################################
### code chunk number 110: DataCleaning.Rnw:2382-2398
###################################################
genoGDS <- GdsGenotypeReader(subj.filt.file)
subjAnnot <- scanAnnot[scanAnnot$scanID %in% getScanID(genoGDS),]
subjAnnot$sex <- as.factor(subjAnnot$sex)
subjAnnot$EV1 <- pca$eigenvect[match(subjAnnot$scanID, pca$sample.id), 1]

genoData <- GenotypeData(genoGDS, scanAnnot=subjAnnot)
chr <- getChromosome(genoData)
assoc.list <- lapply(unique(chr), function(x) {
    ## Y chromsome only includes males, cannot have sex as a covariate
    covar <- ifelse(x == 25, "EV1", c("sex", "EV1"))
    start <- which(chr == x)[1]
    assocRegression(genoData, outcome="status", covar=covar, model.type="logistic", 
                    snpStart=start, snpEnd=start+50)
})
assoc <- do.call(rbind, assoc.list)
close(genoData)


###################################################
### code chunk number 111: DataCleaning.Rnw:2412-2414
###################################################
qqPlot(pval=assoc$Wald.pval,
  truncate=TRUE, main="QQ Plot of Wald Test p-values")


###################################################
### code chunk number 112: DataCleaning.Rnw:2422-2425
###################################################
chrom <- getChromosome(snpAnnot, char=TRUE)
snp.sel <- getSnpID(snpAnnot) %in% assoc$snpID
manhattanPlot(assoc$Wald.pval, chromosome=chrom[snp.sel])


###################################################
### code chunk number 113: DataCleaning.Rnw:2436-2454
###################################################
# Identify SNPs with lowest p-values
snp <- pData(snpAnnot)[snp.sel, c("snpID", "rsID")]
snp$pval <- assoc$Wald.pval
snp <- snp[order(snp$pval),]
snp <- snp[1:9,]
xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
xyGDS <- GdsIntensityReader(xyfile)
xyData <- IntensityData(xyGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genoGDS <- GdsGenotypeReader(genofile)
genoData <- GenotypeData(genoGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
pdf(file="DataCleaning-cluster.pdf")
par(mfrow = c(3,3))
mtxt <- paste("SNP", snp$rsID, "\np =", format(snp$pval, digits=4))
genoClusterPlot(xyData, genoData, snpID=snp$snpID, main.txt=mtxt)
dev.off()
close(xyData)
close(genoData)


###################################################
### code chunk number 114: DataCleaning.Rnw:2461-2462
###################################################
unlink(subj.filt.file)

