#! /usr/bin/env Rscript
cat("Program path:", unlist(strsplit(grep(commandArgs(), pattern = "file=", value = T), split = "="))[2], "\n")
arguments <- commandArgs(trailingOnly = TRUE)
arguments
if(length(arguments) != 9) {stop("This program requires nine arguments. The first one is the destination folder, the second one is a semicolon-delimited CSV file with sample names and grouping variables, the third determines the sample ID variable, the fourth - grouping variable, the fifth - independent factor, the sixth - file with a gene list (one gene per line) to calculate the gene signature, the seventh is the number of threads to use, the eighth says if samples are paired, while the ninth determines the column with pair indicators.")}

csvfile <- arguments[2]
sampleid <- arguments[3]
grouping.cols <- unlist(strsplit(arguments[4], split = "\\,"))
indfactor <- arguments[5]
txt.file <- arguments[6]
threads <- as.numeric(arguments[7])
paired.samples <- arguments[8]
pair.ind <- arguments[9]

library(doMC)
library(foreach)
library(dplyr)

registerDoMC(threads)

csvtable.full <- read.table(file = csvfile, sep = ";", header = T, stringsAsFactors = F)
csvtable.full <- csvtable.full[!grepl(csvtable.full[[indfactor]], pattern = "^([[:space:]]?)+$"),]
csvtable.full <- csvtable.full[!is.na(csvtable.full[[indfactor]]),]
if(paired.samples == "TRUE") {csvtable.full <- csvtable.full[!is.na(csvtable.full[[pair.ind]]),]}
csvtable.full[[sampleid]] <- gsub(csvtable.full[[sampleid]], pattern = "#", replacement = ".")
rownames(csvtable.full) <- csvtable.full[[sampleid]]

save.image <- function(file){save(list=grep(ls(all.names = TRUE, envir = .GlobalEnv), pattern = "^txt.file$", value = T, invert = T), file = file)}

vep.r <- function(groupid) {

workdir <- paste0(arguments[1], "/Summary", "/Grouping_variables:", arguments[4], "/Factor:", indfactor)

runid <- sub(sub(workdir, pattern = "^.*RUNS/", replacement = ""), pattern = "/MAPPINGS.*$", replacement = "")
antype <- sub(sub(workdir, pattern = "\\/Summary\\/.*", replacement = ""), pattern = "^.*\\/", replacement = "")
impacts <- c("HIGH", "(HIGH|MODERATE)")

library(doMC)
library(foreach)
library(dplyr)
library(tidyr)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(pals)
library(openxlsx)
library(data.table)
library(pdftools)
library(ggfortify)
library(limma)

dir.create(workdir, recursive = T)
gene.signature <- gsub(txt.file, pattern = "^.*\\/", replacement = "")
setwd(workdir)
prefix <- paste(runid, antype, groupid, indfactor, paste("paired", paired.samples, sep = ":"), sep = ".")

save.image(paste(prefix, "VEP_analysis.RData", sep = "."))

if(! file.exists(paste(prefix, "VEP_analysis.final.RData", sep = "."))) {
unlink(paste(prefix, "VEP_analysis.xlsx", sep = "."))

registerDoMC(threads)

file.list <- list.files("../../..", pattern = ".*\\.csv$", full.names = T)
s.names <- sub(file.list, pattern = "(\\.\\.\\/\\.\\.\\/\\.\\.\\/)(.*)(\\.sorted(\\.|_).*\\csv$)", replacement = "\\2")
s.names <- gsub(s.names, pattern = "#", replacement = ".")

if(length(file.list) == 0) {stop("The file list is empty.")}

sel.columns <- c("Gene", "SYMBOL", "SYMBOL_SOURCE", "HGVSg", "HGVSc", "HGVSp", "Existing_variation", "Consequence", "SIFT", "PolyPhen", "MAX_AF", "CLIN_SIG", "PUBMED", "ZYG", "IMPACT")

if(!file.exists(paste0("../../", runid, ".", antype, ".df.full.RData"))) {

cat("Generating the first list of hits...\n")
df.list <- foreach(i = file.list) %dopar% {
file.tmp <- fread(file = i, sep = "\t", header = T, select = c(sel.columns))
file.tmp <- file.tmp[file.tmp[["SYMBOL_SOURCE"]] != "Clone_based_ensembl_gene" & grepl(file.tmp[["IMPACT"]], pattern = "(HIGH|MODERATE)"), , drop = F]
file.tmp[["MAX_AF"]][file.tmp[["MAX_AF"]] == "-"] <- -1
file.tmp}
names(df.list) <- s.names

cat("Generating the second list of hits...\n")
df.list <-foreach(s.name = names(df.list)) %dopar% {
file.tmp <- df.list[[s.name]]
if(nrow(file.tmp) > 0) {
foreach(position = seq(1, nrow(file.tmp)), .combine = rbind) %do% {
  if(grepl(file.tmp[["SYMBOL"]][position], pattern = "^-$")) {
    if(file.tmp[["Gene"]][position] == "-") {
    file.tmp[["SYMBOL"]][position] <- make.names(file.tmp[["HGVSg"]][position])} else {
    file.tmp[["SYMBOL"]][position] <- file.tmp[["Gene"]][position]}
    file.tmp[position,]} else {file.tmp[position,]}}
} else {file.tmp[1,]}
}
names(df.list) <- s.names

cat("Generating the third list of hits...\n")
df.full <- foreach(Sample.name = s.names, .combine = rbind) %dopar% {
  df.tmp <- data.frame(Sample.name, df.list[[Sample.name]])
  df.tmp <- subset(df.tmp, subset = !is.na(IMPACT))
   if(nrow(df.tmp) > 0) {
  if(sum(is.na(df.tmp)) > 0) {stop("The table contains missing results.")}
  
  hits <- unique(paste(df.tmp[["SYMBOL"]], df.tmp[["HGVSg"]], df.tmp[["IMPACT"]]))
  
  df.tmp <- foreach(hit = hits, .combine = rbind) %do% {
    hit.vec <- unlist(strsplit(hit, split = " "))
    hit.df <- df.tmp[df.tmp[["SYMBOL"]] == hit.vec[1] & df.tmp[["HGVSg"]] == hit.vec[2] & df.tmp[["IMPACT"]] == hit.vec[3], , drop = F]

    hit.sift <- hit.df[["SIFT"]]
    hit.sift <- sub(hit.sift, pattern = "\\([0-9.]+\\)", replacement = "")
    hit.sift <- names(table(hit.sift)[which.max(table(hit.sift))])
    hit.df[["SIFT"]] <- hit.sift
    
    hit.polyphen <- hit.df[["PolyPhen"]]
    hit.polyphen <- sub(hit.polyphen, pattern = "\\([0-9.]+\\)", replacement = "")
    hit.polyphen <- names(table(hit.polyphen)[which.max(table(hit.polyphen))])
    hit.df[["PolyPhen"]] <- hit.polyphen
    
    hit.clinsig <- hit.df[["CLIN_SIG"]]
    hit.clinsig <- names(table(hit.clinsig)[which.max(table(hit.clinsig))])
    hit.df[["CLIN_SIG"]] <- hit.clinsig
    
    hit.maxaf <- hit.df[["MAX_AF"]]
    hit.maxaf <- names(table(hit.maxaf)[which.max(table(hit.maxaf))])
    hit.df[["MAX_AF"]] <- hit.maxaf
    
    hit.df}
  
    df.tmp <- foreach(hit = hits, .combine = rbind) %do% {
    df.hit <- df.tmp[df.tmp[,c("SYMBOL","HGVSg","IMPACT"), drop = F] %>% unite(col = "United", sep = " ") == hit, , drop = F]
    df.hit[which.min(rowSums(df.hit == "-")), , drop = F]}
  
  df.tmp <- df.tmp[ !grepl(df.tmp[["SIFT"]], pattern = "(tolerated|-)") | !grepl(df.tmp[["PolyPhen"]], pattern = "(benign|-)") | !grepl(df.tmp[["CLIN_SIG"]], pattern = "(benign|-)") | as.numeric(df.tmp[["MAX_AF"]]) < 0.01, , drop = F]
  df.tmp[["MAX_AF"]][df.tmp[["MAX_AF"]] == "-1"] <- "-"
  df.tmp
}}
rm(df.list)

if(!is.null(df.full)) {
col.names <- c("Sample.name", sel.columns)
df.full <- as.matrix(df.full[with(df.full, order(Sample.name, SYMBOL, IMPACT)),])
colnames(df.full) <- col.names
}
save(list = c("df.full"), file = paste0("../../", runid, ".", antype, ".df.full.RData"))
} else {
  cat("Loading a pre-existing RData file:", paste0("../../", runid, ".", antype, ".df.full.RData...\n"))
  load(paste0("../../", runid, ".", antype, ".df.full.RData"))
}
if(!is.null(df.full)) {
  
  if(paired.samples == "TRUE") {csvtable[[sampleid]] <- paste(csvtable[[sampleid]], csvtable[[pair.ind]], sep = "#")}
  rownames(csvtable) <- csvtable[[sampleid]]
  if(length(csvtable[[sampleid]]) != length(csvtable[[indfactor]])) {stop("Some sample ids or group ids are missing.")}
  
  if(paired.samples == "TRUE")
  { anno <- csvtable %>% dplyr::select(indfactor, pair.ind)} else
  { anno <- csvtable %>% dplyr::select(indfactor)}
  anno <- as.data.frame(sapply(anno, factor))
  rownames(anno) <- rownames(csvtable)
  
  if(paired.samples == "TRUE") {
    colnames(anno) <- c(indfactor, pair.ind)
  } else {
    colnames(anno) <- c(indfactor)
  }
  
  anno.sub <- anno
  rownames(anno.sub) <- sub(rownames(anno.sub), pattern = "#.*$", replacement = "")
  df.full <- as.matrix(merge(x = df.full, y = anno.sub, by.x = "Sample.name", by.y = 0))
  if(paired.samples) {df.full <- df.full[, colnames(df.full)[c(1,17,18,2:16)]]} else {
    df.full <- df.full[, colnames(df.full)[c(1,17,2:16)]]}
  if(nrow(df.full) == 0) {df.full <- NULL}}

if(!is.null(df.full)) {

wb <- createWorkbook()

if(nrow(df.full) <= 2**20) {
  sheet.name <- substr(paste("All_mutations", sep = "."), 1, 31)
  if(length(wb$sheet_names) > 0) {if(sheet.name %in% wb$sheet_names) {removeWorksheet(wb, sheet = sheet.name)}}
  addWorksheet(wb = wb, sheetName = sheet.name)
  writeData(x = df.full, wb = wb, sheet = sheet.name) } else {
  All.mut.name <- paste("VEP_analysis", "all_mutations", prefix, "csv" , sep = ".")
  write.table(x = df.full, row.names = F, file = All.mut.name, sep = ";")
}

f.hm <- function(x) {df.tmp <-  x[grepl(x[["IMPACT"]], pattern = impact),]
if(nrow(df.tmp) > 0) {
  return(as.matrix(df.tmp[!duplicated(df.tmp$HGVSg),]))}}

for(impact in impacts) {
impact.name <- paste(sub(sub(unlist(strsplit(impact, split = "\\|")), pattern = "\\(", replacement = ""), pattern = "\\)", replacement = ""), collapse = "_or_")

mut.analysis <<- function(data, anno) {
  
  if(!all(rownames(data) == rownames(anno))) {stop(paste("Sample names:", paste(rownames(data), collapse = ", "), "and annotations:",  paste(rownames(anno), collapse = ", "), "do not match."))}
  if(paired.samples == TRUE) {rownames(data) <- anno.full.names[["full.sample.names"]]}
  data.bool <- data
  data.bool[data.bool == 0] <- FALSE
  data.bool[data.bool > 0] <- TRUE
  
  arg1.value.name <- as.character(sys.call())[2]
  
  data.bool.name <- paste(arg1.value.name, "bool", antype, groupid, make.names(impact), sep = ".")
  assign(data.bool.name, value = data.bool, envir = .GlobalEnv)
  
  if(paired.samples == "TRUE") {
    data.df <- setNames(data.frame(rowSums(data), anno), nm = c("Counts", indfactor, pair.ind))
    data.df <- data.df[order(data.df[[indfactor]], data.df[[pair.ind]]),] 
    data.df[[pair.ind]] <- as.factor(data.df[[pair.ind]]) } else {
      
      data.df <- setNames(data.frame(rowSums(data), anno), nm = c("Counts", indfactor))
      data.df <- data.df[order(data.df[[indfactor]]),] }
  
  plot.annotate <- function(plot) {if(length(unique(data.df[[indfactor]])) > 1) {
    kw.res <- with(data.df, kruskal.test(Counts ~ get(indfactor)))
    
    if(!is.na(kw.res$p.value)) {
    if(kw.res$p.value < 0.001) {kw.res$p.value <- formatC(kw.res$p.value, format = "e", digits = 3)} else {
      kw.res$p.value <- round(kw.res$p.value,4)}}
    
    kw.res.label <- paste(kw.res$method, "p-value =" , kw.res$p.value)
    kw.res.name <- paste("kw.res", antype, impact.name, sep = ".")
    assign(kw.res.name, value = kw.res.label)
    level.names <- levels(as.factor(as.character(data.df[[indfactor]])))
    for(i in seq(1, ncol(combn(level.names, m = 2)))) {sel.levels <- combn(level.names, m = 2)[,i]
    sub.df <- subset(data.df, subset = get(indfactor) %in% sel.levels)
    if(paired.samples == "TRUE") {
      mw.res <- with(sub.df, wilcox.test(Counts ~ get(indfactor), paired = TRUE))} else {
      mw.res <- with(sub.df, wilcox.test(Counts ~ get(indfactor), paired = FALSE))  
      }
    
    if(!is.na(mw.res$p.value)) {
    if(mw.res$p.value < 0.001) {mw.res$p.value <- formatC(mw.res$p.value, format = "e", digits = 3)} else {
      mw.res$p.value <- round(mw.res$p.value,4)}}
    
    sel.levels.name <- paste(sel.levels, collapse = " vs ")
    mw.res.label <- paste0(sel.levels.name, "; ", mw.res$method, " p-value = ", mw.res$p.value)
    mw.res.name <- paste("mw.res", antype, impact.name, paste(sel.levels, collapse = "_"), sep = ".")
    assign(mw.res.name, value = mw.res.label)
    }
    
    top.pos <- 1+(ncol(combn(level.names, m = 2))+1)/20
    plot <- plot + annotate("text", x=(as.numeric(if(plot$labels$x != "Sample IDs") {length(levels(data.df[[plot$labels$x]]))} else {nrow(data.df)})+1)/2, y = max(data.df$Counts)*top.pos, label= get(ls(pattern = paste("kw.res", antype, sep = "."))))
    
    for(i in ls(pattern = paste0("mw.res\\.(SNP|NON_SNP)\\.", impact.name, "\\."))) {
      top.pos <- top.pos - 0.05
      plot <- plot + annotate("text", x=(as.numeric(if(plot$labels$x != "Sample IDs") {length(levels(data.df[[plot$labels$x]]))} else {nrow(data.df)})+1)/2, y = max(data.df$Counts)*top.pos, label= get(i))}
  
    if(class(plot$layers[[1]]$geom)[1] == "GeomBoxplot") {
      if(all(names(table(plot$data[[indfactor]])) == levels(plot$data[[indfactor]]))){plot <- plot + annotate("text", y = max(plot$data[["Counts"]]) * 0.9, x = (length(levels(plot$data[[indfactor]]))+1)/2, label = paste("Group", paste(names(table(plot$data[[indfactor]])), as.vector(table(plot$data[[indfactor]])), sep = ": n. obs. = "), collapse = "\n"))}
    }
  }
    try(print(plot))}
  
  if(length(rownames(data.df))>56) {width <- length(rownames(data.df))/8} else {width <- 7}
  pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis barplot-samples.pdf", sep = "."), width = width)
  
  if(paired.samples == "TRUE") {
        rnames <- factor(rownames(data.df), levels = rownames(data.df))
        ggplot1 <- ggplot(data.df, aes(x = rnames, y = Counts, fill = get(indfactor))) + 
          geom_bar(stat = "identity") + 
          labs(title = paste0("Frequency of ", antype, " alterations ", "(", "impact: ", impact.name, ")"), x = "Sample IDs", fill = indfactor) +
          theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          scale_fill_manual(values = as.vector(cols25()))
        } else {
        ggplot1 <- ggplot(data.df, aes(x = reorder(rownames(data.df), as.numeric(get(indfactor))), y = Counts, fill = get(indfactor))) + 
          geom_bar(stat = "identity") + 
          labs(title = paste0("Frequency of ", antype, " alterations ", "(", "impact: ", impact.name, ")"), x = "Sample IDs", fill = indfactor) +
          theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          scale_fill_manual(values = as.vector(cols25()))}
        
  plot.annotate(plot = ggplot1)
  dev.off()
  
  pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis boxplot-samples.pdf", sep = "."), width = 10)
        ggplot1.2 <- ggplot(data.df, aes(x = get(indfactor), y = Counts, fill = get(indfactor))) + 
          geom_boxplot() + 
          labs(title = paste0("Frequency of ", antype, " alterations ", "(", "impact: ", impact.name, ")"), x = indfactor, fill = indfactor) +
          theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          scale_fill_manual(values = as.vector(cols25()))
  
  plot.annotate(plot = ggplot1.2)
  dev.off()
  
  if(ncol(data) > 0) {
    categories <- unique(anno[[indfactor]])
    all.groups <- foreach(i = categories, .combine = rbind) %dopar% {
      group.tmp <- rownames(anno)[anno[[indfactor]] == i]
      data.sub.tmp <- data[sub(rownames(data), pattern = "#.*$", replacement = "") %in% group.tmp, , drop = F]
      countSums <- colSums(data.sub.tmp)
      countMeans <- colMeans(data.sub.tmp)
      countSds <- colSds(as.matrix(data.sub.tmp))
      countSE <- countSds/sqrt(length(countSds))
      t.dist <- qt((1-0.05)/2 + 0.5, nrow(data.sub.tmp)-1) # Value of the Student's t distribution for alpha = 0.05, tends to be 1.96 if the sample size is big enough.
      count95CI <- countSE * t.dist
      data.frame(countSums, countMeans, countSds, countSE, count95CI, rep(i), colnames(data.sub.tmp))
    }
    
    all.groups <- setNames(all.groups, nm = c("Count_sum", "Count_mean", "Count_SD", "Count_SE", "Count_95CI", indfactor, "Gene_name"))
    all.groups <- all.groups[,c(length(all.groups),length(all.groups)-1, seq(1,length(all.groups)-2))]
    all.groups <- all.groups[all.groups[["Count_sum"]] > 0,]
    all.groups <- all.groups[order(all.groups[[indfactor]]),]
    
    if(nrow(all.groups)>70) {width <- nrow(all.groups)/10} else {width <- 7}
    pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis barplot-genes.pdf", sep = "."), width = width)
    ggplot2 <- ggplot(all.groups, aes(x = Gene_name, y = Count_mean, fill = get(indfactor))) + 
      geom_bar(stat = "identity", position=position_dodge()) + 
      geom_errorbar(aes(y = Count_mean, ymin=Count_mean-Count_95CI, ymax=Count_mean+Count_95CI), position = position_dodge(), color = "darkgrey") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      labs(title = paste0("Frequency of ", antype, " alterations with 95% CI ", "(", "impact: ", impact.name, ")"), fill = indfactor) +
      scale_fill_manual(values = as.vector(cols25()))
    try(print(ggplot2))
    dev.off()
    
    if(! exists("wb")) {
      wb <- createWorkbook()
    }
    sheet.name <- substr(paste("Mut.per.gene.impact",impact.name, sep = "."), 1, 31)
    if(length(wb$sheet_names) > 0) {if(sheet.name %in% wb$sheet_names) {removeWorksheet(wb, sheet = sheet.name)}}
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(x = all.groups, wb = wb, sheet = sheet.name, rowNames = F)
    
    muts.per.gene.name <- paste("muts.per.gene", antype, make.names(impact), sep = ".")
    assign(muts.per.gene.name, value = all.groups)
    
    data.t <- t(data)
    data.t <- data.t[rowSums(data.t) > 0, , drop = F]
    
    if(length(rownames(data.t))>21) {height <- ceiling(length(rownames(data.t))/3)} else {height <- 7}
    if(length(colnames(data.t))>21) {width <- ceiling(length(colnames(data.t))/3)} else {width <- 7}
    
    palette <- rev(stepped())
    palette <- c("#DDDDDD",palette[seq(2,24,2)])
    
    colnames(data.t) <- sub(colnames(data.t), pattern = "#.*$", replacement = "")
    data.t <- data.t[,order(colnames(data.t)), drop = F]
    anno <- anno[order(rownames(anno)), , drop = F]
    if(!all(colnames(data.t) == rownames(anno))) {stop("Annotations do not match sample names.")}
    data.t <- data.t[,order(anno[[indfactor]]), drop = F]
    rownames(data) <- sub(rownames(data), pattern = "#.*$", replacement = "")
    
    pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis.heatmap.pdf", sep = "."), width = width, height = height)
    heatmap <- pheatmap(data.t, cluster_cols = F, cluster_rows = F, color =  palette, annotation_col = anno, cellwidth = 10, cellheight = 10, main = paste0("Heatmap of genetic alterations-", "impact: ", impact.name,"-" ,antype), breaks = c(-1:12), legend_breaks = c(-1:12), legend_labels = append("", append(c(0:11), ">=12")))
    dev.off()
    
    pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis.sampleDist.heatmap.pdf", sep = "."), width = width, height = width)
    data.log2 <- log2(data + 1)
    sampleDists <- dist(data.log2)
    sampleDistsMatrix <- as.matrix(sampleDists)
    tryCatch(expr = {pheatmap(sampleDistsMatrix, 
                              clustering_distance_rows = sampleDists,
                              clustering_distance_cols = sampleDists,
                              annotation_col = anno,
                              cellwidth = 10, cellheight = 10,
                              main = paste0("Sample distances for genetic alterations-", "impact: ", impact.name,"-" ,antype))},
             error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(paste(main="ERROR: All values are identical", "impact", impact.name, sep = "."))}})
    dev.off()
    
    data.prcomp <- prcomp(data.log2)
    ind.length <- length(levels(anno[[indfactor]]))
    pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis.PCA.plot.pdf", sep = "."), width = 10, height = 10)
    if(ncol(data.prcomp$x) > 1) {
      ap1 <- autoplot(data.prcomp, data = anno, size = 3, colour = indfactor) + 
        scale_color_manual(values = if (ind.length > 25) {rainbow(ind.length)} else {as.vector(cols25())}) + coord_fixed() + 
        labs(title = paste("PCA plot", runid, antype, paste0("impact: ", impact.name), sep = "-")) + 
        theme(plot.title = element_text(hjust = 0.5))
      tryCatch(expr = print(ap1),
               error = {function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(paste(main="ERROR: Too low data diversity to draw a PCA plot", "impact", impact.name, sep = "."))}})} else {
                 plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("ERROR: Too few dimensions to draw a PCA plot", "impact", impact.name, sep = "."))
               }
    dev.off()
    
    sheet.name <- substr(paste("Heatmap", "impact", impact.name, sep = "."), 1, 31)
    if(length(wb$sheet_names) > 0) {if(sheet.name %in% wb$sheet_names) {removeWorksheet(wb, sheet = sheet.name)}}
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(x = data.t, wb = wb, sheet = sheet.name, rowNames = T)}
}

New.muts.mx <- df.full[df.full[,"Existing_variation"] == "-", , drop = F]
New.muts.mx <- New.muts.mx[grepl(New.muts.mx[,"IMPACT", drop = F], pattern = impact), , drop = F]

if(nrow(New.muts.mx) <= 2**20) {
sheet.name <- substr(paste("New_mutations", "impact", impact.name, sep = "."), 1, 31)
if(length(wb$sheet_names) > 0) {if(sheet.name %in% wb$sheet_names) {removeWorksheet(wb, sheet = sheet.name)}}
addWorksheet(wb = wb, sheetName = sheet.name)
writeData(x = New.muts.mx, wb = wb, sheet = sheet.name) } else {
New.mut.name <- paste("VEP_analysis", "new_mutations", impact.name, prefix, "csv" , sep = ".")
write.table(x = New.muts.mx, row.names = F, file = New.mut.name, sep = ";")
}

if(paired.samples == "TRUE")
{ anno <- csvtable %>% dplyr::select(indfactor, pair.ind)} else
{ anno <- csvtable %>% dplyr::select(indfactor)}
anno <- as.data.frame(sapply(anno, factor))
rownames(anno) <- rownames(csvtable)

if(paired.samples == "TRUE") {
  colnames(anno) <- c(indfactor, pair.ind)
} else {
  colnames(anno) <- c(indfactor)
}

df.hm.list <- by(data = df.full, INDICES = list(df.full[,"Sample.name"], df.full[,"SYMBOL"]), FUN=f.hm)
names(df.hm.list) <- levels(interaction(df.full[,"Sample.name"], df.full[,"SYMBOL"], sep = "$$"))
df.hm.list <- df.hm.list[sapply(df.hm.list, Negate(is.null))]

geneList <- sort(unique(sub(names(df.hm.list), pattern = "^.*\\$\\$", replacement = "")))

sampleIDs <- s.names
sampleIDs.fin <- NULL
for(i in seq(1,length(sampleIDs))) {
  if(sum(grepl(rownames(anno), pattern = paste0("^", sampleIDs[i], "(#.*)?", "$"))) > 1) {stop(paste(Sample.name, "sample cannot be uniquely identified in the annotation table."))}
  sampleIDs.fin <- append(sampleIDs.fin, rownames(anno)[grep(rownames(anno), pattern = paste0("^", sampleIDs[i], "(#.*)?", "$"))])}
sampleIDs <- sub(sampleIDs.fin, pattern = "#.*$", replacement = "")

res1 <- foreach(gene = geneList, .combine = cbind) %dopar% {foreach(sampleID = sampleIDs, .combine = c) %do% {if(! is.null(df.hm.list[[paste(sampleID, gene, sep = "$$")]])) {nrow(df.hm.list[[paste(sampleID, gene, sep = "$$")]])} else {0}} }
rm(df.hm.list)

res1 <- as.data.frame(res1)
if(nrow(res1) != 0 & ncol(res1) != 0){
rownames(res1) <- sampleIDs
colnames(res1) <- geneList} else 
if(nrow(res1) == 0 & ncol(res1) == 0){
res1 <- data.frame(sampleIDs)
rownames(res1) <- sampleIDs
res1 <- res1[,-1]
} else {stop("The result data frame is invalid.")}

anno.short.sample.names <- anno
anno.short.sample.names[["full.sample.names"]] <- rownames(anno.short.sample.names)
rownames(anno.short.sample.names) <- sub(rownames(anno), pattern = "#.*$", replacement = "")
res1.merged <- merge(x = res1, y = anno.short.sample.names, by.x = 0, by.y = 0)
rownames(res1.merged) <- res1.merged[["Row.names"]]
if(paired.samples == "TRUE") {
  res1.merged <- res1.merged[order(res1.merged[[indfactor]], res1.merged[[pair.ind]]),]} else {
  res1.merged <- res1.merged[order(res1.merged[[indfactor]]),]}

if(paired.samples == "TRUE") {
  res1 <- res1.merged %>% dplyr::select(-c("Row.names", indfactor, pair.ind, "full.sample.names"))
  anno.full.names <- res1.merged %>% dplyr::select(indfactor, pair.ind, "full.sample.names")
  anno <- res1.merged %>% dplyr::select(indfactor, pair.ind)
} else {
  res1 <- res1.merged %>% dplyr::select(-c("Row.names", indfactor, "full.sample.names"))
  anno <- res1.merged %>% dplyr::select(indfactor) }

res1.name <- paste("res1", antype, groupid, make.names(impact), sep = ".")
assign(res1.name, value = res1, envir = .GlobalEnv)

anno.name <- paste("anno", antype, groupid, make.names(impact), sep = ".")
assign(anno.name, value = anno, envir = .GlobalEnv)

mut.analysis(data = res1, anno = anno)
}

pdf_combine(rev(list.files(pattern = "VEP_analysis.*\\.pdf$")), output = paste("VEP_analyses", "all_genes", prefix, "pdf", sep = "."))
unlink(rev(list.files(pattern = "VEP_analysis.*\\.pdf$")))

saveWorkbook(wb, file = paste(prefix, "VEP_analysis.xlsx", sep = "."), overwrite = T)
rm(wb)
} else {
  res1 <- data.frame(s.names)
  rownames(res1) <- s.names
  res1 <- res1[,-1]
  
  for(impact in impacts) {
  res1.name <- paste("res1", antype, groupid, make.names(impact), sep = ".")
  assign(res1.name, value = res1, envir = .GlobalEnv)
  
  res1.bool.name <- paste("res1.bool", antype, groupid, make.names(impact), sep = ".")
  assign(res1.bool.name, value = res1, envir = .GlobalEnv)
  
  anno <- anno[rownames(anno) == s.names, , drop = F]
  anno.name <- paste("anno", antype, groupid, make.names(impact), sep = ".")
  assign(anno.name, value = anno, envir = .GlobalEnv)
  }
  cat(paste("There are no", antype, "mutations (alterations classified as high or moderate in the Ensembl database) in the current data set.\n\nSample IDs:", paste(s.names, collapse = ", "), "\n"), file = paste("VEP_analyses", "all_genes", prefix, "txt", sep = "."))
}
save.image(file = paste(prefix, "VEP_analysis.final.RData", sep = "."))
} else {
  cat("Loading a pre-existing RData file:", paste(prefix, "VEP_analysis.final.RData...\n", sep = "."))
  load(paste(prefix, "VEP_analysis.final.RData", sep = "."))}

library(doMC)
library(foreach)
library(dplyr)
library(tidyr)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(pals)
library(openxlsx)
library(data.table)
library(pdftools)
library(ggfortify)
library(limma)

registerDoMC(threads)

if(txt.file != "NA") {
  con <- file(txt.file)
  gene.list1 <- readLines(con)
  close(con)

  gene.list1.conv <- NULL
  for(i in gene.list1) {if(length(alias2Symbol(i, species = "Hs")) == 0) {gene.list1.conv <- append(gene.list1.conv, values = i)} else {gene.list1.conv <- append(gene.list1.conv, values = alias2Symbol(i, species = "Hs"))}}

  wb <- createWorkbook()
  for(impact in impacts) {
      impact.name <- paste(sub(sub(unlist(strsplit(impact, split = "\\|")), pattern = "\\(", replacement = ""), pattern = "\\)", replacement = ""), collapse = "_or_")
      res1.name <- paste("res1", antype, groupid, make.names(impact), sep = ".")
      res1 <- get(res1.name)
      res1.sel <- res1 %>% dplyr::select(colnames(res1)[colnames(res1) %in% gene.list1.conv])
      anno.name <- paste("anno", antype, groupid, make.names(impact), sep = ".")
      anno <- get(anno.name)
      
      if(ncol(res1.sel) > 1 & nrow(res1.sel) > 0) {
      mut.analysis(data = res1.sel, anno = anno)

      res1.sel.t <- t(res1.sel)
      res1.sel.t.medians <- rowMedians(res1.sel.t)
      names(res1.sel.t.medians) <- rownames(res1.sel.t)
      res1.sel.t.cat <- res1.sel.t
      for(i in seq(1,length(res1.sel.t.medians))) {j <- res1.sel.t.medians[i]
      res1.sel.t.cat[names(j),][res1.sel.t[names(j),] > j] <- 1
      res1.sel.t.cat[names(j),][res1.sel.t[names(j),] <= j] <- 0
      }
      res1.sel.t.cat.t <- t(res1.sel.t.cat)

      if(all(rownames(res1.sel.t.cat.t) == rownames(anno)) == TRUE) {
        res1.sel.t.cat.t.by <- by(data = res1.sel.t.cat.t, INDICES = anno[[indfactor]], FUN = colMeans)} else
        {stop("Annotations do not match sample names.")}
      res1.sel.t.cat.t.by <- res1.sel.t.cat.t.by[sapply(res1.sel.t.cat.t.by, Negate(is.null))]

      res1.sel.t.cat.t.comb <- foreach(name = names(res1.sel.t.cat.t.by), .combine = "cbind") %do% {res1.sel.t.cat.t.by[[name]]}
      res1.sel.t.cat.t.comb <- as.data.frame(res1.sel.t.cat.t.comb)
      colnames(res1.sel.t.cat.t.comb) <- names(res1.sel.t.cat.t.by)

      if(ncol(res1.sel.t.cat.t.comb)>21) {width <- ceiling(ncol(res1.sel.t.cat.t.comb)/3)} else {width <- 7}
      if(nrow(res1.sel.t)>21) {height <- ceiling(nrow(res1.sel.t)/3)} else {height <- 7}

      pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis.gene signature comparison heatmap (median expression-categorized.pdf", sep = "."), height = height, width = width)
      try(heatmap1 <- pheatmap(res1.sel.t.cat.t.comb, main = paste(gene.signature, "gene signature (median expression-categorized", paste("impact:", impact.name), indfactor, "grouped)", sep = "-"), cellwidth = 10, cellheight = 10))
      dev.off()

      res1.sel.t.cat.t.comb[res1.sel.t.cat.t.comb >= 0.5] <- 1
      res1.sel.t.cat.t.comb[res1.sel.t.cat.t.comb < 0.5] <- 0
      res1.sel.t.cat.t.comb.df <- as.data.frame(res1.sel.t.cat.t.comb)

      pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis.gene signature comparison heatmap (median expression-categorized, binarized.pdf", sep = "."), height = height, width = width)
      tryCatch(expr = {heatmap2 <- pheatmap(res1.sel.t.cat.t.comb, main = paste(gene.signature, "gene signature (median expression-categorized, binarized", paste("impact:", impact.name), indfactor, "grouped)", sep = "-"), cellwidth = 10, cellheight = 10)},
               error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(paste(main="ERROR: After binarization all values are identical", "impact", impact.name, sep = "."))}})
      dev.off()

      res1.sel.t.cat.t.comb.df.bool <- res1.sel.t.cat.t.comb.df == res1.sel.t.cat.t.comb.df[[sort(colnames(res1.sel.t.cat.t.comb.df))[1]]]
      res1.sel.t.cat.t.comb.df.bool.t <- t(res1.sel.t.cat.t.comb.df.bool)

      f_fraction <- function(x) {
        return(paste(rowSums(x), "matching out of", length(x), "genes", "(", rowSums(x)/length(x), ")"))}

      fin.summary <- by(res1.sel.t.cat.t.comb.df.bool.t, INDICES = rownames(res1.sel.t.cat.t.comb.df.bool.t), FUN = f_fraction)
      fin.summary <- as.list(fin.summary)
      names(fin.summary) <- gsub(names(fin.summary), pattern = "(.*)", replacement = paste(indfactor, "\\1", sep = ":"))
      sink(file =  paste(runid, antype, groupid, "impact", impact.name, paste(gene.signature, "gene.signature", sep = ":"), "comparison.summary.txt", sep = "."))
      print(fin.summary)
      sink()
      }}
  pdffiles <- rev(list.files(pattern = paste(runid, antype, "impact", "*", "pdf", sep = ".")))
  if(length(pdffiles) > 0) {
  pdf_combine(pdffiles, output = paste("VEP_analyses", paste0("gene.signature:", gene.signature), prefix, "pdf", sep = "."))
  unlink(pdffiles)} else {
  pdf(file = paste("VEP_analyses", paste0("gene.signature:", gene.signature), prefix, "pdf", sep = "."), width = 18)
    plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(paste(main="INFO: At least two genes contained in the", gene.signature, "signature have to harbor", antype, "mutations (impact high or moderate) for the analysis to succeed."))
    dev.off()
    }
  }
}

if(any(grouping.cols != "ALL_SAMPLES")) {
  if(length(grouping.cols) == 1) {
    groupids <- unique(csvtable.full[[grouping.cols]])} else if(length(grouping.cols) == 2) {
    col.pairs <- unique(csvtable.full %>% dplyr::select(grouping.cols))
    groupids <- foreach(i = seq(1,nrow(col.pairs)), .combine = c) %dopar% {paste(as.character(col.pairs[i,]), collapse = "+")}
    } else {stop("This program supports up to two grouping variables only.")}
  
  for(groupid in groupids) {
    if(grepl(groupid, pattern = "\\+")) {col.values <- unlist(strsplit(groupid, split = "\\+"))
    csvtable <- csvtable.full[csvtable.full[[grouping.cols[1]]] == col.values[1] & csvtable.full[[grouping.cols[2]]] == col.values[2], , drop = F]
    } else {
    csvtable <- csvtable.full[csvtable.full[[grouping.cols]] == groupid, , drop = F]}
    if(nrow(csvtable) > 0) {
    vep.r(groupid = groupid)}
  }
} else {
  csvtable <- csvtable.full
  if(nrow(csvtable) > 0) {
  vep.r(groupid = "ALL_SAMPLES")}
}

sessionInfo()
proc.time()

cat("All done.\n")
