#! /usr/bin/env Rscript
cat("Program path:", unlist(strsplit(grep(commandArgs(), pattern = "file=", value = T), split = "="))[2], "\n")

arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) != 10) {stop("This program requires ten arguments.
                                  The first one is the destination folder,
                                  the second one is a semicolon-delimited CSV file with sample names and grouping variables,
                                  the third determines the sample ID variable,
                                  the fourth - grouping variable,
                                  the fifth - independent factor,
                                  the sixth - file with a gene list (one gene per line) to calculate the gene signature,
                                  the seventh is the number of threads to use,
                                  the eighth says if samples are paired,
                                  the ninth determines the column with pair indicators,
                                  while the tenth is a logical indicator determining if the FDR correction should be applied in the TOST analysis of equivalence.")}
arguments
arguments.backup <- arguments

read.args <- function() {
csvfile <<- arguments[2]
sampleid <<- arguments[3]
grouping.cols <<- unlist(strsplit(arguments[4], split = "\\,"))
indfactor <<- arguments[5]
txt.file <<- arguments[6]
threads <<- as.numeric(arguments[7])
paired.samples <<- as.logical(arguments[8])
pair.ident <<- arguments[9]
FDR <<- as.logical(arguments[10])
}
read.args()

genome.ver <- grep(unlist(strsplit(arguments[1], split = "/")), pattern = "\\.genome", value = T)

if(grepl(arguments[1], pattern = "BAMs_wo_dups"))
{dups <- FALSE; csv.pat <- "_sorted\\.no_dups_VEP\\.(SNP|NON_SNP)\\.csv$"} else if(grepl(arguments[1], pattern = "BAMs_w_dups")) 
{dups <- TRUE; csv.pat <- "_sorted_VEP\\.(SNP|NON_SNP)\\.csv$"} else {stop("The duplication status cannot be determined.")}
if(dups) {BAM.type <- "BAM files with duplicates"} else {BAM.type <- "BAM files without duplicates"}

if(grepl(genome.ver, pattern = "hg38|GRCh38")) {genome <- "hg38"} else
if(grepl(genome.ver, pattern = "hg19|GRCh37")) {genome <- "hg19"} else 
  {stop("The genome version is invalid.")}

library(doMC)
library(foreach)
library(dplyr)
library(data.table)
library(quantsmooth)
library(stringi)
library(Hmisc)
library(TOSTER)
library(genefilter)
library(tidyr)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(pals)
library(openxlsx)
library(pdftools)
library(ggfortify)
library(limma)

registerDoMC(threads)

csvtable.full <- fread(file = csvfile, sep = ";", header = T, stringsAsFactors = F, check.names = T)
csvtable.full <- csvtable.full[!grepl(csvtable.full[[indfactor]], pattern = "^([[:space:]]?)+$"),]
csvtable.full <- csvtable.full[!is.na(csvtable.full[[indfactor]]),]
if(paired.samples) {csvtable.full <- csvtable.full[!is.na(csvtable.full[[pair.ident]]),]}

csvtable.full[[sampleid]] <- gsub(csvtable.full[[sampleid]], pattern = "#", replacement = ".")
csvtable.full[[sampleid]] <- sub(csvtable.full[[sampleid]], pattern = "^ *(.*) *$", replacement = "\\1")
rownames(csvtable.full) <- csvtable.full[[sampleid]]
csvtable.full[[indfactor]] <- as.factor(gsub(csvtable.full[[indfactor]], pattern = "^ *(.*) *$", replacement = "\\1"))

save.image <- function(file){save(list=grep(ls(all.names = TRUE, envir = .GlobalEnv), pattern = "^arguments$", value = T, invert = T), file = file)}

vep.r <- function(groupid) {

var.analysis <- function(data, anno, gene.signature, gene.list1.conv = NA, wb) {
    
    if(!all(rownames(data) == sub(rownames(anno), pattern = "#.*$", replacement = ""))) {stop(paste("Sample names:", paste(rownames(data), collapse = ", "), "and annotations:",  paste(rownames(anno), collapse = ", "), "do not match."))}
    rownames(data) <- rownames(anno)
    data.bool <- data
    data.bool[data.bool == 0] <- FALSE
    data.bool[data.bool > 0] <- TRUE
    
    arg1.value.name <- as.character(sys.call())[2]
    
    data.bool.name <- paste(arg1.value.name, "bool", antype, groupid, make.names(impact), sep = ".")
    assign(data.bool.name, value = data.bool, envir = .GlobalEnv)
    
    if(paired.samples) {
      data.df <- setNames(data.frame(rowSums(data), anno), nm = c("Counts", indfactor, pair.ident))
      data.df <- data.df[order(data.df[[indfactor]], data.df[[pair.ident]]),] 
      data.df[[pair.ident]] <- as.factor(data.df[[pair.ident]]) } else {
        
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
      if(paired.samples) {
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
      
      if(class(plot$layers[[1]]$geom)[1] == "GeomBoxplot") {
      plot <- plot + annotate("text", x = (length(levels(plot$data[[indfactor]]))+1)/2, y = max(data.df$Counts)*top.pos, label= get(ls(pattern = paste("kw.res", antype, sep = "."))))
      
      for(i in ls(pattern = paste0("mw.res\\.(SNP|NON_SNP)\\.", impact.name, "\\."))) {
        top.pos <- top.pos - 0.05
        plot <- plot + annotate("text", x = (length(levels(plot$data[[indfactor]]))+1)/2, y = max(data.df$Counts)*top.pos, label= get(i))}

        if(all(names(table(plot$data[[indfactor]])) == levels(plot$data[[indfactor]]))){plot <- plot + annotate("text", y = max(plot$data[["Counts"]]) * 0.9, x = (length(levels(plot$data[[indfactor]]))+1)/2, label = paste("Group", paste(names(table(plot$data[[indfactor]])), as.vector(table(plot$data[[indfactor]])), sep = ": n. obs. = "), collapse = "\n"))}
      } else {
        plot <- plot + annotate("text", x=(as.numeric(if(plot$labels$x != "Sample IDs") {length(levels(data.df[[plot$labels$x]]))} else {nrow(data.df)})+1)/2, y = max(data.df$Counts)*top.pos, label= get(ls(pattern = paste("kw.res", antype, sep = "."))))
        
        for(i in ls(pattern = paste0("mw.res\\.(SNP|NON_SNP)\\.", impact.name, "\\."))) {
          top.pos <- top.pos - 0.05
          plot <- plot + annotate("text", x=(as.numeric(if(plot$labels$x != "Sample IDs") {length(levels(data.df[[plot$labels$x]]))} else {nrow(data.df)})+1)/2, y = max(data.df$Counts)*top.pos, label= get(i))}
      }
    }
      try(print(plot))}
    
    if(length(rownames(data.df))>56) {width <- length(rownames(data.df))/8} else {width <- 7}
    pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis barplot-samples.pdf", sep = "."), width = width)
    
    if(paired.samples) {
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
      labs(title = paste0("Frequency of ", antype, " alterations ", "(", "impact: ", impact.name, ")"), x = "", fill = indfactor) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      scale_fill_manual(values = as.vector(cols25())) +
      geom_jitter(color = "azure4", size = 1) +
      stat_summary(geom = "point", fun = mean, fill = "yellow", shape = 24)
    
    plot.annotate(plot = ggplot1.2)
    dev.off()
    
    if(ncol(data) > 0) {
      categories <- unique(anno[[indfactor]])
      all.groups <- foreach(i = categories, .combine = rbind) %dopar% {
        group.tmp <- rownames(anno)[anno[[indfactor]] == i]
        data.sub.tmp <- data[sub(rownames(data), pattern = "#.*$", replacement = "") %in% sub(group.tmp, pattern = "#.*$", replacement = ""), , drop = F]
        countSums <- colSums(data.sub.tmp)
        countMeans <- colMeans(data.sub.tmp)
        countSds <- colSds(as.matrix(data.sub.tmp))
        countSE <- countSds/sqrt(length(countSds))
        t.dist <- qt((1-0.05)/2 + 0.5, nrow(data.sub.tmp)-1) # Value of the Student's t distribution for alpha = 0.05, tends to be 1.96 if the sample size is big enough.
        count95CI <- countSE * t.dist
        data.frame(countSums, countMeans, countSds, countSE, count95CI, rep(i), colnames(data.sub.tmp))
      }
      
      all.groups <- setNames(all.groups, nm = c("Count_sum", "Count_mean", "Count_SD", "Count_SE", "Count_0.95_CI", indfactor, "Gene_name"))
      all.groups <- all.groups[,c(length(all.groups),length(all.groups)-1, seq(1,length(all.groups)-2))]
      all.groups <- all.groups[all.groups[["Count_sum"]] > 0,]
      all.groups <- all.groups[order(all.groups[[indfactor]]),]
      
      if(nrow(all.groups)>70) {width <- nrow(all.groups)/10} else {width <- 7}
      pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis barplot-genes.pdf", sep = "."), width = width)
      ggplot2 <- ggplot(all.groups, aes(x = Gene_name, y = Count_mean, fill = get(indfactor))) + 
        geom_bar(stat = "identity", position=position_dodge()) + 
        geom_errorbar(aes(y = Count_mean, ymin=Count_mean-Count_0.95_CI, ymax=Count_mean+Count_0.95_CI), position = position_dodge(), color = "darkgrey") +
        theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        labs(title = paste0("Frequency of ", antype, " alterations with 95% CI ", "(", "impact: ", impact.name, ")"), fill = indfactor) +
        scale_fill_manual(values = as.vector(cols25()))
      try(print(ggplot2))
      dev.off()
      
      sheet.name <- substr(paste("Var.per.gene.impact",impact.name, sep = "."), 1, 31)
      if(length(wb$sheet_names) > 0) {if(sheet.name %in% wb$sheet_names) {removeWorksheet(wb, sheet = sheet.name)}}
      addWorksheet(wb = wb, sheetName = sheet.name)
      writeData(x = all.groups, wb = wb, sheet = sheet.name, rowNames = F)
      
      vars.per.gene.name <- paste("vars.per.gene", antype, make.names(impact), sep = ".")
      assign(vars.per.gene.name, value = all.groups)
      
      data.t <- t(data)
      data.t <- data.t[rowSums(data.t) > 0, , drop = F]
      
      if(length(rownames(data.t))>21) {height <- ceiling(length(rownames(data.t))/3)} else {height <- 7}
      if(length(colnames(data.t))>21) {width <- ceiling(length(colnames(data.t))/3)} else {width <- 7}
      
      palette <- rev(stepped())
      palette <- c("#DDDDDD",palette[seq(2,24,2)])
      
      colnames(data.t) <- sub(colnames(data.t), pattern = "#.*$", replacement = "")
      data.t <- data.t[,order(colnames(data.t)), drop = F]
      anno <- anno[order(rownames(anno)), , drop = F]
      
      if(!all(colnames(data.t) == sub(rownames(anno), pattern = "#.*$", replacement = ""))) {stop("Annotations do not match sample names.")}
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
      anno.prcomp <- anno[match(rownames(data.prcomp$x), sub(rownames(anno), pattern = "#.*$", replacement = "")), , drop = F]
      anno.prcomp[[sampleid]] <- rownames(anno.prcomp)
      if(nrow(anno.prcomp) > 25) {
        ind.length <- length(levels(anno.prcomp[[indfactor]]))} else {
          ind.length <- length(anno.prcomp[[sampleid]])
        }
      
      pdf(paste(runid, antype, "impact", impact.name, "VEP_analysis.PCA.plot.pdf", sep = "."), width = 10, height = 10)
      if(ncol(data.prcomp$x) > 1) {
        ap1 <- autoplot(data.prcomp, data = anno.prcomp, size = 3, colour = if(nrow(anno.prcomp) > 25) {indfactor} else {sampleid}, shape = if(nrow(anno.prcomp) > 25) {pty=20} else {indfactor}) + 
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
      writeData(x = data.t, wb = wb, sheet = sheet.name, rowNames = T)
    }
    wb.name <- as.character(sys.call()[["wb"]])
    assign(x = wb.name, value = wb, envir = .GlobalEnv)
}
f_TOST <- function(data, anno, gene.signature) {
  if(ncol(data) > 0) {
    if(length(levels(anno[[indfactor]])) >1) {
      gene.signature.xlsx.name <- paste(prefix, "impact", impact.name, paste("gene.signature", gene.signature, sep = ":"), paste("FDR", FDR, sep = ":"), "comparison_summary.xlsx", sep = ".")
      
      wb3 <- createWorkbook()
      combinations <- combn(levels(anno[[indfactor]]), m = 2)
      for(iter in seq(1,ncol(combinations))) {
        comparison.group <- substr(paste(combinations[,iter], collapse = "_vs_"),1,31)
        addWorksheet(wb = wb3, sheetName = comparison.group)
        anno.sel <- anno[anno[[indfactor]] %in% combinations[, iter], , drop = F]
        anno.sel[[indfactor]] <- factor(anno.sel[[indfactor]])
        data.sel <- data[match(sub(rownames(anno.sel), pattern = "#.*$", replacement = ""), rownames(data)), , drop = F]
        if(! all(rownames(data.sel) == sub(rownames(anno.sel), pattern = "#.*$", replacement = ""))) {stop("Sample names in the data and anno objects do not match.")}
        
        missing.group <- combinations[,iter][! combinations[,iter] %in% anno.sel[[indfactor]]]
        if(length(missing.group) == 0) {
          min.N <- min(table(anno.sel[[indfactor]]))
          min.N.limit <- 2
          if(min.N >= min.N.limit) {
            all.tested.genes <- unique(colnames(data.sel))
            
            f.data.stats <- function(data) {rbind(colMeans(data[,-1], na.rm = T),
                                                  Sds <- colSds(as.matrix(data[,-1]), na.rm = T),
                                                  pooledSd <- sqrt(sum(sapply(Sds, function(x) {x**2}))/2),
                                                  colSums(!is.na(data[,-1])),
                                                  if(paired.samples) {
                                                    sddif <- sd(data[,2] - data[,3])
                                                    cor.res <- cor.test(x = data[,2], data[,3])
                                                    if(is.na(cor.res$estimate)) {cor.res$estimate <- 0} else
                                                      if(cor.res$estimate == 1) {cor.res$estimate <- 0.9999999}
                                                    rbind(rep(sddif,2),
                                                          rep(cor.res$estimate,2)) } )
            }
            ff.10 <- pOverA(0.101,0, na.rm = T)
            
            value.conv <- function(value) {
              if(!is.na(value)) {
                value <- as.numeric(value)
                if(abs(value) < 0.001) {
                  value <- formatC(value, format = "e", digits = 3)} else
                  {
                    value <- round(value,4)}
                return(as.numeric(value))} else {return(value)}}
            
            TOST.res <- foreach(gene.name = all.tested.genes, .combine = rbind) %do% {
              if(! paired.samples) {
                df.tmp <- data.frame(sample = rownames(anno.sel), gene = data.sel[[gene.name]], indfactor = anno.sel[[indfactor]])
                df.tmp <- reshape(data = df.tmp, timevar = "indfactor", direction = "wide", idvar = "sample")
              } else {
                df.tmp <- data.frame(sample = rownames(anno.sel), gene = data.sel[[gene.name]], indfactor = anno.sel[[indfactor]], Pairs = anno.sel[[pair.ident]])
                df.tmp <- reshape(df.tmp[,-1], idvar = "Pairs", timevar = "indfactor", direction = "wide")
              }
              
              if(any(sapply(df.tmp[,-1],ff.10))) {
                df.tmp.stats.two <- f.data.stats(data = df.tmp)
                df.tmp.stats.two[2,][df.tmp.stats.two[2,] == 0] <- 1e-50
                
                sink("/dev/null")
                sink(type = "m")
                
                if(paired.samples) {
                  eqbound <- suppressMessages(round(max(powerTOSTpaired.raw(alpha=0.05, statistical_power=0.8, N = min.N, sdif = df.tmp.stats.two[5,1])),2)) } else {
                    eqbound <- suppressMessages(round(max(powerTOSTtwo.raw(alpha=0.05, statistical_power=0.8, N = min.N, sdpooled = df.tmp.stats.two[3,1])),2)) }
                
                sink()
                
                if(eqbound == 0) {eqbound <- 0.1}
                
                if(paired.samples) {
                  t1 <- tsum_TOST(m1=df.tmp.stats.two[1,1],
                                  m2=df.tmp.stats.two[1,2],
                                  sd1=df.tmp.stats.two[2,1],
                                  sd2=df.tmp.stats.two[2,2],
                                  n1=df.tmp.stats.two[4,1],
                                  n2=df.tmp.stats.two[4,2],
                                  low_eqbound=-eqbound,
                                  high_eqbound=eqbound,
                                  alpha = 0.05,
                                  var.equal=FALSE,
                                  eqbound_type = "raw",
                                  paired = TRUE,
                                  r12 = df.tmp.stats.two[6,1])
                } else {
                  t1 <- tsum_TOST(m1=df.tmp.stats.two[1,1],
                                  m2=df.tmp.stats.two[1,2],
                                  sd1=df.tmp.stats.two[2,1],
                                  sd2=df.tmp.stats.two[2,2],
                                  n1=df.tmp.stats.two[4,1],
                                  n2=df.tmp.stats.two[4,2],
                                  low_eqbound=-eqbound,
                                  high_eqbound=eqbound,
                                  alpha = 0.05,
                                  var.equal=FALSE,
                                  eqbound_type = "raw",
                                  paired = FALSE)
                }
                
                group.names <- gsub(colnames(df.tmp.stats.two), pattern = "gene\\.", replacement = "")
                v1 <- c(t1$method, group.names, df.tmp.stats.two[4,], df.tmp.stats.two[1,], df.tmp.stats.two[2,], -eqbound, eqbound, if(paired.samples) {df.tmp.stats.two[6,1]} else {NA}, t1$TOST["TOST Lower","p.value"], t1$TOST["TOST Upper","p.value"], t1$TOST["t-test","p.value"])
                v1[4:length(v1)] <- sapply(v1[4:length(v1)], value.conv)
                
                if(max(c(t1$TOST["TOST Lower","p.value"], t1$TOST["TOST Upper","p.value"])) < 0.05 & t1$TOST["t-test","p.value"] >= 0.05) {
                  v1 <- c(gene.name, v1, "Yes")
                } else {
                  v1 <- c(gene.name, v1, "No")
                }
                v1
              }
            }
            if(!is.null(TOST.res)) {
              if(is.vector(TOST.res)) {TOST.res <- as.data.frame(t(TOST.res), stringsAsFactors = F)} else
              {TOST.res <- as.data.frame(TOST.res, stringsAsFactors = F)}
              
              rownames(TOST.res) <- NULL
              colnames(TOST.res) <- c("Gene name", "Test type", "Group1 name", "Group2 name", "Group1 N", "Group2 N", "Group1 mean", "Group2 mean", "Group1 SD", "Group2 SD", "Lower equivalence bounds", "Upper equivalence bounds", "Pearson's correlation R2", "TOST Lower p-value", "TOST Upper p-value", "NHST t-test p-value","Equivalence")
              TOST.res <- TOST.res %>% mutate(across(!c(1:4,17), as.numeric))
              if(FDR) {
                TOST.res$`NHST t-test BH-adjusted p-value` <- p.adjust(TOST.res$`NHST t-test p-value`, method = "BH")
                
                TOST.res.no.NA <- TOST.res[rowSums(is.na(TOST.res[,-13])) == 0,, drop = F]
                TOST.res.no.NA[["Difference"]][(TOST.res.no.NA$`TOST Lower p-value` >= 0.05 | TOST.res.no.NA$`TOST Upper p-value` >= 0.05) & TOST.res.no.NA$`NHST t-test BH-adjusted p-value` < 0.05] <- "Yes"
                TOST.res <- merge(x = TOST.res, y = TOST.res.no.NA %>% dplyr::select(c("Gene name", "Difference")), by.x = "Gene name", by.y = "Gene name", all = T)
                TOST.res <- distinct(TOST.res)
                TOST.res$Difference[is.na(TOST.res$Difference)] <- "No"
                TOST.res <- TOST.res[,c(1:16,18,17,19)]
              } else {
                TOST.res.no.NA <- TOST.res[rowSums(is.na(TOST.res[,-13])) == 0,, drop = F]
                TOST.res.no.NA[["Difference"]][(TOST.res.no.NA$`TOST Lower p-value` >= 0.05 | TOST.res.no.NA$`TOST Upper p-value` >= 0.05) & TOST.res.no.NA$`NHST t-test p-value` < 0.05] <- "Yes"
                TOST.res <- merge(x = TOST.res, y = TOST.res.no.NA %>% dplyr::select(c("Gene name", "Difference")), by.x = "Gene name", by.y = "Gene name", all = T)
                TOST.res <- distinct(TOST.res)
                TOST.res$Difference[is.na(TOST.res$Difference)] <- "No"
              }
              TOST.res <- TOST.res[order(as.character(TOST.res[["Equivalence"]]), as.character(TOST.res[["Difference"]]), as.character(TOST.res[["Gene name"]]), decreasing = c(T,T,F), method = "radix"),]
              
              evalGenes.N <- if(gene.signature == "ALL_GENES") {ncol(data.sel)} else {length(gene.list1.conv)}
              equivalent.genes.N <- sum(TOST.res$Equivalence == "Yes")
              non.equivalent.genes.N <- sum(TOST.res$Equivalence == "No")
              different.genes.N <- sum(TOST.res$Difference == "Yes")
              non.different.genes.N <- sum(TOST.res$Difference == "No")
              
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes evaluated to determine the gene signature:", evalGenes.N), startRow = 1)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes with at least one alteration in at least one sample:", length(all.tested.genes)), startRow = 2)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes passing the variant frequency filter (at least one alteration in more than 10% of samples in at least one group):", nrow(TOST.res)), startRow = 3)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes with non-equivalent alteration profiles:", non.equivalent.genes.N), startRow = 5)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes with equivalent alteration profiles:", equivalent.genes.N), startRow = 6)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The percentage of genes with equivalent alteration profiles:", round(equivalent.genes.N/nrow(TOST.res)*100,2), "%"), startRow = 7)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes with non-different alteration profiles:", non.different.genes.N), startRow = 9)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes with different alteration profiles:", different.genes.N), startRow = 10)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The percentage of genes with different alteration profiles:", round(different.genes.N/nrow(TOST.res)*100,2), "%"), startRow = 11)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The table below contains TOST and NHST statistical results for all the genes that passed the variant frequency filtering."), startRow = 13)
              
              writeData(wb = wb3, sheet = comparison.group, x = TOST.res, startRow = 15, keepNA = TRUE, na.string = "NA")
            } else {
              evalGenes.N <- if(gene.signature == "ALL_GENES") {ncol(data.sel)} else {length(gene.list1.conv)}
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes evaluated to determine the gene signature:", evalGenes.N), startRow = 1)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes with at least one alteration in at least one sample:", length(all.tested.genes)), startRow = 2)
              writeData(wb = wb3, sheet = comparison.group, x = paste("The number of genes passing the variant frequency filter (at least one alteration in more than 10% of samples in at least one group):", 0), startRow = 3)
            }
          } else {
            writeData(wb = wb3, sheet = comparison.group, x = paste0("The number of samples in one of the analyzed groups: ", min.N, " is lower than the predefined limit: ", min.N.limit, ". Therefore, the gene signature comparison has not been performed."))
          }
        } else {
          writeData(wb = wb3, sheet = comparison.group, x = paste0("There are no samples in the following analyzed group(s): ", paste(missing.group, collapse = ", "), "."))
        }
      }
      saveWorkbook(wb = wb3, file = gene.signature.xlsx.name, overwrite = T)
    }
  }
}

workdir <- paste0(arguments[1], "/Summary", "/Grouping_variables:", arguments[4], "/Factor:", indfactor)

runid <- sub(sub(workdir, pattern = "^.*RUNS/", replacement = ""), pattern = "/MAPPINGS.*$", replacement = "")
antype <- sub(sub(workdir, pattern = "\\/Summary\\/.*", replacement = ""), pattern = "^.*\\/", replacement = "")
impacts <- c("HIGH", "(HIGH|MODERATE)")

dir.create(workdir, recursive = T)
gene.signature <- gsub(txt.file, pattern = "^.*\\/", replacement = "")
setwd(workdir)
unlink("Warning.no.vcfs.txt")
prefix <- paste(runid, antype, groupid, indfactor, paste("paired", paired.samples, sep = ":"), sep = ".")
save.image(file = paste(prefix, "VEP_analysis.initial.RData", sep = "."))

if(! file.exists(paste(prefix, "VEP_analysis.final.RData", sep = "."))) {

file.list <- list.files("../../..", pattern = csv.pat, full.names = T)
# Exclude empty csv files 
file.list <- file.list[sapply(file.list, file.size) > 0]
s.names <- sub(file.list, pattern = paste0("(\\.\\.\\/\\.\\.\\/\\.\\.\\/)(.*)", paste0("(", csv.pat, ")")), replacement = "\\2")
s.names <- gsub(s.names, pattern = "#", replacement = ".")

if(length(file.list) >0) {

sel.columns <- c("Gene", "SYMBOL", "SYMBOL_SOURCE", "HGVSg", "HGVSc", "HGVSp", "Existing_variation", "Consequence", "SIFT", "PolyPhen", "MAX_AF", "CLIN_SIG", "PUBMED", "ZYG", "IMPACT")

if(!file.exists(paste0("../../", runid, ".", antype, ".df.full.RData"))) {

cat("Generating the first list of hits...\n")
df.list <- foreach(i = file.list) %dopar% {
file.tmp <- fread(file = i, sep = "\t", header = T, select = c(sel.columns))
file.tmp <- file.tmp[file.tmp[["SYMBOL_SOURCE"]] %in% c("EntrezGene", "HGNC") & grepl(file.tmp[["IMPACT"]], pattern = "(HIGH|MODERATE)"), , drop = F]
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
  
  df.tmp <- df.tmp[ !grepl(df.tmp[["SIFT"]], pattern = "(tolerated|-)") | !grepl(df.tmp[["PolyPhen"]], pattern = "(benign|-|unknown)") | !grepl(df.tmp[["CLIN_SIG"]], pattern = "(benign|-|^not_provided$|^protective$)") | as.numeric(df.tmp[["MAX_AF"]]) < 0.01, , drop = F]
  df.tmp[["MAX_AF"]][df.tmp[["MAX_AF"]] == "-1"] <- "-"
  df.tmp
   }}

f.cytoband <- function(x) {
  chromosome <- stri_extract(stri_extract(str = x, regex = "^chr[0-9MXY]+"), regex = "[0-9MXY]+")
  g.loc <- as.integer(stri_extract(stri_extract(str = x, regex = "g\\.[0-9]+"), regex = "[0-9]+"))
  if(chromosome %in% c(1:22, "X", "Y")) {position2Cytoband(chrom = chromosome, position = g.loc, units = genome)} else
  {return ("-")}
}

Cytoband <- foreach(loc = df.full[,"HGVSg"], .combine = c) %dopar% {f.cytoband(loc)}
df.full <- cbind(df.full, Cytoband)[,c(1:7,17,8:16)]

rm(df.list)

if(!is.null(df.full)) {
df.full <- as.matrix(df.full[with(df.full, order(Sample.name, SYMBOL, IMPACT)),])
}
save(list = c("df.full"), file = paste0("../../", runid, ".", antype, ".df.full.RData"))
} else {
  cat("Loading a pre-existing RData file:", paste0("../../", runid, ".", antype, ".df.full.RData...\n"))
  load(paste0("../../", runid, ".", antype, ".df.full.RData"), envir = .GlobalEnv)
  read.args()
}
if(!is.null(df.full)) {
  
  if(paired.samples) {csvtable[[sampleid]] <- paste(csvtable[[sampleid]], csvtable[[pair.ident]], sep = "#")}
  rownames(csvtable) <- csvtable[[sampleid]]
  if(length(csvtable[[sampleid]]) != length(csvtable[[indfactor]])) {stop("Some sample ids or group ids are missing.")}
  
  if(paired.samples)
  { anno <- csvtable %>% dplyr::select(indfactor, pair.ident)} else
  { anno <- csvtable %>% dplyr::select(indfactor)}
  anno <- as.data.frame(sapply(anno, factor))
  rownames(anno) <- rownames(csvtable)
  
  anno.sub <- anno
  rownames(anno.sub) <- sub(rownames(anno.sub), pattern = "#.*$", replacement = "")
  df.full <- as.matrix(merge(x = df.full, y = anno.sub, by.x = "Sample.name", by.y = 0))
  if(paired.samples) {df.full <- df.full[, colnames(df.full)[c(1,18,19,2:17)], drop = F]} else {
    df.full <- df.full[, colnames(df.full)[c(1,18,2:17)], drop = F]}
  if(nrow(df.full) == 0) {df.full <- NULL}
  }

if(!is.null(df.full)) {

wb <<- createWorkbook()

if(nrow(df.full) <= 2**20) {
  sheet.name <- substr(paste("All_variants", sep = "."), 1, 31)
  if(length(wb$sheet_names) > 0) {if(sheet.name %in% wb$sheet_names) {removeWorksheet(wb, sheet = sheet.name)}}
  addWorksheet(wb = wb, sheetName = sheet.name)
  writeData(x = df.full, wb = wb, sheet = sheet.name) } else {
  All.var.name <- paste("VEP_analysis", "all_variants", prefix, "csv" , sep = ".")
  write.table(x = df.full, row.names = F, file = All.var.name, sep = ";")
}

f.hm <- function(x) {df.tmp <-  x[grepl(x[["IMPACT"]], pattern = impact),]
if(nrow(df.tmp) > 0) {
  return(as.matrix(df.tmp[!duplicated(df.tmp$HGVSg),]))}}

for(impact in impacts) {
impact.name <- paste(sub(sub(unlist(strsplit(impact, split = "\\|")), pattern = "\\(", replacement = ""), pattern = "\\)", replacement = ""), collapse = "_or_")
New.vars.mx <- df.full[df.full[,"Existing_variation"] == "-", , drop = F]
New.vars.mx <- New.vars.mx[grepl(New.vars.mx[,"IMPACT", drop = F], pattern = impact), , drop = F]

if(nrow(New.vars.mx) <= 2**20) {
sheet.name <- substr(paste("New_variants", "impact", impact.name, sep = "."), 1, 31)
if(length(wb$sheet_names) > 0) {if(sheet.name %in% wb$sheet_names) {removeWorksheet(wb, sheet = sheet.name)}}
addWorksheet(wb = wb, sheetName = sheet.name)
writeData(x = New.vars.mx, wb = wb, sheet = sheet.name) } else {
New.var.name <- paste("VEP_analysis", "new_variants", impact.name, prefix, "csv" , sep = ".")
write.table(x = New.vars.mx, row.names = F, file = New.var.name, sep = ";")
}

df.hm.list <- by(data = df.full, INDICES = list(df.full[,"Sample.name"], df.full[,"SYMBOL"]), FUN=f.hm)
names(df.hm.list) <- levels(interaction(df.full[,"Sample.name"], df.full[,"SYMBOL"], sep = "$$"))
df.hm.list <- df.hm.list[sapply(df.hm.list, Negate(is.null))]

geneList <- sort(unique(sub(names(df.hm.list), pattern = "^.*\\$\\$", replacement = "")))

sampleIDs <- s.names
for(i in seq(1,length(sampleIDs))) {
  if(sum(grepl(rownames(anno), pattern = paste0("^", sampleIDs[i], "(#.*)?", "$"))) > 1) {stop(paste(sampleIDs[i], "sample cannot be uniquely identified in the annotation table."))}}

res1 <- foreach(gene = geneList, .combine = cbind) %:% foreach(sampleID = sampleIDs, .combine = c) %dopar% {if(! is.null(df.hm.list[[paste(sampleID, gene, sep = "$$")]])) {nrow(df.hm.list[[paste(sampleID, gene, sep = "$$")]])} else {0}}
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
if(paired.samples) {
  res1.merged <- res1.merged[order(res1.merged[[indfactor]], res1.merged[[pair.ident]]),]} else {
  res1.merged <- res1.merged[order(res1.merged[[indfactor]]),]}

if(paired.samples) {
  res1 <- res1.merged %>% dplyr::select(-c("Row.names", indfactor, pair.ident, "full.sample.names"))
  anno <- res1.merged %>% dplyr::select(indfactor, pair.ident, "full.sample.names")
  rownames(anno) <- anno[["full.sample.names"]]
  anno <- anno %>% dplyr::select(-c("full.sample.names"))
} else {
  res1 <- res1.merged %>% dplyr::select(-c("Row.names", indfactor, "full.sample.names"))
  anno <- res1.merged %>% dplyr::select(indfactor) }

New.vars.name <- paste("new.vars", antype, groupid, make.names(impact), sep = ".")
assign(New.vars.name, value = New.vars.mx, envir = .GlobalEnv)

res1.name <- paste("res1", antype, groupid, make.names(impact), sep = ".")
assign(res1.name, value = res1, envir = .GlobalEnv)

anno.name <- paste("anno", antype, groupid, make.names(impact), sep = ".")
assign(anno.name, value = anno, envir = .GlobalEnv)

var.analysis(data = res1, anno = anno, gene.signature = "ALL_GENES", wb = wb)
f_TOST(data = res1, anno = anno, gene.signature = "ALL_GENES")
}

pdf_combine(rev(list.files(pattern = "VEP_analysis.*\\.pdf$")), output = paste("VEP_analyses", "ALL_GENES", prefix, "pdf", sep = "."))
unlink(rev(list.files(pattern = "VEP_analysis.*\\.pdf$")))

saveWorkbook(wb = wb, file = paste(prefix, paste("gene.signature", "ALL_GENES", sep = ":"), "VEP_analysis.xlsx", sep = "."), overwrite = T)
} else {
  
  if(paired.samples) {csvtable[[sampleid]] <- paste(csvtable[[sampleid]], csvtable[[pair.ident]], sep = "#")}
  rownames(csvtable) <- csvtable[[sampleid]]
  if(length(csvtable[[sampleid]]) != length(csvtable[[indfactor]])) {stop("Some sample ids or group ids are missing.")}
  
  if(paired.samples)
  { anno <- csvtable %>% dplyr::select(indfactor, pair.ident)} else
  { anno <- csvtable %>% dplyr::select(indfactor)}
  anno <- as.data.frame(sapply(anno, factor))
  rownames(anno) <- rownames(csvtable)
  
  for(impact in impacts) {
  
  anno.name <- paste("anno", antype, groupid, make.names(impact), sep = ".")
  assign(anno.name, value = anno, envir = .GlobalEnv)
  
  res1 <- data.frame(s.names)
  rownames(res1) <- s.names
  res1 <- res1[,-1]
  res1 <- res1[match(rownames(anno), rownames(res1)), , drop = F]
  
  res1.name <- paste("res1", antype, groupid, make.names(impact), sep = ".")
  assign(res1.name, value = res1, envir = .GlobalEnv)
  
  res1.bool.name <- paste("res1.bool", antype, groupid, make.names(impact), sep = ".")
  assign(res1.bool.name, value = res1, envir = .GlobalEnv)
  
  New.vars.name <- paste("new.vars", antype, groupid, make.names(impact), sep = ".")
  New.vars.mx <- t(matrix(rep("",18)))[-1,]
  colnames(New.vars.mx) <- c("Sample.name", "Group", "Gene", "SYMBOL", "SYMBOL_SOURCE", "HGVSg", "HGVSc", "HGVSp", "Cytoband", "Existing_variation", "Consequence", "SIFT", "PolyPhen", "MAX_AF", "CLIN_SIG", "PUBMED", "ZYG", "IMPACT")
  assign(New.vars.name, value = New.vars.mx, envir = .GlobalEnv)
  }
  cat(paste("There are no", antype, "variants (alterations classified as high or moderate in the Ensembl database) in the current data set.\n\nSample IDs:", paste(rownames(res1), collapse = ", "), "\n"), file = paste("VEP_analyses", "ALL_GENES", prefix, "txt", sep = "."))
}
save.image(file = paste(prefix, "VEP_analysis.final.RData", sep = "."))
} else {
  cat("There are no matching CSV files with VEP results to be processed.\nPlease, check if the corresponding", BAM.type, "and VCF files are present.\n", file = "Warning.no.vcfs.txt")
}
} else {
  cat("Loading a pre-existing RData file:", paste(prefix, "VEP_analysis.final.RData...\n", sep = "."))
  load(paste(prefix, "VEP_analysis.final.RData", sep = "."), envir = .GlobalEnv)
  read.args()

for(impact in impacts) {
  impact.name <- paste(sub(sub(unlist(strsplit(impact, split = "\\|")), pattern = "\\(", replacement = ""), pattern = "\\)", replacement = ""), collapse = "_or_")

  res1.name <- paste("res1", antype, groupid, make.names(impact), sep = ".")
  
  anno.name <- paste("anno", antype, groupid, make.names(impact), sep = ".")
  
  f_TOST(data = get(res1.name), anno = get(anno.name), gene.signature = "ALL_GENES")
}
}

if(! file.exists("Warning.no.vcfs.txt")) {

if(txt.file != "NA") {
  con <- file(txt.file)
  gene.list1 <- unique(gsub(readLines(con), pattern = "\\s", replacement = ""))
  close(con)

  gene.list1.conv <- NULL
  for(i in gene.list1) {if(length(alias2Symbol(i, species = "Hs")) == 0) {gene.list1.conv <- append(gene.list1.conv, values = i)} else {gene.list1.conv <- append(gene.list1.conv, values = alias2Symbol(i, species = "Hs"))}}
  gene.list1.conv <- unique(gene.list1.conv)
  
  cat("Loading a pre-existing RData file:", paste0("../../", runid, ".", antype, ".df.full.RData...\n"))
  load(paste0("../../", runid, ".", antype, ".df.full.RData"), envir = .GlobalEnv)
  read.args()
  
  if(!is.null(df.full)) {
    
    if(paired.samples) {csvtable[[sampleid]] <- paste(csvtable[[sampleid]], csvtable[[pair.ident]], sep = "#")}
    rownames(csvtable) <- csvtable[[sampleid]]
    if(length(csvtable[[sampleid]]) != length(csvtable[[indfactor]])) {stop("Some sample ids or group ids are missing.")}
    
    if(paired.samples)
    { anno <- csvtable %>% dplyr::select(indfactor, pair.ident)} else
    { anno <- csvtable %>% dplyr::select(indfactor)}
    anno <- as.data.frame(sapply(anno, factor))
    rownames(anno) <- rownames(csvtable)
    
    anno.sub <- anno
    rownames(anno.sub) <- sub(rownames(anno.sub), pattern = "#.*$", replacement = "")
    df.full <- as.matrix(merge(x = df.full, y = anno.sub, by.x = "Sample.name", by.y = 0))
    if(paired.samples) {df.full <- df.full[, colnames(df.full)[c(1,18,19,2:17)], drop = F]} else {
      df.full <- df.full[, colnames(df.full)[c(1,18,2:17)], drop = F]}
    df.full <- df.full[df.full[,"SYMBOL"] %in% gene.list1.conv, , drop = F]
    if(nrow(df.full) == 0) {df.full <- NULL}
  }
  
  if(!is.null(df.full)) {
    
    wb2 <<- createWorkbook()
    
    if(nrow(df.full) <= 2**20) {
      sheet.name <- substr(paste("All_variants", sep = "."), 1, 31)
      if(length(wb2$sheet_names) > 0) {if(sheet.name %in% wb2$sheet_names) {removeWorksheet(wb2, sheet = sheet.name)}}
      addWorksheet(wb = wb2, sheetName = sheet.name)
      writeData(x = df.full, wb = wb2, sheet = sheet.name) } else {
        All.var.name <- paste("VEP_analysis", "all_variants", paste("gene_signature", gene.signature, sep = ":"), prefix, "csv" , sep = ".")
        write.table(x = df.full, row.names = F, file = All.var.name, sep = ";")
      }
  
  for(impact in impacts) {
    
    impact.name <- paste(sub(sub(unlist(strsplit(impact, split = "\\|")), pattern = "\\(", replacement = ""), pattern = "\\)", replacement = ""), collapse = "_or_")
    New.vars.mx <- df.full[df.full[,"Existing_variation"] == "-", , drop = F]
    New.vars.mx <- New.vars.mx[grepl(New.vars.mx[,"IMPACT", drop = F], pattern = impact), , drop = F]
    
    if(nrow(New.vars.mx) <= 2**20) {
      sheet.name <- substr(paste("New_variants", "impact", impact.name, sep = "."), 1, 31)
      if(length(wb2$sheet_names) > 0) {if(sheet.name %in% wb2$sheet_names) {removeWorksheet(wb2, sheet = sheet.name)}}
      addWorksheet(wb = wb2, sheetName = sheet.name)
      writeData(x = New.vars.mx, wb = wb2, sheet = sheet.name) } else {
        New.var.name <- paste("VEP_analysis", "new_variants", impact.name, paste("gene_signature", gene.signature, sep = ":"), prefix, "csv" , sep = ".")
        write.table(x = New.vars.mx, row.names = F, file = New.var.name, sep = ";")
      }
    
      res1.name <- paste("res1", antype, groupid, make.names(impact), sep = ".")
      res1 <- get(res1.name)
      res1.sel <- res1 %>% dplyr::select(colnames(res1)[colnames(res1) %in% gene.list1.conv])
      anno.name <- paste("anno", antype, groupid, make.names(impact), sep = ".")
      anno <- get(anno.name)
      
      if(ncol(res1.sel) > 0 & nrow(res1.sel) > 0) {
      var.analysis(data = res1.sel, anno = anno, gene.signature = gene.signature, gene.list1.conv = gene.list1.conv, wb = wb2)
      f_TOST(data = res1.sel, anno = anno, gene.signature = gene.signature)
        } else {
        sheet.name <- substr(paste("Var.per.gene.impact",impact.name, sep = "."), 1, 31)
        if(length(wb2$sheet_names) > 0) {if(sheet.name %in% wb2$sheet_names) {removeWorksheet(wb2, sheet = sheet.name)}}
        addWorksheet(wb = wb2, sheetName = sheet.name)
        writeData(wb = wb2, sheet = sheet.name, x = c("No genes with genetic alterations have been found."))
      }}
  pdffiles <- rev(list.files(pattern = paste(runid, antype, "impact", "*", "pdf", sep = ".")))
  if(length(pdffiles) > 0) {
  pdf_combine(pdffiles, output = paste("VEP_analyses", paste0("gene.signature:", gene.signature), prefix, "pdf", sep = "."))
  unlink(pdffiles)} else {
  pdf(file = paste("VEP_analyses", paste0("gene.signature:", gene.signature), prefix, "pdf", sep = "."), width = 18)
    plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(paste(main="INFO: At least two genes contained in the", gene.signature, "signature have to harbor", antype, "variants (impact high or moderate) for the analysis to succeed."))
    dev.off()
  }
  saveWorkbook(wb = wb2, file = paste(prefix, paste("gene.signature", gene.signature, sep = ":"), "VEP_analysis.xlsx", sep = "."), overwrite = T)
  }}
save.image(file = paste(prefix, "VEP_analysis.final.RData", sep = "."))
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
date()

cat("All done.\n")
