library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(VennDiagram)

##Utility functions for processing ACSNI-derive predictions

#Get path and name of the required input files.
#ACSNI-derive predictions. Gene name should be the 1st column
pfilename <- "/path/to/ACSNI-derive/predictions/<gene>_finalpredict.csv"
#Full dataset. Gene name should be the 1st column
pfdilename <- "/path/to/ACSNI-derive/predictions/<gene>_fulldataset.csv"
#Ground-Truth data to compare predictions against. Usually these are DE genes from an independent study. Gene name should be the 1st column and there should also be a column with the name "adj.P.Val" containing adjusted pvalues.
gfilename <- "/path/to/ground_truth_data/<DEgenes>.csv"

ground_truth_label <- 'Ground-Truth'

#Function to design permutation test.
permutation_test <- function(n_boot, target, space, size, observed) {
  # Function to do permutation test
  # n_boot : number of bootstrap
  # target : Character vector to overlap
  # space : Character vector to sample from
  # size : the number of predicted in the other side
  # observed : size observed
  #
  p = 0
  out <- c()

  for (i in 1:n_boot) {
    sample <- sample(space, size, replace = FALSE)
    lap <- length(sample[sample %in% target])

    if (lap >= observed) {
      p = p + 1
    }

    out[i] <- lap
  }
  CI <- quantile(out, c(.025, .975))
  print(paste("95% CI [", CI[1], ",", CI[2], "]"))

  #
  p_value  = paste("p_value = ", p/n_boot)
  #

  hist(out, breaks = 10, col = "lightblue",
       xlab = "Overlap with random gene set",
       main = p_value,
       sub = paste("95% CI [", CI[1], ",", CI[2], "]"))
}
#Run permutation test.
run_permut <- function(pr,pf,g,n_boot=100000,pval=0.05) {
  ps <-  read.csv(file = pr)
  pfd <- read.csv(file = pf)
  gf <- read.csv(file = g)
  gf2 <- gf[gf[,1] %in% pfd[,1],]
  poverlp <- intersect(as.character(gf2[,1]), as.character(ps[,1]))
  poverlp <- as.data.frame(poverlp)
  rownames(poverlp) <- poverlp$poverlp
  gfde <- gf2[gf2$adj.P.Val <= pval,]
  de_overlap <- gfde[as.character(gfde[,1]) %in% as.character(poverlp$poverlp),]
  observed <- nrow(de_overlap)
  target <- as.factor(gfde[,1])
  space <- as.factor(gf2[,1])
  size <- nrow(poverlp)
  set.seed(23)
  out <- permutation_test(n_boot = n_boot, target = target, space = space, size = size, observed = observed)
  #return(out)
  get_inputs <- list(length(target),size,observed)
  return(get_inputs)
}
out <- run_permut(pr=pfilename,pf=pfdilename,g=gfilename)

#Draw VennDiagram
grid.newpage()
draw.pairwise.venn(area1 = out[[1]], area2 = out[[2]], cross.area = out[[3]], category = c(ground_truth_label, "ACSNI"), lty = rep("blank",2), fill = c("light blue", "orange"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE) #3,3 #area1 = target, area2 = size, cross.area=observed

#Function to get input for Cluster Profiler
get_inputs_ont <- function(pf,g,pval=0.05){
  pfd <- read.csv(file = pf)
  gf <- read.csv(file = g)
  names(pfd) <- c('gene',"pfd")
  pfd_1 <- pfd[pfd$pfd == 'True',]
  pfd_1$group = "pfd"
  gfde <- gf[gf$adj.P.Val < pval,]
  gfde$pfd <- 'DE'
  gfde$group <- 'DE'
  names(gfde)[1] <- 'gene'
  gfde <- gfde[,c("gene","pfd","group")]
  mpfdg <- rbind(pfd_1,gfde)
  return(mpfdg)
}
ont <- get_inputs_ont(pf=pfdilename,g=gfilename)

#Run Cluster Profiler to compare ontologies between ACSNI-derive predictions and ground-truth. You can change the keyType based on the identifier (for eg., ENSEMBL,ENTREZID,kegg,etc)
DF1a <- clusterProfiler::compareCluster(gene~group,
                                        keyType  = "SYMBOL",
                                        data = ont,
                                        fun = "enrichGO",
                                        ont = "BP",
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff  = 0.05,
                                        OrgDb = "org.Hs.eg.db")

#Plot significantly enriched ontologies
dotplot(DF1a, showCategory = 25, font.size = 5,  includeAll = TRUE)
#Simplify
DF1asim <- clusterProfiler::simplify(DF1a, cutoff = 0.60, measure = 'Wang', semData = NULL)
#Plot simplified version
dotplot(DF1asim, showCategory = 25, font.size = 5,  includeAll = TRUE)
