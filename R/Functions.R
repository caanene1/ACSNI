####### Library
library(ggpubr)
library(ComplexHeatmap)
library(ComplexHeatmap)
library(clusterProfiler)
library(circlize)
library(dendextend)
library(gplots)
library(VennDiagram)
######################## Start of Functions ####################################
### Geometric mean function
geo_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=TRUE) / length(x))}
# Function to count zeros and classify
predict_z <- function(x, y){
  ## Function to process the prediction
  #x_num <- x[c(1, ncol(x)-1)]
  x_num <- x[c("name", "Sum_stat")]
  print(names(x_num))
  #boot <- mean(x[[ncol(x)]]) * y / 100
  #x_num[2] <- (x_num[[2]] / boot) * 100
  x_num$P <- ifelse(x_num$Sum_stat >= y, "P", "B")
  print(summary(as.factor(x_num$P)))
  return(x_num) }
# Function to call Benchmark methods
predict_z_b <- function(x, y, cut){
  ## Function to process the prediction for benchmark
  x_num <- x[c("name", y)]
  x_num$P <- ifelse(abs(x_num[[y]]) >= cut, "P", "B")
  print(summary(as.factor(x_num$P)))
  return(x_num) }
# Function to convert ASCNI DERIVE output
predict_z_d <- function(x, y){
  ## Function to process the prediction from derive
  x_num <- x
  names(x_num) <- c("name", "stat")
  x_num$P <- ifelse(x_num$stat == "True", "P", "B")
  print(summary(as.factor(x_num$P)))
  return(x_num) }
# Function to table groups
tablate_group <- function(x) {
  x <- unname(summary(as.factor(x)))
  return(x)
}
## Function to calculate DE stats
stat_DE <- function(x=x, y=datasets, f=0.1, fc=0.5) {
  for (i in y) {
    x$FDR <- as.numeric(gsub(".*_", "", x[[i]]))
    x$FC <- abs(as.numeric(gsub("_.*", "", x[[i]])))
    x1 <- x[!is.na(x$FDR), ]
    x1$DE <- ifelse(x1$FDR <= f & x1$FC >= fc, "DE", "nDE")
    #
    tbl <- table(x1$P, x1$DE)
    print(tbl)
    Chiq <- chisq.test(tbl)
    #
    x_squared <- round(unname(Chiq[["statistic"]]), 2)
    p_value <- unname(Chiq[["p.value"]])
    p_value <- formatC(p_value, digits = 3, format = "e")
    #
    propT <- prop.table(tbl, 1)
    print(i)
    print(propT[2:2, 1:1] / propT[1:1, 1:1])
    print(p_value)
    print("________________")
    #
    mosaicplot(propT, shade = F,
               las = 1, off = 15,
               main = i,
               sub = paste("p =", p_value),
               ylab = "Genetic Pertubation",
               xlab = "Prediction",
               type = "pearson",
               color = c("darkred", "gray")) }
}
## Function to calculate DE per weight
stat_w_DE <- function(x=x, df=df, y=datasets, z="KLF6_GSE115763", f=0.1, fc=1.5) {
  for (i in y) {
    df$PP <- ifelse(df$P == "P" & df[[i]] == "P", "P", "B")
    x$FDR <- as.numeric(gsub(".*_", "", x[[z]]))
    x$FC <- abs(as.numeric(gsub("_.*", "", x[[z]])))
    x1 <- x[!is.na(x$FDR), ]
    x1$DE <- ifelse(x1$FDR <= f & x1$FC >= fc, "DE", "nDE")
    x1 <- merge(df, x1, by.x = "name", by.y = "Symbol")
    #
    tbl <- table(x1$PP, x1$DE)
    Chiq <- chisq.test(tbl)
    #
    x_squared <- round(unname(Chiq[["statistic"]]), 2)
    p_value <- unname(Chiq[["p.value"]])
    p_value <- formatC(p_value, digits = 3, format = "e")
    #
    propT <- prop.table(tbl, 1)
    print(i)
    print(propT[2:2, 1:1] / propT[1:1, 1:1])
    print(p_value) } }
## Function to perform permutation test
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

  hist(out, breaks = 10, col = "blue",
       xlab = "Overlap with random gene set",
       main = p_value,
       sub = paste("95% CI [", CI[1], ",", CI[2], "]"))
}
## Function to calculate ChIP enrichment for global
stat_chip_global <- function(x, y, chip){
  x_al <- merge(x, chip, by.x = "name", by.y = "gene", all.x = T)
  x_alb <- x_al
  for (i in y){
    x_alb[[i]] <- ifelse(is.na(x_alb[[i]]) | x_alb[[i]] <= 0, "P", "B")
  }

  Res <- data.frame(TF="base", ACSNI=0, p_density=0, ACSNI_Macs=0, p_Macs=0)
  for (i in y) {
    x_al1 <- x_al[!is.na(x_al[i]), ]
    p <- mean(x_al1[x_al1$P == "P",][[i]], na.rm = T)
    b <- mean(x_al1[x_al1$P == "B",][[i]], na.rm = T)
    p_density <- wilcox.test(x_al1[x_al1$P == "P",][[i]], x_al1[x_al1$P == "B",][[i]], na.rm = T)[["p.value"]]
    # Enrichment
    tbl <- table(x_alb[["P"]], x_alb[[i]])
    Chiq <- chisq.test(tbl)
    p_value <- unname(Chiq[["p.value"]])
    p_Macs <- formatC(p_value, digits = 3, format = "e")
    propT <- prop.table(tbl, 1)
    propT <- propT[2:2, 1:1] / propT[1:1, 1:1]
    #
    res <- data.frame(TF=i, ACSNI= (p/b), ACSNI_Macs=propT,
                      p_density = p_density,
                      p_Macs=p_Macs)
    Res <- rbind(Res, res)
  }
  Res <- Res[-1,]
}
## Function calculate the Enrichment for ChIP by layer
stat_chip_layer <- function(x, y, chip, boot){
  x_al <- merge(x, chip, by = "name", all.x = T)
  # For Enrichment type like case one
  x_alb <- x_al
  for (i in y){
    x_alb[[i]] <- ifelse(is.na(x_alb[[i]]) | x_alb[[i]] <= 0, "P", "B")
  }
  #
  Res <- data.frame(TF="Base", layer = "Base", Enrich=0, p_density=0, Enrich_Macs=0, p_Macs=0)
  #
  for (l in names(x[2:(ncol(x)-1)])) {
    x_al$P <- ifelse(x_al[[l]] == "P" & x_al$Sum_stat >= boot, "P", "B")
    x_alb$P <- ifelse(x_al[[l]] == "P" & x_al$Sum_stat >= boot, "P", "B")
    for (i in y) {
      x_al1 <- x_al[!is.na(x_al[i]), ]
      p <- mean(x_al1[x_al1$P == "P",][[i]], na.rm = T)
      b <- mean(x_al1[x_al1$P == "B",][[i]], na.rm = T)
      p_density <- wilcox.test(x_al1[x_al1$P == "P",][[i]], x_al1[x_al1$P == "B",][[i]], na.rm = T)[["p.value"]]
      # Enrichment
      tbl <- table(x_alb[["P"]], x_alb[[i]])
      Chiq <- chisq.test(tbl)
      p_value <- unname(Chiq[["p.value"]])
      p_Macs <- formatC(p_value, digits = 3, format = "e")
      propT <- prop.table(tbl, 1)
      propT <- propT[2:2, 1:1] / propT[1:1, 1:1]
      #
      res <- data.frame(TF=i, layer = l, Enrich=(p/b), p_density=p_density,
                        Enrich_Macs=propT, p_Macs=p_Macs)
      Res <- rbind(Res, res) }
  }
  Res <- Res[-1,]
  return(Res)
}
## Function to process CHIP density and combine for HUVCE
process_Chip_density <- function(path){
  old_path <- getwd()
  setwd(path)

  files <- dir(pattern = ".tsv")
  # read function
  rd <- function(x){ read.delim(file = x) }
  # sub function
  sb <- function(x){ sub(".1.tsv","", x) }

  # sub function
  sbg <- function(x){
    gsub("^.*\\.","",x)
  }
  #
  fil <- "JUN.1.tsv"
  get_HAEC <- rd(fil)
  names(get_HAEC) <- sbg(names(get_HAEC))
  get_HAEC <- get_HAEC[,names(get_HAEC) %in% c("Target_genes", "HAEC")]
  nu <- ncol(get_HAEC) - 1
  get_HAEC[["HAEC"]] <- rowSums(get_HAEC[-1]) / nu
  get_HAEC <- get_HAEC[1:2]
  names(get_HAEC)[2] <- sb(fil)
  names(get_HAEC) <- c("gene", sb(fil))

  #
  for (i in c("JUNB.1.tsv", "IRF1.1.tsv", "ERG.1.tsv", "EP300.1.tsv",
              "RELA.1.tsv", "NFE2L2.1.tsv", "CEBPD.1.tsv")) {
    print(i)
    fil_a <- rd(i)
    names(fil_a) <- sbg(names(fil_a))
    fil_a <- fil_a[,names(fil_a) %in% c("Target_genes", "HAEC")]
    nu <- ncol(fil_a) - 1
    fil_a[["HAEC"]] <- rowSums(fil_a[-1]) / nu
    fil_a <- fil_a[1:2]
    names(fil_a) <- c("gene", sb(i))
    get_HAEC <- merge(get_HAEC, fil_a, by = "gene", all = TRUE)
  }
  setwd(old_path)
  return(get_HAEC)
}
## Function to perform robust overlap
robust_overP <- function(x=PP, y=PP2, b=100){
  P1 <- x[c(1,3)]
  names(P1)[2] <- "P1"
  P2 <- y[c(1,3)]
  names(P2)[2] <- "P2"
  df <- merge(P1, P2, by = "name")
  df$OP <- ifelse(df$P1 == "P" & df$P2 == "P", "P", "B")
  res <- data.frame(ID = "real", overlap=nrow(df[df$OP == "P",]))

  df_new <- df
  for(i in 1:b){
    df_new$P2 <- sample(df$P2)
    df_new$OP <- ifelse(df_new$P1 == "P" & df_new$P2 == "P", "P", "B")
    res_new <- data.frame(ID = "random", overlap=nrow(df_new[df_new$OP == "P",]))
    res <- rbind(res, res_new)
  }

  return(res)
}
###################### End of Functions ########################################


##################### Simulation studies functions #############################
##
## Function to stimulate geneset expression pattern
sim_geneSet_pat <- function(exp=AA, geneSet=mTOR) {
  ExpressionData <- exp[exp$gene %in% geneSet, ]
  ExpressionData$gene <- paste("Ingene_", ExpressionData$gene, sep = "")
  return(ExpressionData)
}

## Function to get the prior matrix
get_sim_prior <- function(y){
  simSet <- data.frame("gene" = y, "simSet" = 1)
  return(simSet)
}

## Function to shuffle expression mat
shuffle_exp <- function(x){
  # Note: No seed is intentional to avoid unwanted patterns
  x_class <- sapply(x, class)
  num_cols <- names(x[x_class %in% "numeric"])
  if (length(num_cols) == 0) {
    stop("Expression matrix must have numeric columns")
  }

  for (i in num_cols) {
    cat("Shuffling", i, "\n")
    x[i] <- sample(x[[i]])
  }
  return(x)
}

## Function to whole transcriptome expression pattern
sim_transcriptome <- function(exp=AA, geneSet=mTOR){
  ExpressionData <- exp[!exp$gene %in% geneSet,]
  ExpressionData <- shuffle_exp(ExpressionData)
  ExpressionData$gene <- paste("trans_", ExpressionData$gene, sep = "")
  return(ExpressionData)
}

## Function to generate inverse of the the gene set pattern positive control
sim_invert_geneSet_pat <- function(frame){
  new_frame <- frame
  num <- nrow(frame)

  #Reverse the columns
  for (i in 1:num){
    original <- as.numeric(frame[i:i, -1])
    rev_vec <- rev(original)
    new_frame[i:i, -1] <- rev_vec
  }
  new_frame$gene <- paste("invert_", new_frame$gene, sep = "")
  return(new_frame)
}

## Function to add lognormal noise to data frame
sim_add_lnorm_noise <- function(frame, factor=0.1){
  new_frame <- frame
  num <- nrow(frame)
  n_col <- ncol(frame) - 1
  ## Generate random noise

  #Reverse the columns
  for (i in 1:num){
    original <- as.numeric(frame[i:i, -1])
    m <- mean(original)
    # Generate noise using factor
    k <- rlnorm(original, meanlog = log(m+1), sdlog = log((factor*m)+1))
    noise_vec <- original + k
    new_frame[i:i, -1] <- noise_vec
  }
  new_frame$gene <- paste("noise", factor, "_", new_frame$gene,   sep = "")
  return(new_frame)
}

## Function to combine simulation data
get_and_save <- function(x=AA, y=prior, namesave="path") {
  df1 <- sim_geneSet_pat(exp=x, geneSet=y)
  prior <- get_sim_prior(df1)
  df2 <- sim_transcriptome(exp=x, geneSet=y)
  # Invert
  df3 <- sim_invert_geneSet_pat(df1)
  # Noise
  df4 <- sim_add_lnorm_noise(frame = df1, factor = 0.05)
  df5 <- sim_add_lnorm_noise(frame = df1, factor = 0.06)
  df6 <- sim_add_lnorm_noise(frame = df1, factor = 0.07)
  df7 <- sim_add_lnorm_noise(frame = df1, factor = 0.08)
  df8 <- sim_add_lnorm_noise(frame = df1, factor = 0.09)
  df9 <- sim_add_lnorm_noise(frame = df1, factor = 0.1)
  df10 <- sim_add_lnorm_noise(frame = df1, factor = 0.2)
  df11 <- sim_add_lnorm_noise(frame = df1, factor = 0.3)
  df12 <- sim_add_lnorm_noise(frame = df1, factor = 0.4)
  df13 <- sim_add_lnorm_noise(frame = df1, factor = 0.5)

  ## Do shuffled
  df14 <- shuffle_exp(df1)
  df14$gene <- paste("shuffled_", df14$gene, sep = "")

  ### Save the data
  df <- do.call(rbind, list(df1, df2, df3, df4, df5, df6, df7, df8,
                            df9, df10, df11, df12, df13, df14))

  ## Write table
  write.csv(prior, file = paste(namesave, "_P1.csv", sep = ""), row.names = F)
  write.csv(df, file = paste(namesave, "_P2.csv", sep = ""), row.names = F)
}

############# Load data for simualtion
BiocManager::install("qusage")
setG <- qusage::read.gmt("~/Desktop/Install ACSNI-Linux/c2.cp.pid.v7.2.symbols.gmt")

# Filter list
newlist <- list()
for (i in names(setG)) {
  if (length(setG[[i]]) > 25 & length(setG[[i]]) < 200) {
    newlist[i] <- setG[i]
  }
}

# Get sample
setGSample <- sample(newlist, 50)

## Get the expressions
Breast <- read.csv("~/Desktop/Install ACSNI-Linux/Exp/Breast.csv")
Liver <- read.csv("~/Desktop/Install ACSNI-Linux/Exp/Liver.csv")
Pancreas <- read.csv("~/Desktop/Install ACSNI-Linux/Exp/Pancreas.csv")
Prostate <- read.csv("~/Desktop/Install ACSNI-Linux/Exp/Prostate.csv")
Stomach <- read.csv("~/Desktop/Install ACSNI-Linux/Exp/Stomach.csv")

### Usage of the code for the reported simulations
for (i in names(setGSample)) {
  print(i)
  nnm1 <- paste("Breast_", i, sep = "")
  nnm2 <- paste("Liver_", i, sep = "")
  nnm3 <- paste("Pancreas_", i, sep = "")
  nnm4 <- paste("Prostate_", i, sep = "")
  nnm5 <- paste("Stomach_", i, sep = "")


  get_and_save(x=Breast, y = setGSample[[i]], namesave = nnm1)
  get_and_save(x=Liver, y = setGSample[[i]], namesave = nnm2)
  get_and_save(x=Pancreas, y = setGSample[[i]], namesave = nnm3)
  get_and_save(x=Prostate, y = setGSample[[i]], namesave = nnm4)
  get_and_save(x=Stomach, y = setGSample[[i]], namesave = nnm5)
}

## Function to calculate the rates in a simulation
calc_simulate_rate <- function(path, files) {
  Result <- data.frame("ID"="base", "Sum"=0, "B" = 0,
                       "P"=0, "Rate"=0, "mad"=0, "alpha"=0, "percentage"=0,
                       "totalP" = 0, "Size" = 0,  "pathway" = 0, "reg"=0)

    for (f in files){
    P <- read.csv(paste(path, "/", f, sep = ""))
    P$Group <- gsub("(_).*", "\\1", P$gene)
    set_group <- unique(P$Group)
    P$PP <- ifelse(P$Sum_stat >= 2, "P", "B")
    totalP <- nrow(P[P$PP == "P",])
    geneSetSize <- nrow(P[P$Group == "Ingene_", ])

    temp <- data.frame("ID"="base", "Sum"=0, "B" = 0,
                         "P"=0, "Rate"=0)
    #
    for(i in set_group){
      n_df <- P[P$Group == i, ]
      summ <- nrow(n_df)
      B_l <- nrow(n_df[n_df$PP == "B",])
      P_l <- nrow(n_df[n_df$PP == "P",])
      rate <- P_l/summ
      res <- data.frame("ID"=i, "Sum"=summ, "B" = B_l, "P"=P_l,
                        "Rate"=rate)
      temp <- rbind(temp, res)
    }

    temp <- temp[-1,]

    temp$mad <- unique(P$mad)
    temp$alpha <- unique(P$alpha)
    temp$percentage <- unique(P$percentage)
    temp$totalP <- totalP
    temp$Size <- geneSetSize
    temp$pathway <- unique(P$pathway)
    temp$reg <- unique(P$status)


    Result <- rbind(Result, temp)
    }
  Result <- Result[-1, ]
  return(Result)
}

### Run data
output <- calc_simulate_rate(path, files)
####

## Function to create or add data to simulation database
save_simulation <- function(z=res_db){
  if(file.exists("db_sim.RData")){
    load("db_sim.RData")
    db_sim <- rbind(db_sim, z)
    save(db_sim, file = "db_sim.RData")
    print("added to simulation database")
  } else {
    db_sim <- z
    save(db_sim, file = "db_sim.RData")
    print("created simulation database with x")
  } }
## Function to process the data folder
load_process_save <- function(pat){
  filess <- dir(path = pat, pattern = ".csv")
  # Run all file
  for(i in filess){
    f_path <- paste(pat, i, sep = "")
    dat <- read.csv(f_path)
    typeU <- gsub(".*_","",i)
    typeU <- substr(typeU, 1, nchar(typeU)-4)
    setU <- gsub("_.*","",i)
    save_simulation(calc_simulate_rate(x=dat, y=typeU, z=setU)) }

}
######### Robustness functions
## Function to compute concordance across sample splits
## This is for the robustness analysis
get_concordance <- function(x1, x2, x3){
  names_set <- names(x1[-1])
  x1[-1] <- ifelse(x1[-1] >= 1, "P", "B")
  x2[-1] <- ifelse(x2[-1] >= 1, "P", "B")
  x3[-1] <- ifelse(x3[-1] >= 1, "P", "B")

  # Set Frame
  Results <- data.frame("ID"="base", "Split1_size"=0,
                        "Split2_size"=0,
                        "Shuffled_size"=0,
                        "Overlap_n1_n2"=0,
                        "Overlap_n1_sf"=0)

  #i <- names_set[1]
  for (i in names_set){
    x1p <- x1[x1[[i]] == "P",][["gene"]]
    x2p <- x2[x2[[i]] == "P",][["gene"]]
    x3p <- x3[x3[[i]] == "P",][["gene"]]
    overS1_S2 <- length(x1p[x1p %in% x2p])
    overS1_sf <- length(x1p[x1p %in% x3p])

    res <- data.frame("ID"=i, "Split1_size"=length(x1p),
                      "Split2_size"=length(x2p),
                      "Shuffled_size"=length(x3p),
                      "Overlap_n1_n2"=overS1_S2,
                      "Overlap_n1_sf"=overS1_sf)
    Results <- rbind(Results, res)}

  # Add Jaccard index
  Results <- Results[-1,]
  Results$cover_real <- Results$Overlap_n1_n2*100/Results$Split2_size
  Results$cover_shuffeled <- Results$Overlap_n1_sf*100/Results$Split2_size

  Sum_real <- Results$Split1_size + Results$Split2_size
  Sum_shuffeled <- Results$Split1_size + Results$Shuffled_size
  ##
  Union_real <- Sum_real - Results$Overlap_n1_n2
  Union_shuffeled <- Sum_shuffeled - Results$Overlap_n1_sf

  Results$Jaccard_real <- Results$Overlap_n1_n2 / Union_real
  Results$Jaccard_shuffeled <- Results$Overlap_n1_sf / Union_shuffeled

  return(Results)
}
## Function to convert wide to long for Jaccard index
get_long <- function(x,y=c("cover", "jaccard")){
  if(y == "jaccard"){
    x1 <- x[c("ID", "Jaccard_real", "Jaccard_shuffeled")]
    x1 <- reshape2::melt(x1, id.vars=c("ID"))
    names(x1)[2:3] <- c("Group", "Jaccard")
  } else {
    x1 <- x[c("ID", "cover_real", "cover_shuffeled")]
    x1 <- reshape2::melt(x1, id.vars=c("ID"))
    names(x1)[2:3] <- c("Group", "Cover")
  }
  return(x1) }
##############################################################################

##############################################################################
## Benchmarking methods for pathway activity scores
get_benchmark <- function(expr, set){
  # Function to get benchmark method
  # expr: Expression matrix as per gsva() spec
  # set: Single gene set as per ACSNI use


  method=c("gsva", "ssgsea", "zscore", "plage")
  #
  n = 0
  for (i in method){
    print(i)
    n = n + 1
    if (n == 1){
      res <- GSVA::gsva(as.matrix(X), set, method = i, verbose=TRUE)
      res <- as.data.frame(res)
      rownames(res) <- i
    } else {
      re <- GSVA::gsva(as.matrix(X), set, method = i, verbose=TRUE)
      re <- as.data.frame(re)
      rownames(re) <- i
      res <- rbind(res, re)
    }
  }

  ### Get network

  X_net <- rbind(res, X)

  out <- data.frame("gene" = rownames(X_net))
  for (i in 1:4){
    vec <- c()
    for (z in 1:nrow(X_net)){
      cor <- cor(as.numeric(X_net[i,]), as.numeric(X_net[z, ]))
      vec <- c(vec, cor)
    }
    out[method[i]] <- vec
  }
  #
  out <- out[-c(1:4), ]
  return(out)
}

{
  ## Example data preparation using AA and ATF2 case
  AA <- read.csv("~/Documents/ACSNI/2.Data/ATF4-case/AA_.csv")
  ATF2 <- read.csv("~/Documents/ACSNI/2.Data/ATF4-case/ATF2.csv")

  X <- AA
  rownames(X) <- X$gene
  X <- X[2:433]

  set <- mTOR[mTOR$PID_MTOR_4PATHWAY == 1,]["gene"]
  set <- ATF2[ATF2$PID_ATF2_PATHWAY == 1,]["gene"]
}


##############################################################################
{
  # Particular calculations
  TCGA <- read.csv("~/Desktop/mTOR_TCGA/TCGA_PID_MTOR_4PATHWAY-Y69LHUN/resultsdbsTCGA_.ptl.csv")
  GSE <- read.csv("~/Desktop/mTOR_GSE/GSE_PID_MTOR_4PATHWAY-9429F82/resultsdbsGSE_.ptl.csv")
dat <- GSE

val = c("KLF6_GSE115763", "EPAS1_GSE115389")
res <- predict_z(x=dat, y=2)
db_res <- merge(res, db, by.x = "name", by.y = "Symbol")
out <- stat_DE(x=db_res, y=val, f=0.1, fc=0.5)

method <- names(dat[2:5])
for (m in method){
  print(m)
  res <- predict_z_b(x=dat, y=m, cut=0.4)
  db_res <- merge(res, db, by.x = "name", by.y = "Symbol")
  out <- stat_DE(x=db_res, y=val, f=0.1, fc=1.5) }
}

{
 tab <- process_Chip_density("/Users/chineduanene/Documents/ACSNI/2.Data/ATF4-case/CHIP")
  AA <- read.csv("~/Desktop/ATF2_AA/AA_PID_ATF2_PATHWAY-U050TS5/resultsdbsAA_.ptl.csv")
  dat <- AA
  res <- predict_z(x=dat, y=2)
  out <- stat_chip_global(x=res, y=names(tab[2:9]), chip=tab)

  db_res <- merge(res, db, by.x = "name", by.y = "Symbol")
  out <- stat_DE(x=db_res, y=val, f=0.1, fc=0.5)

  method <- names(dat[2:5])
  n = 0
  for (m in method) {
    n = n+1
    res <- predict_z_b(x=dat, y=m, cut=0.4)

     if(n == 1) {
       out <- stat_chip_global(x=res, y=names(tab[2:9]), chip=tab)
       out$method <- m
     } else {
       ou <- stat_chip_global(x=res, y=names(tab[2:9]), chip=tab)
       ou$method <- m
       out <- rbind(out, ou) }
  }

  write.csv(out, "AA_ATF2.csv", row.names = F)

}

### Get the tables
###
{
  # Benchmark plot
  plot <- plot[plot$Model == "Real",]
  plot <- plot[plot$Method == "ACSNI",]


  ggpubr::ggbarplot(plot, x = "Gene", y = "Ratio",
            fill = "P",
            add = "mean_se",
            width = 0.6,
            size = 0.4,
            lab.size = 1,
            ylab = "Enrichment (P/B)",
            xlab = "Gene",
            palette = c("orange", rev(blues9)),
            sort.by.groups = T,
            position = position_dodge(0.9),
            ggtheme = theme_minimal(),
            lab.pos = "out"
            #x.text.angle = 90
            )


  ggpubr::ggbarplot(plot, x = "Group", y = "Coverage",
                    fill = "Group",
                    width = 0.6,
                    size = 0.4,
                    lab.size = 1,
                    ylab = "ACSNI/DE overlap (%)",
                    xlab = "Method",
                    palette = c("darkred", "grey"),
                    sort.by.groups = T,
                    position = position_dodge(0.9),
                    ggtheme = theme_minimal(),
                    lab.pos = "out"
                    #x.text.angle = 90
  )

}



#### Process data for comparing the different types of reduction
# Load the full P data
vect <- names(P)

AE <- P[1]
PCA <- P[1]
NMF <- P[1]

for(i in vect){
  if (startsWith(i, "AE_")){
    AE <- cbind(AE, P[i])
  } else if (startsWith(i, "PCA_")){
    PCA <- cbind(PCA, P[i])
  } else if (startsWith(i, "NMF_")){
    NMF <- cbind(NMF, P[i])
  }
}

# Get dimensions
d_AE <- ncol(AE[endsWith(names(AE[-1]), ".1")])
d_PCA <- ncol(PCA[endsWith(names(PCA[-1]), ".1")])
d_NMF <- ncol(NMF[endsWith(names(NMF[-1]), ".1")])

# get number of bootstraps
b_AE <- ncol(AE[-1]) / d_AE
b_PCA <- ncol(PCA[-1]) / d_PCA
b_NMF <- ncol(NMF[-1]) / d_NMF

## Process in turn
PP <- lapply(1:b_NMF, function(x){
  en <- d_NMF * x + 1
  st <- en - d_NMF + 1
  pr <- rowSums(NMF[st:en] == "P")
  pr <- ifelse(pr >= 1, 1, 0)
})
PP <- do.call(cbind, PP)
##

## Add each sum
AE$PP <- rowSums(PP)
PCA$PP <- rowSums(PP)
NMF$PP <- rowSums(PP)
##
ae <- predict_z_b(x=AE, y="PP", cut=2)
pca <- predict_z_b(x=PCA, y="PP", cut=2)
nmf <- predict_z_b(x=NMF, y="PP", cut=2)
allre <- predict_z_b(x=P, y = "Sum_stat", cut = 2)

##
val = c("KLF6_GSE115763", "EPAS1_GSE115389")
db_res <- merge(allre, db, by.x = "gene", by.y = "Symbol")
out <- stat_DE(x=db_res, y=val, f=0.1, fc=0.5)

###
tab <- process_Chip_density("/Users/chineduanene/Documents/ACSNI/2.Data/ATF4-case/CHIP")
out <- stat_chip_global(x=allre, y=names(tab[2:9]), chip=tab)
out$Method <- "all"

f_out <- out
f_out <- rbind(f_out, out)
write.csv(f_out, file = "ATF_DIMREMETH.csv", row.names = F)

# Get the overlap too
pca <- pca[c(1, 3)]
names(ae)[2] <- "AE"
overlpa <- merge(ae, pca, by = "gene")
overlpa <- merge(overlpa, nmf, by = "gene")
overlpa$suu <- rowSums(overlpa[-1] == "P")
overlpa$ae_pca <- rowSums(overlpa[2:3] == "P")
overlpa$ae_nmf <- rowSums(overlpa[c(2,4)] == "P")
overlpa$pca_nmf <- rowSums(overlpa[3:4] == "P")

# number of predicted MTOR case
ae = 1150
pca = 1
nmf = 16
all = 1166
# diagram
all = 0
ae_pca = 0
ae_nmf = 0
pca_nmf = 1


## number of predicted for mTOR case
ae = 699
pca = 237
nmf = 157
all = 912
# STAT FILE IS f_out
# diagram
all = 44
ae_pca = 69
ae_nmf = 55
pca_nmf = 101

### Generate plots
grid.newpage()
draw.triple.venn(area1 = 699, area2 = 237, area3 = 157, n12 = 69, n13 = 55, n23 = 101, n123 = 44,
                 category = c("AE", "PCA", "NMF") ,
                   cat.col = c("black", "black", "black"),
                # lty = "blank",
                 fill = NULL,
                 scaled = F)

plot <- plot[plot$Method != "all",]

# Benchmark plot
ggpubr::ggbarplot(plot, x = "TF", y = "Enrichment",
                  fill = "P",
                  #add = "mean_se",
                  width = 0.6,
                  size = 0.4,
                  lab.size = 1,
                  ylab = "Enrichment (P/B)",
                  xlab = "Gene",
                  palette = c("orange", rev(blues9)),
                  sort.by.groups = T,
                  position = position_dodge(0.9),
                  ggtheme = theme_minimal(),
                  lab.pos = "out",
                  x.text.angle = 90,
                  order = c("JUNB", "CEBPD", "NFE2L2",
                            "JUN", "IRF1", "RELA", "ERG", "EP300")
)


ggpubr::ggbarplot(plot, x = "Signal", y = "ACSNI",
                  fill = "Method",
                  #add = "mean_se",
                  width = 0.6,
                  size = 0.4,
                  lab.size = 1,
                  ylab = "Enrichment (P/B)",
                  xlab = "Gene",
                  palette = c("orange", rev(blues9)),
                  sort.by.groups = T,
                  position = position_dodge(0.9),
                  ggtheme = theme_minimal(),
                  lab.pos = "out",
                  x.text.angle = 90
)

