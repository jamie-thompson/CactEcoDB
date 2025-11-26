# CactEcoDB DR statistic code

library(ape)
library(dplyr)
library(phytools)

setwd("")

tree1 <- read.tree("CactEcoDB_phylogeny_V1.nwk")
tree2 <- read.tree("CactEcoDB_phylogeny_V2_de_Vos_backbone.tre")

DR_statistic <- function(tree, return.mean = FALSE){
  rootnode <- length(tree$tip.label) + 1
  sprates <- numeric(length(tree$tip.label))
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- tree$edge.length[tree$edge[,2] == node]
      node <- tree$edge[,1][tree$edge[,2] == node]			
      qx <- qx + el* (1 / 2^(index-1))			
      index <- index + 1
    }
    sprates[i] <- 1/qx
  }
  if (return.mean){
    return(mean(sprates))		
  }else{
    names(sprates) <- tree$tip.label
    return(sprates)
  }
}

dr_1 = DR_statistic(tree1)
dr_1 = data.frame(species = names(dr_1), dr1 = dr_1)
dr_2 = DR_statistic(tree2)
dr_2 = data.frame(species = names(dr_2), dr2 = dr_2)
dr_all = full_join(dr_1, dr_2)
write.csv(dr_all,"CactEcoDB DR statistic both trees.txt")

sink("DR statistic congruence.txt")
"Correlation"
cor.test(log(dr_all$dr1), log(dr_all$dr2))
"Linear regression"
summary(lm(log(dr_all$dr2) ~ log(dr_all$dr1)))
sink()
