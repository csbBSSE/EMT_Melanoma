#This script contains the code to generate plots for "Mapping phenotypic heterogeneity in melanoma onto the epithelial-hybrid-mesenchymal axis"
#R version 4.1.1 

library(ggplot2)
library(ggpubr)
library(ggExtra)

#Figure 1B-E/S1A-D
#Volcano plots
corr <- read.csv("Correlations/Correlations_PI_KS_GS_EM.csv")

EM <- c("KS","GS","E","M")
names(corr) <- gsub("KS_E","E",names(corr))
names(corr) <- gsub("KS_M","M",names(corr))
PI <- c("Ver","Hoek")
count <- data.frame()
c=1

for (i in EM) {
  for (j in PI) {
    corr_sub <- corr[grepl(i, names(corr))]  #Correlation between score i and j
    corr_sub <- corr_sub[grepl(j,names(corr_sub))]
    corr_sub <- data.frame('R'= c(corr_sub[,1],corr_sub[,3]),'p'= c(corr_sub[,2],corr_sub[,4]),
                           'Score'=c(rep(substr(names(corr_sub)[1],3,5), nrow(corr_sub)),rep(substr(names(corr_sub)[3],3,5), nrow(corr_sub)) ))
    corr_sub$Score <- factor(corr_sub$Score, levels = c("Pro","Inv"))
    corr_sub$p[corr_sub$p<(exp(-23))] <- exp(-23)
    count[c,1:4] <- c(sum(corr_sub$R > 0.36 & corr_sub$p < 0.05 & corr_sub$Score=="Pro"),  #Calculate number of significant points that are significantly 
                      sum(corr_sub$R < 0.36 & corr_sub$p < 0.05 & corr_sub$Score=="Pro"),  #positively and negatively correlated with Proliferative and 
                      sum(corr_sub$R > 0.36 & corr_sub$p < 0.05 & corr_sub$Score=="Inv"),  #invasive scores
                      sum(corr_sub$R < 0.36 & corr_sub$p < 0.05 & corr_sub$Score=="Inv"))
    rownames(count)[c] <- paste0(i,"_",j)
    c=c+1
    if(j=="Ver"){
      j="Verfaillie"
    }
    print(ggplot(corr_sub, aes(x = R, y = -log10(p), col= Score))+
            geom_point(size=3)+
            geom_hline(yintercept= 1.3, col = 'red') +
            geom_vline(xintercept = 0.3,linetype = 'dotted', col = 'blue')+
            geom_vline(xintercept = -0.3,linetype = 'dotted', col = 'blue')+
            ggtitle(paste0("Correlation between ",i, " and ",j," scores"))+
            ylab("-log10(p)")+
            xlab("R (Spearman's correlation coefficient)")+
            theme_classic()+
            theme(text =  element_text(size = 20)))
    ggsave(paste("Figures/Fig.1/Volcano_plot_",i,"_",j,".png"))
  }
  names(count) <- c("Positive_Proliferative", "Negative_Proliferative", "Positive_Invasive","Negative_Invasive")
}
write.csv(count, "Correlations/corr_counts.csv")  #Saving counts for significant points in each plot

#Box plots
for (i in list(c("KS","GS"), c("E","M"))) {
  corr_sub <- corr[,grepl(i[1],names(corr))|grepl(i[2],names(corr))] #Selecting correlations for each set of EMT scores
  for (j in c("Ver","Hoek")) {
    corr_sub1 <- corr_sub[,grepl(j,names(corr_sub))] #Selecting PI score 
    for (k in 1:4) {
      corr_sub1[corr_sub1[,(2*k)]>0.05, ((2*k)-1):(2*k)] <- NA #Removing points with p>0.05
    }
    corr_sub1 <- corr_sub1[,grepl("R",names(corr_sub1))] #Extract correlation coefficients
    corr_sub1$Dataset <- rownames(corr_sub1) 
    corr_sub1 <- reshape2::melt(corr_sub1)
    corr_sub1 <- corr_sub1[complete.cases(corr_sub1),] #Remove points with p>0.05
    corr_sub1$score <- as.factor(substr(corr_sub1$variable,3,5)) #Extracting information of PI score for each point
    n <- nchar(as.character(corr_sub1$variable[1]))
    corr_sub1$EM <- as.factor(substr(corr_sub1$variable, n-1,n)) #Extracting information on EM score for each point
    corr_sub1$EM <- gsub("_","",corr_sub1$EM)
    corr_sub1$score <- factor(corr_sub1$score, levels = c("Pro","Inv"), labels = c("Proliferative","Invasive"))
    print(ggplot(corr_sub1, aes(y = value, x = score))+
            geom_boxplot(aes(fill=score)) + facet_grid(.~EM)+
            theme_classic()+
            theme(text = element_text(size=16))+
            border()+
            rremove("legend")+
            ylab("R"))
    ggsave(paste0("Figures/Fig.1/box_",j,"_",i[1],"_",i[2],".png"))
  }
}

#Figure 2 : 2D EMT and Moving window average
#Datasets too be used for plot
ds <- c("CCLE", "GSE7127","GSE4843", "GSE65904","GSE81383", "GSE72056")
#Moving window function
mov_avg <- function(E, prol){
  win <- 0.6*max(E) + 0.4*min(E) #60% initial window
  end <- seq(win,max(E),length.out =41) #Range on end points of window
  diff <- end[2]- end[1] #Window displacement
  score <- vector()
  for (m in end) {
    s <- min(E) + ((which(end == m)-1)*diff) #rRange of start points of window
    score[which(end==m)] <- mean(prol[s<E&E<m]) #Moving average
  }
  return(score)
}

#Two-dimensional plots
for (i in ds) {
  scores <- read.csv(paste0("Compiled Scores/",i,".csv"), row.names = 1)[c(1,2,5,6)] #Extract Verfaillie PI scores and EM scores
  fig_p <- ggplot(data = scores)+   #2D-EMT plot with Verfaillie proliferative score
    geom_point(aes(x = KS_E, y = KS_M, color= Pro_Ver), size =4) +  
    scale_color_continuous(name ="Proliferative score", low = "red", high = "blue", guide = "colourbar", aesthetics = "colour") + 
    theme_classic()+
    theme(text=element_text(size = 22))+
    xlab("E score")+
    ylab("M score")
  fig_i <- ggplot(data = scores)+  #2D-EMT plot with Verfaillie invasive score
    geom_point(aes(x = KS_E, y = KS_M, color= Inv_Ver), size =4) +  
    scale_color_continuous(name ="Invasive score", low = "red", high = "blue", guide = "colourbar", aesthetics = "colour") + 
    theme_classic()+
    theme(text=element_text(size = 22))+
    xlab("E score")+
    ylab("M score")
  
  #Moving window
  PI <- read.csv(paste0("Compiled Scores/",i,".csv"), row.names = 1)[c(1,2)] #Extract Verfaillie PI scores 
  EM <- read.csv(paste0("Compiled Scores/",i,".csv"), row.names = 1)[c(5,6)] #Extract Verfaillie EM scores
  mov_scores <- data.frame()
  for (j in 1:2) {
    for (k in 1:2) {
      mov_scores[1:41,(((j-1)*2)+k)] <- mov_avg(EM[,j],PI[,k])  #Calculate moving window average for kth score among jth axis
      names(mov_scores)[(((j-1)*2)+k)] <- paste0(names(EM)[j],"_", names(PI)[k])
    }
  }
  mov_scores <- cbind(mov_scores,Window = c(1:41))
  mov_scores <- reshape2::melt(mov_scores, id.vars = "Window")
  names(mov_scores)[2:3] <- c("Comparison", "Average")
  fig_mov <- ggplot(mov_scores, aes(x = Window, y=Average, color = Comparison))+  #Plot moving window average
    geom_line(size=2)+
    scale_color_brewer(palette="Paired")+theme_classic() + 
    theme(text=element_text(size=22))+rremove("legend")
  
  plot <- ggarrange(fig_p,fig_i, fig_mov, nrow = 1, ncol=3) #Combine all 3 plots
  annotate_figure(plot, top = text_grob(paste0(i),size = 24, face = "bold"))
  ggsave(paste0("Figures/Fig.2/EM_",i,".png"), height = 5.5, width = 20)
  
  fig_e <- ggplot(data = scores)+  #E score along 2D-PI plane
    geom_point(aes(x = Pro_Ver, y = Inv_Ver, color= KS_E), size =4) +  
    scale_color_continuous(name ="E score", low = "red", high = "blue", guide = "colourbar", aesthetics = "colour") + 
    theme_classic()+
    theme(text=element_text(size = 16))+
    xlab("Proliferative score")+
    ylab("Invasive score")
  fig_m <- ggplot(data = scores)+ #M score along 2D-PI plane
    geom_point(aes(x = Pro_Ver, y = Inv_Ver, color= KS_M), size =4) +  
    scale_color_continuous(name ="M score", low = "red", high = "blue", guide = "colourbar", aesthetics = "colour") + 
    theme_classic()+
    theme(text=element_text(size = 16))+
    xlab("Proliferative score")+
    ylab("Invasive score")
  plot <- ggarrange(fig_e,fig_m) #Combine plots
  annotate_figure(plot, top = text_grob(paste0(i),size = 24, face = "bold"))
  ggsave(paste0("Figures/Fig.2/PI_",i,".png"), height = 3.5, width = 10)
}


#Figure 3: Use Classification_Mscore.py to identify number of samples in each regime of M score
#Conditional probabilities were calculated as described in methods section

#Figure 4: M score v/s Invasive score plots depicting non-monotonic increase in M score
#Datasets used for plots
ds <- c("GSE19234", "GSE158607","GSE80829" ,"GSE7127", "GSE65904","GSE101434")
#Function to assign phenotypes:
#X is a logical indicator of if a sample lies in top 10% of a score
assignment <- function(x,df){ 
  s <- rowSums(x)
  pos <- which(s>1) #Select row if more than 1 assignment, to ensure only one assignment per sample
  for (i in pos) {
    columns <- which(x[i,])
    df_sub <- scale(df)  #Scale all scores for comparison
    ind <- names(which.max(df_sub[i,columns])) #Identify phenotype with highest scaled score
    x[i,] <- FALSE  #Assign all phenotypes as FASLE
    x[i,ind] <- TRUE #Assign identified phenotype as TRUE
  }
  return(x)
}

for (i in ds) {
  scores <- read.csv(paste0("Compiled Scores/",i,".csv"),row.names = 1) #Extract Mscore and invasive scores
  tsoi <-  read.csv(paste0("Tsoi geneset scores/",i,".csv"), row.names = 1)[,c(1,3,5,7)] #Extract ssGSEA scores for Melanocytic, Transitory, NCSC, Undifferentiated genesets
  names(tsoi)[3] <-"NCSC" #Rename NCSC phenotype
  col <- apply(tsoi, 2, function(x){m <- rep(FALSE, length(x)); #Set all phenotypes as FALSE for all samples
  m[sort(x, index.return=TRUE, decreasing=TRUE)$ix[1:(length(x)/10)]] <- TRUE; return(m)}) #Identify samples with top 10% score
  assign <- assignment(col,tsoi)
  assign <- apply(assign, 1, function(x){names(which(x))}) #Assign label to each cell
  assign <- unlist(lapply(assign, function(x) if(identical(x, character(0))) "UNX" else x)) #Assign label to unlabelled cells
  scores <- scores[!assign=="UNX",] #Remove unlabelled cells
  assign <- factor(assign[assign!="UNX"], levels = c("Melanocytic","Transitory", "NCSC", "Undifferentiated"))
  p <-ggplot(data = scores)+ #Plot M score v/s Invasive score with phenotype assignment
    geom_point(aes(x = Inv_Ver, y = KS_M, color= assign), size = 4) + 
    theme_classic()+
    theme(text=element_text(size = 16))+
    labs(color="Phenotype")+
    ylab("M score")+
    xlab("Invasive score")
  
  p2 <- ggMarginal(p, groupColour = TRUE, groupFill = TRUE, margins = "x") #Plot kernel density estimates for each phenotype
  annotate_figure(p2, top = text_grob(i,size = 18, face = "bold"))
  ggsave(paste0("Figures/Fig.4/",i, ".png"))
}

#Figure 5: Metabolic scores
#For calculating scores, use scripts under Codes/Metabolism
#For generating plots use Codes/Metabolism/Volcano_plots.py




