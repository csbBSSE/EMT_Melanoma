#Created:22.06.22
#This script contains the code to generate 2D EMT scatter plots for TCGA GBM and SARC
#for "Mapping phenotypic heterogeneity in melanoma onto the epithelial-hybrid-mesenchymal axis"
#R version 4.1.1 


require(ggplot2)
require(ggpubr)
ds <- c("GBM", "SARC")
for (i in ds) {
    df <- read.csv(paste0("Compiled Scores/",i,".csv"), row.names = 1) #Read each score sheet
    print(ggplot(df, aes(x=Epithelial, y=Mesenchymal, color=Phenotypes)) +  #Plotting 2D EMT scatterplot
        geom_point() + 
        geom_rug() + 
        theme_minimal() + 
        labs(title="Phenotypes across \nEpithelial and Mesenchymal Axes")+
        scale_color_brewer(palette="Dark2") +  
        theme(axis.title = element_text(size = 20),plot.title = element_text(size = 20),
              legend.title = element_text(size = 20),axis.text = element_text(size = 20),
              legend.text = element_text(size = 10))  + 
        border() )
    
    
}
