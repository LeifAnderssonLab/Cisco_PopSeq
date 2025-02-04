# 3_PCA_admix

## 3.1 run PCAngsd

```bash
#!/bin/bash -l
#SBATCH -A naiss2024-5-277
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 04-00:00:00
#SBATCH -J PCA_del_feg_sto
#SBATCH -e PCA_del_feg_sto_%A_%a.err
#SBATCH -o PCA_del_feg_sto_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@imbim.uu.se

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load PCAngsd/1.11

# Qiaoling Deng, 2024-10-8
######################################################################################

cd /proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/PCA_analysis/2_pca_add_mt

IN_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/PCA_analysis/2_pca_add_mt

#AdmixRun PCAngsd


###round 2
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix2  --admix --admix_K 2 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix3  --admix --admix_K 3 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix4  --admix --admix_K 4 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix5  --admix --admix_K 5 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix6  --admix --admix_K 6 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix7  --admix --admix_K 7 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix8  --admix --admix_K 8 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix9  --admix --admix_K 9 -t 8
pcangsd --beagle $IN_DIR/PCA_del_feg_sto_MinMAF_01_with_mt.beagle.gz  --out ./2_PCA_del_feg_sto_MinMAF_01_with_mt_Admix10  --admix --admix_K 10 -t 8

```

## 3.2 visualization in R --PCA

```R
##set work dir
setwd("D:/1_Cisco_remapping/7_remapping/2_ANGSD/1_PCA")
##load file for PCA
C <- as.matrix(read.table("PCA_brackish_LulA19_20_MinMAF_01_output_Admix.cov"))
e <- eigen(C)

eigen.data <- eigen(C)
eigenvectors <- as.data.frame(eigen.data$vectors)
eigenvalues <-  eigen.data$values

#eigenvectors$Sample <- bamList$SampleID
#eigenvectors$Morph_short <- bamList$Morph_short
#eigenvectors$Lake <-  bamList$Lake

#eigenvectors <- eigenvectors[!eigenvectors$Sample %in% Myv_rm_list, ]



#Get vars
pca.eigenval.sum = sum(eigen.data$values)
varPC1 <- (eigen.data$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (eigen.data$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (eigen.data$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (eigen.data$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
#We can then extract the eigenvectors from the pca object and format them into a dataframe for plotting, e.g. using ggplot().
eigenvalues <-  e$values
##load popinfo file
pop<-read.table("PCA_Brackish_LulA19_20.txt")
head(pop)
# extract pop
pop$group <- sub("^(.*?)_.*$", "\\1", pop$V1)
pop$sample <- sub("^([^_]*_[^_]*).*", "\\1", pop$V1)

values <- unique(pop$group)

library(ggplot2)
PCA_all_northern <- as.data.frame(e$vectors[, 1:2])

pca_d <- ggplot(PCA_all_northern, aes(x = V1, y = V2,color = pop$group)) +
  geom_point() +
  xlab(paste0("PC1", ": ", round(varPC1,1),"% variance")) +
  ylab(paste0("PC2", ": ", round(varPC2,1),"% variance")) +
  labs(color = " ") +
  #ggtitle("PCA Brackish cisco populations")+
  #ggtitle(" ")+
  geom_point(size = 3) +
  theme_bw() + theme(panel.grid = element_blank())+
  theme(text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size =18,face = "bold"),
        axis.text = element_text(size = 16,face = "bold"),
        panel.border = element_rect(colour = "black", size = 2),
        legend.position = "none") +  
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) + 
  scale_color_manual(values = c(#"FegA" = "#a6cee3", 
    #"FegS" = "#1f78b4", 
    "KalA" = "#e31a1c", 
    "KalK" = "#d9d9d9",
    "LulA" = "#fb9a99", 
    "PitA" = "#b2df8a", 
    "PitK" = "#fdbf6f", 
    #"Mal" = "#ff7f00",
    #"Sto" = "#6a3d9a", 
    "UleO" = "#cab2d6"
    #"VanO" = "#33a02c", 
    #"VanV" = "#b15928"
    ))#+#+
  #labs(title = "PCA_all_northern")
pca_d 
library(export)
graph2pdf(file="PCA4_Brackish_LulA19_20.pdf",width=7,height=6)

#unique(pop$group)
################################################################Done
library(patchwork)
pca_a/ pca_b / pca_c/pca_d
```

