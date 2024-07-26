library("limma")
library("affy")
library("dplyr")
library("openxlsx")
library("ggplot2")
library("xlsx")
# ----------------------------------------------------------------------

de_analysis <- function(grps,d,d.matrix,i) {
  f = factor(grps,levels=c("control","dis"))
  design = model.matrix(~ 0 + f)
  colnames(design) = c("control","dis")
  d.fit = lmFit(d.matrix,design)
  contrast.matrix = makeContrasts(dis-control,levels=design)
  d.fit.con = contrasts.fit(d.fit,contrast.matrix)
  d.fit.eb = eBayes(d.fit.con)
  
  options(digits=2)
  tab = topTable(d.fit.eb,number=200,adjust.method="BH")
  topgenes = tab[tab[, "adj.P.Val"] < 0.05, ]
  fgenes = topgenes[topgenes[, "logFC"] >= 1 | topgenes[, "logFC"] <= -1,]
  fgenes$genes <- rownames(fgenes)
  rownames(fgenes) <- NULL
  
  topups = topgenes[topgenes[, "logFC"] >= 1, ]
  topups$genes <- rownames(topups)
  rownames(topups) <- NULL
  #write.csv(topups, paste0("topups", i,".csv"))
  topdowns = topgenes[topgenes[, "logFC"] <= -1, ]
  #write.csv(topdowns, paste0("topdowns", i,".csv"))

  return(fgenes)
}


p1 = "D:/RStudio/Comorbidity Network Analysis/GSE2171/" #HIV
d1 = ReadAffy(celfile.path=p1)
d1@phenoData@data[ ,1] = c("control1","control2","control3",
                   "control4","control5","control6",
                   "control7","control8","control9",
                   "control10","control11","control12",
                   "dis1","dis2","dis3", "dis4","dis5","dis6",
                   "dis7","dis8","dis9","dis10")
# RMA Normalization
d1.rma = affy::rma(d1)
d1.matrix = exprs(d1.rma)

# DE Analysis
d1@phenoData@data[ ,1] = c("control","control","control",
                   "control","control","control",
                   "control","control","control",
                   "control","control","control",
                   "dis","dis","dis", "dis","dis","dis",
                   "dis","dis","dis","dis")

colnames(d1@phenoData@data)[1]= "source"
t1 = de_analysis(d1@phenoData@data$source,d1,d1.matrix,1)

# -----------------------------------------------------------------------

p2 = "D:/RStudio/Comorbidity Network Analysis/GSE8823/" #COPD
d2 = ReadAffy(celfile.path=p2)
d2@phenoData@data[ ,1] = c("control1","control2","control3",
                           "control4","control5","control6",
                           "control7","control8","control9",
                           "control10","control11","dis1",
                           "dis2","dis3","dis4","dis5","dis6",
                           "dis7","dis8","dis9","dis10","dis11",
                           "dis12","dis13")
# RMA Normalization
d2.rma = affy::rma(d2)
d2.matrix = exprs(d2.rma)

# DE Analysis
d2@phenoData@data[ ,1] = c("control","control","control",
                           "control","control","control",
                           "control","control","control",
                           "control","control","dis",
                           "dis","dis","dis","dis","dis",
                           "dis","dis","dis","dis","dis",
                           "dis","dis")

colnames(d2@phenoData@data)[1]= "source"
t2 = de_analysis(d2@phenoData@data$source,d2,d2.matrix,2)

# -----------------------------------------------------------------------

p3 = "D:/RStudio/Comorbidity Network Analysis/GSE25724/" #Diabetes
d3 = ReadAffy(celfile.path=p3)
d3@phenoData@data[ ,1] = c("control1","control2","control3",
                          "control4","control5","control6",
                          "control7","dis1","dis2","dis3",
                          "dis4","dis5","dis6")
# RMA Normalization
d3.rma = affy::rma(d3)
d3.matrix = exprs(d3.rma)

# DE Analysis
d3@phenoData@data[ ,1] = c("control","control","control",
                           "control","control","control",
                           "control","dis","dis","dis",
                           "dis","dis","dis")

colnames(d3@phenoData@data)[1]= "source"
t3 = de_analysis(d3@phenoData@data$source,d3,d3.matrix,3)


# -----------------------------------------------------------------------

p4 = "D:/RStudio/Comorbidity Network Analysis/GSE54992/" #TB
d4 = ReadAffy(celfile.path=p4)
d4@phenoData@data[ ,1] = c("dis1","dis2","dis3", "dis4","dis5","dis6",
                           "dis7","dis8","dis9","dis10", "control1",
                           "control2","control3","control4","control5",
                           "control6")
# RMA Normalization
d4.rma = affy::rma(d4)
d4.matrix = exprs(d4.rma)

# DE Analysis
d4@phenoData@data[ ,1] = c("dis","dis","dis", "dis","dis","dis","dis",
                           "dis","dis","dis","control","control","control",
                          "control","control","control")

colnames(d4@phenoData@data)[1]= "source"
t4 = de_analysis(d4@phenoData@data$source,d4,d4.matrix,4)

# -----------------------------------------------------------------------

p5 = "D:/RStudio/Comorbidity Network Analysis/GSE1402/" #RA
d5 = ReadAffy(celfile.path=p5)
d5@phenoData@data[ ,1] = c("dis1","dis2","dis3", "dis4","dis5","dis6",
                            "dis7","dis8","dis9","dis10","dis11","dis12",
                            "dis13", "dis14","dis15","dis16","dis17","dis18",
                            "dis19","dis20","dis21","dis22","dis23", "dis24",
                            "dis25","dis26","control1","control2","control3",
                            "control4","control5","control6","control7",
                            "control8","control9","control10","control11",
                            "dis26","dis27","dis28", "dis29","dis30","dis31",
                            "dis32","dis33","dis34","dis35","dis36","dis37",
                            "dis38", "dis39","dis40","dis41","dis42","dis43",
                            "dis44","dis45")

# RMA Normalization
d5.rma = affy::rma(d5)
d5.matrix = exprs(d5.rma)

# DE Analysis
d5@phenoData@data[ ,1] = c("dis","dis","dis", "dis","dis","dis",
                            "dis","dis","dis","dis","dis","dis",
                            "dis", "dis","dis","dis","dis","dis",
                            "dis","dis","dis","dis","dis", "dis",
                            "dis","dis","control","control","control",
                            "control","control","control","control",
                            "control","control","control","control",
                            "dis","dis","dis", "dis","dis","dis",
                            "dis","dis","dis","dis","dis","dis",
                            "dis", "dis","dis","dis","dis","dis",
                            "dis","dis")

colnames(d5@phenoData@data)[1]= "source"
t5 = de_analysis(d5@phenoData@data$source,d5,d5.matrix,5)

# -----------------------------------------------------------------------

p6 = "D:/RStudio/Comorbidity Network Analysis/GSE10072/" #Lung Cancer
d6 = ReadAffy(celfile.path=p6)
d6@phenoData@data[ ,1] = c("dis1","control1","dis2","control2","dis3","dis4","dis5","control3","dis6",
                            "control4","control5","dis7","dis8","control6","dis9","control7","dis10","dis11",
                            "control8","control9","dis12",'control10',"dis13","dis14","control11","dis15","control12",
                            "dis16","control13","dis17","control14","dis18","dis19","control15","dis20","control16",
                            "dis21","control17","dis22","dis23","control18","dis24","control19","dis25","control20",
                            "dis26","control21","dis27","control22","dis28","dis29","control23","control24","dis30",
                            "control25","dis31","control26","dis32","control27","dis33","control28","dis34","dis35",
                            "dis36","control29","dis37","control30","dis38","control31","dis39","control32","dis40",
                            "dis41","dis42","control33","dis43","dis44","control34","control35","dis45","dis46",
                            "control36","dis47","control37","dis48","control38","control39","control40","control41",
                            "dis49","control42","dis50","control43","dis51","control44","dis52","dis53","dis54",
                            "control45","dis55","control46","dis56","control47","dis57","dis58","control48","control49")
# RMA Normalization
d6.rma = affy::rma(d6)
d6.matrix = exprs(d6.rma)

# DE Analysis
d6@phenoData@data[ ,1] = c("dis","control","dis","control","dis","dis","dis","control","dis",
                            "control","control","dis","dis","control","dis","control","dis","dis",
                            "control","control","dis",'control',"dis","dis","control","dis","control",
                            "dis","control","dis","control","dis","dis","control","dis","control",
                            "dis","control","dis","dis","control","dis","control","dis","control",
                            "dis","control","dis","control","dis","dis","control","control","dis",
                            "control","dis","control","dis","control","dis","control","dis","dis",
                            "dis","control","dis","control","dis","control","dis","control","dis",
                            "dis","dis","control","dis","dis","control","control","dis","dis",
                            "control","dis","control","dis","control","control","control","control",
                            "dis","control","dis","control","dis","control","dis","dis","dis",
                            "control","dis","control","dis","control","dis","dis","control","control")

colnames(d6@phenoData@data)[1]= "source"
t6 = de_analysis(d6@phenoData@data$source,d6,d6.matrix,6)

# take gene names from these and convert from DAVID

write.csv(t1, paste0("HIV",".csv"), row.names = FALSE)
write.csv(t2, paste0("COPD",".csv"), row.names = FALSE)
write.csv(t3, paste0("Diabetes",".csv"), row.names = FALSE)
write.csv(t4, paste0("TB",".csv"), row.names = FALSE)
write.csv(t5, paste0("RA",".csv"), row.names = FALSE)
write.csv(t6, paste0("Lung Cancer",".csv"), row.names = FALSE)


# -----------------------------------------------------------------------

# gProfileR2 --> convert affy gene ids to gene names
# library("gprofiler2")
#
# t1c = gconvert(t1$genes, organism = "hsapiens", target = "HGNC")
# t1c = subset(t1c %>% group_by(input_number) %>% filter(row_number()==1), select = c(input, name)) 
# colnames(t1c)[which(names(t1c) == "input")] <- "genes"
# 
# t2c = gconvert(t2$genes, organism = "hsapiens", target = "HGNC")
# t2c = subset(t2c %>% group_by(input_number) %>% filter(row_number()==1), select = c(input, name))
# colnames(t2c)[which(names(t2c) == "input")] <- "genes"
# 
# t3c = gconvert(t3$genes, organism = "hsapiens", target = "HGNC")
# t3c = subset(t3c %>% group_by(input_number) %>% filter(row_number()==1), select = c(input, name))
# colnames(t3c)[which(names(t3c) == "input")] <- "genes"
# 
# t4c = gconvert(t4$genes, organism = "hsapiens", target = "HGNC")
# t4c = subset(t4c %>% group_by(input_number) %>% filter(row_number()==1), select = c(input, name))
# colnames(t4c)[which(names(t4c) == "input")] <- "genes"
# 
# t5c = gconvert(t5$genes, organism = "hsapiens", target = "HGNC")
# t5c = subset(t5c %>% group_by(input_number) %>% filter(row_number()==1), select = c(input, name))
# colnames(t5c)[which(names(t5c) == "input")] <- "genes"
# 
# t6c = gconvert(t6$genes, organism = "hsapiens", target = "HGNC")
# t6c = subset(t6c %>% group_by(input_number) %>% filter(row_number()==1), select = c(input, name))
# colnames(t6c)[which(names(t6c) == "input")] <- "genes"

# -----------------------------------------------------------------------

# read the converted gene names from DAVID and merge with parameter dataframe

t1c = read.csv("Converted from DAVID/HIV.csv")
t1c = subset(t1c %>% group_by(genes) %>% filter(row_number()==1), select = c(genes, name))

t2c = read.csv("Converted from DAVID/COPD.csv")
t2c = subset(t2c %>% group_by(genes) %>% filter(row_number()==1), select = c(genes, name))

t3c = read.csv("Converted from DAVID/Diabetes.csv")
t3c = subset(t3c %>% group_by(genes) %>% filter(row_number()==1), select = c(genes, name))

t4c = read.csv("Converted from DAVID/TB.csv")
t4c = subset(t4c %>% group_by(genes) %>% filter(row_number()==1), select = c(genes, name))

t5c = read.csv("Converted from DAVID/RA.csv")
t5c = subset(t5c %>% group_by(genes) %>% filter(row_number()==1), select = c(genes, name))

t6c = read.csv("Converted from DAVID/Lung Cancer.csv")
t6c = subset(t6c %>% group_by(genes) %>% filter(row_number()==1), select = c(genes, name))


t1 <- merge(t1, t1c, by = c("genes"))
t1$type <- ifelse(t1$logFC >= 1, 'upregulated', 'downregulated')

t2 <- merge(t2, t2c, by = c("genes"))
t2$type <- ifelse(t2$logFC >= 1, 'upregulated', 'downregulated')

t3 <- merge(t3, t3c, by = c("genes"))
t3$type <- ifelse(t3$logFC >= 1, 'upregulated', 'downregulated')

t4 <- merge(t4, t4c, by = c("genes"))
t4$type <- ifelse(t4$logFC >= 1, 'upregulated', 'downregulated')

t5 <- merge(t5, t5c, by = c("genes"))
t5$type <- ifelse(t5$logFC >= 1, 'upregulated', 'downregulated')

t6 <- merge(t6, t6c, by = c("genes"))
t6$type <- ifelse(t6$logFC >= 1, 'upregulated', 'downregulated')


# -----------------------------------------------------------------------

# list of all DEGS in each dataset for common genes and venn diagram

write.csv(t1, paste0("HIV",".csv"), row.names = FALSE)
write.csv(t2, paste0("COPD",".csv"), row.names = FALSE)
write.csv(t3, paste0("Diabetes",".csv"), row.names = FALSE)
write.csv(t4, paste0("TB",".csv"), row.names = FALSE)
write.csv(t5, paste0("RA",".csv"), row.names = FALSE)
write.csv(t6, paste0("Lung Cancer",".csv"), row.names = FALSE)

# -------------------------------------------------------------------------

# graph for up and downregulated genes

t1$disease <- 'HIV'
t2$disease <- 'COPD'
t3$disease <- 'Diabetes'
t4$disease <- 'TB'
t5$disease <- 'RA'
t6$disease <- 'Lung Cancer'

brgph = rbind(t1,t2,t3,t4,t5)
ggplot(brgph, aes(factor(disease), fill = factor(type))) + 
  geom_bar(position = position_dodge(0.9)) +
  labs(y = "Frequency", x = "Comorbidities of Tuberculosis", fill = "Type") +
  ggtitle("Comparison of Differentially Expressed Genes")

# write.xlsx(brgph, file = 'updown.xlsx')

# -------------------------------------------------------------------------

# list of all distinct DEGs from 14,24,34,54,64 for STRING network analysis

t14 = rbind(t1,t4)
deg14 = distinct(t14, name)
t24 = rbind(t2,t4)
deg24 = distinct(t24, name)
t34 = rbind(t3,t4)
deg34 = distinct(t34, name)
t54 = rbind(t5,t4)
deg54 = distinct(t54, name)
t64 = rbind(t6,t4)
deg64 = distinct(t64, name)

distinct <- list(deg14,deg24,deg34,deg64,deg54)
write.xlsx(distinct, file = 'Distinct DEGs in Each Dataset Combination.xlsx')

# -------------------------------------------------------------------------

# common genes in each dataset and distinct genes respectively for Cytoscape
# use net1 files in python to mark the node tables

dis1 <- t1$name
lung1 <- t4$name
s16 <- base::intersect(t1$name, t4$name)
dis1 <- dis1[! dis1 %in% s16]
lung1 <- lung1[! lung1 %in% s16]
dn1 <- list('dis' = dis1, 'common' = s16, 'tb' = lung1)
# write.xlsx(dn1, file = 'net1 - hiv.xlsx')


dis2 <- t2$name
lung2 <- t4$name
s26 <- base::intersect(t2$name, t4$name)
dis2 <- dis2[! dis2 %in% s26]
lung2 <- lung2[! lung2 %in% s26]
dn2 <- list('dis' = dis2, 'common' = s26, 'tb' = lung2)
# write.xlsx(dn2, file = 'net2 - copd.xlsx')


dis3 <- t3$name
lung3 <- t4$name
s36 <- base::intersect(t3$name, t4$name)
dis3 <- dis3[! dis3 %in% s36]
lung3 <- lung3[! lung3 %in% s36]
dn3 <- list('dis' = dis3, 'common' = s36, 'tb' = lung3)
# write.xlsx(dn3, file = 'net3 - diabetes.xlsx')


dis4 <- t6$name
lung4 <- t4$name
s46 <- base::intersect(t6$name, t4$name)
dis4 <- dis4[! dis4 %in% s46]
lung4 <- lung4[! lung4 %in% s46]
dn4 <- list('dis' = dis4, 'common' = s46, 'tb' = lung4)
# write.xlsx(dn4, file = 'net4 - lung cancer.xlsx')


dis5 <- t5$name
lung5 <- t4$name
s56 <- base::intersect(t5$name, t4$name)
dis5 <- dis5[! dis5 %in% s56]
lung5 <- lung5[! lung5 %in% s56]
dn5 <- list('dis' = dis5, 'common' = s56, 'tb' = lung5)
# write.xlsx(dn5, file = 'net5 - ra.xlsx')

# ----------------------------------------------------------------------------------

# compute MCI for common genes

# mci <- data.frame (Disease = c("HIV", "COPD", "Diabetes", "Lung Cancer", "RA"),
#                    Total = c(nrow(t1), nrow(t2), nrow(t3), nrow(t6), nrow(t5)),
#                    Common = c(length(s16), length(s26), length(s36), length(s46), length(s56)))
# ggplot(mci) + geom_col(aes(x = Disease, y = Total, fill="Total number of genes"), size = 1, width=0.4) + 
#   ylim(-10, 200)+ scale_fill_manual(name = NULL, values = c("Total number of genes" = "steelblue")) + 
#   geom_line(aes(x = Disease, y = Common, color="MCI"), size = 1.5, group = 1) + scale_y_continuous(trans= 'log10')  +
#   scale_color_manual(name = NULL, values = c("MCI" = "red")) + labs(y = "Frequency") + ggtitle("Molecular Comorbidity Index")




mci <- data.frame (Disease = c("HIV", "COPD", "Diabetes"),
                   Total = c(nrow(t1), nrow(t2), nrow(t3)),
                   Common = c(length(s16), length(s26), length(s36)))
ggplot(mci) + geom_col(aes(x = Disease, y = Total, fill="Total number of genes"), size = 1, width=0.4) + 
  ylim(-10, 200)+ scale_fill_manual(name = NULL, values = c("Total number of genes" = "steelblue")) + 
  geom_line(aes(x = Disease, y = Common, color="MCI"), size = 1.5, group = 1) + scale_y_continuous(trans= 'log10')  +
  scale_color_manual(name = NULL, values = c("MCI" = "red")) + labs(y = "Frequency") + ggtitle("Molecular Comorbidity Index")