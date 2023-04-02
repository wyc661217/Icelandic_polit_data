#analyze the Icelandic pilot data
rm(list=ls())
setwd("~/Desktop/Iceland/reads_statis/")
library(readxl)
library(ggplot2)

#QC report 1
##functions for extracting CGG ID
CGGID_extract = function(df,colNO){
  
  for (i in 1:dim(df)[1]) {
    df[i,colNO] = paste0("CGG3_",strsplit(strsplit(df[i,colNO], split="-")[[1]][4], split="_")[[1]][1])
  }
  
  df = df[-which(df[,colNO] == "CGG3_NA"),]
  df = df[order(df[,colNO]),]
  return(df)
  
}

##functions for merging the 4 lanes CGG ID
Mgerge_4lanes = function(df){
  
  df$raw_reads = as.numeric(df$raw_reads) *2
  df$fastp_input = as.numeric(df$fastp_input) * 1000000
  df$fastp_after_filtering = as.numeric(df$fastp_after_filtering) * 1000000
  df$fastp_passed_filters = as.numeric(df$fastp_passed_filters) * 1000000
  df$fastp_too_short = as.numeric(gsub("%", "", df$fastp_too_short))
  df$fastp_low.complexity = as.numeric(gsub("%", "", df$fastp_low.complexity))
  
  df_new = as.data.frame(matrix(nrow = dim(df)[1]/4, ncol = dim(df)[2]))
  colnames(df_new) = colnames(df)
  
  for (i in seq(1,dim(df)[1],4)) {
    j = (i-1)/4 +1
    df_new[j,1] = df[i,1]
    df_new[j,2] = sum(df[c(i:(i+3)),2])
    df_new[j,3] = sum(df[c(i:(i+3)),3])
    df_new[j,4] = sum(df[c(i:(i+3)),4])
    df_new[j,5] = sum(df[c(i:(i+3)),5])
    df_new[j,6] = sum(df[c(i:(i+3)),6])
    df_new[j,7] = sum(df[c(i:(i+3)),4] * df[c(i:(i+3)),7]) / sum(df[c(i:(i+3)),4])
    df_new[j,8] = sum(df[c(i:(i+3)),4] * df[c(i:(i+3)),8]) / sum(df[c(i:(i+3)),4])
  }
  
  df_new = CGGID_extract(df_new,1)
  return(df_new)
  
}

#functon for merge two flowcells
Mgerge_2flowcells=function(df1,df2,x,y,z){
  
  if (all(df1[,1] == df2[,1])) {
    
    df_new = df1
    
    for (i in x) {
      df_new[,i] = df1[,i] + df2[,i]
    }
    
    for (i in y) {
      df_new[,i] = (df1[,i]*df1[,z] + df2[,i]*df2[,z]) / (df1[,z] + df2[,z])
    }
    
    return(df_new)
  }
}


qc1_28_raw = read.csv("QC_report1_028.csv")
qc1_28 = Mgerge_4lanes(qc1_28_raw)
qc1_29_raw = read.csv("QC_report1_029.csv")
qc1_29 = Mgerge_4lanes(qc1_29_raw[-1,])
rm(list=(c("qc1_28_raw","qc1_29_raw")))
qc1 = Mgerge_2flowcells(qc1_28,qc1_29,c(2:6),c(7,8),4)

sum(qc1$raw_reads)
mean(qc1$fastp_too_short)


#clean fq report
cleaned_28 = read.table("~/Desktop/Iceland/reads_statis/lib028.cleaned.statreport.txt", header=TRUE, quote="", comment.char="")
cleaned_28 = CGGID_extract(df=cleaned_28,colNO=1)
cleaned_28 = cleaned_28[,c(1,4,7)]
cleaned_28$num_seqs = as.numeric(gsub(",", "", cleaned_28$num_seqs))

cleaned_29 = read.table("~/Desktop/Iceland/reads_statis/lib029.cleaned.statreport.txt", header=TRUE, quote="", comment.char="")
cleaned_29 = CGGID_extract(df=cleaned_29,colNO=1)
cleaned_29 = cleaned_29[,c(1,4,7)]
cleaned_29$num_seqs = as.numeric(gsub(",", "", cleaned_29$num_seqs))

all(cleaned_28$file == cleaned_29$file)
cleaned = cleaned_28
cleaned$num_seqs = cleaned_28$num_seqs + cleaned_29$num_seqs
cleaned$avg_len = (cleaned_28$num_seqs*cleaned_28$avg_len + cleaned_29$num_seqs*cleaned_29$avg_len)/cleaned$num_seqs


#mapped reads
mapped_28 = read.csv("classified_read_028.txt",header = F)
mapped_28 = mapped_28[,1:2]
mapped_28 = mapped_28[order(mapped_28$V1),]

mapped_29 = read.csv("classified_read_029.txt",header = F)
mapped_29 = mapped_29[,1:2]
mapped_29 = mapped_29[-which(!mapped_29$V1 %in% mapped_28$V1),]
mapped_29 = mapped_29[order(mapped_29$V1),]

all(mapped_28$V1 == mapped_29$V1)
mapped = mapped_28
mapped$V2 = mapped_28$V2 + mapped_29$V2
mapped$V2 = mapped$V2 / 1000000
sum(mapped$V2)


#two flowcelled merged
d30 = read.table("2flowcells_merged.complexfilter30.txt",header = T)
d30 = d30[-seq(2,dim(d30)[1],2),c(1,4,5)]
d30$num_seqs = as.numeric(gsub(",", "", d30$num_seqs))
d30$sum_len = as.numeric(gsub(",", "", d30$sum_len))

d3 = read.table("2flowcells_merged.complexfilter3.txt",header = T)
d3 = d3[-seq(2,dim(d3)[1],2),c(1,4,5)]
d3$num_seqs = as.numeric(gsub(",", "", d3$num_seqs))
d3$sum_len = as.numeric(gsub(",", "", d3$sum_len))

dm = read.table("2flowcells_merged.statreport.txt",header = T)
dm = dm[-seq(2,dim(dm)[1],2),c(1,4,5)]
dm$num_seqs = as.numeric(gsub(",", "", dm$num_seqs))
dm$sum_len = as.numeric(gsub(",", "", dm$sum_len))

duplicated(dm$file)
dm = dm[order(dm$file),]
dm.1 = cleaned[which(cleaned$file %in% dm$file),]
dm.2 = cleaned_28[which(cleaned_28$file %in% dm$file),]

all(dm$file == dm.1$file)
all(dm$file == dm.2$file)

dm = cbind(dm,dm.1,dm.2)
dm = dm[,c(1,2,5,8)]
dm$new_rate = (dm$num_seqs - dm$num_seqs.2) / (dm$num_seqs.1 - dm$num_seqs.2) *100



#sample metadata
mdata1 = read_excel("~/Desktop/Iceland/sample/metadata/Tjørnen_Samp_221011.xlsx")
mdata1 = mdata1[-c(65,129),c(1,4,17,19)]
colnames(mdata1) = c("sample","depth","age","site")
mdata1$depth = round(mdata1$depth)

m_age = read.delim("~/Desktop/Iceland/sample/metadata/TJON23_85_ages.txt")
m_age = m_age[,c(1,5)]

for (i in 1:dim(mdata1)[1]) {
  mdata1$age[i] = m_age$mean[which(m_age$depth == mdata1$depth[i])]
}

mdata1$age = as.numeric(mdata1$age)
mdata1 = mdata1[which(mdata1$sample %in% cleaned_28$file),]
mdata1 = as.data.frame(mdata1)
mdata1[17,c(2,3)] = c(147,694)
mdata1[32,c(2,3)] = c(161,772)
which(duplicated(mdata1$age))
mdata1$site = "Tjornen"

mdata2 = read_excel("~/Desktop/Iceland/sample/metadata/Wesley_CGG_Database_Sediment_Template_220809_PSO_WRF.xlsx")
mdata2 = mdata2[1:52,c(1,2,17,19)]
colnames(mdata2) = c("sample","depth","age","site")
mdata2$depth = sub("ISL22-20.7_","",mdata2$depth)
mdata2$depth = sub("ISL22-20.2A_","",mdata2$depth)
mdata2$depth = sub("ISL22-20.6_","",mdata2$depth)
mdata2$age[1:20] = "<1 ka"
mdata2$age[21:30] = "1-2 ka"
mdata2$age[31:52] = "~10 ka"
mdata2$site = "Borgarvatnet"

mdata = rbind(mdata1,mdata2)
all(cleaned_28$file %in% mdata$sample)
mdata = mdata[order(mdata$sample),]
rm(list="m_age")


#correlate the lab factors
lab_lib = read_excel("lib_building.xlsx")
lab_lib_28 = lab_lib[lab_lib$`Library Plate ID` == "eDNALib028",]
lab_lib_28 = lab_lib_28[order(lab_lib_28$`Sample ID (CGG No)`),]
lab_lib_29 = lab_lib[lab_lib$`Library Plate ID` == "eDNALib029",]
lab_lib_29 = lab_lib_29[order(lab_lib_29$`Sample ID (CGG No)`),]

lab_ext = read_excel("extraction.xlsx")
lab_ext = lab_ext[order(lab_ext$`Sample ID (CGG No)`),]
lab_ext = lab_ext[1:92,]

all(lab_lib_28$`Sample ID (CGG No)` == lab_lib_29$`Sample ID (CGG No)`)
all(lab_ext$`Sample ID (CGG No)` == lab_lib_29$`Sample ID (CGG No)`)



#kranken taxa profile
taxa_list = read.delim("~/Desktop/Iceland/taxa_list/taxa_list/niis_genus_and_species_list.with_taxaID.tsv", header=FALSE)
taxa_list$taxa = paste(taxa_list$V2,taxa_list$V1,taxa_list$V3,sep = ":")
taxa_list = taxa_list[,c(2,4)]

list_28 = dir("~/Desktop/Iceland/kraken_results/niis_restuls_028/")
taxa_list_28 = taxa_list

for (i in 1:length(list_28)) {
  
  taxa_list_28$new = 0
  df1 = read.table(paste0("~/Desktop/Iceland/kraken_results/niis_restuls_028/",list_28[i]), quote="\"", comment.char="")
  
  for (j in 1:dim(df1)[1]) {
    
    taxa_list_28$new[which(taxa_list_28$V2 == df1$V1[j])] =  df1[j,2]
    
  }
  colnames(taxa_list_28)[i+2] = gsub(".niis.kmer_freq.txt","",list_28[i])
}


list_29 = dir("~/Desktop/Iceland/kraken_results/niis_restuls_029/")
taxa_list_29 = taxa_list
for (i in 1:length(list_29)) {
  
  taxa_list_29$new = 0
  df1 = read.table(paste0("~/Desktop/Iceland/kraken_results/niis_restuls_029/",list_29[i]), quote="\"", comment.char="")
  
  for (j in 1:dim(df1)[1]) {
    
    taxa_list_29$new[which(taxa_list_29$V2 == df1$V1[j])] =  df1[j,2]
    
  }
  colnames(taxa_list_29)[i+2] = gsub(".niis.kmer_freq.txt","",list_29[i])
}

colnames(taxa_list_28) == colnames(taxa_list_29)
taxa_list_29 = taxa_list_29[,-7]
all(colnames(taxa_list_28) == colnames(taxa_list_29))

taxa_list = taxa_list_28
for (i in 3:93) {
  taxa_list[,i] = taxa_list_28[,i] + taxa_list_29[,i]
}

taxa_list = taxa_list[,-1]
taxa_list_1 = taxa_list[,c(1,which(colnames(taxa_list) %in% mdata1$sample))]
taxa_list_2 = taxa_list[,c(1,which(colnames(taxa_list) %in% mdata2$sample))]
taxa_list_1 = taxa_list_1[which(rowSums(taxa_list_1[,-1]) != 0),]
taxa_list_2 = taxa_list_2[which(rowSums(taxa_list_2[,-1]) != 0),]

write.table(taxa_list_1,"lake1.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(taxa_list_2,"lake2.txt",sep = "\t",col.names = T,row.names = F,quote = F)








#figure1
cleaned$file == qc1$sample
fd1 = cbind(qc1,cleaned)
mean(fd1$fastp_too_short)
mean(fd1$fastp_low.complexity)
fd1$fastp_too_short = fd1$fastp_too_short * fd1$fastp_input / 100
fd1$fastp_too_short = fd1$fastp_input - fd1$fastp_too_short
fd1$fastp_low.complexity = round(fd1$fastp_low.complexity * fd1$fastp_input / 100)
fd1 = fd1[,c(1,2,3,7,6,10)]
fd1[,-1] = fd1[,-1]/1000000

fd1.1 = data.frame(
  Step = c("Raw reads","Adepter removed","Adepter residue removed","Low complexity 1 controlled","Final QCed","Mapped"),
  'Reads number' = c(mean(fd1[,2]),mean(fd1[,3]),mean(fd1[,4]),mean(fd1[,5]),mean(fd1[,6]),mean(mapped$V2)),
  error = c(qnorm(0.95) * sd(fd1[,2]) / sqrt(length(fd1[,2])),
            qnorm(0.95) * sd(fd1[,3]) / sqrt(length(fd1[,3])),
            qnorm(0.95) * sd(fd1[,4]) / sqrt(length(fd1[,4])),
            qnorm(0.95) * sd(fd1[,5]) / sqrt(length(fd1[,5])),
            qnorm(0.95) * sd(fd1[,6]) / sqrt(length(fd1[,6])),
            qnorm(0.95) * sd(mapped$V2) / sqrt(length(mapped$V2))))


fd1.1$Step = factor(fd1.1$Step,ordered = T,c("Raw reads","Adepter removed","Adepter residue removed","Low complexity 1 controlled","Final QCed","Mapped"))

fp1=ggplot(fd1.1, aes(x = Step, y = Reads.number)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Reads.number - error, ymax = Reads.number + error), width = 0.4, position = position_dodge(0.9))+
                  labs(x = element_blank(), y = "Reads number (million)") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))

ggsave(plot=fp1, height=6, width=10, dpi=300,
       filename="figure1.pdf",
       useDingbats=FALSE)

                  

#figure2



#figure3
all(qc1_28$sample == cleaned_28$file)
all(lab_lib_28$`Sample ID (CGG No)` == cleaned_28$file)
all(lab_lib_28$`Sample ID (CGG No)` == mdata$sample)
fd3.28 = cbind(qc1_28[,c(1,2)],cleaned_28[,c(2,3)],lab_ext$`Sample Mass (g)`,lab_lib_28[,c(2,3,4)],mdata[,c(2,3,4)])
fd3.28$raw_reads = fd3.28$raw_reads / 1000000
fd3.28$num_seqs = fd3.28$num_seqs  / 1000000
fd3.28$percentage_passed_QC = fd3.28$num_seqs / fd3.28$raw_reads

all(qc1_29$sample == cleaned_29$file)
all(lab_lib_29$`Sample ID (CGG No)` == cleaned_29$file)
all(lab_lib_29$`Sample ID (CGG No)` == mdata$sample)
fd3.29 = cbind(qc1_29[,c(1,2)],cleaned_29[,c(2,3)],lab_ext$`Sample Mass (g)`,lab_lib_29[,c(2,3,4)],mdata[,c(2,3,4)])
fd3.29$raw_reads = fd3.29$raw_reads / 1000000
fd3.29$num_seqs = fd3.29$num_seqs  / 1000000
fd3.29$percentage_passed_QC = fd3.29$num_seqs / fd3.29$raw_reads

fd3 = rbind(fd3.28,fd3.29)
colnames(fd3)[5] = "Sample Mass"

fp3.1=ggplot(fd3, aes(x = fd3$`Library Concentration (nM)`, y = fd3$percentage_passed_QC)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)+
  labs(x = "Library Concentration (nM)", y = "Percentage passed QC") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))

fp3.2=ggplot(fd3, aes(x = fd3$Ct, y = fd3$percentage_passed_QC)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)+
  labs(x = "qPCR Ct value", y = "Percentage passed QC") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))


fp3.3=ggplot(fd3, aes(x = fd3$Ct, y = fd3$raw_reads)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)+
  labs(x = "qPCR Ct value", y = "Raw reads number (million)") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))

fp3.4=ggplot(fd3, aes(x = fd3$`Library Concentration (nM)`, y = fd3$raw_reads)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)+
  labs(x = "Library Concentration (nM)", y = "Raw reads number (million)") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))


ggsave(plot=fp3.1, height=6, width=9, dpi=300,
       filename="figure3.1.pdf",
       useDingbats=FALSE)
ggsave(plot=fp3.2, height=6, width=9, dpi=300,
       filename="figure3.2.pdf",
       useDingbats=FALSE)

ggsave(plot=fp3.3, height=6, width=9, dpi=300,
       filename="figure3.3.pdf",
       useDingbats=FALSE)
ggsave(plot=fp3.4, height=6, width=9, dpi=300,
       filename="figure3.4.pdf",
       useDingbats=FALSE)

#figure4
fd4.1 = fd3[which(fd3$site == "Tjornen"),]
fd4.1$depth = as.numeric(fd4.1$depth )
fp4.1=ggplot(fd4.1, aes(x = depth, y = Ct)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)+
  labs(x = "Sample depth (cm)", y = "qPCR Ct value") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))

fd4.2 = fd3[which(fd3$site == "Borgarvatnet"),]
fd4.2$depth = as.numeric(fd4.2$depth )
fp4.2=ggplot(fd4.2, aes(x = depth, y = Ct)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)+
  labs(x = "Sample depth (cm)", y = "qPCR Ct value") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))

ggsave(plot=fp4.1, height=6, width=9, dpi=300,
       filename="figure4.1.pdf",
       useDingbats=FALSE)
ggsave(plot=fp4.2, height=6, width=9, dpi=300,
       filename="figure4.2.pdf",
       useDingbats=FALSE)


#figure 5
fd5.28 = fd3.28
fd5.28 = fd5.28[which(fd5.28$sample %in% mapped_28$V1),]
all(fd5.28$sample == mapped_28$V1)
fd5.28 = cbind(fd5.28,mapped_28)
fd5.28$mapped_rate = fd5.28$V2 / fd5.28$num_seqs / 1000000

fd5.29 = fd3.29
fd5.29 = fd5.29[which(fd5.29$sample %in% mapped_29$V1),]
all(fd5.29$sample == mapped_29$V1)
fd5.29 = cbind(fd5.29,mapped_29)
fd5.29$mapped_rate = fd5.29$V2 / fd5.29$num_seqs / 1000000

fd5 = rbind(fd5.28,fd5.29)

fp5=ggplot(fd5, aes(x = avg_len, y = mapped_rate)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)+
  labs(x = "QCed reads mean lengh", y = "Percentage of mapped reads") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))

ggsave(plot=fp5, height=6, width=9, dpi=300,
       filename="figure5.pdf",
       useDingbats=FALSE)



#figure6
fd6 = dm[,c(1,5)]
fd6$file = 1:40

fp6=ggplot(fd6, aes(x = file, y = new_rate)) + 
  geom_point() +
  labs(x = "Sample", y = "Percentage of unique new read") +
  ylim(0,100) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, colour = "black" ,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        axis.line.x = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.5, linetype = "solid"),
        axis.ticks.length = unit(0.2, "cm"))

ggsave(plot=fp6, height=6, width=9, dpi=300,
       filename="figure6.pdf",
       useDingbats=FALSE)




