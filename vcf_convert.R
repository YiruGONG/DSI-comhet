# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("thierrygosselin/radiator")
# library(radiator)
library(tidyverse)

data = read.csv('2022-03-28_16-11-18_MCD289_genotypes.csv')

idx = which(colnames(data) %in% c("AA", "AC", "AF", "AFR_AF", "AMR_AF", "AN", "DP", 
                            "EAS_AF", "EUR_AF", "HRun", "NS", "SAS_AF", "VT", "EFF"))
## How to integrate the columns
data[243:244,idx]

## integrate phasing
data = data[which(data$Covered.Case==3),]

trio = data %>% 
  select(Variant.ID,Rs.Number,Sample.Name) %>% 
  mutate(value=1) %>% 
  pivot_wider(
    names_from = Sample.Name,
    values_fill = 0,
    values_from = value
  ) %>% 
  filter(GHART291b==1) %>% 
  mutate(GHART291b=paste(GHART1506,GHART1542a,sep="|") )
rownames(trio) = trio$Variant.ID

child = data[which(data$Sample.Name=="GHART291b"),]
child$GHART291b = trio[child$Variant.ID,]$GHART291b

info = strsplit(child$Variant.ID,"-")
info = Reduce(rbind,info)
info = as.data.frame(info)
colnames(info) = c('CHROM',"POS","REF","ALT")
rownames(info) = NULL

vcf = cbind(info,child) %>% 
  select(-Variant.ID) %>% 
  mutate(INFO=".",
         CHROM=paste0('chr',CHROM)) %>% 
  rename('ID'=Rs.Number,'QUAL'=Qual,"#CHROM"=CHROM) %>% 
  select(`#CHROM`,POS,ID,REF,ALT,QUAL,FILTER,INFO,everything())
vcf[is.na(vcf$ID),]$ID = '.'

view = select(vcf, GT, GHART291b)

write.table(vcf,"chr22.vcf",sep = '\t',quote = F,row.names = F)

vcf2 = vcf[,c(1:8,224)]
write.table(vcf2,"simple_chr22.vcf",sep = '\t',quote = F,row.names = F)



# 
# 
# ##old
# info = strsplit(child$Variant.ID,"-")
# info = Reduce(rbind,info)
# info = as.data.frame(info)
# colnames(info) = c('CHROM',"POS","REF","ALT")
# rownames(info) = NULL
# 
# vcf = cbind(info[,1:2],ID = child$Rs.Number, info[3:4],
#             QUAL = child$Qual, FILTER = child$FILTER,
#             INFO =".",sample=child$Sample.Name,value=1)
# vcf[is.na(vcf$ID),]$ID = '.'
# vcf$CHROM = paste0('chr',vcf$CHROM)
# 
# vcf2 = vcf %>% 
#   pivot_wider(
#     names_from = sample,
#     values_fill = 0,
#     values_from = value
#   ) %>% 
#   group_by(CHROM,POS,REF,ALT) %>% 
#   summarize(
#     ID = paste(unique(ID),collapse = ","),
#     # REF = paste(unique(REF),collapse = ","),
#     ALT = paste(unique(ALT),collapse = ","),
#     QUAL = ".",
#     FILTER = paste(FILTER,collapse = ","),
#     INFO = ".",
#     GHART1506 = sum(GHART1506),
#     GHART1542a = sum(GHART1542a),
#     GHART291b = sum(GHART291b)
#   ) %>% 
#   relocate(CHROM,POS,ID,everything())
# 
# vcf2[,9:11] = apply(vcf2[,9:11],2,as.character)
# vcf2[vcf2==1] = "1/1"
# vcf2[vcf2==0] = "0/0"
# 
# colnames(vcf2)[1] = "#CHROM"
# write.table(vcf2,"chr22.vcf",sep = '\t',quote = F,row.names = F)
# 
# ##Calculation
# r1 = rownames(unique(info))
# r2 = rownames(unique(info[,1:2]))
# setdiff(r1,r2)
