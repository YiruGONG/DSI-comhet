# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("thierrygosselin/radiator")
# library(radiator)
library(tidyverse)

data = read.csv('2022-03-28_16-11-18_MCD289_genotypes.csv')

idx = which(colnames(data) %in% c("AA", "AC", "AF", "AFR_AF", "AMR_AF", "AN", "DP", 
                            "EAS_AF", "EUR_AF", "HRun", "NS", "SAS_AF", "VT", "EFF"))

## integrate phasing
# data = data[which(data$Covered.Case==3),]


###### generate mother/father info ########
trio = data %>% 
  select(Variant.ID,Rs.Number,Sample.Name) %>% 
  mutate(value=1) %>% 
  pivot_wider(
    names_from = Sample.Name,
    values_fill = 0,
    values_from = value
  ) %>% 
  filter(GHART291b==1) %>% 
  rename("mother"=GHART1506,"father"=GHART1542a) %>% 
  mutate(GHART291b=paste(mother,father,sep="|") )
trio = as.data.frame(trio)
rownames(trio) = trio$Variant.ID


###### Add CH tags ########
## other tags tbd!!!
child = data[which(data$Sample.Name=="GHART291b"),]
child = cbind(child[,1:27],trio[child$Variant.ID,3:5],mut_flag=NA ,child[,28:219])
idx = child %>% 
  select(UpToDate.Gene.Name,GT,mother,father) %>% 
  filter(GT=='het' | mother+father==1) %>% 
  group_by(UpToDate.Gene.Name) %>% 
  count() %>% 
  filter(n>=2) %>% 
  select(UpToDate.Gene.Name)
idx = idx$UpToDate.Gene.Name
child[which(child$UpToDate.Gene.Name %in% idx),'mut_flag'] = 'CH'

view = select(child, GT, GHART291b)
unique(view)
# ## add tags
# child[which(child$GT=="het" & child$GHART291b=="0|0"),"mut_flag"]="denovo"
# child[which(child$GT=="het" & (child$GHART291b=="1|0" | child$GHART291b== "0|1" )),"mut_flag"]="CH"
# child[which(child$GT=="het" & (child$GHART291b=="1|1")),"mut_flag"]="other"
# child[which(child$GT=="hom" & child$GHART291b=="0|0"),"mut_flag"]="other"
# child[which(child$GT=="hom" & (child$GHART291b=="1|0" | child$GHART291b== "0|1" )),"mut_flag"]="other"
# child[which(child$GT=="hom" & (child$GHART291b=="1|1")),"mut_flag"]="homo"

write.csv(child,"child_chr22.csv",row.names = F,na="")


###### identify CH ########
CH = child %>% 
  filter(Covered.Case==3) %>% 
  filter(GT=="het", 
         GHART291b %in% c("1|0","0|1")) %>% 
  # separate(GHART291b,into = c('mother','father'),sep = '\\|') %>% 
  mutate( gene = strsplit(UpToDate.Gene.Name,"\\|") ) %>% 
  unnest(gene) %>% 
  relocate(gene, .before = UpToDate.Gene.Name) %>% 
  group_by(gene) %>% 
  summarise(mother_sum = sum(as.numeric(mother)),
            father_sum = sum(as.numeric(father)) ) %>% 
  filter(mother_sum!=0, father_sum!=0)
##Problem: one row only corresponding to one gene?

write.csv(CH,'chr22_CH.csv',row.names = F)


####### Compare with VarCount #######
vc = read.table('CH_CountFile.txt',header = T)
vc = vc[,-4] %>% filter(COUNT != 0) %>% select(-COUNT)
setdiff(CH$gene,vc$CODING_GENE)
intersect(CH$gene,vc$CODING_GENE)
setdiff(vc$CODING_GENE,CH$gene)


####### Generate VCF #######
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

# view = select(vcf, GT, GHART291b)

write.table(vcf,"chr22.vcf",sep = '\t',quote = F,row.names = F)

vcf2 = cbind(vcf[,c(1:8)],GHART291b=vcf$GHART291b) 
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
