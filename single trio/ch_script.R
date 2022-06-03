library(tidyverse)

data = read.csv('2022-03-28_16-11-18_MCD289_genotypes.csv')
data = data %>% filter(Covered.Case==3)

###Step 1: filter gene with >1 variants in child, and coverage=3
gene = data %>% 
  select(Variant.ID, UpToDate.Gene.Name,Sample.Name) %>% 
  filter(Sample.Name=='GHART291b') %>%
  group_by(UpToDate.Gene.Name) %>% 
  count() %>% 
  filter(n>=2) %>% 
  select(UpToDate.Gene.Name)
gene = gene$UpToDate.Gene.Name
two_var = data %>% filter(UpToDate.Gene.Name %in% gene)

###Step 2: convert data into child variant and parent info
trio = two_var %>% 
  select(Variant.ID,UpToDate.Gene.Name,GT,Sample.Name) %>% 
  pivot_wider(
    names_from = Sample.Name,
    values_fill = "0",
    values_from = GT
  ) %>% 
  rename("child"=GHART291b,"mother"=GHART1506,"father"=GHART1542a) %>% 
  filter(child!="0") 
  # mutate(GHART291b=paste(mother,father,sep="|") )

###Step 3: add mut_flag and phasing info
trio = trio %>% 
  mutate(mother_n=0,father_n=0, phasing="",mut_flag="",
         mother_n = replace(mother_n,mother!='0',1),
         father_n = replace(father_n,father!='0',1),
         sum = mother_n+father_n) %>% 
  mutate(mut_flag=replace(mut_flag,sum==0,'denovo'),
         mut_flag=replace(mut_flag,child=='hom','hom'),
         mut_flag=replace(mut_flag,(child=='het' & mother=='hom' & father=='hom'),'strange'),
         phasing=case_when(
           child=='hom' ~ '1|1',
           (child=='het' & sum==1) ~ paste(mother_n,father_n,sep = '|'),
           (child=='het' & sum==0) ~ '.|.', #'denovo'
           (child=='het' & sum==2 & mother=='hom' & father=='hom') ~ '.|.', #strange
           (child=='het' & sum==2 & (mother=='hom' | father=='hom')) ~ paste(as.numeric(mother=="hom"),as.numeric(father=="hom"),sep = '|'),
           (child=='het' & sum==2 & (mother=='het' & father=='het')) ~ '.|.' #check the other variant
         )
  ) %>% 
  select(-mother_n, -father_n, -sum)

###Step 4: add CH tag
CH = trio %>% 
  filter(mut_flag !='hom',mut_flag !='denovo') %>% 
  separate(phasing,into = c('mother_n','father_n'),sep = '\\|') %>%
  select(UpToDate.Gene.Name,mother_n,father_n) %>% 
  pivot_longer(
    c(mother_n,father_n),
    names_to = 'parent',
    values_to = 'haplo'
  ) %>% 
  count(UpToDate.Gene.Name,parent,haplo) %>% 
  filter(haplo!='0') %>% 
  pivot_wider(
    names_from = parent,
    values_from = n,
    values_fill = 0
  )

ch = NULL
for (gene in unique(CH$UpToDate.Gene.Name)){
  tmp = filter(CH,UpToDate.Gene.Name==gene)
  if (nrow(tmp)==2){ ### when one het phased and one het unknown phasing --> CH
    ch = c(ch,gene)
  } else if (tmp$haplo=="1"){ ### all successfully phased with var on both allele
    if (tmp$father_n >0 & tmp$mother_n >0) ch=c(ch,gene)
  } ### both het are unknown phased --> same
}
CH = CH %>% filter(UpToDate.Gene.Name %in% ch)

trio = trio %>% mutate(mut_flag=replace(mut_flag,UpToDate.Gene.Name %in% ch, 'CH'))
trio = as.data.frame(trio)
rownames(trio) = trio$Variant.ID

###Step 5: generate new file
child = two_var %>% filter(Sample.Name=="GHART291b")
child = cbind(child[,1:27],trio[child$Variant.ID,4:7],child[,28:219])
rownames(child) = NULL
write.csv(child,'chr22_ch.csv',row.names = F)

