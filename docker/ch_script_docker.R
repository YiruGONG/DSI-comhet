library(tidyverse)

startTime <- Sys.time()

data = read.csv("/data/chr22_RHD_335trios.csv")

child = "1"
mother = "2"
father = "3"
dp_cut = 10
dp_col = colnames(data)[str_starts(colnames(data),'DP_')]
gene_col = "SYMBOL"

GT_child = paste0('GT_',child)
GT_mother = paste0('GT_',mother)
GT_father = paste0('GT_',father)

data = data %>%
  ###Step 0: select the first gene name from SYMBOL
  mutate(gene=sapply(strsplit(!!as.symbol(gene_col),","),function(x) x[1])) %>% 
  relocate(gene, .before = !!as.symbol(gene_col)) %>% 
  ###Step 1: filter coverage=3
  # filter(DP_1>dp_cut, DP_2>dp_cut, DP_3>dp_cut) %>%
  filter(if_all(.cols = all_of(dp_col),
                .fns = ~ .x >= dp_cut)) %>% 
  ###Step 2: filter variables in child
  filter(!!as.symbol(GT_child) !="'0/0'")

###Step 3: obtain results in each trio
# switch to parallel calculation if needed
# EG. sample = 'BS_1F9EFPM7'

#### Main function
CH_identifier = function(subset){
  # subset = subset %>%
  #   mutate( gene = strsplit(!!as.symbol(gene_col),",") ) %>%
  #   unnest(gene) %>%
  #   relocate(gene, .before = !!as.symbol(gene_col))
  
  ###Step 3: filter gene with >1 variants in child
  gene_idx = subset %>% 
    group_by(gene) %>% 
    count() %>% 
    filter(n>=2) %>% 
    select(gene)
  gene_idx = gene_idx$gene
  two_var = subset %>% filter(gene %in% gene_idx)
  
  ###Step 4: add mut_flag and phasing info
  two_var = two_var %>% 
    mutate(phasing="",mut_flag="",
           child_het = (!!as.symbol(GT_child) %in% c("'0/1'","'1/0'") ),
           mother='ref', father= 'ref',
           mother = case_when(
             !!as.symbol(GT_mother)=="'1/1'" ~ "hom",
             !!as.symbol(GT_mother) %in% c("'0/1'","'1/0'") ~ "het",
             !!as.symbol(GT_mother)=="'0/0'" ~ "ref"
           ),
           father = case_when(
             !!as.symbol(GT_father)=="'1/1'" ~ "hom",
             !!as.symbol(GT_father) %in% c("'0/1'","'1/0'") ~ "het",
             !!as.symbol(GT_father)=="'0/0'" ~ "ref"
           ),
           mother_n=0,father_n=0,
           mother_n = replace(mother_n,mother!='ref',1),
           father_n = replace(father_n,father!='ref',1),
           parent_var_sum = mother_n+father_n
    ) %>% 
    relocate(phasing,mut_flag,.before = !!as.symbol(paste0('PGT_',child))) %>% 
    mutate(mut_flag=replace(mut_flag,parent_var_sum==0,'denovo'),
           mut_flag=replace(mut_flag,!!as.symbol(GT_child)=="'1/1'",'hom'),
           mut_flag=replace(mut_flag,(child_het & mother=='hom' & father=='hom'),'strange'),
           phasing = case_when(
             !!as.symbol(GT_child)=="'1/1'" ~ '1|1',
             (child_het & parent_var_sum==1) ~ paste(mother_n,father_n,sep = '|'),
             (child_het & parent_var_sum==0) ~ '.|.', #'denovo'
             (child_het & parent_var_sum==2 & mother=='hom' & father=='hom') ~ '.|.', #strange
             (child_het & parent_var_sum==2 & (mother=='hom' | father=='hom')) ~ paste(as.numeric(mother=="hom"),as.numeric(father=="hom"),sep = '|'),
             (child_het & parent_var_sum==2 & (mother=='het' & father=='het')) ~ '.|.' #check the other variant
           )
    ) %>% 
    select(-child_het,-mother,-father,-mother_n,-father_n,-parent_var_sum)
  
  ###Step 5: add CH tag
  CH = two_var %>% 
    filter(mut_flag !='hom',mut_flag !='denovo') %>% 
    separate(phasing,into = c('mother_n','father_n'),sep = '\\|') %>%
    select(gene,mother_n,father_n) %>% 
    pivot_longer(
      c(mother_n,father_n),
      names_to = 'parent',
      values_to = 'haplo'
    ) %>% 
    count(gene,parent,haplo) %>% 
    pivot_wider(
      names_from = parent,
      values_from = n,
      values_fill = 0
    ) %>% 
    filter(haplo!='0')
  
  ch = NULL
  pot_ch = NULL
  for (gene in unique(CH$gene)){
    tmp = CH[which(CH$gene==gene),]
    if (nrow(tmp)==2){ ### when one het phased and one het unknown phasing --> potential CH
      pot_ch = c(pot_ch,gene)
    } else if (tmp$haplo=="1"){ ### all successfully phased with var on both allele
      if (tmp$father_n >0 & tmp$mother_n >0) ch=c(ch,gene)
    } ### both het are unknown phased --> same
  }
  # CH = CH %>% filter(gene %in% ch)
  two_var = two_var %>% mutate(mut_flag=replace(mut_flag,gene %in% ch, 'CH'),
                               mut_flag=replace(mut_flag,gene %in% pot_ch, 'potential CH'))
  
  return(two_var)
}

### The actual execution
output = NULL
for ( sample in unique(data[[paste0("SAMPLE_",child)]]) ){
  subset = data %>% filter(!!as.symbol(paste0('SAMPLE_',child)) ==sample)
  out = CH_identifier(subset)
  output = rbind(output,out)
}

write.csv(output,'/data/chr22_335_flagged.csv',row.names = F)
endTime <- Sys.time()
print(endTime - startTime)
### 26.7 secs for chr22, 335 trios


#################### Phase 2: add another comparative output file
var_info = c("VAR","GT_1","GT_2","GT_3") ##add the desired column info here

#output = read.csv('chr22_335_flagged.csv')
ch = output %>% 
  filter(mut_flag %in% c("CH","potential CH")) %>% 
  select(SAMPLE_1,gene,mut_flag,VAR,starts_with('GT_'),phasing) %>% 
  nest(info=all_of(var_info)) %>% #c(VAR,starts_with('GT_'))
  group_by(SAMPLE_1,gene,mut_flag) %>% 
  pivot_wider(
    names_from = phasing,
    values_from = info
  ) %>% 
  rename(mother=`1|0`,father=`0|1`,unknown=`.|.`) %>% 
  unnest(c(mother, father, unknown),names_sep = '_')

write.csv(ch,"/data/comparative_335.csv",row.names = F, na="")

# tag: category: CH, potential CH

