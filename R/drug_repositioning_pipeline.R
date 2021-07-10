#script to perform all steps in the drug repositioning framework for
#chronic inflammatory diseases, following the same rationale used
#for neuropsychiatric diseases in: 
#"Drug repositioning for psychiatric and neurological disorders
#through a network medicine approach"
#https://doi.org/10.1038/s41398-020-0827-5

#load required packages
library(data.table)
library(dplyr)
library(igraph)
options(stringsAsFactors = F)

#make dis_gene_edges and dis_drug_edges (STEP A)
path <- "data/ID_till2018"
df <- list.files(path = path,pattern = ".csv",full.names = T) %>%
  lapply(read.csv,check.names = T) %>% 
  bind_rows %>%
  as.data.frame()

write.csv(df,file = "tables/table_S2_genes_and_drugs_CIDs_WDD.csv",
          row.names = F)

#get drugs that have known connections with inflammatory diseases
dis_drug_edges <- df %>%
  dplyr::filter(Source.type=="DRUG" & Target.type=="CONDITION") %>%
  dplyr::filter(Confidence >= 50 | Confidence == 1) %>%
  dplyr::filter(Documents >= 2) %>%
  dplyr::filter(!duplicated(paste0(Source.Display.Name,Target.Display.Name))) %>%
  dplyr::select(3,6,8,9,10)

colnames(dis_drug_edges) <- c("Source","Target","Confidence","Documents","Doc_IDs")

#get genes with known connections with inflammatory diseases
dis_gene_edges <- df %>%
  dplyr::filter(Source.type=="GENE" & Target.type=="CONDITION") %>%
  dplyr::filter(Confidence >= 50 | Confidence == 1) %>%
  dplyr::filter(Documents >= 2) %>%
  dplyr::filter(!duplicated(paste0(Source.Display.Name,Target.Display.Name))) %>%
  dplyr::select(3,6,8,9,10)

colnames(dis_gene_edges) <- c("Source","Target","Confidence","Documents","Doc_IDs")

genes <- unique(dis_gene_edges$Source) #1351 genes
drugs <- unique(dis_drug_edges$Source) #711 drugs

#calculate the degree of genes and keep only genes connected to ONE disease
#(Step B)
g <- graph_from_data_frame(dis_gene_edges %>%
                             dplyr::select(Source,Target),
                           directed = F)

degree <- data.frame(Id=names(degree(g)),degree=degree(g))

mods <- cluster_louvain(g)

exclusive_genes <- degree %>% #827 exclusive genes
  filter(degree==1)

exclusive_genes_diseases <- dis_gene_edges %>%
  filter(Source %in% exclusive_genes$Id)

write.csv(exclusive_genes_diseases,"tables/table_S3_exclusive_genes.csv",row.names = F)

#load all files for genes searched in WDD (big drug-gene interaction table)
drugs_all <- fread("data/drugs_all_50percent_2documents.csv")

#get all drugs that interact with the exclusive 827 genes 
#(Step C)
drugs_exclusive <- drugs_all %>%
  filter(Target %in% exclusive_genes$Id)

length(unique(drugs_exclusive$Source)) #2367 drugs
length(unique(drugs_exclusive$Target)) #704 genes

write.csv(drugs_exclusive,file="tables/table_S4_drugs_affecting_exclusive_CID_genes.csv",
          row.names = F,quote = T)

#remove drugs that are directly connected to ICDs
#(Step D1)
drugs_exclusive_b <- drugs_exclusive %>%
  filter(!Source %in% drugs) 

length(unique(drugs_exclusive_b$Source)) #1833 drugs
length(unique(drugs_exclusive_b$Target)) #587 genes

#remove drugs that target more than one gene
#(Step D2)
dupli_drugs <- unique(drugs_exclusive_b$Source[duplicated(drugs_exclusive_b$Source)])
length(unique(dupli_drugs)) #944 drugs that affect more than one gene

drugs_exclusive_c <- drugs_exclusive_b %>%
  filter(!Source %in% dupli_drugs)

length(unique(drugs_exclusive_c$Source))  #889 drugs
length(unique(drugs_exclusive_c$Target)) #276 genes

###drugs_exclusive_c is the source for the table_WDD_mistakes_curation###
table_S5_potential_repositioning_drugs <- merge(drugs_exclusive_c,dis_gene_edges,
                      by.x="Target",by.y="Source") %>%
  dplyr::select(2,1,7,5,10) %>%
  rename(Drug=Source,Gene=Target,ICD=Target.y,
         gene_ICD_ref=Doc_IDs.y,drug_gene_ref=Doc_IDs.x) %>%
  dplyr::select(Drug,Gene,ICD,gene_ICD_ref,drug_gene_ref)

#table to be used for curation of potential repositioning candidates
write.csv(table_S5_potential_repositioning_drugs,
          file = "tables/table_S5_potential_repositioning_candidates.csv",
          row.names = F,quote = T)

###select 225 drugs to be manually curated###
#we curated ~25% of the relationships to detected drugs with 
#potential for repositioning to be used for CIDs
# table_S6_selected_drugs_to_curate <- table_S5_potential_repositioning_drugs %>%
#   dplyr::sample_n(size=225,replace = F)
# 
# write.csv(table_S6_selected_drugs_to_curate,
#           file="tables/table_S6_WDD_mistakes_curation.csv",
#           row.names = F)
#THESE LINES ARE COMMENTED SINCE RUNNING THEM AGAIN WILL GENERATE A NEW
#SAMPLE OF 225 DRUGS. NEXT, WE PROVIDE THE ORIGINAL 225 CURATED DRUGS
#THAT WERE INCLUDED IN THE MANUSCRIPT WITH THE CURATION PERFORMED

table_S6_selected_drugs_to_curate <- fread("data/table_S6_WDD_mistakes_curation.csv")

length(unique(table_S6_selected_drugs_to_curate$Gene)) #112 exclusive genes
length(unique(table_S6_selected_drugs_to_curate$Drug)) #225 drugs
length(unique(table_S6_selected_drugs_to_curate$ICD)) #18 CIDs

#number of false gene-disease associations
table_S6_selected_drugs_to_curate %>%
  filter(!duplicated(Gene)) %>%
  pull(gene_ICD_TRUE_FALSE) %>%
  table() #29 False, 83 True

#number of false drug-gene associations
table_S6_selected_drugs_to_curate %>%
  filter(gene_ICD_TRUE_FALSE=="TRUE") %>%
  filter(!duplicated(paste0(Drug,Gene))) %>%
  pull(drug_gene_TRUE_FALSE) %>%
  table() #50 False, 117 True

#Drug-gene-disease table without WDD mistakes
drug_gene_CID_to_curate_for_potential <- table_S6_selected_drugs_to_curate %>%
  filter(gene_ICD_TRUE_FALSE=="TRUE",
         drug_gene_TRUE_FALSE=="TRUE")

length(unique(drug_gene_CID_to_curate_for_potential$Drug)) #117 drugs
length(unique(drug_gene_CID_to_curate_for_potential$Gene)) #60 genes
length(unique(drug_gene_CID_to_curate_for_potential$ICD)) #12 CIDs

#load curated drug-gene-disease with YES/NO for actual repo potential 
#this table was created outside R following the manual curation steps
#described in the manuscript
table_S7_curated_repositioning_table <- readxl::read_excel("tables/table_S7_curated_potential.xlsx")

#Potential after manual curation (drug-gene and gene-disease mechanism)
table(table_S7_curated_repositioning_table$`Drug repo potential`)
#NO potential = 28
#YES potential = 90

only_true_potential <- table_S7_curated_repositioning_table %>%
  filter(`Drug repo potential`=="YES")

length(unique(only_true_potential$Gene)) #44 genes
length(unique(only_true_potential$ICD)) #11 CIDs

#FDA approval all
table(table_S7_curated_repositioning_table$`FDA approval`)
#NO (discontinued) = 4
#Not found = 74
#YES = 40

#FDA approval potential
table(only_true_potential$`FDA approval`)
#NO (discontinued) = 3
#Not found = 55
#YES = 32

#EMA approval all
table(table_S7_curated_repositioning_table$`EMA approval`)
#NO (discontinued) = 13
#Not found = 65
#YES = 40

#EMA approval potential
table(only_true_potential$`EMA approval`)
#NO (discontinued) = 11
#Not found = 49
#YES = 30

#Litereature evidence pre 2018
table(table_S7_curated_repositioning_table$`Literature evidence (pre 2018)`)
#NO 78
#YES 40

#Litereature evidence post 2018
table(table_S7_curated_repositioning_table$`Literature evidence (post 2018)`)
#NO 83
#YES 30

#Clinical trial
table(table_S7_curated_repositioning_table$`Clinical trials`)
#NO (discontinued) = 108
#YES = 10


#No previouse literature evidence
no_literature <- only_true_potential %>%
  filter(`Literature evidence (pre 2018)`=="NO",
         `Literature evidence (post 2018)`=="NO",
         `Clinical trials`=="NO")

no_literature_mabs <- no_literature %>%
  filter(grepl("mab",tolower(Drug)))

write.csv(no_literature,"data/drugs_without_literature.csv",row.names = F)
