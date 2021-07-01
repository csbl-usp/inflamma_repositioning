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

write.csv(exclusive_genes_diseases,"tables/table_S1_exclusive_genes.csv",row.names = F)

#load all files for genes searched in WDD (big drug-gene interaction table)
drugs_all <- fread("data/drugs_all_50percent_2documents.csv")

#get all drugs that interact with the exclusive 827 genes 
#(Step C)
drugs_exclusive <- drugs_all %>%
  filter(Target %in% exclusive_genes$Id)

length(unique(drugs_exclusive$Source)) #2367 drugs
length(unique(drugs_exclusive$Target)) #704 genes

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

###drugs_exclusive_c is the source for the new table_S2_potential###
#with this table, we curated ~25% of the relationships to detected
#drugs with further potential for repositioning to be used in chronic
#inflammatory diseases
table_S2_new <- merge(drugs_exclusive_c,dis_gene_edges,
                      by.x="Target",by.y="Source") %>%
  dplyr::select(2,1,7,5,10) %>%
  rename(Drug=Source,Gene=Target,ICD=Target.y,
         gene_ICD_ref=Doc_IDs.y,drug_gene_ref=Doc_IDs.x) %>%
  dplyr::select(Drug,Gene,ICD,gene_ICD_ref,drug_gene_ref)

#table to be used for curation of potential repositioning candidates
write.csv(table_S2_new,file="tables/table_S2_potential_drugs_all.csv",row.names = F)

###select 200 drugs to be manually curated (keep drugs already curated)###
#these drugs had been curated previously following the method described
#in the paper
table_S2_selected <- fread("data/selected_S2.csv") 

table_S2_selected_new_keep <- table_S2_new %>%
  dplyr::filter(Drug %in% table_S2_selected$Drug) %>%
  mutate(keep_new="keep")

table_S2_selected_new_rest <- table_S2_new %>%
  dplyr::filter(!Drug %in% table_S2_selected_new_keep$Drug) %>%
  dplyr::sample_n(size=105,replace = F) %>%
  mutate(keep_new="new")

table_S2_selected_new_final <- rbind(table_S2_selected_new_keep,table_S2_selected_new_rest)  
#save this table as table S3 - selected drugs (~25% of the 889 potential drugs)
write.csv(table_S2_selected_new_final,
          file="tables/table_S3_potential_drugs_selected.csv",
          row.names = F)

