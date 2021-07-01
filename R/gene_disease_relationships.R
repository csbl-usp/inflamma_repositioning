#script to analyze gene-disease relationships of chronic inflammatory diseases

#load required packages
library(data.table)
library(dplyr)
library(igraph)
library(clusterProfiler)
library(ggplot2)
library(ggraph)
library(ForceAtlas2)
library(resolution)
options(stringsAsFactors = F)

#load Reactome pathways gmt
reactome_gmt <- readRDS("data/all_level_reactome.RDS")

#make dis_gene_edges and dis_drug_edges (STEP A)
path <- "data/ID_till2018"
df <- list.files(path = path,pattern = ".csv",full.names = T) %>%
  lapply(read.csv,check.names = T) %>% 
  bind_rows %>%
  as.data.frame()

#get genes with known connections with inflammatory diseases
dis_gene_edges <- df %>%
  dplyr::filter(Source.type=="GENE" & Target.type=="CONDITION") %>%
  dplyr::filter(Confidence >= 50 | Confidence == 1) %>%
  dplyr::filter(Documents >= 2) %>%
  dplyr::filter(!duplicated(paste0(Source.Display.Name,Target.Display.Name))) %>%
  dplyr::select(3,6,8,9,10)

colnames(dis_gene_edges) <- c("Source","Target","Confidence","Documents","Doc_IDs")

genes <- unique(dis_gene_edges$Source) #1351 genes

#convert genes not in gmt file to adequate format
genes_not_in_reactome <- genes[!genes %in% reactome_gmt$gene]

#convert genes to official gene symbol (genes that are not in correct format)
#this table was created mannually to match WDD genes to official gene symbols elsewhere
load("data/gene_convertion_table.RData")
genes_not_in_reactome_df <- gene_convertion_table %>%
  filter(alias %in% genes_not_in_reactome)

dis_gene_edges_converted <- dis_gene_edges %>%
  dplyr::left_join(genes_not_in_reactome_df,by=c("Source"="alias")) %>%
  dplyr::mutate(converted=ifelse(is.na(converted),Source,converted),
         Source=converted) %>%
  dplyr::select(-converted)


#run Fisher test between diseases to detect cluster of similar diseases
diseases <- unique(dis_gene_edges_converted$Target)
diseases_gmt <- dis_gene_edges_converted %>%
  dplyr::select(2,1) %>%
  dplyr::rename(disease=1,gene=2)

disease_genes <- split(diseases_gmt$gene,diseases_gmt$disease)

enrichment <- lapply(disease_genes,FUN = enricher,
                     pAdjustMethod = "BH",
                     universe = unique(diseases_gmt$gene),
                     TERM2GENE = diseases_gmt)

for(i in 1:length(enrichment)){
  enr <- enrichment[[i]]
  dis <- names(enrichment)[i]
  if(class(enr)=="NULL"){
    next()
  }else{
    enr_res <- enr@result
    enr_res <- enr_res %>%
      filter(ID != dis) %>%
      mutate(disease=dis)
  }
  if(i==1){
    enr_res_all <- enr_res
  }else{
    enr_res_all <- rbind(enr_res_all,enr_res)
  }
}

#find modules of diseases that significantly share genes
dis_gene_network <- dis_gene_edges_converted %>%
  dplyr::select(Source,Target) %>%
  dplyr::rename(from=1,to=2)

g1 <- graph_from_data_frame(dis_gene_network,directed = F)
V(g1)$class <- c(rep("gene",1349),
                rep("disease",27))

V(g1)$degree <- degree(g1)

dis_degree <- data.frame(node=V(g1)$name,
                         class=V(g1)$class,
                         degree=as.numeric(degree(g1)))

x <- resolution::cluster_resolution_RandomOrderFULL(graph = g1,t = 0.9)

V(g1)$module <- as.character(cluster_louvain(g1)$membership)
V(g1)$labels <- c(rep(NA,1349),
                  dis_degree$node[1350:1376])

col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
pal <- sample(col_vector,size = 15,replace = F)

layout <- layout.forceatlas2(graph = g1,iterations = 1000,plotstep = 0)
p <- ggraph(g1,layout = layout)+
  geom_edge_fan(color="lightgrey",alpha=0.5)+
  geom_node_point(aes(color=module,size=degree))+
  geom_node_text(aes(label=labels),size=2)+
  scale_color_manual(values = pal)+
  theme_void()
p

network <- enr_res_all %>%
  dplyr::select(disease,ID,p.adjust) %>%
  dplyr::filter(p.adjust < 0.01) %>%
  dplyr::mutate(weight=-log(p.adjust)) %>%
  dplyr::select(-p.adjust) %>%
  dplyr::rename(from=1,to=2)

g <- graph_from_data_frame(network,directed = F)
E(g)$weight <- network$weight
V(g)$degree <- degree(g)
p <- ggraph(graph = g,layout = "auto")+
  geom_node_point(aes(size=degree))+
  geom_edge_link(aes(width=weight))
p
