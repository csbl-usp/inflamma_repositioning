#script to analyze gene-disease relationships of chronic inflammatory diseases

#load required packages
library(data.table)
library(dplyr)
library(tidyverse)
library(igraph)
library(reshape2)
library(clusterProfiler)
library(pheatmap)
library(ggplot2)
library(ggraph)
library(ForceAtlas2)
library(resolution)
options(stringsAsFactors = F)
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

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
  dplyr::select(-converted) %>%
  dplyr::filter(!duplicated(.))

#run Fisher test between diseases to detect cluster of similar diseases
diseases <- unique(dis_gene_edges_converted$Target) #27 diseases
diseases_gmt <- dis_gene_edges_converted %>%
  dplyr::select(2,1) %>%
  dplyr::rename(disease=1,gene=2)

converted_genes <- unique(dis_gene_edges_converted$Source) #1349 genes

disease_genes <- split(diseases_gmt$gene,diseases_gmt$disease)

#find modules of diseases that significantly share genes in Gephi (outside R)
dis_gene_network <- dis_gene_edges_converted %>%
  dplyr::select(Source,Target) %>%
  dplyr::rename(from=1,to=2)

#export to gephi for outside module detection
dis_gene_edges <- dis_gene_network %>%
  rename(Source=1,Target=2)

dis_gene_nodes <- data.frame(Id=unique(c(dis_gene_edges$Source,
                                         dis_gene_edges$Target)),
                             Label=c(rep("",length(unique(dis_gene_edges$Source))),
                                     unique(dis_gene_edges$Target)),
                             Class=c(rep("gene",length(unique(dis_gene_edges$Source))),
                                     rep("disease",length(unique(dis_gene_edges$Target)))))

write.csv(dis_gene_edges,file = "data/dis_gene_edges.csv",
          row.names = F,quote = F)

write.csv(dis_gene_nodes,file = "data/dis_gene_nodes.csv",
          row.names = F,quote = F)

#load module assignment made in Gephi
#Module detection in Gephi was performed using the "Modularity" function
#with a resolution of 1.5 to find five modules of diseases and genes.
#Figure 01A was created using Gephi's "ForceAtlas2" layout.
#Nodes were collored according to modularity class. 
#Node size was proportional to degree calculated with the
#"Avarage Degree" function in Gephi.
gephi_modules <- fread("data/gephi_modules.csv")

disease_modules <- gephi_modules %>%
  filter(class=="disease") %>%
  select(Id,Degree,modularity_class)

gene_modules <- gephi_modules %>%
  filter(class=="gene") %>%
  select(Id,Degree,modularity_class)

#plot genes per disease plot
p <- disease_modules %>%
  mutate(Id=firstup(Id)) %>%
  arrange(desc(Degree)) %>%
  mutate(Id=factor(Id,levels = rev(Id))) %>%
  ggplot(aes(y=Id,x=Degree))+
  geom_segment(aes(x=0, xend=Degree, y=Id, yend=Id),
               color="pink",size=1)+
  geom_point(aes(size=Degree))+
  scale_x_continuous(name = "#Genes")+
  scale_y_discrete(name="Disease")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15))
pdf("figures/genes_per_disease.pdf",width = 6.75,height = 6.95)
print(p)
dev.off()

#plot diseases per gene plot
p <- gene_modules %>%
  top_n(wt = Degree,n = 20) %>%
  arrange(desc(Degree)) %>%
  mutate(Id=factor(Id,levels = rev(Id))) %>%
  ggplot(aes(y=Id,x=Degree))+
  geom_segment(aes(x=0, xend=Degree, y=Id, yend=Id),
               color="skyblue",size=1)+
  geom_point(aes(size=Degree))+
  scale_x_continuous(name = "#Diseases")+
  scale_y_discrete(name="Gene")+
  theme_minimal()+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15))
pdf("figures/diseases_per_gene.pdf",width = 6.75,height = 6.95)
print(p)
dev.off()

#use clusterProfiler enricher function to calculate dis-dis enrichment
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

#plot dis-dis network using enrichment results as edge weights
network <- enr_res_all %>%
  dplyr::select(disease,ID,p.adjust) %>%
  dplyr::filter(p.adjust < 0.01) %>% #keep only significant dis-dis relations
  dplyr::mutate(weight=-log(p.adjust)) %>%
  dplyr::select(-p.adjust) %>%
  dplyr::rename(from=1,to=2)

#plot using ggraph
g <- graph_from_data_frame(network,directed = F)
E(g)$weight <- network$weight

diseases_in_network <- V(g)$name

degrees <- disease_modules %>%
  filter(Id %in% diseases_in_network)

degrees <- degrees[match(diseases_in_network,degrees$Id),]
  
V(g)$degree <- degrees$Degree
V(g)$module <- degrees$modularity_class

lay <- ForceAtlas2::layout.forceatlas2(g,directed = F,
                                       plotstep = 0,
                                       iterations = 1000)

#color in palette obtained from Figure 01A - modules in Gephi dis-gene
#network
pal <- c("#8cb900","#d97dd8","#23966f","#00c7ff","#ff7045")
p <- ggraph(graph = g,layout = lay)+
  geom_edge_link(aes(color=weight,width=weight))+
  scale_edge_width(range = c(.2,4))+
  scale_edge_color_gradient(low = "pink",high = "red")+
  geom_node_point(aes(size=degree,color=as.character(module)))+
  geom_node_text(aes(label=V(g)$name),size=2.5)+
  scale_size_continuous(range = c(2,12))+
  scale_color_manual(values = pal)+
  theme_void()
pdf(file = "figures/dis_dis_network_modules.pdf",
    width = 7.9,height = 5.8)
print(p)
dev.off()
