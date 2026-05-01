
dir<-"/Users/srujansingh/Documents/Personal/GeoMx"
setwd(dir)

library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(msigdbr) #For HALLMARK Pathway
library(enrichplot) #For ClusterProfiler()
library(biomaRt)

#Loading the results from the differential expression analysis performed earlier:
results2<-readRDS("results2.rds")
head(results2)

res_glom<-results2[(results2$Subset=="glomerulus")&(results2$'Pr(>|t|)' <0.05),]
res_tub<-results2[(results2$Subset=="tubule")&(results2$'Pr(>|t|)' <0.05),]

df_glom<-res_glom[rev(order(res_glom$Estimate)),]
df_tub<-res_tub[rev(order(res_tub$Estimate)),]

df_glom_HALL<-df_glom$Estimate
names(df_glom_HALL)<-df_glom$Gene
head(df_glom_HALL)

gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

df_glom.HALL.GSEA <- GSEA(
  geneList = df_glom_HALL,
  TERM2GENE = background,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

df_glom.HALL.GSEA.Plot<-df_glom.HALL.GSEA
gseaplot(df_glom.HALL.GSEA.Plot, geneSetID = 1)
saveRDS(df_glom.HALL.GSEA.Plot, file="EnrichmentPlots/Glomerulus_Normal_vs_DKD_HALL.GSEA.rds")

if(sign(max(df_glom.HALL.GSEA$NES))==1) 
{ df_glom.HALL.GSEA<-df_glom.HALL.GSEA[rev(order(df_glom.HALL.GSEA$NES)),]
} else { df_glom.HALL.GSEA<-df_glom.HALL.GSEA[order(df_glom.HALL.GSEA$NES),] }

df_glom.HALL.GSEA$core_enrichment<-gsub("/",",",df_glom.HALL.GSEA$core_enrichment)
head(df_glom.HALL.GSEA$Description) #Visualize the enriched pathways.

df_glom_vis<-data.frame(Pathway=df_glom.HALL.GSEA$Description,
                        NES=df_glom.HALL.GSEA$NES)

df_glom_vis$Pathway<-gsub("HALLMARK_","", df_glom_vis$Pathway)
df_glom_vis$Enrichment<-"Positive"
df_glom_vis[df_glom_vis$NES<0,]$Enrichment<-"Negative"
head(df_glom_vis)

df_glom_vis$Pathway<-factor(df_glom_vis$Pathway, levels=rev(df_glom_vis$Pathway))
#pdf("NES_GeoMx_Glom.pdf", width=8, height=6)

ggplot(df_glom_vis, aes(x=Pathway, y=NES, fill=Enrichment)) + 
  geom_bar(stat="identity") +
  coord_flip() +
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + 
  coord_flip() + theme_minimal()+
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.x = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.text.y = element_text(size = 14, color="black"),
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14))
#dev.off()

#VISUALIZATION OF ENRICHED PATHWAYS:
#pdf("MYC.pdf", width=7, height=3)
ID1="HALLMARK_MYC_TARGETS_V1"
gseaplot2(df_glom.HALL.GSEA.Plot, geneSetID = c(ID1),  subplots=1:2, title="MYC_TARGETS_V1_PATHWAY")
#dev.off()

#pdf("EMT.pdf", width=7, height=3)
ID2="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
gseaplot2(df_glom.HALL.GSEA.Plot, geneSetID = c(ID2),  subplots=1:2, title="EPITHELIAL_MESENCHYMAL_TRANSITION_PATHWAY")
d#ev.off()

#pdf("Unfolded_Protein.pdf", width=7, height=3)
ID3="HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
gseaplot2(df_glom.HALL.GSEA.Plot, geneSetID = c(ID3),  subplots=1:2, title="UNFOLDED_PROTEIN_RESPONSE_PATHWAY")
#dev.off()

#pdf("Glycolysis.pdf", width=7, height=3)
ID4="HALLMARK_GLYCOLYSIS"
gseaplot2(df_glom.HALL.GSEA.Plot, geneSetID = c(ID4),  subplots=1:2, title="GLYCOLYSIS_PATHWAY")
#dev.off()

#pdf("NFKB.pdf", width=7, height=3)
ID5="HALLMARK_TNFA_SIGNALING_VIA_NFKB"
gseaplot2(df_glom.HALL.GSEA.Plot, geneSetID = c(ID5),  subplots=1:2, title="TNFA_SIGNALING_VIA_NFKB_PATHWAY")
#dev.off()

#pdf("All_Pathways.pdf", width=8, height=6)
gseaplot2(df_glom.HALL.GSEA.Plot, geneSetID = c(ID1, ID2, ID3, ID4, ID5),  subplots=1:2, title="ALL_ENRICHED_PATHWAYS")
#dev.off()
