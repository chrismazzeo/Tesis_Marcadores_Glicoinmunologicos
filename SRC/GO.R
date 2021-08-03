#TOP GO

library(readr)
geneList = ""
geneList <- read_delim("Results/MOUSE/GSE43338/DEA/colorectal_tumor_colitis_associated_VS_colorectal_control_epithelium_colitis_associated/DEA/RES_SIGNIFICANT_GLYCO_GENES.csv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
dim(geneList)
geneList = na.omit(geneList)
dim(geneList)
geneList = geneList[order(-abs(geneList$log2FoldChange)),]
#geneList = ONTHO_GO_GENE_LIST(tempResults,padj = project.settings.DEA.padj, log2FC_UP = 1, log2FC_DOWN = 1)
dim(geneList)
HELPER_SaveCSV(geneList, "./go.csv")


library(GOplot)
data(EC)
head(EC$david)

circ <- circle_dat(EC$david, EC$genelist)
GOBubble(circ, labels = 3)
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)



top_GO_data <- new("topGOdata",
                   description = "Simple session",
                   ontology = "BP",
                   allGenes = geneList,
                   nodeSize = 10, 
                   annot = annFUN.db, 
                   affyLib = paste0(annotation(affyData),".db"))



#testeamos los terminos GO sobrerepresentados dentro de los genes con expresion diferencial
result_top_GO_classic =  runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

#The algorithm starts processing the nodes / GO categories from the highest (bottommost) level and then iteratively moves to nodes from a lower level. If a node is scored as significant, all of its genes are marked as removed in all ancestor nodes. This way, the “elim” algorithm aims at finding the most specific node for every gene.
result_top_GO_elim =     runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")


res_top_GO <- GenTable(top_GO_data, 
                       Fisher.classic = result_top_GO_classic,
                       Fisher.elim = result_top_GO_elim,
                       orderBy = "Fisher.elim" , ranksOf = "Fisher.classic", topNodes = 100)

#summary
genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = paste0(annotation(affyData),".db"), geneCutOff = 1000)


#Los genes que estan en una term y son significativos tienen  pvalue = 2
res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"]," - "),
        collapse = "")
})
res_top_GO$Term=gsub(",", " - ", res_top_GO$Term)
HELPER_SAVE_DATA_FRAME(res_top_GO, paste0(outputDirGroup,"GO/EnrichedTerms.csv"))
pdf(paste0(outputDirGroup,"GO/EnrichedTerms.pdf"), width = 10, height = 20)
graphicsOver =  res_top_GO %>% 
  top_n(n = 100, wt=-`Rank in Fisher.classic`) %>% 
  mutate(hitsPerc=Significant*100/Annotated) %>% 
  ggplot(aes(x=hitsPerc, 
             y=Term, 
             colour=`Rank in Fisher.classic`, 
             size=Significant)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="Rank", size="Count")

print(graphicsOver)
dev.off()
