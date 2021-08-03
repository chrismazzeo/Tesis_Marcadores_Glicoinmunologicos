path = "./results/human/gdc/LGALS1_KO-6KO/htseq/DEA/LGALS1_KO_VS_AOM_DSS/DEA/RES_SIGNIFICANT_GLYCO_GENES.csv"
glycoGenes = HELPER_LoadTSV(path)



tempGlyco =glycoGenes
tempGlyco = tempGlyco[order(tempGlyco$log2FoldChange, decreasing = TRUE),]

tempGlyco.all = as.matrix(tempGlyco$log2FoldChange)
rownames(tempGlyco.all) = tempGlyco$external_gene_name
breaksList = seq(-8, 8, by = 2)

grid.arrange(rectGrob(), rectGrob())
max.row = 55
max.col = 3
max.pages = ceiling(length(tempGlyco.all)/max.row/max.col)

col.number = ceiling(length(tempGlyco.all)/max.row)


for (page in 1:max.pages){
  plot_list=list()
  start =  max.row*max.col*(page-1)+page
  offset = start +max.row*max.col
  if (offset> length(tempGlyco.all)){
    offset = length(tempGlyco.all) 
  }
  
  col.number = ceiling((offset-start)/max.row)
  for (i in 1:col.number){
      base = start + max.row*(i-1)
      if ((base + max.row) > length(tempGlyco.all)){
        end = length(tempGlyco.all)
      }else{
        end = base + max.row
      }
  
        temp = as.matrix( tempGlyco.all[base:end,])
      
       
        x = pheatmap(temp, 
                     color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaksList)),
      
                     breaks = breaksList,
                     cluster_rows = FALSE,cluster_cols = FALSE, cellwidth = 40, cellheight =14,display_numbers = TRUE,fontsize_number = 10,fontsize_row = 13,
                     na_col = "black", silent = TRUE)
        
        plot_list[[i]] = x[[4]]
      
  }

  g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=col.number))
  ggsave(paste0("./LGALS1_KO_VS_AOM_DSS","_",page,".pdf"), g)
    
}    
    
   

