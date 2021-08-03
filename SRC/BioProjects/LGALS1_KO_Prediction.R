
res.pca = prcomp(vsd.glyco)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             label = "none",
             col.ind = group,
             repel = FALSE)
fviz_pca_var(
  res.pca,     # Avoid text overlapping,
  label = "var",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  col.var = "contrib"
)

fviz_pca(res.pca,   axes = c(1,2)   ,  # Avoid text overlapping,
         label = "none",
         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
         col.ind = group,
         select.var = list(contrib = 90)
)

fviz_pca(res.pca,  axes = c(1, 3) ,  # Avoid text overlapping,
         label = "none",
         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
         col.ind = group,
         select.var = list(contrib = 90)
)

fviz_pca(res.pca,  axes = c(2, 3) ,  # Avoid text overlapping,
         label = "none",
         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
         col.ind = group,
         select.var = list(contrib = 90)
)


#