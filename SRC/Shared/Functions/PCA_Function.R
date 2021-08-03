#ojo la matrix tienen que tener los objetos en las columnas
#escalo o no escalo? uso vsd o rawcounts?******************
htseq.rawMatrix.PCA = PCA(t(htseq.rawCountMatrix), scale.unit = FALSE, graph = FALSE)
htseq.vsd.PCA = PCA(t(assay(htseq.vsd)), scale.unit = FALSE, graph = FALSE)

eig.val <- get_eigenvalue(htseq.rawMatrix.PCA) 
eig.val

eig.val <- get_eigenvalue(htseq.vsd.PCA) 
eig.val

data(decathlon2)
X = decathlon2[, 1:10]
#PAra poder usar esto las variables tienen que estar correlacionadas
#Cada PC es ortogonal al resto
#obtenemos reduccion de dimensiones
res.pca = PCA(X, scale.unit = TRUE, ncp = 5, graph = FALSE)  #scale unit TRUE para que las variables con mayor numero u otra unidad tengan el mismo peso
print ("Eigen vaules -->Amount of variation retaind by each Principal Component") 
# Para ver con cuantos PC nos quedamos
eig.val <- get_eigenvalue(res.pca) 
#los eigen vector nos dice cuanta informacion retiene cada PC, si un eig>1 es porque retuvo mas informacion que una de las variables originales
# nos podemos quedar con PC que tenga eigen values > 1 o ver cuanto suman los primeros y si nos parece bien utilziamos eso
eig.val
#otra manera para elegires graficar los eigen values ordenaos de mayor a menor
# vemos cual tiene la mayor pendientes (al que le llega), y de ahi nos quedamos como el proximo  3 en este caso 
# o descartar los que son mas parecidos en este caso nos quedamos con los primeros 5
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50)) 

#extraemos las informacion de las variables del PCA
var <- get_pca_var(res.pca)
var

# Coordinates
print ("Coordinates for the variables for scatter plot ")
head(var$coord)

print ("Correlations between variables and dimensions")
head(var$cor)
print ("The correlation between a variable and a principal component (PC) is used as the coordinates of the variable on the PC. The representation of variables differs from the plot of the observations: The observations are represented by their projections, but the variables are represented by their correlations (Abdi and Williams 2010).")
print ("Positively correlated variables are grouped together.
       Negatively correlated variables are positioned on opposite sides of the plot origin (opposed quadrants).
       The distance between variables and the origin measures the quality of the variables on the factor map. Variables that are away from the origin are well represented on the factor map.
       variables with low cos2 values will be colored in “white”
       variables with mid cos2 values will be colored in “blue”
       variables with high cos2 values will be colored in red")

fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE )

# Cos2: quality on the factore map
print ("Cos2 for the variables  represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.")
print ("A high cos2 indicates a good representation of the variable on the principal component. In this case the variable is positioned close to the circumference of the correlation circle.
       
       A low cos2 indicates that the variable is not perfectly represented by the PCs. In this case the variable is close to the center of the circle.
       For a given variable, the sum of the cos2 on all the principal components is equal to one.")
head(var$cos2)
corrplot(var$cos2, is.corr=FALSE, title = "Cos2")
# 
print ("Contributions  (in percentage) of each variable to  principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
       The contributions of variables in accounting for the variability in a given principal component are expressed in percentage.
       
       Variables that are correlated with PC1 (i.e., Dim.1) and PC2 (i.e., Dim.2) are the most important in explaining the variability in the data set.
       Variables that do not correlated with any PC or correlated with the last dimensions are variables with low contribution and might be removed to simplify the overall analysis.
       
       ")
head(var$contrib,10)
#si tengo pocas
corrplot(var$contrib, is.corr=FALSE, title = "Contribution") 

#si tengo muchas variables
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)

#contributions to PC1 y PC2
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)