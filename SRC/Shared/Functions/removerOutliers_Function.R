
#esta funcion no esta probada
#si es para genes no se si hay que aplicarla, le saca variabilidad
removeSamplesOutliers = function(data, group){
  data = countMatrix
  group = samples$`Sample Type`
  i = 2
  
  
  grupos = levels(factor(group))
  
  #creamos una nueva matriz para ir guardando los objetos que no son ouliers
  final.data = data.frame()

  
  for (i in 1:length(grupos)){
    tempBd = data[group == grupos[i],]
    tempGroup =  group[group == grupos[i]]
    #calculamos la distancia con mahalanobis para ver outliers
    mahalanobisDistance<-mahalanobis(tempBd,colMeans(tempBd), cov(tempBd))
    #trazamos una linea por debajo p<=0.95, degreeof fredoom
    cutoff = qchisq(0.95, ncol(data), lower.tail = TRUE)

    #vemos los outliers
    outliers = which(mahalanobisDistance > cutoff)
    #los removemos
    subgrupo = subgrupo[-outliers,]
    bd.sinOutliers = rbind(subgrupo)
    
    
    print (paste0("*******Outliers del  GRUPO: ",grupos[i]))
    print (ouliers)
  }
  cat("\n")
  print("********************")
  print("Outliers a remover")
  print (paste0("Cantidad: " ,length(outlierList)))
  print(outlierList)
  
  #removemos los outliers de la bd
  data.noOutliers = data[!rownames(data) %in% outlierList,]
  
  #removemos los outliers de los grupos
  data.noOutliers.group = group[!rownames(as.data.frame(group)) %in% outlierList]
  
 
  return (outlierList)
  
}



removerOutliersFunc = function(data, group = NULL, graficos = TRUE){
  #calcuamos la disntacia de mahlanobis de los objetos
  if (is.null(group)){
    #dev.new()

      centroide = colMeans(data)
      
      mahalanobisDistance<-mahalanobis(data,centroide, cov(data))
      cutoff = qchisq(0.95, ncol(data), lower.tail = TRUE)
      if (graficos){
        plot(mahalanobisDistance,main ="Distancia de Mahalanobis")
        #trazamos una linea por debajo p<=0.95, degreeof fredoom
        abline(cutoff,0, col="red")
      }
      print("********************\n")
      print("Outliers a remover\n")
      outliersIndex = (which(mahalanobisDistance > cutoff))
      print(outliersIndex)
      data.noOutliers = data[!rownames(data) %in% outliersIndex,]
      
      #dev.new()
      par(mfrow = c(2,ncol(data)))
      #vemos la distribucion de las variables con los outliers antes de removerlos
      if (graficos){
          for (i in 1:ncol(data)){
              boxplot(data[i],main =colnames(data[i]),outpch=25, outcol="red", outbg="red")
          }
          
          #vemos como queda sin los outliers
          for (i in 1:ncol(data)){
              boxplot(data.noOutliers[i],main = "Sin outliers",outpch=25, outcol="red", outbg="red")
          }
      }
      return (outliersIndex)
  }
  else{
    grupos = levels(factor(group))
    centroides = NULL
    offset = 0
    outlierList  = c()
    par(mfrow = c(1,length(grupos)))
    for (i in 1:length(grupos)){
       subgrupo = subset(data, group == grupos[i])
       print (paste0("******* GRUPO: ",grupos[i]))
       centroides = rbind(centroides,colMeans(subgrupo))
       
       # dev.new()
       #distancia de mahalanobis para ver los outliers
       mahalanobisDistance<-mahalanobis(subgrupo,colMeans(subgrupo), cov(subgrupo))
       
       
       #trazamos una linea por debajo p<=0.95, degreeof fredoom
       cutoff = qchisq(0.95, ncol(data), lower.tail = TRUE)
       if (graficos){
          plot(mahalanobisDistance,main =paste0("Distancia de Mahalanobis -",grupos[i]))
          abline(cutoff,0, col="red")
       }
     
       #guardamos los outliers para cada grupo
       outliers = (which(mahalanobisDistance > cutoff))
       outlierList = c(outlierList,outliers+offset)
      #manejamos un offset por grupo
       offset =offset + nrow(subgrupo)
    }
    cat("\n")
    print("********************")
    print("Outliers a remover")
    print (paste0("Cantidad: " ,length(outlierList)))
    print(outlierList)
    
    #removemos los outliers de la bd
    data.noOutliers = data[!rownames(data) %in% outlierList,]
    
    #removemos los outliers de los grupos
    data.noOutliers.group = group[!rownames(as.data.frame(group)) %in% outlierList]
    
 
  #  grupos = levels(factor(data.noOutliers.group))
   # for (i in 1:length(grupos)){
     
    #  subgrupo = subset(data.noOutliers, data.noOutliers.group == grupos[i])
    #   mahalanobisDistance<-mahalanobis(subgrupo,colMeans(subgrupo), cov(subgrupo))
     #  plot(mahalanobisDistance,main =paste0("Distancia de Mahalanobis -  Sin outliers -",grupos[i]))
        #trazamos una linea por debajo p<=0.95, degreeof fredoom
      #cutoff = qchisq(0.95, ncol(data.noOutliers), lower.tail = TRUE)
      # abline(cutoff,0, col="red")
    #}
  
    #par(mfrow = c(2,ncol(data)))
    #for (i in 1:length(grupos)){
       # dev.new()
     #  subgrupo = subset(data, group == grupos[i])
       
       #vemos la distribucion de las variables con los outliers antes de removerlos
      # for (j in 1:ncol(subgrupo)){
       #  boxplot(subgrupo[j],main = paste0(colnames(subgrupo[j]),"-",grupos[i]),outpch=25, outcol="red", outbg="red")
    
       #}
       #vemos como queda sin los outliers
       # for (j in 1:ncol(subgrupo)){
        #  boxplot(data.noOutliers[j],main = "Sin outliers",outpch=25, outcol="red", outbg="red")
        #}
     #}
    grupos = c(grupos,"Total")
    centroides = rbind(centroides,colMeans(data))
    
  
    print(cbind("Grupo:",grupos, "Centroide",centroides))
    return (outlierList)
    
  }
  par(mfrow = c(1,1))
}
