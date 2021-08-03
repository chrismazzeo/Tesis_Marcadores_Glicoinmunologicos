pieGraphic = function(var, title="Title", path = NULL){
  
  var = na.omit(var)
  var.table = table(var)
  if (!is.null(path)){
    pdf(path) 
  }
  labels = paste0(names(var.table)," (n = ",var.table," ,", round(var.table/length(var),2)*100,"%)")
  pie(var.table,main = title,labels = labels,   col = settings.color[1:length(names(var.table))])
  if (!is.null(path)){
    dev.off() 
  }
  #legend("topright",legend=toupper(names(gender.counts)),bty="n", fill = colors)
}
barplotGraphic = function (var,breaks = NULL,labels = NULL,beside,  title="Title", xTitle="", yTitle="", path = NULL){
  var = na.omit(var)
  
  
  
  if (!is.null(breaks)){
    breaks = cut(as.numeric(as.character(var[,2])),breaks = breaks)
    var.table = table(var[,1], breaks)
  }
  else{
    var.table = table(var)
  }
  
  if (!is.null(labels)){
    colnames(var.table) = labels
  }
  
  if (!is.null(path)){
    pdf(path) 
  }
  
  barplot(var.table,  beside = beside,main =title ,xlab = xTitle, ylab = yTitle, names.arg=colnames(var.table), col=settings.color[1:length(rownames(var.table))],legend.text = rownames(var.table))
  
  if (!is.null(path)){
    dev.off() 
  }
}