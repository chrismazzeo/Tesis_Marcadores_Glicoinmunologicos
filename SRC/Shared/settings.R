set.seed(123)
#*******************************************************************************************************
if (project.type != "Clinical"){
    if (is.null(project.type) || is.null(project.resourcePath) || is.null(project.specie) || is.null(project.name) || is.null(project.resultsPath)){
      print ("**** SETUP PROJECT SETTINGS FIRS ****")
    }
}
#******************************************************************************************************************************
#key para exportar en ORCA
Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1IjoiY2hyaXNtYXp6ZW8iLCJhIjoiY2puZ2NsM3dxMDFycTNwbzUwa3ZrdXp1byJ9.epVfBOfjT3GSkZG5OqsXQA')
#******************************************************************************************************************************
setttings.coresNumber = detectCores(TRUE)-1
#******************************************************************************************************************************
#Correlation
#******************************************************************************************************************************
settings.correlation.cutoff = 0.6
#******************************************************************************************************************************
#Kmeans
#******************************************************************************************************************************
settings.kmeans.numberOfinteractions = 100
settings.kmeans.numbersOfRandomCenters = 20
#******************************************************************************************************************************
#heatmap Settings
#******************************************************************************************************************************
settings.heatmap.maxColumn = 130
settings.heatmap.naColor = "#black"
settings.heatmap.de.color=c("blue","black","red")
settings.heatmap.correlation.color=c("green","white","brown")
settings.heatmp.distance.colorRamp = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
settings.heatmap.cellWidth = 6
settings.heatmap.cellHeight = 3

#******************************************************************************************************************************
#Color Settings
#******************************************************************************************************************************
settings.graphics.theme = theme_classic()
settings.graphics.borderColor = scale_color_igv(palette = c("default", "alternating"))
settings.graphics.colors = rep(c(pal_igv("default")(51),pal_ucscgb()(26)),times = 3)
settings.graphics.fillColor = scale_fill_manual(settings.graphics.colors)



settings.color = c(  
  
  "chartreuse1","deeppink3","cyan3","yellow", "darkblue","darkorchid1","darkgoldenrod1",
  "chartreuse4","hotpink","cadetblue1","darkorchid4","coral3", "brown4","azure4","darkgoldenrod4",
  "forestgreen","deeppink4","cadetblue3","chocolate1","coral1","brown3","cornsilk4",
  "chartreuse3", "deeppink1","cyan1","darkviolet", "orangered", "chocolate3","azure3","blue",
  
  "darkolivegreen1","hotpink4","aquamarine3","brown1","darkgoldenrod3",
  "darkolivegreen3","cyan4","aquamarine4",  
  "darkolivegreen4", "cadetblue4","chocolate4","cornsilk1","black",
  "chartreuse1","deeppink3","cyan3","yellow", "darkblue","darkorchid1","darkgoldenrod1",
  "chartreuse4","hotpink","cadetblue1","darkorchid4","coral3", "brown4","azure4","darkgoldenrod4",
  "forestgreen","deeppink4","cadetblue3","chocolate1","coral1","brown3","cornsilk4",
  "chartreuse3", "deeppink1","cyan1","darkviolet", "orangered", "chocolate3","azure3","blue",
  
  "darkolivegreen1","hotpink4","aquamarine3","brown1","darkgoldenrod3",
  "darkolivegreen3","cyan4","aquamarine4",  
  "darkolivegreen4", "cadetblue4","chocolate4","cornsilk1","black",
  "chartreuse1","deeppink3","cyan3","yellow", "darkblue","darkorchid1","darkgoldenrod1",
  "chartreuse4","hotpink","cadetblue1","darkorchid4","coral3", "brown4","azure4","darkgoldenrod4",
  "forestgreen","deeppink4","cadetblue3","chocolate1","coral1","brown3","cornsilk4",
  "chartreuse3", "deeppink1","cyan1","darkviolet", "orangered", "chocolate3","azure3","blue",
  
  "darkolivegreen1","hotpink4","aquamarine3","brown1","darkgoldenrod3",
  "darkolivegreen3","cyan4","aquamarine4",  
  "darkolivegreen4", "cadetblue4","chocolate4","cornsilk1","black"
  
)

#******************************************************************************************************************************
#package description
settings.fileName.packageFileName = paste0(project.resultsPath,"sessionInfo.txt")
#******************************************************************************************************************************
if (project.type == "MICROARRAY"){
  source("./src/shared/microarray/Microarray.r")
}
if (project.type == "RNASEQ"){
  source("./src/shared/rnaseq/rnaseq.r")
}
if (project.type == "TCGA"){
  source("./src/shared/rnaseq/rnaseq.r")
}

