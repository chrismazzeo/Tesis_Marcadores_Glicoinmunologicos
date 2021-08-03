#http://www.sthda.com/english/wiki/survival-analysis-basics
#http://www.sthda.com/english/wiki/cox-proportional-hazards-model
# #************************************************************************************************************************************************************************
# La tabla se arma, con los dias hasta un cierto evento por ejemplos, dias hasta la muerte, y de los que no se murieron, son se ponen los dias hasta last follow up
# 
# Remember, the Cox regression analyzes the continuous variable over the whole range of its distribution, 
# where the log-rank test on the Kaplan-Meier plot can change depending on how you categorize your continuous variable. 
# They’re answering a similar question in a different way: the regression model is asking, 
# “what is the effect of age on survival?”, while the log-rank test and the KM plot is asking,
# “are there differences in survival between those less than 70 and those greater than 70 years old?”.
# 
# kaplan con 1 variable categorica con 2 levels, 
# para variable continuas hay pasarlas a categoricas con la media o mediana o con surv_cutpoint() and surv_categorize()  
# 
# 
# cox varias variabls categoricas
# y automaticmente toma las variables continuas 
# )

#kaplan Meier 
  #They describe the survival according to one factor under investigation, but ignore the impact of any others.
  # hago de 1 variable cualitativa , 
  #las variables cuantitativas las tengo para sar discretas:
    #package para ver densidades (recomendado)
    #media, mediana, .25 qunatile -75quantile, kmeans con variable dummy 
#coxph
  #todas las variables juntas
  #cuantitativa o cualititavia
# HR = 1: No effect
# HR < 1: Reduction in the hazard
# HR > 1: Increase in Hazard

#******************************************************************************************************************************
#Survival   kaplan-meier
#******************************************************************************************************************************
GET_SURVIVAL_TIME = function (timeToLastFollow, vitalStatus, daysToDeath, otherVars){
  t = data.frame(Time = as.numeric(as.character(timeToLastFollow)),VitalStatus = as.character(vitalStatus), otherVars = otherVars,stringsAsFactors = FALSE)
  
  #asignamos a los que murieron los dias
  t$Time[t$VitalStatus == "dead"] = daysToDeath[t$VitalStatus == "dead"] 
  
  #removemos los que no tenga tiempo o tiempo 0 y los que  no tenga el vital status 
  t = t[!(is.na(t$Time) | t$Time <= 0 | t$VitalStatus == "--"),]
  t$VitalStatus[t$VitalStatus == "alive"] = 0
  t$VitalStatus[t$VitalStatus == "dead"] = 1
  t$VitalStatus = as.numeric(t$VitalStatus)
  colnames(t) = c("Time", "VitalStatus", colnames(otherVars))
  return (t)
  
}
GET_SURVIVAL_GROUP = function(var, mode = "MEDIAN", numberOfGroup){
#TODO Quantile

 if (numberOfGroup >2 | mode == "KMEANS"){
    data = data.frame(var,mean(var)) #agregamos una variable dummy
    cluster = kmeans(data, centers = numberOfGroup, iter.max = 100, nstart = 10)
    return (cluster$cluster)
  }
  if (mode == "MEDIAN"){
     median = median(var)
     data = data.frame(var, cluster = 1)
     data$cluster[var>median] = 2
     return (data$cluster)
  }
  
}

GET_KAPLAN_MEIER_PLOT = function (survivalData, outputDir){

   dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)
  
    for (i in 3:dim(survivalData)[2]){
      
      try({
            print (paste0(i, "    ",colnames(survivalData)[i]))
            tempData = data.frame (Time = survivalData$Time, VitalStatus = survivalData$VitalStatus, Variable = survivalData[i])
            colnames(tempData)[3] = colnames(survivalData)[i]
            tempData = tempData[!(tempData[,3] == "not reported" | tempData[,3] == "--" | is.na(tempData[,3])) ,]
            
           
             varContinuos = FALSE
            if (dim(tempData)[1] > 0){
              legendTitle = gsub("_", " ",colnames(tempData)[3])
              legendTitle = paste0(toupper(substring(legendTitle,1,1)),substring(legendTitle,2))
              tempData.continues = tempData
              #categorizamos si cuantitativa
                if (class(tempData[,3]) == "integer" | class(tempData[,3]) == "numeric"){
                  try({  
                       res.cut <- surv_cutpoint(tempData, time = "Time", event = "VitalStatus",
                                                variables = c(colnames(tempData)[3]))
                      
                        summary(res.cut)
                        categories = surv_categorize(res.cut)[3]
                        data.boxplot = data.frame(value = tempData[,3], categories)
                        colnames(data.boxplot) = c("value" ,"group")
                        
                        tempData <- surv_categorize(res.cut)
                        varContinuos = TRUE
                  })
                }
                
                surObj = Surv(time = tempData$Time, event =tempData$VitalStatus)
          
                fit.surv <- survfit(formula(paste("surObj ~ ",colnames(tempData)[3])), data = tempData )
                fit.surv$call$data = tempData
                 fit.surv$call$formula = formula(paste("surObj ~ ",colnames(tempData)[3]))
                
                pvalue = surv_pvalue(fit.surv, data = tempData)
                
             
                if  (pvalue$pval <0.05){
                  pdf(paste0(outputDir,colnames(tempData)[3],".pdf"), width = 13, height = 10)
                   if (varContinuos){
                    print(plot(res.cut, main = legendTitle, palette = "npg"))
                    print( boxplot (value ~ group , data = data.boxplot,  col=c("#f2aaa1","#abdde9"), main =colnames(tempData)[3]))
                    
                   }
                    try({
                     
                        tempData[,3] = as.factor(tempData[,3])
                        levels(tempData[,3]) =  paste0(toupper(substring(  levels(tempData[,3]),1,1)),substring(  levels(tempData[,3]),2))
                       
                        legendLabs = levels(tempData[,3])
                
                        plot = ggsurvplot(fit.surv, data = tempData,  conf.in = FALSE, pval = TRUE,  surv.plot.height = 0.2, 
                                          risk.table = "abs_pct",
                                          risk.table.y.text = FALSE,
                                          surv.median.line = "hv",
                                          fun =  "pct",
                                          title = "Overal Survival",
                                          xlab = "Time in days",
                                          ylab = "Surviving (%)",
                                          legend.title = legendTitle,
                                          legend.labs = legendLabs,
                                          font.legend = c(12,"plain", "black"),
                                          legend = "right",
                                          tables.height = 0.35,
                                          tables.theme =settings.graphics.theme,
                                          
                                          # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                                          # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                                          
                                          ggtheme = settings.graphics.theme # Change ggplot2 
                        )
                        print(plot)
                        
                     
                    })
                  # #la tomamos como continua en vez de pasarla a categorica
                  # if (varContinuos | length(levels(tempData[,3])) >2){
                  #      fit.coxph <- coxph(formula(paste("surObj ~ ",colnames(tempData.continues )[3])) , data = tempData.continues)
                  #      plot2 = ggforest(fit.coxph, data = tempData.continues)
                  #      print(plot2)
                  # }
                  
                  dev.off()
                }
           
            }
      })
  }
  
}

GET_OVERALL_SURVIVAL = function(survivalData, outputDir){
  dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)
  surObj = Surv(time = survivalData$Time, event =survivalData$VitalStatus)
  fit.surv <- survfit(surObj ~  1 , data = survivalData )
  fit.surv$call$data = survivalData
  fit.surv$call$formula = formula(surObj ~  1)

  pdf(paste0(outputDir,"overallSurvival.pdf"), width = 13, height = 10)
  plot = ggsurvplot(fit.surv, data = survivalData,  conf.in = FALSE, pval = FALSE,  surv.plot.height = 0.2, 
                    risk.table = "abs_pct",
                    risk.table.y.text = FALSE,
                    surv.median.line = "hv",
                    fun =  "pct",
                    title = "Overal Survival",
                    xlab = "Time in days",
                    ylab = "Surviving (%)",

                    font.legend = c(12,"plain", "black"),
                    legend = "right",
                    tables.height = 0.35,
                    tables.theme =settings.graphics.theme,
                    
                    # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                    # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                    
                    ggtheme = settings.graphics.theme # Change ggplot2 
  )
  print(plot)
  dev.off()
 
  
}

GET_COXPH_PLOT = function(survivalData, outputDir){
  
  dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)
  
 vars = colnames(survivalData)[3:dim(survivalData)[2]]
 # res.cut <- surv_cutpoint(survivalData, time = "Time", event = "VitalStatus",
  #                         variables = c(vars), minprop = 0.0001)
  
#  summary(res.cut)
 # survivalData <- surv_categorize(res.cut)
  
  vars =  paste(vars, collapse = " + ")
  
  surObj = Surv(time = survivalData$Time, event =survivalData$VitalStatus)
  data.formula = formula(paste("surObj ~ ",vars))

  fit.coxph <- coxph(data.formula , data = survivalData)
  pdf(paste0(outputDir,"Coph_allVars.pdf"), width = 13, height = 10)
  plot2 = ggforest(fit.coxph, data = survivalData,fontsize = 0.5)
  print(plot2)
  dev.off()
}

#para hacer agrupados
# survivalData.filtered = data.frame(TimeToLastFollow = survivalData$Time,
#                                    VitalStatus = survivalData$VitalStatus,
#                                    gender = survivalData$gender,
#                                    tumor_stage = survivalData$tumor_stage,
#                                    msi_status =  survivalData$msi_status,
#                                    lymphatic_invasion = survivalData$lymphatic_invasion,
#                                    loss_mismatch_repair_protein = survivalData$loss_expression_of_mismatch_repair_proteins_by_ihc)
# dim(survivalData.filtered)
# survivalData.filtered = na.omit(survivalData.filtered)              
# dim(survivalData.filtered)
# survivalData.filtered = survivalData.filtered[survivalData.filtered$msi_status != "indeterminate" ,]
# dim(survivalData.filtered)
# 
# surObj = Surv(survivalData.filtered$Time, survivalData.filtered$VitalStatus)
# 
# fit.surv <- survfit(surObj ~  tumor_stage + msi_status  + lymphatic_invasion , data = survivalData.filtered )
# 
# 
# 
# 
# ggsurv = ggsurvplot(fit.surv,   pval = TRUE,  surv.plot.height = 0.2,
#                     risk.table = TRUE, risk.table.col="strata", 
#                     legend = "right",
#                     ggtheme = theme_bw())
# 
# curv_facet <- ggsurv$plot + facet_grid(tumor_stage ~ msi_status)
# curv_facet





