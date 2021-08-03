  #Immunoinfiltraciones
  #*******************************************************************************************************
  #*********************
  immuCCPath = paste0(project.resultsPath,"ImmuCC/")
  dir.create(file.path(immuCCPath), showWarnings = FALSE)

  #ImmuCC RNASeq
  #Nota immuneCCQuantileNormalization hay que guardarlo como csv y subirlo a la web  http://218.4.234.74:3200/immune/# para procesarlo y descargar el resultado y luego analizarlo
  #*********************
  immuCCQuantileNormalization = IMMUNE_ImmuCC_RNASEQ(rawExpressionMatrix = htseq.rawCountMatrix, outputPath = paste0(immuCCPath,"immuneCCQuantileNormalization.csv"))
  #*******************
  #SVR
  immuCCResult.svr = IMMUNE_ImmuCC_RNASEQ_STATS(filePath = project.immucc.svrPath,
                            outputDir = immuCCPath ,
                            title =  "ImmuCC SVR",
                            samplesName = samplesList$SAMPLE,
                            group = samplesList$GROUP,
                            decovolutionalMethod = "SVR"
  )

  #LLSR
  immuCCResult.LLSR = IMMUNE_ImmuCC_RNASEQ_STATS(filePath = project.immucc.llsrPath,
                             outputDir = immuCCPath ,
                             title =  "ImmuCC LLSR",
                             samplesName = samplesList$SAMPLE,
                             group = samplesList$GROUP,
                             decovolutionalMethod = "LLSR"
  )
  

  #*******************************************************************************************************




