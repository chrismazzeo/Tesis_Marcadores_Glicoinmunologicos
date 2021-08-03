
  #************************************************************************************
  #  Merge the individual htseq result files into one
  htseq_merge <- function(){
      # Arguments:

      options(stringsAsFactors=F)
      files <- list.files()
      files <- files[grep("txt", files)]

      for (file in files) {
          cat(file, "\n")
          temp <- read.table(file, row.names=1, header=F, sep="")
          if (ncol(temp)==0) temp <- read.table(file, row.names=1, header=F, sep="\t")
          if (file==files[1]) {
              counts <- temp
              col <- ncol(temp)
              row <- nrow(temp)
          } else {
              counts <- cbind(counts, temp[rownames(counts), ])
              col <- c(col, ncol(temp))
              row <- c(row, nrow(temp))
          }
      }
      colnames(counts) <- gsub("\\.txt", "", files)
      counts
  }

  #************************************************************************************
  # Merge immune receptor genes into one
  receptor_merge <- function(expression, gene.list=receptor.ensemble.merge){
      # Arguments:
      #    expression:input expression matrix  
      #    gene.list: list of immune receptor genes

      gene.total <- as.character(unlist(gene.list))
      nonreceptor.expression <- expression[setdiff(rownames(expression), gene.total), ]
      family <- names(gene.list)

      receptor.expression <- c()
      for (gene in family) {
          gene.temp <- as.character(unlist(gene.list[gene]))
          gene.temp <- intersect(gene.temp, rownames(expression))
          cat("The family number of ", gene, " is ",length(gene.temp), "\n")
          if (length(gene.temp) > 1) {
              expression.temp <- apply(expression[gene.temp, ], 2, sum)
          } else {
              expression.temp <- expression[gene.temp, ]
          }
          receptor.expression <- rbind(receptor.expression, expression.temp)          
      }
      rownames(receptor.expression) <- family
      expression <- rbind(nonreceptor.expression, receptor.expression)
      expression
  }

  #************************************************************************************
  # quantile normalization
  quartile <- function(filename, p) {
      # Arguments:
      #    filename:input file name  
      #    p: quantile

      options(stringsAsFactors=F)
      if (is.data.frame(filename)|is.matrix(filename)) {
          data <- filename
      } else {
          if (file.exists(filename)) {
              if (grep("csv", filename)==1)      { data <- read.csv(filename, row.names=1)} 
              else if (grep("txt", filename)==1) { data <- read.table(filename, row.names=1,header=T)}
              else { break }
          }
      }

      n <- ncol(data)
      result <- matrix(nrow=nrow(data), ncol=0)
      for (i in seq(n)) {
          expres <- as.numeric(data[, i])
          value <- expres[expres!=0]
         
          value <- order(value, decreasing=T)
          value.quantile <- quantile(value, p)
          scale <- as.numeric(value.quantile)/1000
          expres <- ceiling(expres/scale)
          result <- cbind(result, expres)
      }

      colnames(result) <- colnames(data)
      rownames(result) <- rownames(data)
      return(result)
  }



  #################################################################################
  getRNASEQInfiltrations = function (expressionMatrix, receptorPath){
     load(receptorPath)

   counts.merge <- receptor_merge(expressionMatrix, gene.list=receptor.ensemble.merge)

    p <- 0.75
   counts.quatile <- quartile(counts.merge, p)
   return (as.data.frame(counts.quatile))
}
