
#******
# geneLength must be in kb
#******
rawCountsTorpkm<- function(counts, geneLength) {
  rate <- colSums(counts)/ 10^6
  a = counts / rep(rate,each=nrow(counts))
  rpkm = a/rep(geneLength,times=ncol(a))

  return (rpkm)
}
rawCountsTotpm<- function(counts, geneLength) {
  a = counts / rep(geneLength,times=ncol(a)) 
  rate = colSums(a)/10^6
  tpm = a/rep(rate,each=nrow(a))

  return (tpm)
}

FPKMtoTPM = function (counts){
  b = colSums(counts)
  tpm = counts/rep(b,each=nrow(counts))*10^6
  
  return (tpm)
}


#********************************************************
#test cases

# geneLength = matrix(c(2,4,1,10), ncol=1, byrow = TRUE)
# matrixCounts = matrix(c(10,20,5,0,12,25,8,0,30,60,15,1),nrow = 4, ncol = 3)

# a = rawCountsTorpkm(matrixCounts,geneLength)
# a
# rawCountsTotpm(matrixCounts,geneLength)
# 
# rpkmToTpm(a)

