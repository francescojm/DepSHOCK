
source('DepSHOCK_DE_analysis_library.R')

fn<-dir('../DepSHOCK_data/DiffExp/')
fn<-setdiff(fn,"MV411-FLT3_DE.Robj")

consensusUP<-lapply(fn,function(x){
  print(x)
  load(paste('../DepSHOCK_data/DiffExp/',x,sep=''))
  DE<-DS_split_DE_by_perturbation(DE_obj = DE_obj)
  
  NT_DE<-DE$sgAAVS1
  
  GenesUP<-unlist(names(which(rowSums(NT_DE$adjp<0.05 & NT_DE$lfc> 0.5)>0)))
})

alg<-sort(unique(unlist(consensusUP)))

membUP<-do.call(cbind,lapply(consensusUP,function(x){
  is.element(alg,x)
}))

rownames(membUP)<-alg

SUMMUP<-summary(as.factor(rowSums(membUP)))

consensusDOWN<-lapply(fn,function(x){
  print(x)
  load(paste('../DepSHOCK_data/DiffExp/',x,sep=''))
  DE<-DS_split_DE_by_perturbation(DE_obj = DE_obj)
  
  NT_DE<-DE$sgAAVS1
  
  GenesDOWN<-unlist(names(which(rowSums(NT_DE$adjp<0.05 & NT_DE$lfc< -0.5)>0)))
})

alg<-sort(unique(unlist(consensusDOWN)))

membDOWN<-do.call(cbind,lapply(consensusDOWN,function(x){
  is.element(alg,x)
}))

rownames(membDOWN)<-alg

SUMDOWN<-summary(as.factor(rowSums(membDOWN)))


