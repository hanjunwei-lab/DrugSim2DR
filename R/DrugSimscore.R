
go_p_score_row<-function(genecount,descore,path_size){
  score_row<-rep(0,path_size)
  for(j in 1:length(genecount)){
    if(genecount[j]!=""){
      gene<-unlist(strsplit(genecount[j], split = ","))
      location<-match(gene,descore[,1])
      if(length(which(is.na(location)))>0){
        de_score1<-0
      }else{
        location<-location[which(location!="NA")]
        de1<-descore[location,2]
        de_score1<-median(as.numeric(de1))
      }

      if (!is.na(de_score1)) {
        score_row[j]<-de_score1
      }
    }
  }
  return(score_row)
}

#'@title DrugSimscore
#'@description The function "DrugSimscore" is used in calculating the drug functional similarity score.
#'@param DE A matrix with one column of zscore.
#'@param nperm Number of random permutations (default: 0).
#'@return A dataframe with four columns those are drug1, drug2, drug1 name, drug2 name, functional similarity score and FDR.
#'@importFrom stats p.adjust
#'@importFrom stats median
#'@importFrom stats pnorm
#'@importFrom stats sd
#'@importFrom reshape2 melt
#'@usage DrugSimscore(DE,nperm = 0)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-Gettest("GEP")
#' label<-Gettest("label")
#' # Run the function
#' DEscore<-CalDEscore(GEP,label)
#' # Run the function
#' \donttest{drug_drug<-DrugSimscore(DE=DEscore,nperm = 0)}

DrugSimscore<-function(DE,nperm = 0){
  old <- options()
  on.exit(options(old))
  options(stringsAsFactors = F)
  drugname<-Gettest("drugname")
  jaccard<-Gettest("Jaccard_ind")
  colGO<-Gettest("commongenes_ind")
  Go<-Gettest("GO_MF")
  drug_genes<-Gettest("Drugs_ind")
  Go_SubPath_gene<-colGO
  go_path_gene<-as.matrix(Go_SubPath_gene)
  path_size<-length(drug_genes[,1])
  go_size<-length(Go[,1])
  DE<-cbind(names(DE[,1]),abs(DE[,1]))
  median_score<-matrix(0,nrow=path_size,ncol=go_size)
  for(k in 1:path_size){
    con_gene<-go_path_gene[k,]
    row<-go_p_score_row(con_gene,DE,go_size)
    median_score[k,]<-row
  }
  rownames(median_score) = drug_genes[,1]
  colnames(median_score) = Go[,1]
  jaccard<-t(jaccard)
  go_drug<-median_score*jaccard
  edge<-as.matrix(go_drug)
  edget<-t(edge)
  drug_drug<-edge%*%edget
  ##需要r包reshape2
  drug_drug[lower.tri(drug_drug)] = 0
  drug_drug<-reshape2::melt(drug_drug, na.rm=TRUE)
  a<-drug_drug
  a<-drug_drug[which(drug_drug$value != 0),]
  a<-a[a[,1]!=a[,2],]
  colnames(a)<-c("drug1","drug2","similarity")
  colnames(drugname)[1] = "drug1"
  a<-merge(a,drugname,by = "drug1")
  colnames(a)<-c("drug1","drug2","similarity","drug1_name")
  colnames(drugname)[1] = "drug2"
  a<-merge(a,drugname,by = "drug2")
  colnames(a)<-c("drug1","drug2","similarity","drug2_name","drug1_name")
  drug_drug1<-a[,c(1,2,5,4,3)]
  #####扰动
  if (nperm == 0) {
    drug_drug1<-drug_drug1[order(drug_drug1[,5],decreasing = TRUE),]
    return(drug_drug1)
  }else{
    DElist<-list()
    DElist<-lapply(1:nperm, function(x){
      new_rowname<-rownames(DE)
      new_rowname<-sample(new_rowname,length(new_rowname),replace = FALSE)
      rownames(DE)<-new_rowname
      DE[,1]<-new_rowname
      return(DE)
    })
    result<-list()
    result<-lapply(DElist,function(x){
      median_score<-matrix(0,nrow=path_size,ncol=go_size)
      for(k in 1:path_size){
        con_gene<-go_path_gene[k,]
        row<-go_p_score_row(con_gene,x,go_size)
        median_score[k,]<-row
      }
      rownames(median_score) = drug_genes[,1]
      colnames(median_score) = Go[,1]
      go_drug<-median_score*jaccard
      edge<-as.matrix(go_drug)
      edget<-t(edge)
      drug_drug<-edge%*%edget
      ##需要r包reshape2
      drug_drug[lower.tri(drug_drug)] = 0
      drug_drug<-reshape2::melt(drug_drug, na.rm=TRUE)
      a<-drug_drug[which(drug_drug$value != 0),]
      a = a[a[,1] != a[,2],]
      return(a)
    })
    for(i in 1:nperm){
      a<-cbind(a,result[[i]][,3])
    }
    #####计算p值
    perm_rank<-a[,-c(1,2,3,4,5)]
    orig_rank<-a[,3]
    if (nperm == 1) {
      pval=pnorm(as.numeric(orig_rank),mean = mean(as.numeric(perm_rank)),sd=sd(as.numeric(perm_rank)),lower.tail = FALSE)
      p_padjust<-p.adjust(pval,method = "fdr")
      drug_drug<-cbind(a[,c(1,2,3)],p_padjust)
      colnames(drug_drug) = c("drug1","drug2","similarity","FDR")
      colnames(drugname)[1] = "drug1"
      drug_drug<-merge(drug_drug,drugname,by = "drug1")
      colnames(drug_drug)<-c("drug1","drug2","similarity","FDR","drug1_name")
      colnames(drugname)[1] = "drug2"
      drug_drug<-merge(drug_drug,drugname,by = "drug2")
      colnames(drug_drug)<-c("drug1","drug2","similarity","FDR","drug2_name","drug1_name")
      drug_drug<-drug_drug[order(drug_drug[,3],decreasing = TRUE),]
      drug_drug<-drug_drug[,c(1,2,6,5,3,4)]
      return(drug_drug)
    }else{
      pval=pnorm(as.numeric(orig_rank),mean = mean(apply(perm_rank,2,as.numeric)),sd=sd(apply(perm_rank,2,as.numeric)),lower.tail = FALSE)
      p_padjust<-p.adjust(pval,method = "fdr")
      drug_drug<-cbind(a[,c(1,2,3)],p_padjust)
      colnames(drug_drug) = c("drug1","drug2","similarity","FDR")
      colnames(drugname)[1] = "drug1"
      drug_drug<-merge(drug_drug,drugname,by = "drug1")
      colnames(drug_drug)<-c("drug1","drug2","similarity","FDR","drug1_name")
      colnames(drugname)[1] = "drug2"
      drug_drug<-merge(drug_drug,drugname,by = "drug2")
      colnames(drug_drug)<-c("drug1","drug2","similarity","FDR","drug2_name","drug1_name")
      drug_drug<-drug_drug[order(drug_drug[,3],decreasing = TRUE),]
      drug_drug<-drug_drug[,c(1,2,6,5,3,4)]
      return(drug_drug)
    }

  }

}

