random_network<-function(kegg_random,r,p){
  adj.final1<-as.matrix(kegg_random)
  graph1 = graph.adjacency(adj.final1,mode=c("undirected"), weighted=TRUE,add.rownames=TRUE)
  temp1 = page.rank(graph1, vids=V(graph1), directed=FALSE, damping=r, weights=NULL, options = list(niter = 10^10, eps = p))
  rank1 = temp1$vector
  rank2 = as.matrix(rank1)
  return(rank2)
}



#'@title DrugReposition
#'@description The function "DrugReposition" is used in drug repositioning by calculating the eigenvector centrality of drugs.
#'@param DE A matrix with one column of zscore.
#'@param nperm Number of random permutations (default: 1000).
#'@param r Restart the probability of the random-walk algorithm (default: 0.9).
#'@param p For each node, if the difference in centrality score between iterations changes less than this value, the algorithm considers the calculation complete (default: 10^-10).
#'@return A dataframe with seven columns those are drugbankid, centralscore, p.value,fdr,number of targets, drug targets,drugname.
#'@importFrom igraph graph.adjacency
#'@importFrom igraph V
#'@importFrom igraph page.rank
#'@importFrom stats p.adjust
#'@importFrom stats median
#'@importFrom fastmatch fmatch
#'@importFrom utils setTxtProgressBar
#'@importFrom utils txtProgressBar
#'@usage DrugReposition(DE,nperm = 1000,r = 0.9,p = 10^-10)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-Gettest("GEP")
#' label<-Gettest("label")
#' # Run the function
#' DEscore<-CalDEscore(GEP,label)
#' # Run the function
#' \donttest{drug_centrality<-DrugReposition(DE=DEscore,nperm = 1000,r = 0.9,p = 10^-10)}

DrugReposition<-function(DE,nperm = 1000,r = 0.9,p = 10^-10){
  library(igraph)
  go_p_score_row<-function(genecount,descore,path_size){
    score_row<-rep(0,path_size)
    aaa<-descore[,1]
    bbb<-descore[,2]

    for(j in 1:length(genecount)){
      if(genecount[j]!=""){
        gene<-unlist(strsplit(genecount[j], split = ","))
        location<-fmatch(gene,aaa)
        if(length(which(is.na(location)))>0){
          de_score1<-0
        }else{
          location<-location[which(location!="NA")]
          de1<-bbb[location]
          de_score1<-median(as.numeric(de1))
        }

        if (!is.na(de_score1)) {
          score_row[j]<-de_score1
        }
      }
      #print(j)
    }
    return(score_row)
  }
  old <- options()
  on.exit(options(old))
  options(stringsAsFactors = FALSE)
  drugname<-Gettest("drugname")
  jaccard<-Gettest("Jaccard")
  colGO<-Gettest("commongenes")
  Go<-Gettest("GO_MF")
  drug_genes<-Gettest("Drugs")
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
  adj.final<-as.matrix(drug_drug)
  diag(adj.final)<-0
  graph = graph.adjacency(adj.final,mode=c("undirected"), weighted=TRUE,add.rownames=TRUE)
  temp = page.rank(graph, vids=V(graph), directed=FALSE, damping=r, weights=NULL,options = list(niter = 10^10, eps = p))
  rank = temp$vector
  rank1 = as.matrix(rank)
  #####扰动
  subpathway<-as.character(colnames(adj.final))
  if (nperm == 0) {
    p <- rep(1, length(subpathway[1]))
    fdr <- p
    allresult <- cbind(subpathway, rank1)
    allresult <- cbind(allresult, p)
    allresult <- cbind(allresult, fdr)
    allresult[, 1] <- as.character(allresult[, 1])
    changdu<-c()
    for (i in 1:length(adj.final[,1])) {
      m = length(unlist(strsplit(drug_genes[i,3],split = ",")))
      changdu<-c(changdu,m)

    }
    result1<-cbind(allresult,changdu)
    result1<-cbind(result1,drug_genes[,3])
    colnames(result1) = c("DBID","Centrality","P.value","FDR","Target_length","Target")

    result1<-as.data.frame(result1)
    result1[,2]<-as.numeric(as.character(result1[,2]))
    result1[,2]<-round(result1[,2],5)
    result1[,3]<-as.numeric(result1[,3])
    result1[,4]<-as.numeric(result1[,4])
    result1[,5]<-as.numeric(result1[,5])

    part1<-result1[grep("_",result1[,1]),]
    part2<-result1[-grep("_",result1[,1]),]
    part1_str<-as.data.frame(matrix(nrow = 1,ncol = 6))
    for (i in 1:length(part1[,1])) {
      a<-length(unlist(strsplit(part1[i,1],split = "_")))
      b<-matrix(nrow = a,ncol = length(result1[1,]))
      b<-as.data.frame(b)
      b[1:a,] = part1[i,]
      b[,1] = unlist(strsplit(part1[i,1],split = "_"))
      part1_str<-rbind(part1_str,b)
    }
    colnames(part1_str) = c("DBID","Centrality","P.value","FDR","Target_length","Target")
    result2<-rbind(part2,part1_str[-1,])
    colnames(drugname)[1] = "DBID"
    result2<-merge(result2,drugname,by = "DBID")
    result2<-result2[,c(1,7,2,3,4,6,5)]
    result2<-result2[order(result2[,3],decreasing = TRUE),]
    return(result2)
  }else{
    iter<-nperm
    Centrality_Scores<-matrix(nrow=path_size,ncol=iter+1)
    real.centra<-rank1
    Centrality_Scores[,1]<-real.centra

    print("Note: Completing 1000 perturbations may be a time-consuming task.")
    pb <- txtProgressBar(style=3)

    for (i in 1:iter) {

      DEnames <- sample(rownames(DE), replace = FALSE)
      rownames(DE) <- DEnames
      DE[,1]<-DEnames
      median_score <- matrix(0, nrow = path_size, ncol = go_size)
      for (k in 1:path_size) {
        con_gene <- go_path_gene[k, ]
        row <- go_p_score_row(con_gene, DE, go_size)
        median_score[k, ] <- row
        #print(k)
      }
      rownames(median_score) = drug_genes[, 1]
      colnames(median_score) = Go[, 1]
      # jaccard <- t(jaccard)
      go_drug <- median_score * jaccard
      edge <- as.matrix(go_drug)
      edget <- t(edge)
      drug_drug <- edge %*% edget
      adj.final <- as.matrix(drug_drug)
      diag(adj.final) <- 0
      graph = graph.adjacency(adj.final, mode = c("undirected"),
                              weighted = TRUE, add.rownames = TRUE)
      temp = page.rank(graph, vids = V(graph), directed = FALSE,
                       damping = r, weights = NULL, options = list(niter = 10^10,
                                                                   eps = p))
      rank = temp$vector
      rank1 = as.matrix(rank)

      Centrality_Scores[, i + 1] <- rank1
      setTxtProgressBar(pb, i/iter)
    }
    close(pb)
    adj = as.matrix(Centrality_Scores)
    perm_rank = adj[,2:(iter+1)]
    perm_rank<-as.matrix(perm_rank)
    orig_rank = adj[,1]
    pval<-c()

    for (i in 1:length(orig_rank)) {

      pvall <- sum(perm_rank[i, ] >= orig_rank[i])/nperm
      pval <- c(pval,pvall)

      print(i)
    }

    p_padjust <- p.adjust(pval, method = "fdr")
    pa <- as.numeric(p_padjust)
    fdr <- round(pa, 3)
    allresult<-cbind.data.frame(subpathway,Centrality_Scores[,1])
    allresult<-cbind.data.frame(allresult,pval)
    allresult<-cbind.data.frame(allresult,fdr)
    allresult[,1]<-as.character(allresult[,1])
    changdu<-c()
    for (i in 1:length(drug_drug[,1])) {
      m = length(unlist(strsplit(drug_genes[i,3],split = ",")))
      changdu<-c(changdu,m)

    }
    result1<-cbind(allresult,changdu)
    result1<-cbind(result1,drug_genes[,3])
    colnames(result1) = c("DBID","Centrality","P.value","FDR","Target_length","Target")

    result1<-as.data.frame(result1)
    result1[,2]<-as.numeric(as.character(result1[,2]))
    result1[,2]<-round(result1[,2],5)
    result1[,3]<-as.numeric(result1[,3])
    result1[,4]<-as.numeric(result1[,4])
    result1[,5]<-as.numeric(result1[,5])

    part1<-result1[grep("_",result1[,1]),]
    part2<-result1[-grep("_",result1[,1]),]
    part1_str<-as.data.frame(matrix(nrow = 1,ncol = 6))
    for (i in 1:length(part1[,1])) {
      a<-length(unlist(strsplit(part1[i,1],split = "_")))
      b<-matrix(nrow = a,ncol = length(result1[1,]))
      b<-as.data.frame(b)
      b[1:a,] = part1[i,]
      b[,1] = unlist(strsplit(part1[i,1],split = "_"))
      part1_str<-rbind(part1_str,b)
    }
    colnames(part1_str) = c("DBID","Centrality","P.value","FDR","Target_length","Target")
    result2<-rbind(part2,part1_str[-1,])
    colnames(drugname)[1] = "DBID"
    result2<-merge(result2,drugname,by = "DBID")
    result2<-result2[,c(1,7,2,3,4,6,5)]
    result2<-result2[order(result2[,4],decreasing = F),]
    return(result2)
}
}

