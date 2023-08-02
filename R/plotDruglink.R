
#'@title plotDruglink
#'@description The function "plotDruglink" is used to plot a bipartite network of drugs and shared molecular functions.
#'@param drug1 The drugbank ID of drug1.
#'@param drug2 The drugbank ID of drug2.
#'@param i Specifies the number of outputs molecular functions, which is 5 by default.
#'@param color_MF Defines the color of MF nodes in the network.
#'@param color_drug Defines the color of drug nodes in the network.
#'@param layout_type layout_type used to set the appropriate arrangement, there is an option to choose from "circle","dh",and "sugiyama".
#'@return A bipartite network of drugs and shared molecular functions.
#'@importFrom igraph graph_from_data_frame
#'@importFrom igraph layout_with_dh
#'@usage plotDruglink(drug1,drug2,i = 5,color_MF = "#43AAEF",color_drug = "#F7525B",
#'layout_type = "circle")
#'@export
#'@examples
#' # Set drug1
#' drug1<-"DB02721"
#' # Set drug2
#' drug2<-"DB01213"
#' # Run the function
#' library(igraph)
#' plotDruglink(drug1,drug2,i = 5)



plotDruglink<-function(drug1,drug2,i = 5,color_MF = "#43AAEF",color_drug = "#F7525B",layout_type = "circle"){
  library(igraph)
  drugname<-Gettest("drugname")
  col<-Gettest("commongenes_ind")
  if (layout_type == "circle") {
    a = layout.circle
  }else if (layout_type == "dh") {
    a = layout_with_dh
  }else if (layout_type == "sugiyama") {
    a = layout_with_sugiyama
  }else{stop("No valid output format was specified")}
  drug1_mf<-colnames(col)[col[grep(drug1,rownames(col)),] != ""]
  drug2_mf<-colnames(col)[col[grep(drug2,rownames(col)),] != ""]
  if(length(grep(drug1,rownames(col))) == 0||length(grep(drug2,rownames(col))) == 0){
    stop("The drug was not identified")
  }else{
    if (length(intersect(drug1_mf,drug2_mf)) == 0)  {
      stop("The two drugs are not shared GOMF")
    }else{
      if (i >length(intersect(drug1_mf,drug2_mf))) {
        drug1<-drugname[which(drugname[,1] == drug1),2]
        drug2<-drugname[which(drugname[,1] == drug2),2]
        actors <- data.frame(name=c(intersect(drug1_mf,drug2_mf),drug1,drug2))
        charact<-c(rep("MF",length(intersect(drug1_mf,drug2_mf))),rep("DRUG",2))
        actors<-cbind(actors,charact)
        relations <- data.frame(from=c(rep(drug1,length(actors[,1])-2),rep(drug2,length(actors[,1])-2)),
                                to=c(rep(intersect(drug1_mf,drug2_mf),2)))
        g <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
        my_color<-c(rep(color_MF,length(intersect(drug1_mf,drug2_mf))),rep(color_drug,2))
        #coul  <- brewer.pal(3, "Set1")
        #my_color <- coul[as.numeric(as.factor(V(g)$charact))]
        plot(g, layout = a,vertex.color=my_color)
      }else{
        drug1<-drugname[which(drugname[,1] == drug1),2]
        drug2<-drugname[which(drugname[,1] == drug2),2]
        actors <- data.frame(name=c(intersect(drug1_mf,drug2_mf)[1:i],drug1,drug2))
        charact<-c(rep("MF",i),rep("DRUG",2))
        actors<-cbind(actors,charact)
        relations <- data.frame(from=c(rep(drug1,length(actors[,1])-2),rep(drug2,length(actors[,1])-2)),
                                to=c(rep(intersect(drug1_mf,drug2_mf)[1:i],2)))
        g <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
        #coul  <- brewer.pal(3, "Set1")
        #my_color <- coul[as.numeric(as.factor(V(g)$charact))]
        my_color<-c(rep(color_MF,i),rep(color_drug,2))
        plot(g, layout = a,vertex.color=my_color)
      }
    }
  }



}

