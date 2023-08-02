#'@title plotTargetheatmap
#'@description The function "plotTargetheatmap" is used to plot a heat map of drug targets expression.
#'@param drugid The drugbank ID of a drug.
#'@param ExpData A gene expression profile of interest (rows are genes, columns are samples).
#'@param label A character vector consists of "0" and "1" which represent sample class in the gene expression profile. "0" means normal sample and "1" means disease sample.
#'@param significance This parameter controls whether the p-value of differential expression is displayed.
#'@param cluster.rows Logical value that represents whether row clustering is used.
#'@param cluster.cols Logical value that represents whether col clustering is used.
#'@param bk This parameter adjusts the range of values displayed by the color bar.
#'@param show.rownames This parameter controls whether row names are displayed.
#'@param show.colnames This parameter controls whether column names are displayed.
#'@param ann_colors Vector of colors used to define groups.
#'@param col Vector of colors used in the heatmap.
#'@return A heat map of drug targets expression.
#'@importFrom pheatmap pheatmap
#'@importFrom grDevices colorRampPalette
#'@importFrom tidyr unite
#'@importFrom stats wilcox.test
#'@usage plotTargetheatmap(drugid,ExpData,label,significance=FALSE,
#'cluster.rows=FALSE,cluster.cols=FALSE,bk=c(-2.4,2.3),show.rownames=TRUE,
#'show.colnames=FALSE,ann_colors=c("#FFAA2C","#2CBADA"),col=c("#2A95FF","#FF1C1C"))
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-Gettest("GEP")
#' label<-Gettest("label")
#' # Run the function
#' plotTargetheatmap("DB00780",GEP,label)


plotTargetheatmap<-function(drugid,ExpData,label,significance=FALSE,cluster.rows=FALSE,cluster.cols=FALSE,bk=c(-2.4,2.3),
                          show.rownames=TRUE,show.colnames=FALSE,ann_colors=c("#FFAA2C","#2CBADA"),
                          col=c("#2A95FF","#FF1C1C")){
  old <- options()
  on.exit(options(old))
  GEP<-ExpData
  drug.target<-Gettest('Drugs_ind')
  drugname<-Gettest('drugname')
  drug.target.frame<-drug.target
  sig.genes<-drug.target.frame[which(drug.target.frame[,1]==drugid),"Gene"]
  sig.genes<-unlist(strsplit(sig.genes,","))
  sig.genes<-intersect(sig.genes,rownames(GEP))
  GEP<-GEP[sig.genes,]
  colann<-data.frame(sample=colnames(GEP),group=label)
  colann.disease<-colann[which(colann[,2]=='1'),]
  colann.normal<-colann[which(colann[,2]=='0'),]
  colann<-rbind(colann.disease,colann.normal)
  rownames(colann)<-colann[,1]
  colann.1<-as.data.frame(colann[,-1])
  rownames(colann.1)<-rownames(colann)
  colnames(colann.1)<-'Group'
  GEP<-GEP[,match(rownames(colann.1),colnames(GEP))]
  colann.1[,1]<-as.character(colann.1[,1])
  colann.1[which(colann.1[,1]=="1"),1]<-"Disease"
  colann.1[which(colann.1[,1]=="0"),1]<-"Normal"

  if(significance==TRUE){
    p<-c()
    options(digits = 3,scipen = 0)
    for (i in 1:length(GEP[,1])) {
      p1<-wilcox.test(as.numeric(GEP[i,which(colann.1[,1]=="Disease")]),
                      as.numeric(GEP[i,which(colann.1[,1]=="Normal")]),alternative = "two.sided",paired = FALSE)[[3]]
      p1<-format(p1,scientific = TRUE)
      p1<-paste("(p=",p1,")",sep = "")
      p<-c(p,p1)
    }
    GEP<-cbind.data.frame(rownames(GEP),p,GEP)
    colnames(GEP)[1]<-"gene"
    GEP<-unite(GEP,"genes_p",gene,p,sep = " ",remove = TRUE)
    rownames(GEP)<-GEP[,1]
    GEP<-GEP[,-1]
  }
  drugname<-drugname[which(drugname[,1]==drugid),2]
  breaks<-c(seq(bk[1],-0.1,by=0.1),seq(0,bk[2],by=0.1))
  ann_colors=list(Group=c(Disease=ann_colors[1],Normal=ann_colors[2]))
  pheatmap(GEP,
           scale = 'row',
           color = colorRampPalette(c(col[1], "#FFFFFF", col[2]))(50),
           breaks = breaks,
           cluster_rows=cluster.rows,
           cluster_cols=cluster.cols,
           annotation_col =colann.1,
           annotation_colors = ann_colors,
           show_rownames=show.rownames,
           show_colnames=show.colnames,
           main=paste("Heatmap of the targets of",paste(drugname,paste("(",drugid,")",sep = ""),sep = " "),sep = " ")
  )

}

