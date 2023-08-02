
#' @title Gettest
#' @description Get the example data
#' @param exampleData  A character,should be one of"Jaccard","commongenes","GO_MF","Drugs","Drugbankid_CID","drugname","GEP","label"
#' @return data
#' @export

Gettest<-function(exampleData){

  utils::data("myenv",package="DrugSim2DR")

  if (exampleData=="Jaccard")
  {
    dataset<- get("Jaccard",envir=myenv)
    return(dataset)
  }
  if (exampleData=="commongenes")
  {
    dataset<- get("commongenes",envir=myenv)
    return(dataset)
  }
  if (exampleData=="GO_MF")
  {
    dataset<- get("GO_MF",envir=myenv)
    return(dataset)
  }
  if (exampleData=="Drugs")
  {
    dataset<- get("Drugs",envir=myenv)
    return(dataset)
  }
  if (exampleData=="Drugbankid_CID")
  {
    dataset<- get("Drugbankid_CID",envir=myenv)
    return(dataset)
  }
  if (exampleData=="GEP")
  {
    dataset<- get("GEP",envir=myenv)
    return(dataset)
  }
  if (exampleData=="label")
  {
    dataset<- get("label",envir=myenv)
    return(dataset)
  }
  if (exampleData=="drug_centrality")
  {
    dataset<- get("drug_centrality",envir=myenv)
    return(dataset)
  }
  if (exampleData=="Drugs_ind")
  {
    dataset<- get("Drugs_ind",envir=myenv)
    return(dataset)
  }
  if (exampleData=="commongenes_ind")
  {
    dataset<- get("commongenes_ind",envir=myenv)
    return(dataset)
  }
  if (exampleData=="Jaccard_ind")
  {
    dataset<- get("Jaccard_ind",envir=myenv)
    return(dataset)
  }
  if (exampleData=="drug_drug")
  {
    dataset<- get("drug_drug",envir=myenv)
    return(dataset)
  }
  if (exampleData=="drugname")
  {
    dataset<- get("drugname",envir=myenv)
    return(dataset)
  }
}
