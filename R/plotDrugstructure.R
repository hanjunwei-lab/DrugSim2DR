#' @title plotDrugstructure
#' @description The function "plotDrugstructure" can plot the chemical structure of a drug.
#' @param drugid A drugbank ID.
#' @return A chemical structure of specific drug
#' @importFrom rvest html_text
#' @importFrom ChemmineR read.SDFset
#' @importFrom rvest read_html
#' @importFrom sp plot
#' @usage plotDrugstructure(drugid = "")
#' @export
#' @examples
#' # Load depend package
#' library(ChemmineR)
#' library(rvest)
#' # Obtain molecular formula and visualize it.
#' \donttest{plotDrugstructure(drugid ="DB00780")}


plotDrugstructure<-function(drugid=NULL) {
  Drugbankid_CID<-Gettest("Drugbankid_CID")
  drugname<-Gettest("drugname")
  Drugs_CID<-Drugbankid_CID
  drugCid <- Drugs_CID[which(Drugs_CID[, 1] == drugid),2]
  drug_url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/",
                    drugCid, "/record/SDF/?record_type=2d&response_type=display",
                    sep = "")

  cw <- try(read_html(drug_url))
  if ("try-error" %in% class(cw)) {
    stop("Please ensure smooth network connection")
  }
  drugnr <- html_text(cw)
  drugnr <- strsplit(drugnr, "\n")
  drugnr <- unlist(drugnr)
  sdfset <- read.SDFset(drugnr)


  drugname<-drugname[which(drugname[,1]==drugid),2]
  sdfset@ID <- paste(drugname,paste("(",drugid,")",sep = ""),sep = " ")

  sp::plot(sdfset)

}

