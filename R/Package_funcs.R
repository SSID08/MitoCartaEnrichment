#' Perform fisher the fisher test per pathway in the Mitocarta. This is a private function
#' @param x List of genes/ENTREZ Ids in the pathway
#' @param y list of genes/ENTREZ Ids from the set in question/of interest
#' @param z Background/Gene universe set
#' @import stats
#' @return Fisher test results and list of intersection between DE and pathway genes
#' @noRd
mito_enrichment_function=function(x,y,z){
  #x=list of genes in pathway/ontology
  #y=list of genes of interest (DE proteins for instance)
  #z=background set
  x=unlist(stats::na.omit(x))
  y=stats::na.omit(y)
  z=stats::na.omit(z)
  k=(setdiff(z,y)) # not DE proteins
  DE_in_pathway=length(intersect(x,y))
  no_DE_in_pathway=length(intersect(k,x))
  DE_not_in_pathway=length(setdiff(y,x))
  no_DE_not_in_pathway=length(setdiff(k,x))
  contingency_matrix=matrix(c(DE_in_pathway,no_DE_in_pathway,DE_not_in_pathway,no_DE_not_in_pathway),nrow = 2,byrow = T)
  test=stats::fisher.test(contingency_matrix,alternative = 'greater')
  gene_ratio=DE_in_pathway/length(x)
  ret_val=list(c(intersect(x,y)),gene_ratio,test$estimate,test$p.value)
  #names(ret_val)=c('Genes of interest in Pathway','Odds Ratio','P-value')
  return(ret_val)
}

#' Run enrichment across the different pathways based on the the input gene set of interest and
#' background set
#' @import rlang
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import stringr
#' @import dplyr
#' @import tidyverse
#' @param x List of ENTREZ Ids of interest
#' @param y List of ENTREZ Ids of background set
#' @param q_threshold FDR threshold for retaining enrichment terms. Default is 0.05.
#' To retain results across all pathway terms, use 1.
#' @return Dataframe of significant enrichment results at FDR< 0.05
#' @export
run_enrichment=function(x,y,q_threshold=0.05){
  test_results=lapply(X = pathways$ENTREZ,FUN = mito_enrichment_function,y=x,z=y)
  test_results=data.frame(do.call(rbind, test_results))
  test_results$Pathway=pathways$MitoPathway
  test_results=test_results%>%dplyr::filter(X3!=0)
  colnames(test_results)=c('ENTREZ IDs','Gene-ratio','Odds-ratio','P-values','Pathway')
  test_results=test_results%>%dplyr::select(Pathway,`ENTREZ IDs`,`Odds-ratio`,`P-values`,`Gene-ratio`)
  test_results$`Odds-ratio`=as.numeric(test_results$`Odds-ratio`)
  test_results$`P-values`=as.numeric(test_results$`P-values`)
  test_results$`Gene-ratio`=as.numeric(test_results$`Gene-ratio`)
  test_results$adj.P=stats::p.adjust(test_results$`P-values`,method = 'BH')
  test_results_final=test_results%>%dplyr::filter(adj.P<q_threshold)
  test_results_final$Gene_Names=apply(test_results_final,FUN = function(x){
    Gene_Symbols=AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=x["ENTREZ IDs"][[1]],keytype ="ENTREZID",column="SYMBOL")
    names(Gene_Symbols)=NULL
    #print(Gene_Symbols)
    return(Gene_Symbols)
  },simplify = T,MARGIN = 1)
  test_results_final$Gene_Names=sapply(test_results_final$Gene_Names,FUN = paste,collapse=',')
  test_results_final$`ENTREZ IDs`=sapply(test_results_final$`ENTREZ IDs`,FUN = paste,collapse=',')
  #test_results_final=test_results_final%>%select(-Gene_namees_split)
  test_results_final=test_results_final%>%dplyr::mutate(Gene_count = stringr::str_count(Gene_Names, pattern = ',') + 1)
  test_results_final=test_results_final%>%dplyr::mutate_at(c("Odds-ratio","Gene-ratio","P-values"),~round(.x,3))
  return(test_results_final)
}


#' Convert input biological IDs to ENTREZ Ids
#' @param x vector of input IDs
#' @param from string to denote input ID type
#' @description For a list of acceptable input IDs and acceptable strings; use valid_ids()
#' @return Vector of input IDs converted to ENTREZ IDs
#' @export
#'
ID_to_ENTREZ = function(x, from) {
  tryCatch(
    {
      out_vec <- AnnotationDbi::mapIds(x = org.Hs.eg.db::org.Hs.eg.db, keys = x, keytype = from, column = 'ENTREZID')
      return(out_vec)
    },
    error = function(e) {
      message("Error: Unable to map IDs to ENTREZ IDs -", e$message)
      return(NULL)
    }
  )
}

#' Check valid input IDs and corresponding strings for converting to ENTREZ
#' @description Use one of the valid input IDs returned by this function to output ENTREZ IDs
#' for enrichment analysis
#' @export
valid_input=function(){
  return(AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db))
}

#' Plot enrichment results
#' @import ggplot2
#' @param enrichment_results The dataframe of the enrichment analysis results from run_enrichment
#' @param n The number of top 'n' terms to include in the plot
#' @param text_size The size of axis tick labels on both the X and Y axis
#'
#' @description Plot results of enrichment analysis. By default, the ten most significant terms are plotted
#' @export
barplot_enrichment=function(enrichment_results,n=NULL,text_size=10){
  enrichment_results=enrichment_results%>%dplyr::arrange(adj.P)
  if (is.null(n)){
    n=dim(enrichment_results)[1]
  }
  else if(n>dim(enrichment_results)[1]){
    n=dim(enrichment_results)[1]
  }
  else{
    n=n
  }

  plot=ggplot(data=enrichment_results[c(1:n),])+geom_bar(mapping = aes(x=reorder(Pathway,`Gene-ratio`),y=`Gene-ratio`,fill=adj.P),stat = 'identity')+
    scale_fill_gradient(low='blue',high='red',name='Adjusted P value')+
    xlab('Pathway')+ylab('Gene Ratio')+coord_flip()+theme_minimal()+
    theme(legend.title = element_text(size=10),legend.text = element_text(size=8),axis.title = element_text(size=text_size+2),axis.text = element_text(size=text_size))
  return((plot))
}
