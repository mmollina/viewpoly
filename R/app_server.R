#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # Your application server logic
  # Upload size
  options(shiny.maxRequestSize=500*1024^2)
  
  ## Start modules
  datas <- callModule(mod_upload_server,
                      "upload_ui_1", 
                      parent_session=session)
  
  # QTL view
  callModule(mod_qtl_view_server,
             "qtl_view_ui_1",
             loadExample = datas$loadExample,
             loadMap_custom = datas$loadMap_custom,
             loadMap_mappoly = datas$loadMap_mappoly,
             loadQTL_custom = datas$loadQTL_custom,
             loadQTL_qtlpoly = datas$loadQTL_qtlpoly,
             loadQTL_diaQTL = datas$loadQTL_diaQTL,
             loadQTL_polyqtlR = datas$loadQTL_polyqtlR,
             parent_session=session)
  
  # Genes view
  callModule(mod_genes_view_server,
             "genes_view_ui_1", 
             loadExample = datas$loadExample,
             loadMap_custom = datas$loadMap_custom,
             loadMap_mappoly = datas$loadMap_mappoly,
             loadQTL_custom = datas$loadQTL_custom,
             loadQTL_qtlpoly = datas$loadQTL_qtlpoly,
             loadQTL_diaQTL = datas$loadQTL_diaQTL,
             loadQTL_polyqtlR = datas$loadQTL_polyqtlR,
             loadJBrowse_fasta = datas$loadJBrowse_fasta, 
             loadJBrowse_gff3 = datas$loadJBrowse_gff3, 
             loadJBrowse_vcf = datas$loadJBrowse_vcf, 
             parent_session=session)
  
  # Map view
  callModule(mod_map_view_server,
             "map_view_ui_1", 
             loadExample = datas$loadExample,
             loadMap_custom = datas$loadMap_custom,
             loadMap_mappoly = datas$loadMap_mappoly,
             loadQTL_custom = datas$loadQTL_custom,
             loadQTL_qtlpoly = datas$loadQTL_qtlpoly,
             loadQTL_diaQTL = datas$loadQTL_diaQTL,
             loadQTL_polyqtlR = datas$loadQTL_polyqtlR,
             parent_session=session)
}
