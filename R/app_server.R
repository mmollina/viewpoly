#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # Your application server logic
  # Upload size
  options(shiny.maxRequestSize=5000*1024^2)
  
  ## Start modules
  datas <- callModule(mod_upload_server,
                      "upload_ui_1", 
                      parent_session=session)
  
  # QTL view
  callModule(mod_qtl_view_server,
             "qtl_view_ui_1",
             loadMap = datas$loadMap,
             loadQTL = datas$loadQTL,
             parent_session=session)
  
  # Genes view
  callModule(mod_genes_view_server,
             "genes_view_ui_1", 
             loadMap = datas$loadMap,
             loadQTL = datas$loadQTL,
             loadJBrowse_fasta = datas$loadJBrowse_fasta, 
             loadJBrowse_gff3 = datas$loadJBrowse_gff3, 
             loadJBrowse_vcf = datas$loadJBrowse_vcf, 
             loadJBrowse_align = datas$loadJBrowse_align, 
             loadJBrowse_wig = datas$loadJBrowse_wig, 
             parent_session=session)
  
  # Map view
  callModule(mod_map_view_server,
             "map_view_ui_1", 
             loadMap = datas$loadMap,
             loadQTL = datas$loadQTL,
             parent_session=session)
}
