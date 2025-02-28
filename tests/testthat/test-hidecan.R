# test_that("Tests HIDECAN",{
#   
#   library(hidecan)
#   gwaspoly_file <- system.file("extdata/gwaspoly_res_thr.rda", package = "hidecan")
#   gwaspoly_file <- "~/Documents/Sweetpotato_Phill_GBS_DArT/GWAS_results/NormalizedBLUE_fb_DArTag_viewpoly.RData"
# 
#   input <- list()
#   input$gwaspoly$datapath <- c(gwaspoly_file, gwaspoly_file2, gwaspoly_file4)
#   input$gwaspoly$datapath <- gwaspoly_file
#   
#   for(i in 1:length(input$gwaspoly$datapath)){
#     temp <- load(input$gwaspoly$datapath[i])
#     gwaspoly_temp <- get(temp)
#     gwaspoly_list <- hidecan::GWAS_data_from_gwaspoly(gwaspoly_temp)
#     
#     if(i == 1) gwaspoly <- gwaspoly_list else {
#       if(!all(gwaspoly$chrom_length == gwaspoly_list$chrom_length)) {
#         # If same chromosome but different chromosomes length - keep the maximum
#         if(all(gwaspoly$chrom_length$chromosome == gwaspoly_list$chrom_length$chromosome)){
#           idx <- which(gwaspoly$chrom_length$length < gwaspoly_list$chrom_length$length)
#           if(length(idx) > 0)
#             gwaspoly$chrom_length$length[idx] <-  gwaspoly_list$chrom_length$length[idx]
#         } else stop("Not same reference genome used")
#       }
#       gwaspoly$gwas_data_list <- c(gwaspoly$gwas_data_list, gwaspoly_list$gwas_data_list)
#       gwaspoly$gwas_data_thr_list <- c(gwaspoly$gwas_data_thr_list, gwaspoly_list$gwas_data_thr_list)
#     }
#   }
#   
#   # Merging gwaspoly objects
#   
#   custom_files <- get_example_data()
#   loadHidecan <- list("GWASpoly" = gwaspoly,
#                       "GWAS" = list(GWAS_data(custom_files[["GWAS"]])),
#                       "DE" = list(DE_data(custom_files[["DE"]])),
#                       "CAN" = list(CAN_data(custom_files[["CAN"]])))
#   
#   x <- loadHidecan[["GWASpoly"]]$gwas_data_thr_list
#   
#   csv_names_gwas <- names(loadHidecan[["GWAS"]])
#   
#   chrom_length <- combine_chrom_length(
#     loadHidecan[["GWASpoly"]][["gwas_data_list"]]
#   )
#   
#   hidecan_data <- list(x, chrom_length)
#   
#   
#   ## Function to create a name for each dataset to use when choosing which
#   ## dataset should be plotted
#   make_names_hidecan_data <- function(hidecan_list){
#     
#     data_type_labels <- c("GWAS_data_thr" = "GWAS data",
#                           "DE_data_thr" = "DE data",
#                           "CAN_data_thr" = "Candidate genes list")
#     
#     labels <- sapply(hidecan_list, function(x){class(x)[[1]]})
#     
#     labels <- paste0(
#       data_type_labels[labels],
#       " (",
#       names(hidecan_list),
#       ")"
#     )
#     
#     labels <- sub(" ( )", "", labels, fixed = TRUE)
#     
#     labels
#   }
#   
#   track_choices <- as.list(make_names_hidecan_data(hidecan_data[[1]]))
#   names(track_choices) <- make_names_hidecan_data(hidecan_data[[1]])
#   
#   x <- hidecan_data[[1]]
#   #x <- x[match(input$tracks, make_names_hidecan_data(x))]
#   
#   #x <- lapply(x, function(y) y[which(y$chromosome %in% input$chrom),])
#   
#   chrom_length <- hidecan_data[[2]]
#   #chrom_length <- chrom_length[match(input$chrom, chrom_length$chromosome),]
#   
#   x <- x[-which(sapply(x, nrow) == 0)]
#   
#   str(x)
#   if(all(sapply(x, nrow) == 0)) stop("No QTL found")
#   print(chrom_length)
#   
#   p <- create_hidecan_plot(x,
#                             chrom_length,
#                            colour_genes_by_score = TRUE,
#                            remove_empty_chrom = TRUE,
#                            title = NULL,
#                            subtitle = NULL,
#                            n_rows = NULL,
#                            n_cols = 2,
#                            legend_position = "none",
#                            point_size = 3,
#                            label_size = 3.5,
#                            label_padding = 0.15)
#   p
#   
#   input$ncols <- 2
#   n.chr <- length(unique(p$data$chromosome))
#   ## Also use the number of tracks on the y axis
#   n.ytracks <- length(unique(p$data$dataset))
#   
#   size <- (n.ytracks * n.chr/input$ncols)*80
#   
#   size
#   
# })