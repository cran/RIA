# #' @title Calculates GLCM-based statistics
# #'
# #' @description Calculates GLCM-based statistics for given GLCM matrix.
# #'
# #' @param RIA_data_in \emph{RIA_image}.
# #'
# #' @param use_type string, can be \emph{"single"} which runs the function on a single image,
# #' which is determined using \emph{"use_orig"} or \emph{"use_slot"}. \emph{"glcm"}
# #' takes all datasets in the \emph{RIA_image$glcm} slot and runs the analysis on them.
# #'
# #' @param use_orig logical, indicating to use image present in \emph{RIA_data$orig}.
# #' If FALSE, the modified image will be used stored in \emph{RIA_data$modif}. However, GLCM matrices
# #' are usually not present in either slots, therefore giving the slot name using \emph{use_slot} is advised.
# #'
# #' @param use_slot string, name of slot where data wished to be used is. Use if the desired image
# #' is not in the \emph{data$orig} or \emph{data$modif} slot of the \emph{RIA_image}. For example,
# #' the desired dataset is in \emph{RIA_image$glcm$ep_4_111}, then \emph{use_slot} should be
# #' \emph{glcm$ep_4_111}. The results are automatically saved. If the results are not saved to
# #' the desired slot, then please use \emph{save_name} parameter. If the string contains "-" characters
# #' use "`" before the last slot name, for example: \emph{glcm$`ep_4_-1-1-1`}
# #'
# #' @param save_name string, indicating the name of subslot of \emph{$glcm} to save results to.
# #' If left empty, then it will be automatically determined.
# #'
# #' @param verbose_in logical, indicating whether to print detailed information.
# #' Most prints can also be suppressed using the \code{\link{suppressMessages}} function.
# #'
# #' @return \emph{RIA_image} containing the statistical information.
# #'
# #' @examples \dontrun{
# #' #Discretize loaded image and then calculate GLCM statistics
# #' RIA_image <- discretize(RIA_image, bins_in = 8, equal_prob = TRUE)
# #' RIA_image <- glcm(RIA_image, use_orig = FALSE, use_slot = "discretized$ep_8",
# #' off_right = 0, off_down = 1, off_z = 0)
# #' RIA_image <- glcm_stat(RIA_image, use_orig = FALSE, use_slot = "glcm$ep_8_010")
# #' 
# #' #Batch calculation of GLCM-based statistics on all calculated GLCMs
# #' RIA_image <- glcm_stat(RIA_image, use_type = "discretized")
# #' }
# #'
# #' @references Márton KOLOSSVÁRY et al.
# #' Radiomic Features Are Superior to Conventional Quantitative Computed Tomographic
# #' Metrics to Identify Coronary Plaques With Napkin-Ring Sign
# #' Circulation: Cardiovascular Imaging (2017).
# #' DOI: 10.1161/circimaging.117.006843
# #' \url{https://pubmed.ncbi.nlm.nih.gov/29233836/}
# #' 
# #' Márton KOLOSSVÁRY et al.
# #' Cardiac Computed Tomography Radiomics: A Comprehensive Review on Radiomic Techniques.
# #' Journal of Thoracic Imaging (2018).
# #' DOI: 10.1097/RTI.0000000000000268
# #' \url{https://pubmed.ncbi.nlm.nih.gov/28346329/}
# #' @encoding UTF-8



glcm_stat <- function(RIA_data_in, use_type = "single", use_orig = FALSE, use_slot = "glcm$es_8_111", save_name = NULL, verbose_in = TRUE)
{
  
  
  data_in_orig <- check_data_in(RIA_data_in, use_type = use_type, use_slot = use_slot, verbose_in = verbose_in)
  
  if(any(class(data_in_orig) != "list")) data_in_orig <- list(data_in_orig)
  list_names <- names(data_in_orig)
  if(!is.null(save_name) & (length(data_in_orig) != length(save_name))) {stop(paste0("PLEASE PROVIDE THE SAME NUMBER OF NAMES AS THERE ARE IMAGES!\n",
                                                                                     "NUMBER OF NAMES:  ", length(save_name), "\n",
                                                                                     "NUMBER OF IMAGES: ", length(data_in), "\n"))
  }
  
  for (i in 1: length(data_in_orig))
  {
    data_in <- data_in_orig[[i]]
    data_NA <- as.vector(data_in)
    data_NA <- data_NA[!is.na(data_NA)]
    if(length(data_NA) == 0) {message("WARNING: SUPPLIED RIA_image DOES NOT CONTAIN ANY DATA!!!")}
    
    
    Contrast            <- contrast(data_in)
    Contrast_s          <- contrast(data_in, type_in = "squared")
    Contrast_e          <- contrast(data_in, type_in = "entropy")
    
    Homogeneity2        <- homogeneity2(data_in)
    Homogeneity2_s      <- homogeneity2(data_in, type_in = "squared")
    Homogeneity2_e      <- homogeneity2(data_in, type_in = "entropy")
    
    Homogeneity2_nd        <- homogeneity2(data_in, diag = FALSE)
    Homogeneity2_s_nd      <- homogeneity2(data_in, type_in = "squared", diag = FALSE)
    Homogeneity2_e_nd      <- homogeneity2(data_in, type_in = "entropy", diag = FALSE)
    
    Dissimilarity       <- dissimilarity(data_in)
    Dissimilarity_s     <- dissimilarity(data_in, type_in = "squared")
    Dissimilarity_e     <- dissimilarity(data_in, type_in = "entropy")
    
    Homogeneity1        <- homogeneity1(data_in)
    Homogeneity1_s      <- homogeneity1(data_in, type_in = "squared")
    Homogeneity1_e      <- homogeneity1(data_in, type_in = "entropy")
    
    Homogeneity1_nd        <- homogeneity1(data_in, diag = FALSE)
    Homogeneity1_s_nd      <- homogeneity1(data_in, type_in = "squared", diag = FALSE)
    Homogeneity1_e_nd      <- homogeneity1(data_in, type_in = "entropy", diag = FALSE)
    
    DMN                 <- dmn(data_in)
    DMN_s               <- dmn(data_in, type_in = "squared")
    DMN_e               <- dmn(data_in, type_in = "entropy")
    
    IDMN                <- idmn(data_in)
    IDMN_s              <- idmn(data_in, type_in = "squared")
    IDMN_e              <- idmn(data_in, type_in = "entropy")
    
    IDMN_nd                 <- idmn(data_in, diag = FALSE)
    IDMN_s_nd               <- idmn(data_in, type_in = "squared", diag = FALSE)
    IDMN_e_nd               <- idmn(data_in, type_in = "entropy", diag = FALSE)
    
    DN                  <- dn(data_in)
    DN_s                <- dn(data_in, type_in = "squared")
    DN_e                <- dn(data_in, type_in = "entropy")
    
    IDN                 <- idn(data_in)
    IDN_s               <- idn(data_in, type_in = "squared")
    IDN_e               <- idn(data_in, type_in = "entropy")
    
    IDN_nd                  <- idn(data_in, diag = FALSE)
    IDN_s_nd                <- idn(data_in, type_in = "squared", diag = FALSE)
    IDN_e_nd                <- idn(data_in, type_in = "entropy", diag = FALSE)
    
    Autocorrelation     <- autocorrelation(data_in)
    Autocorrelation_s   <- autocorrelation(data_in, type_in = "squared")
    Autocorrelation_e   <- autocorrelation(data_in, type_in = "entropy")
    
    Autocorrelation_nd     <- autocorrelation(data_in, diag = FALSE)
    Autocorrelation_s_nd   <- autocorrelation(data_in, type_in = "squared", diag = FALSE)
    Autocorrelation_e_nd   <- autocorrelation(data_in, type_in = "entropy", diag = FALSE)
    
    Inv_autocorrelation   <- inv_autocorrelation(data_in)
    Inv_autocorrelation_s <- inv_autocorrelation(data_in, type_in = "squared")
    Inv_autocorrelation_e <- inv_autocorrelation(data_in, type_in = "entropy")
    
    Inv_autocorrelation_nd   <- inv_autocorrelation(data_in, diag = FALSE)
    Inv_autocorrelation_s_nd <- inv_autocorrelation(data_in, type_in = "squared", diag = FALSE)
    Inv_autocorrelation_e_nd <- inv_autocorrelation(data_in, type_in = "entropy", diag = FALSE)
    
    Gauss               <- gauss(data_in)
    Gauss_s             <- gauss(data_in, type_in = "squared")
    Gauss_e             <- gauss(data_in, type_in = "entropy")
    
    Gauss_nd               <- gauss(data_in, diag = FALSE)
    Gauss_s_nd             <- gauss(data_in, type_in = "squared", diag = FALSE)
    Gauss_e_nd             <- gauss(data_in, type_in = "entropy", diag = FALSE)
    
    Gauss_lp               <- gauss(data_in, loc = 1)
    Gauss_lp_s             <- gauss(data_in, type_in = "squared", loc = 1)
    Gauss_lp_e             <- gauss(data_in, type_in = "entropy", loc = 1)
    
    Gauss_lp_nd               <- gauss(data_in, diag = FALSE, loc = 1)
    Gauss_lp_s_nd             <- gauss(data_in, type_in = "squared", diag = FALSE, loc = 1)
    Gauss_lp_e_nd             <- gauss(data_in, type_in = "entropy", diag = FALSE, loc = 1)
    
    Gauss_lf               <- gauss(data_in, loc = 2)
    Gauss_lf_s             <- gauss(data_in, type_in = "squared", loc = 2)
    Gauss_lf_e             <- gauss(data_in, type_in = "entropy", loc = 2)
    
    Gauss_lf_nd               <- gauss(data_in, diag = FALSE, loc = 2)
    Gauss_lf_s_nd             <- gauss(data_in, type_in = "squared", diag = FALSE, loc = 2)
    Gauss_lf_e_nd             <- gauss(data_in, type_in = "entropy", diag = FALSE, loc = 2)
    
    Gauss_rf               <- gauss(data_in, loc = 4)
    Gauss_rf_s             <- gauss(data_in, type_in = "squared", loc = 4)
    Gauss_rf_e             <- gauss(data_in, type_in = "entropy", loc = 4)
    
    Gauss_rf_nd               <- gauss(data_in, diag = FALSE, loc = 4)
    Gauss_rf_s_nd             <- gauss(data_in, type_in = "squared", diag = FALSE, loc = 4)
    Gauss_rf_e_nd             <- gauss(data_in, type_in = "entropy", diag = FALSE, loc = 4)
    
    Gauss_rp               <- gauss(data_in, loc = 5)
    Gauss_rp_s             <- gauss(data_in, type_in = "squared", loc = 5)
    Gauss_rp_e             <- gauss(data_in, type_in = "entropy", loc = 5)
    
    Gauss_rp_nd               <- gauss(data_in, diag = FALSE, loc = 5)
    Gauss_rp_s_nd             <- gauss(data_in, type_in = "squared", diag = FALSE, loc = 5)
    Gauss_rp_e_nd             <- gauss(data_in, type_in = "entropy", diag = FALSE, loc = 5)
    
    Inv_Gauss           <- gauss(data_in, inverse = TRUE)
    Inv_Gauss_s         <- gauss(data_in, inverse = TRUE, type_in = "squared")
    Inv_Gauss_e         <- gauss(data_in, inverse = TRUE, type_in = "entropy")
    
    Inv_Gauss_nd           <- gauss(data_in, inverse = TRUE, diag = FALSE)
    Inv_Gauss_s_nd         <- gauss(data_in, inverse = TRUE, type_in = "squared", diag = FALSE)
    Inv_Gauss_e_nd         <- gauss(data_in, inverse = TRUE, type_in = "entropy", diag = FALSE)
    
    Inv_Gauss_lp           <- gauss(data_in, inverse = TRUE, loc = 1)
    Inv_Gauss_lp_s         <- gauss(data_in, inverse = TRUE, type_in = "squared", loc = 1)
    Inv_Gauss_lp_e         <- gauss(data_in, inverse = TRUE, type_in = "entropy", loc = 1)
    
    Inv_Gauss_lp_nd           <- gauss(data_in, inverse = TRUE, diag = FALSE, loc = 1)
    Inv_Gauss_lp_s_nd         <- gauss(data_in, inverse = TRUE, type_in = "squared", diag = FALSE, loc = 1)
    Inv_Gauss_lp_e_nd         <- gauss(data_in, inverse = TRUE, type_in = "entropy", diag = FALSE, loc = 1)
    
    Inv_Gauss_lf           <- gauss(data_in, inverse = TRUE, loc = 2)
    Inv_Gauss_lf_s         <- gauss(data_in, inverse = TRUE, type_in = "squared", loc = 2)
    Inv_Gauss_lf_e         <- gauss(data_in, inverse = TRUE, type_in = "entropy", loc = 2)
    
    Inv_Gauss_lf_nd           <- gauss(data_in, inverse = TRUE, diag = FALSE, loc = 2)
    Inv_Gauss_lf_s_nd         <- gauss(data_in, inverse = TRUE, type_in = "squared", diag = FALSE, loc = 2)
    Inv_Gauss_lf_e_nd         <- gauss(data_in, inverse = TRUE, type_in = "entropy", diag = FALSE, loc = 2)
    
    Inv_Gauss_rf           <- gauss(data_in, inverse = TRUE, loc = 4)
    Inv_Gauss_rf_s         <- gauss(data_in, inverse = TRUE, type_in = "squared", loc = 4)
    Inv_Gauss_rf_e         <- gauss(data_in, inverse = TRUE, type_in = "entropy", loc = 4)
    
    Inv_Gauss_rf_nd           <- gauss(data_in, inverse = TRUE, diag = FALSE, loc = 4)
    Inv_Gauss_rf_s_nd         <- gauss(data_in, inverse = TRUE, type_in = "squared", diag = FALSE, loc = 4)
    Inv_Gauss_rf_e_nd         <- gauss(data_in, inverse = TRUE, type_in = "entropy", diag = FALSE, loc = 4)
    
    Inv_Gauss_rp           <- gauss(data_in, inverse = TRUE, loc = 5)
    Inv_Gauss_rp_s         <- gauss(data_in, inverse = TRUE, type_in = "squared", loc = 5)
    Inv_Gauss_rp_e         <- gauss(data_in, inverse = TRUE, type_in = "entropy", loc = 5)
    
    Inv_Gauss_rp_nd           <- gauss(data_in, inverse = TRUE, diag = FALSE, loc = 5)
    Inv_Gauss_rp_s_nd         <- gauss(data_in, inverse = TRUE, type_in = "squared", diag = FALSE, loc = 5)
    Inv_Gauss_rp_e_nd         <- gauss(data_in, inverse = TRUE, type_in = "entropy", diag = FALSE, loc = 5)
    
    Gauss_2f            <- gauss2f(data_in)
    Gauss_2f_s          <- gauss2f(data_in, type_in = "squared")
    Gauss_2f_e          <- gauss2f(data_in, type_in = "entropy")
    
    Gauss_2f_nd            <- gauss2f(data_in, diag = FALSE)
    Gauss_2f_s_nd          <- gauss2f(data_in, type_in = "squared", diag = FALSE)
    Gauss_2f_e_nd          <- gauss2f(data_in, type_in = "entropy", diag = FALSE)
    
    Inv_Gauss_2f        <- gauss2f(data_in, inverse = TRUE)
    Inv_Gauss_2f_s      <- gauss2f(data_in, inverse = TRUE, type_in = "squared")
    Inv_Gauss_2f_e      <- gauss2f(data_in, inverse = TRUE, type_in = "entropy")
    
    Inv_Gauss_2f_nd        <- gauss2f(data_in, inverse = TRUE, diag = FALSE)
    Inv_Gauss_2f_s_nd      <- gauss2f(data_in, inverse = TRUE, type_in = "squared", diag = FALSE)
    Inv_Gauss_2f_e_nd      <- gauss2f(data_in, inverse = TRUE, type_in = "entropy", diag = FALSE)
    
    Gauss_2p            <- gauss2p(data_in)
    Gauss_2p_s          <- gauss2p(data_in, type_in = "squared")
    Gauss_2p_e          <- gauss2p(data_in, type_in = "entropy")
    
    Gauss_2p_nd            <- gauss2p(data_in, diag = FALSE)
    Gauss_2p_s_nd          <- gauss2p(data_in, type_in = "squared", diag = FALSE)
    Gauss_2p_e_nd          <- gauss2p(data_in, type_in = "entropy", diag = FALSE)
    
    Inv_Gauss_2p        <- gauss2p(data_in, inverse = TRUE)
    Inv_Gauss_2p_s      <- gauss2p(data_in, inverse = TRUE, type_in = "squared")
    Inv_Gauss_2p_e      <- gauss2p(data_in, inverse = TRUE, type_in = "entropy")
    
    Inv_Gauss_2p_nd        <- gauss2p(data_in, inverse = TRUE, diag = FALSE)
    Inv_Gauss_2p_s_nd      <- gauss2p(data_in, inverse = TRUE, type_in = "squared", diag = FALSE)
    Inv_Gauss_2p_e_nd      <- gauss2p(data_in, inverse = TRUE, type_in = "entropy", diag = FALSE)
    
    Cluster_p            <- cluster(data_in,  pow = 4, inverse = FALSE, type_in = "single", diag = TRUE)
    Cluster_p_s          <- cluster(data_in,  pow = 4, inverse = FALSE, type_in = "squared", diag = TRUE)
    Cluster_p_e          <- cluster(data_in,  pow = 4, inverse = FALSE, type_in = "entropy", diag = TRUE)
    
    Cluster_p_nd            <- cluster(data_in,  pow = 4, inverse = FALSE, type_in = "single", diag = FALSE)
    Cluster_p_s_nd          <- cluster(data_in,  pow = 4, inverse = FALSE, type_in = "squared", diag = FALSE)
    Cluster_p_e_nd          <- cluster(data_in,  pow = 4, inverse = FALSE, type_in = "entropy", diag = FALSE)
    
    Inv_Cluster_p        <- cluster(data_in,  pow = 4, inverse = TRUE, type_in = "single", diag = TRUE)
    Inv_Cluster_p_s      <- cluster(data_in,  pow = 4, inverse = TRUE, type_in = "squared", diag = TRUE)
    Inv_Cluster_p_e      <- cluster(data_in,  pow = 4, inverse = TRUE, type_in = "entropy", diag = TRUE)
    
    Inv_Cluster_p_nd        <- cluster(data_in,  pow = 4, inverse = TRUE, type_in = "single", diag = FALSE)
    Inv_Cluster_p_s_nd      <- cluster(data_in,  pow = 4, inverse = TRUE, type_in = "squared", diag = FALSE)
    Inv_Cluster_p_e_nd      <- cluster(data_in,  pow = 4, inverse = TRUE, type_in = "entropy", diag = FALSE)
    
    Cluster_s            <- cluster(data_in,  pow = 3, inverse = FALSE, type_in = "single", diag = TRUE)
    Cluster_s_s          <- cluster(data_in,  pow = 3, inverse = FALSE, type_in = "squared", diag = TRUE)
    Cluster_s_e          <- cluster(data_in,  pow = 3, inverse = FALSE, type_in = "entropy", diag = TRUE)
    
    Cluster_s_nd            <- cluster(data_in,  pow = 3, inverse = FALSE, type_in = "single", diag = FALSE)
    Cluster_s_s_nd          <- cluster(data_in,  pow = 3, inverse = FALSE, type_in = "squared", diag = FALSE)
    Cluster_s_e_nd          <- cluster(data_in,  pow = 3, inverse = FALSE, type_in = "entropy", diag = FALSE)
    
    Inv_Cluster_s        <- cluster(data_in,  pow = 3, inverse = TRUE, type_in = "single", diag = TRUE)
    Inv_Cluster_s_s      <- cluster(data_in,  pow = 3, inverse = TRUE, type_in = "squared", diag = TRUE)
    Inv_Cluster_s_e      <- cluster(data_in,  pow = 3, inverse = TRUE, type_in = "entropy", diag = TRUE)
    
    Inv_Cluster_s_nd        <- cluster(data_in,  pow = 3, inverse = TRUE, type_in = "single", diag = FALSE)
    Inv_Cluster_s_s_nd      <- cluster(data_in,  pow = 3, inverse = TRUE, type_in = "squared", diag = FALSE)
    Inv_Cluster_s_e_nd      <- cluster(data_in,  pow = 3, inverse = TRUE, type_in = "entropy", diag = FALSE)
    
    Cluster_t            <- cluster(data_in,  pow = 2, inverse = FALSE, type_in = "single", diag = TRUE)
    Cluster_t_s          <- cluster(data_in,  pow = 2, inverse = FALSE, type_in = "squared", diag = TRUE)
    Cluster_t_e          <- cluster(data_in,  pow = 2, inverse = FALSE, type_in = "entropy", diag = TRUE)
    
    Cluster_t_nd            <- cluster(data_in,  pow = 2, inverse = FALSE, type_in = "single", diag = FALSE)
    Cluster_t_s_nd          <- cluster(data_in,  pow = 2, inverse = FALSE, type_in = "squared", diag = FALSE)
    Cluster_t_e_nd          <- cluster(data_in,  pow = 2, inverse = FALSE, type_in = "entropy", diag = FALSE)
    
    Inv_Cluster_t        <- cluster(data_in,  pow = 2, inverse = TRUE, type_in = "single", diag = TRUE)
    Inv_Cluster_t_s      <- cluster(data_in,  pow = 2, inverse = TRUE, type_in = "squared", diag = TRUE)
    Inv_Cluster_t_e      <- cluster(data_in,  pow = 2, inverse = TRUE, type_in = "entropy", diag = TRUE)
    
    Inv_Cluster_t_nd        <- cluster(data_in,  pow = 2, inverse = TRUE, type_in = "single", diag = FALSE)
    Inv_Cluster_t_s_nd      <- cluster(data_in,  pow = 2, inverse = TRUE, type_in = "squared", diag = FALSE)
    Inv_Cluster_t_e_nd      <- cluster(data_in,  pow = 2, inverse = TRUE, type_in = "entropy", diag = FALSE)
    
    Cluster_d            <- cluster(data_in,  pow = 1, inverse = FALSE, type_in = "single", diag = TRUE)
    Cluster_d_s          <- cluster(data_in,  pow = 1, inverse = FALSE, type_in = "squared", diag = TRUE)
    Cluster_d_e          <- cluster(data_in,  pow = 1, inverse = FALSE, type_in = "entropy", diag = TRUE)
    
    Cluster_d_nd            <- cluster(data_in,  pow = 1, inverse = FALSE, type_in = "single", diag = FALSE)
    Cluster_d_s_nd          <- cluster(data_in,  pow = 1, inverse = FALSE, type_in = "squared", diag = FALSE)
    Cluster_d_e_nd          <- cluster(data_in,  pow = 1, inverse = FALSE, type_in = "entropy", diag = FALSE)
    
    Inv_Cluster_d        <- cluster(data_in,  pow = 1, inverse = TRUE, type_in = "single", diag = TRUE)
    Inv_Cluster_d_s      <- cluster(data_in,  pow = 1, inverse = TRUE, type_in = "squared", diag = TRUE)
    Inv_Cluster_d_e      <- cluster(data_in,  pow = 1, inverse = TRUE, type_in = "entropy", diag = TRUE)
    
    Inv_Cluster_d_nd        <- cluster(data_in,  pow = 1, inverse = TRUE, type_in = "single", diag = FALSE)
    Inv_Cluster_d_s_nd      <- cluster(data_in,  pow = 1, inverse = TRUE, type_in = "squared", diag = FALSE)
    Inv_Cluster_d_e_nd      <- cluster(data_in,  pow = 1, inverse = TRUE, type_in = "entropy", diag = FALSE)
    
    Average               <- avg(data_in, type_in = "single")
    Average_s             <- avg(data_in, type_in = "squared")
    Average_e             <- avg(data_in, type_in = "entropy")
    
    Variances     <- variance(data_in, type_in = "single")
    Variances_s   <- variance(data_in, type_in = "squared")
    Variances_e   <- variance(data_in, type_in = "entropy")
    
    Correlation           <- correlation(data_in, type_in = "single")
    Correlation_s         <- correlation(data_in, type_in = "squared")
    Correlation_e         <- correlation(data_in, type_in = "entropy")
    
    Sum_average  <- sum_f(data_in, type_in = "single")
    Sum_energy   <- sum_f(data_in, type_in = "squared")
    Sum_entropy  <- sum_f(data_in, type_in = "entropy")
    Sum_variance <- sum_f(data_in, type_in = "variance")
    
    Dif_average  <- dif_f(data_in, type_in = "single")
    Dif_energy   <- dif_f(data_in, type_in = "squared")
    Dif_entropy  <- dif_f(data_in, type_in = "entropy")
    Dif_variance <- dif_f(data_in, type_in = "variance")
    
    Inv_sum_average  <- sum_f(data_in, inverse = TRUE, type_in = "single")
    Inv_sum_energy   <- sum_f(data_in, inverse = TRUE, type_in = "squared")
    Inv_sum_variance <- sum_f(data_in, inverse = TRUE, type_in = "variance")
    
    Inv_dif_average  <- dif_f(data_in, inverse = TRUE, type_in = "single")
    Inv_dif_energy   <- dif_f(data_in, inverse = TRUE, type_in = "squared")
    Inv_dif_variance <- dif_f(data_in, inverse = TRUE, type_in = "variance")
    
    
    IMC1                <- imc1(data_in)
    IMC2                <- imc2(data_in)
    
    
    data <- as.vector(data_in)
    data <- data[!is.na(data)]
    
    Energy       <- energy(data)
    Entropy      <- entropy(data, 2)
    
    Mean         <- base::mean(data)
    Median       <- stats::median(data)
    Mode         <- mode(data)[1]
    Geo_mean     <- geo_mean(data)
    Geo_mean2    <- geo_mean2(data)
    Geo_mean3    <- geo_mean3(data)
    Har_mean     <- har_mean(data)
    Trim_mean_5  <- base::mean(data, trim = 0.025)
    Trim_mean_10 <- base::mean(data, trim = 0.05)
    Trim_mean_20 <- base::mean(data, trim = 0.1)
    IQ_mean      <- base::mean(data, trim = 0.25)
    Tri_mean     <- (as.numeric(stats::quantile(data, 0.25)) +2*as.numeric(stats::quantile(data, 0.50)) + as.numeric(stats::quantile(data, 0.25)))/4
    Mn_AD_mn     <- mn_AD_mn(data)
    Mn_AD_md     <- mn_AD_md(data)
    Md_AD_mn     <- md_AD_mn(data)
    Md_AD_md     <- md_AD_md(data)
    MAD          <- stats::mad(data)
    Max_AD_mn    <- max_AD_mn(data)
    Max_AD_md    <- max_AD_md(data)
    RMS          <- rms(data)
    Min          <- base::min(data)
    Max          <- base::max(data)
    Quartiles    <- stats::quantile(data, seq(0.25, 0.75, 0.50))
    IQR          <- stats::IQR(data)
    Low_notch    <- as.numeric(Quartiles[1])-1.5*IQR
    High_notch   <- as.numeric(Quartiles[1])+1.5*IQR
    Range        <- abs(abs(base::range(data)[2] - base::range(data)[1]))
    Deciles      <- stats::quantile(data, seq(0.1, 0.9, 0.1))
    
    Variance     <- ifelse(length(data)>1, stats::var(data), 0)
    SD           <- ifelse(length(data)>1, stats::sd(data), 0)
    Skew         <- ifelse(length(data)>1, skew(data), 0)
    Kurtosis     <- ifelse(length(data)>1, kurtosis(data), 0)
    
    Uniformity   <- uniformity(data)
    
    
    
    
    metrics <- list(
      Contrast            <- Contrast,
      Contrast_s          <- Contrast_s,
      Contrast_e          <- Contrast_e,
      Homogeneity2        <- Homogeneity2,
      Homogeneity2_s      <- Homogeneity2_s,
      Homogeneity2_e      <- Homogeneity2_e,
      Homogeneity2_nd        <- Homogeneity2_nd,
      Homogeneity2_s_nd      <- Homogeneity2_s_nd,
      Homogeneity2_e_nd      <- Homogeneity2_e_nd,
      Dissimilarity       <- Dissimilarity,
      Dissimilarity_s     <- Dissimilarity_s,
      Dissimilarity_e     <- Dissimilarity_e,
      Homogeneity1        <- Homogeneity1,
      Homogeneity1_s      <- Homogeneity1_s,
      Homogeneity1_e      <- Homogeneity1_e,
      Homogeneity1_nd        <- Homogeneity1_nd,
      Homogeneity1_s_nd      <- Homogeneity1_s_nd,
      Homogeneity1_e_nd      <- Homogeneity1_e_nd,
      DMN                 <- DMN,
      DMN_s               <- DMN_s,
      DMN_e               <- DMN_e,
      IDMN                <- IDMN,
      IDMN_s              <- IDMN_s,
      IDMN_e              <- IDMN_e,
      IDMN_nd                 <- IDMN_nd,
      IDMN_s_nd               <- IDMN_s_nd,
      IDMN_e_nd               <- IDMN_e_nd,
      DN                  <- DN,
      DN_s                <- DN_s,
      DN_e                <- DN_e,
      IDN                 <- IDN,
      IDN_s               <- IDN_s,
      IDN_e               <- IDN_e,
      IDN_nd                  <- IDN_nd,
      IDN_s_nd                <- IDN_s_nd,
      IDN_e_nd                <- IDN_e_nd,
      Autocorrelation     <- Autocorrelation,
      Autocorrelation_s   <- Autocorrelation_s,
      Autocorrelation_e   <- Autocorrelation_e,
      Autocorrelation_nd     <- Autocorrelation_nd,
      Autocorrelation_s_nd   <- Autocorrelation_s_nd,
      Autocorrelation_e_nd   <- Autocorrelation_e_nd,
      Inv_autocorrelation   <- Inv_autocorrelation,
      Inv_autocorrelation_s <- Inv_autocorrelation_s,
      Inv_autocorrelation_e <- Inv_autocorrelation_e,
      Inv_autocorrelation_nd   <- Inv_autocorrelation_nd,
      Inv_autocorrelation_s_nd <- Inv_autocorrelation_s_nd,
      Inv_autocorrelation_e_nd <- Inv_autocorrelation_e_nd,
      Gauss               <- Gauss,
      Gauss_s             <- Gauss_s,
      Gauss_e             <- Gauss_e,
      Gauss_nd               <- Gauss_nd,
      Gauss_s_nd             <- Gauss_s_nd,
      Gauss_e_nd             <- Gauss_e_nd,
      Gauss_lp               <- Gauss_lp,
      Gauss_lp_s             <- Gauss_lp_s,
      Gauss_lp_e             <- Gauss_lp_e,
      Gauss_lp_nd               <- Gauss_lp_nd,
      Gauss_lp_s_nd             <- Gauss_lp_s_nd,
      Gauss_lp_e_nd             <- Gauss_lp_e_nd,
      Gauss_lf               <- Gauss_lf,
      Gauss_lf_s             <- Gauss_lf_s,
      Gauss_lf_e             <- Gauss_lf_e,
      Gauss_lf_nd               <- Gauss_lf_nd,
      Gauss_lf_s_nd             <- Gauss_lf_s_nd,
      Gauss_lf_e_nd             <- Gauss_lf_e_nd,
      Gauss_rf               <- Gauss_rf,
      Gauss_rf_s             <- Gauss_rf_s,
      Gauss_rf_e             <- Gauss_rf_e,
      Gauss_rf_nd               <- Gauss_rf_nd,
      Gauss_rf_s_nd             <- Gauss_rf_s_nd,
      Gauss_rf_e_nd             <- Gauss_rf_e_nd,
      Gauss_rp               <- Gauss_rp,
      Gauss_rp_s             <- Gauss_rp_s,
      Gauss_rp_e             <- Gauss_rp_e,
      Gauss_rp_nd               <- Gauss_rp_nd,
      Gauss_rp_s_nd             <- Gauss_rp_s_nd,
      Gauss_rp_e_nd             <- Gauss_rp_e_nd,
      Inv_Gauss           <- Inv_Gauss,
      Inv_Gauss_s         <- Inv_Gauss_s,
      Inv_Gauss_e         <- Inv_Gauss_e,
      Inv_Gauss_nd           <- Inv_Gauss_nd,
      Inv_Gauss_s_nd         <- Inv_Gauss_s_nd,
      Inv_Gauss_e_nd         <- Inv_Gauss_e_nd,
      Inv_Gauss_lp           <- Inv_Gauss_lp,
      Inv_Gauss_lp_s         <- Inv_Gauss_lp_s,
      Inv_Gauss_lp_e         <- Inv_Gauss_lp_e,
      Inv_Gauss_lp_nd           <- Inv_Gauss_lp_nd,
      Inv_Gauss_lp_s_nd         <- Inv_Gauss_lp_s_nd,
      Inv_Gauss_lp_e_nd         <- Inv_Gauss_lp_e_nd,
      Inv_Gauss_lf           <- Inv_Gauss_lf,
      Inv_Gauss_lf_s         <- Inv_Gauss_lf_s,
      Inv_Gauss_lf_e         <- Inv_Gauss_lf_e,
      Inv_Gauss_lf_nd           <- Inv_Gauss_lf_nd,
      Inv_Gauss_lf_s_nd         <- Inv_Gauss_lf_s_nd,
      Inv_Gauss_lf_e_nd         <- Inv_Gauss_lf_e_nd,
      Inv_Gauss_rf           <- Inv_Gauss_rf,
      Inv_Gauss_rf_s         <- Inv_Gauss_rf_s,
      Inv_Gauss_rf_e         <- Inv_Gauss_rf_e,
      Inv_Gauss_rf_nd           <- Inv_Gauss_rf_nd,
      Inv_Gauss_rf_s_nd         <- Inv_Gauss_rf_s_nd,
      Inv_Gauss_rf_e_nd         <- Inv_Gauss_rf_e_nd,
      Inv_Gauss_rp           <- Inv_Gauss_rp,
      Inv_Gauss_rp_s         <- Inv_Gauss_rp_s,
      Inv_Gauss_rp_e         <- Inv_Gauss_rp_e,
      Inv_Gauss_rp_nd           <- Inv_Gauss_rp_nd,
      Inv_Gauss_rp_s_nd         <- Inv_Gauss_rp_s_nd,
      Inv_Gauss_rp_e_nd         <- Inv_Gauss_rp_e_nd,
      Gauss_2f            <- Gauss_2f,
      Gauss_2f_s          <- Gauss_2f_s,
      Gauss_2f_e          <- Gauss_2f_e,
      Gauss_2f_nd            <- Gauss_2f_nd,
      Gauss_2f_s_nd          <- Gauss_2f_s_nd,
      Gauss_2f_e_nd          <- Gauss_2f_e_nd,
      Inv_Gauss_2f        <- Inv_Gauss_2f,
      Inv_Gauss_2f_s      <- Inv_Gauss_2f_s,
      Inv_Gauss_2f_e      <- Inv_Gauss_2f_e,
      Inv_Gauss_2f_nd        <- Inv_Gauss_2f_nd,
      Inv_Gauss_2f_s_nd      <- Inv_Gauss_2f_s_nd,
      Inv_Gauss_2f_e_nd      <- Inv_Gauss_2f_e_nd,
      Gauss_2p            <- Gauss_2p,
      Gauss_2p_s          <- Gauss_2p_s,
      Gauss_2p_e          <- Gauss_2p_e,
      Gauss_2p_nd            <- Gauss_2p_nd,
      Gauss_2p_s_nd          <- Gauss_2p_s_nd,
      Gauss_2p_e_nd          <- Gauss_2p_e_nd,
      Inv_Gauss_2p        <- Inv_Gauss_2p,
      Inv_Gauss_2p_s      <- Inv_Gauss_2p_s,
      Inv_Gauss_2p_e      <- Inv_Gauss_2p_e,
      Inv_Gauss_2p_nd        <- Inv_Gauss_2p_nd,
      Inv_Gauss_2p_s_nd      <- Inv_Gauss_2p_s_nd,
      Inv_Gauss_2p_e_nd      <- Inv_Gauss_2p_e_nd,
      Cluster_p            <- Cluster_p,
      Cluster_p_s          <- Cluster_p_s,
      Cluster_p_e          <- Cluster_p_e,
      Cluster_p_nd            <- Cluster_p_nd,
      Cluster_p_s_nd          <- Cluster_p_s_nd,
      Cluster_p_e_nd          <- Cluster_p_e_nd,
      Inv_Cluster_p        <- Inv_Cluster_p,
      Inv_Cluster_p_s      <- Inv_Cluster_p_s,
      Inv_Cluster_p_e      <- Inv_Cluster_p_e,
      Inv_Cluster_p_nd        <- Inv_Cluster_p_nd,
      Inv_Cluster_p_s_nd      <- Inv_Cluster_p_s_nd,
      Inv_Cluster_p_e_nd      <- Inv_Cluster_p_e_nd,
      Cluster_s            <- Cluster_s,
      Cluster_s_s          <- Cluster_s_s,
      Cluster_s_e          <- Cluster_s_e,
      Cluster_s_nd            <- Cluster_s_nd,
      Cluster_s_s_nd          <- Cluster_s_s_nd,
      Cluster_s_e_nd          <- Cluster_s_e_nd,
      Inv_Cluster_s        <- Inv_Cluster_s,
      Inv_Cluster_s_s      <- Inv_Cluster_s_s,
      Inv_Cluster_s_e      <- Inv_Cluster_s_e,
      Inv_Cluster_s_nd        <- Inv_Cluster_s_nd,
      Inv_Cluster_s_s_nd      <- Inv_Cluster_s_s_nd,
      Inv_Cluster_s_e_nd      <- Inv_Cluster_s_e_nd,
      Cluster_t            <- Cluster_t,
      Cluster_t_s          <- Cluster_t_s,
      Cluster_t_e          <- Cluster_t_e,
      Cluster_t_nd            <- Cluster_t_nd,
      Cluster_t_s_nd          <- Cluster_t_s_nd,
      Cluster_t_e_nd          <- Cluster_t_e_nd,
      Inv_Cluster_t        <- Inv_Cluster_t,
      Inv_Cluster_t_s      <- Inv_Cluster_t_s,
      Inv_Cluster_t_e      <- Inv_Cluster_t_e,
      Inv_Cluster_t_nd        <- Inv_Cluster_t_nd,
      Inv_Cluster_t_s_nd      <- Inv_Cluster_t_s_nd,
      Inv_Cluster_t_e_nd      <- Inv_Cluster_t_e_nd,
      Cluster_d            <- Cluster_d,
      Cluster_d_s          <- Cluster_d_s,
      Cluster_d_e          <- Cluster_d_e,
      Cluster_d_nd            <- Cluster_d_nd,
      Cluster_d_s_nd          <- Cluster_d_s_nd,
      Cluster_d_e_nd          <- Cluster_d_e_nd,
      Inv_Cluster_d        <- Inv_Cluster_d,
      Inv_Cluster_d_s      <- Inv_Cluster_d_s,
      Inv_Cluster_d_e      <- Inv_Cluster_d_e,
      Inv_Cluster_d_nd        <- Inv_Cluster_d_nd,
      Inv_Cluster_d_s_nd      <- Inv_Cluster_d_s_nd,
      Inv_Cluster_d_e_nd      <- Inv_Cluster_d_e_nd,
      Average               <- Average,
      Average_s             <- Average_s,
      Average_e             <- Average_e,
      Variances     <- Variances,
      Variances_s   <- Variances_s,
      Variances_e   <- Variances_e,
      Correlation           <- Correlation,
      Correlation_s         <- Correlation_s,
      Correlation_e         <- Correlation_e,
      Sum_average  <- Sum_average,
      Sum_energy   <- Sum_energy,
      Sum_entropy  <- Sum_entropy,
      Sum_variance <- Sum_variance,
      Dif_average  <- Dif_average,
      Dif_energy   <- Dif_energy,
      Dif_entropy  <- Dif_entropy,
      Dif_variance <- Dif_variance,
      Inv_sum_average  <- Inv_sum_average,
      Inv_sum_energy   <- Inv_sum_energy,
      Inv_sum_variance <- Inv_sum_variance,
      Inv_dif_average  <- Inv_dif_average,
      Inv_dif_energy   <- Inv_dif_energy,
      Inv_dif_variance <- Inv_dif_variance,
      IMC1                <- IMC1,
      IMC2                <- IMC2,
      Energy       <- Energy,
      Entropy      <- Entropy,
      
      Mean         <- Mean,
      Median       <- Median,
      Mode         <- Mode,
      Geo_mean     <- Geo_mean,
      Geo_mean2    <- Geo_mean2,
      Geo_mean3    <- Geo_mean3,
      Har_mean     <- Har_mean,
      Trim_mean_5  <- Trim_mean_5,
      Trim_mean_10 <- Trim_mean_10,
      Trim_mean_20 <- Trim_mean_20,
      IQ_mean      <- IQ_mean,
      Tri_mean     <- Tri_mean,
      Mn_AD_mn     <- Mn_AD_mn,
      Mn_AD_md     <- Mn_AD_md,
      Md_AD_mn     <- Md_AD_mn,
      Md_AD_md     <- Md_AD_md,
      MAD          <- MAD,
      Max_AD_mn    <- Max_AD_mn,
      Max_AD_md    <- Max_AD_md,
      RMS          <- RMS,
      Min          <- Min,
      Max          <- Max,
      Quartiles    <- Quartiles,
      IQR          <- IQR,
      Low_notch    <- Low_notch,
      High_notch   <- High_notch,
      Range        <- Range,
      Deciles      <- Deciles,
      Variance     <- Variance,
      SD           <- SD,
      Skew         <- Skew,
      Kurtosis     <- Kurtosis,
      Uniformity   <- Uniformity
    )
    
    
    stat_names <- c(    "Contrast",
                        "Contrast_s",
                        "Contrast_e",
                        "Homogeneity2",
                        "Homogeneity2_s",
                        "Homogeneity2_e",
                        "Homogeneity2_nd",
                        "Homogeneity2_s_nd",
                        "Homogeneity2_e_nd",
                        "Dissimilarity",
                        "Dissimilarity_s",
                        "Dissimilarity_e",
                        "Homogeneity1",
                        "Homogeneity1_s",
                        "Homogeneity1_e",
                        "Homogeneity1_nd",
                        "Homogeneity1_s_nd",
                        "Homogeneity1_e_nd",
                        "DMN",
                        "DMN_s",
                        "DMN_e",
                        "IDMN",
                        "IDMN_s",
                        "IDMN_e",
                        "IDMN_nd",
                        "IDMN_s_nd",
                        "IDMN_e_nd",
                        "DN",
                        "DN_s",
                        "DN_e",
                        "IDN",
                        "IDN_s",
                        "IDN_e",
                        "IDN_nd",
                        "IDN_s_nd",
                        "IDN_e_nd",
                        "Autocorrelation",
                        "Autocorrelation_s",
                        "Autocorrelation_e",
                        "Autocorrelation_nd",
                        "Autocorrelation_s_nd",
                        "Autocorrelation_e_nd",
                        "Inv_autocorrelation",
                        "Inv_autocorrelation_s",
                        "Inv_autocorrelation_e",
                        "Inv_autocorrelation_nd",
                        "Inv_autocorrelation_s_nd",
                        "Inv_autocorrelation_e_nd",
                        "Gauss",
                        "Gauss_s",
                        "Gauss_e",
                        "Gauss_nd",
                        "Gauss_s_nd",
                        "Gauss_e_nd",
                        "Gauss_lp",
                        "Gauss_lp_s",
                        "Gauss_lp_e",
                        "Gauss_lp_nd",
                        "Gauss_lp_s_nd",
                        "Gauss_lp_e_nd",
                        "Gauss_lf",
                        "Gauss_lf_s",
                        "Gauss_lf_e",
                        "Gauss_lf_nd",
                        "Gauss_lf_s_nd",
                        "Gauss_lf_e_nd",
                        "Gauss_rf",
                        "Gauss_rf_s",
                        "Gauss_rf_e",
                        "Gauss_rf_nd",
                        "Gauss_rf_s_nd",
                        "Gauss_rf_e_nd",
                        "Gauss_rp",
                        "Gauss_rp_s",
                        "Gauss_rp_e",
                        "Gauss_rp_nd",
                        "Gauss_rp_s_nd",
                        "Gauss_rp_e_nd",
                        "Inv_Gauss",
                        "Inv_Gauss_s",
                        "Inv_Gauss_e",
                        "Inv_Gauss_nd",
                        "Inv_Gauss_s_nd",
                        "Inv_Gauss_e_nd",
                        "Inv_Gauss_lp",
                        "Inv_Gauss_lp_s",
                        "Inv_Gauss_lp_e",
                        "Inv_Gauss_lp_nd",
                        "Inv_Gauss_lp_s_nd",
                        "Inv_Gauss_lp_e_nd",
                        "Inv_Gauss_lf",
                        "Inv_Gauss_lf_s",
                        "Inv_Gauss_lf_e",
                        "Inv_Gauss_lf_nd",
                        "Inv_Gauss_lf_s_nd",
                        "Inv_Gauss_lf_e_nd",
                        "Inv_Gauss_rf",
                        "Inv_Gauss_rf_s",
                        "Inv_Gauss_rf_e",
                        "Inv_Gauss_rf_nd",
                        "Inv_Gauss_rf_s_nd",
                        "Inv_Gauss_rf_e_nd",
                        "Inv_Gauss_rp",
                        "Inv_Gauss_rp_s",
                        "Inv_Gauss_rp_e",
                        "Inv_Gauss_rp_nd",
                        "Inv_Gauss_rp_s_nd",
                        "Inv_Gauss_rp_e_nd",
                        "Gauss_2f",
                        "Gauss_2f_s",
                        "Gauss_2f_e",
                        "Gauss_2f_nd",
                        "Gauss_2f_s_nd",
                        "Gauss_2f_e_nd",
                        "Inv_Gauss_2f",
                        "Inv_Gauss_2f_s",
                        "Inv_Gauss_2f_e",
                        "Inv_Gauss_2f_nd",
                        "Inv_Gauss_2f_s_nd",
                        "Inv_Gauss_2f_e_nd",
                        "Gauss_2p",
                        "Gauss_2p_s",
                        "Gauss_2p_e",
                        "Gauss_2p_nd",
                        "Gauss_2p_s_nd",
                        "Gauss_2p_e_nd",
                        "Inv_Gauss_2p",
                        "Inv_Gauss_2p_s",
                        "Inv_Gauss_2p_e",
                        "Inv_Gauss_2p_nd",
                        "Inv_Gauss_2p_s_nd",
                        "Inv_Gauss_2p_e_nd",
                        "Cluster_p",
                        "Cluster_p_s",
                        "Cluster_p_e",
                        "Cluster_p_nd",
                        "Cluster_p_s_nd",
                        "Cluster_p_e_nd",
                        "Inv_Cluster_p",
                        "Inv_Cluster_p_s",
                        "Inv_Cluster_p_e",
                        "Inv_Cluster_p_nd",
                        "Inv_Cluster_p_s_nd",
                        "Inv_Cluster_p_e_nd",
                        "Cluster_s",
                        "Cluster_s_s",
                        "Cluster_s_e",
                        "Cluster_s_nd",
                        "Cluster_s_s_nd",
                        "Cluster_s_e_nd",
                        "Inv_Cluster_s",
                        "Inv_Cluster_s_s",
                        "Inv_Cluster_s_e",
                        "Inv_Cluster_s_nd",
                        "Inv_Cluster_s_s_nd",
                        "Inv_Cluster_s_e_nd",
                        "Cluster_t",
                        "Cluster_t_s",
                        "Cluster_t_e",
                        "Cluster_t_nd",
                        "Cluster_t_s_nd",
                        "Cluster_t_e_nd",
                        "Inv_Cluster_t",
                        "Inv_Cluster_t_s",
                        "Inv_Cluster_t_e",
                        "Inv_Cluster_t_nd",
                        "Inv_Cluster_t_s_nd",
                        "Inv_Cluster_t_e_nd",
                        "Cluster_d",
                        "Cluster_d_s",
                        "Cluster_d_e",
                        "Cluster_d_nd",
                        "Cluster_d_s_nd",
                        "Cluster_d_e_nd",
                        "Inv_Cluster_d",
                        "Inv_Cluster_d_s",
                        "Inv_Cluster_d_e",
                        "Inv_Cluster_d_nd",
                        "Inv_Cluster_d_s_nd",
                        "Inv_Cluster_d_e_nd",
                        "Average",
                        "Average_s",
                        "Average_e",
                        "Variances",
                        "Variances_s",
                        "Variances_e",
                        "Correlation",
                        "Correlation_s",
                        "Correlation_e",
                        "Sum_average",
                        "Sum_energy",
                        "Sum_entropy",
                        "Sum_variance",
                        "Dif_average",
                        "Dif_energy",
                        "Dif_entropy",
                        "Dif_variance",
                        "Inv_sum_average",
                        "Inv_sum_energy",
                        "Inv_sum_variance",
                        "Inv_dif_average",
                        "Inv_dif_energy",
                        "Inv_dif_variance",
                        "IMC1",
                        "IMC2",
                        "Energy",
                        "Entropy",
                        
                        "Mean",
                        "Median",
                        "Mode",
                        "Geo_mean",
                        "Geo_mean2",
                        "Geo_mean3",
                        "Har_mean",
                        "Trim_mean_5",
                        "Trim_mean_10",
                        "Trim_mean_20",
                        "IQ_mean",
                        "Tri_mean",
                        "Mn_AD_mn",
                        "Mn_AD_md",
                        "Md_AD_mn",
                        "Md_AD_md",
                        "MAD",
                        "Max_AD_mn",
                        "Max_AD_md",
                        "RMS",
                        "Min",
                        "Max",
                        "Quartiles",
                        "IQR",
                        "Low_notch",
                        "High_notch",
                        "Range",
                        "Deciles",
                        "Variance",
                        "SD",
                        "Skew",
                        "Kurtosis",
                        "Uniformity")
    
    
    names(metrics) <- stat_names
    
    
    
    if(use_type == "single") {
      if(any(class(RIA_data_in) == "RIA_image"))
      {
        if(is.null(save_name)) {
          txt <- automatic_name(RIA_data_in, use_orig, use_slot)
          RIA_data_in$stat_glcm[[txt]] <- metrics
          
        }
        if(!is.null(save_name)) {RIA_data_in$stat_glcm[[save_name]] <- metrics
        }
      }
    }
    
    if(use_type == "glcm") {
      if(any(class(RIA_data_in) == "RIA_image"))
      {
        if(is.null(save_name[i])) {
          txt <- list_names[i]
          RIA_data_in$stat_glcm[[txt]] <- metrics
        }
        if(!is.null(save_name[i])) {RIA_data_in$stat_glcm[[save_name[i]]] <- metrics
        }
      }
    }
    
    if(is.null(save_name)) {txt_name <- txt
    } else {txt_name <- save_name}
    if(verbose_in) {message(paste0("GLCM STATISTICS WAS SUCCESSFULLY ADDED TO '", txt_name, "' SLOT OF RIA_image$stat_glcm\n"))}
    
    
  }
  
  if(any(class(RIA_data_in) == "RIA_image") ) return(RIA_data_in)
  else return(metrics)
}

