# Connecting helpers functions

#' @include helpers.R


#' @export get_batch_loci
get_batch_loci <- function(loci, proj_path, group_list, genome) {
  
  # proj_path <- '/media/HEAP-EPI/Claes-Jensen/'
  
  # loci_gs <- loci_gs[loci_gs$plot==1,] # Только те у кого есть pathologic blood
  
  llply(
    1: length(loci), 
    function(locus_i, loci, proj_path) {
      # locus_i <- 1
      locus <- loci[locus_i]
      
      # chr <-  loci_df$chr[loci_df$Locus==locus]
      chr <- as.character(seqnames(locus))
      start_g <- as.numeric(start(locus))
      end_g <- as.numeric(end(locus))

      if (locus$bs_type == 'OT') {
        
        refG <- Biostrings::getSeq(genome, chr, start_g, end_g)
        refBS <- Biostrings::chartr('C','T', refG)
        refBS <- as.character(refBS)
        refG <- as.character(refG)
        G_position <- paste0(chr,':',seq(start_g, end_g))
        coef <- c('A'=1,'G'=1,'T'=1)
        
      }
      else if (locus$bs_type == 'CTOT') {
        
        refG <- Biostrings::getSeq(genome, chr, start_g, end_g)
        refBS <- Biostrings::chartr('C','T', refG)
        refBS <- as.character(Biostrings::reverseComplement(refBS))
        refG <- as.character(Biostrings::reverseComplement(refG))
        G_position <- paste0(chr,':',seq(end_g, start_g))
        coef <- c('A'=1,'C'=1,'T'=1)
        
      }
      else {
        stop(sprintf("Unexpected bs_type='%s', only 'OT' or 'CTOT' allowed.", locus$bs_type))
      }
      
      locus_name <- names(locus)
      locus_path <- file.path(proj_path, locus_name)
      
      get_locus(locus, locus_path, group_list, refG, refBS, coef, G_position)  

    }, loci, proj_path)
}


get_locus <- function(locus, locus_path, group_list, refG, refBS, coef, G_position) {
  
  #' @description Main function to read ab1 files, calculates corrected amplification values for nucleotides and methylation values. 
  #' @param locus GRanges The genomic locus of sequence, "TUBGstart" for example
  #' @param path String path to files with three folders: 'nometh', 'group1', 'group2'
  #' @param refG String Genomic sequence 
  #' @param refBS String Bisulfite treated genomic sequence 
  #' @param coef Numeric Nucleotide coefficients for optimization purposes, For "OT" it should be c('A'=1,'G'=1,'T'=1), whereas for CTOT  is c('A'=1,'C'=1,'T'=1). Feel Free to experiment with numbers for coefficients.
  #' @param G_position String For OT its genomic positions from start to end, for CTOT its positions from end to start. Example: paste0(chr,':',seq(end, start)) for CTOT.
  
  locus_name <- names(locus)
  chr_g <- as.character(seqnames(locus))
  start_g <- as.numeric(start(locus))
  end_g <- as.numeric(end(locus))
  names(group_list) <- group_list
  
  # пути до файлов nobis и blood
  # llply_with_names <- function(l, fun, ...){
  #   res <- llply(l, fun, ...)
  #   names(res) <- names(l)
  # }
  
  groups_path <- llply(
    group_list,
    function(group_name){
      with(
        list(file_list = list.files(file.path(locus_path, group_name), full.names = T, pattern = '\\.ab1$')),
        file_list[order(as.numeric(str_match(file_list, '(\\d+)-')[,2]))]  
      )
    })

  groups_seq <- llply(
    group_list,
    function(group_name, groups_path, locus_name){
      cat(sprintf("Getting abif files for group '%s' of locus '%s'.\n", group_name, locus_name))
      get_seq_list(groups_path[[group_name]], group_name)
    }, groups_path, locus_name)
  

  cat(paste0('Getting amplification positions for nobis to locus: ', locus_name, '\n'))
  
  groups_exp_peaks <- llply(
    group_list,
    function(group_name, refBS, locus){
      cat(sprintf("Getting amplification positions for '%s' of locus '%s'.\n", group_name, locus_name))
      exp_peaks_am(groups_seq[[group_name]], refBS)
    }, refBS, locus)
  
  # exp_peaks_nobis <- exp_peaks_am(seq_list_nobis, refBS)
  
  cat(sprintf("Collecting SDs, Supports and Direction for 'nometh' to locus '%s'.", locus_name, '\n'))
  
  nometh_SD <- apply(groups_exp_peaks[[1]], 1, sd, na.rm = TRUE)
  
  # TODO: Mask peaks with SD greater than 1.5 SD of 'nometh' SD.
  
  nometh_SD_outliers <- c(1)
  while( sum(nometh_SD_outliers, na.rm = TRUE) > 0){
    nometh_SD_sd <- sd(nometh_SD, na.rm = TRUE)
    nometh_SD_outliers <- nometh_SD > (5 * nometh_SD_sd)
    nometh_SD <- ifelse(nometh_SD_outliers, NA, nometh_SD)
    cat(sprintf("sd of 'nometh' SD %2f\n", nometh_SD_sd))
  }

  cpg_mask <- mask_cpg(refG, refBS)
  cpg <- ifelse(cpg_mask==1, 0, 1)

  # чтобы CG не участовали в оптимизации
  nometh_Sup <- pmin(!is.na(nometh_SD), cpg_mask)
  
  nometh_Dir <- ifelse(nometh_Sup == 1, apply(groups_exp_peaks[[1]], 1, median, na.rm=TRUE), NA)
  
  groups_corrected <- llply(
    group_list,
    function(group_name, groups_exp_peaks, nometh_Dir, refBS, coef){
      optCorr(peaks = groups_exp_peaks[[group_name]], nometh_Dir = nometh_Dir, refBS = refBS, coef = coef)
    }, groups_exp_peaks, nometh_Dir, refBS, coef)
  

  # colnames(nobis_corrected) <- paste('nobisD', seq(1, ncol(exp_peaks_nobis)), sep='')
  # colnames(blood_corrected) <- colnames(exp_peaks_blood)
  # colnames(bloodP_corrected) <- colnames(exp_peaks_bloodP)
  
  nometh_Dscale <- 100/apply(groups_corrected[[1]], 1, median, na.rm = TRUE)
  
  groups_meth <- llply(
    group_list,
    function(group_name, groups_corrected, nometh_Dscale){
      # group_name <- group_list[[2]]
      100 - (groups_corrected[[group_name]] * t(nometh_Dscale))
    }, groups_corrected, nometh_Dscale)
  
  # meth_blood_corr <- 100 - (blood_corrected * t(nometh_Dscale))
  # meth_bloodP_corr <- 100 - (bloodP_corrected * t(nometh_Dscale))

  nometh_DirD <- apply(groups_corrected[[1]], 1, median, na.rm=TRUE)

  cat(paste0('Getting corrected amplitudes for bloodN and nobis to locus: ', locus_name, '\n'))
  
  groups_SD_corr <- llply(
    group_list,
    function(group_name, groups_corrected){
      apply(groups_corrected[[group_name]], 1, sd, na.rm=TRUE)
    }, groups_corrected)
  
  # nobisSD_corr <- apply(nobis_corrected, 1, sd, na.rm=TRUE)
  # bloodN_SD_corr <- apply(blood_corrected, 1, sd, na.rm=TRUE)
  
  seq_df <- data.frame(
    'Position' = seq(1, length(unlist(strsplit(refG,'')))), 
    'refG' = unlist(strsplit(refG,'')),
    'refBS' = unlist(strsplit(refBS, '')),
    'G_position' = G_position,
    'CpG' = cpg)
  
  output <- list(
    'Locus' = locus_name, 'refG' = refG, 'refBS'=refBS, 
    'groups_seq' = groups_seq,
    # 'nobis_seq_list' = seq_list_nobis, 'blood_seq_list' = seq_list_blood, 
    'groups_peaks'= groups_exp_peaks,
    # 'peaks_nobis' = exp_peaks_nobis, 'peaks_blood' = exp_peaks_blood, 
    'nometh_SD' = nometh_SD, 'nometh_Sup' = nometh_Sup, 'nometh_Dir' = nometh_Dir, 
    'groups_corrected' = groups_corrected,
    # 'peaks_nobis_corr' = nobis_corrected, 'peaks_blood_corrected' = blood_corrected, 
    'groups_SD_corr' = groups_SD_corr,
    # 'nobisSD_corr' = nobisSD_corr, 'bloodN_SD_corr' = bloodN_SD_corr, 
    'nometh_DirD' = nometh_DirD, 'nometh_Dscale' = nometh_Dscale, 
    'seq_df' = seq_df, 
    'groups_meth' = groups_meth
    # 'methBlood' = meth_blood_corr, 'methBloodP' = meth_bloodP_corr
    )
  
  
  class(output) <- 'Locus'
  
  cat(paste0('DONE for locus: ', locus_name, '\n'))
  
  
  #table_google <- load_sheet(ss_id = ss_id, sheet = locus)
  #table_google <- as.data.frame(table_google)
  #names_locus <- colnames(table_google)
  

  return(output)
}


#Function to write Locus variable to google sheets

write_locus_to_sheet <- function(ss_id, locus_obj) {
  
  # locus_obj <- batch_locus[[1]]
  
  
  locus <- locus_obj$Locus
  raw <- range_read_cells(ss=ss_id, sheet=locus)
  cell_values <- c()
  for (cells in 1:nrow(raw )) {
    value <- raw$cell[[cells]]$formattedValue
    cell_values <- rbind(cell_values, value)
  }
  raw$name_of_cell <- cell_values
  
  # refG and refBS
  
  
  cat(paste0('Writing refG and refBS to locus: ', locus, '\n'))
  
  l_ply(colnames(locus_obj$seq_df), function(name, raw, seq, ss_id, locus) {
    
    loc <- raw$loc[grep(paste0('^',name,'$'), raw$name_of_cell)]
    start_data <- paste0(unlist(strsplit(loc,''))[1],4)
    range_write(ss=ss_id, seq[,name, drop=F], sheet = locus, range=start_data, col_names = F, reformat = F)
    
    
  }, raw, locus_obj$seq_df, ss_id, locus)
  
  
  
  cat(paste0('Writing nobis data before correction to locus: ', locus,'\n'))
  
  #nobis до коррекции
  
  l_ply(colnames(locus_obj$peaks_nobis), function(name, raw, seq, ss_id, locus) {
    
    # name <- colnames(exp_am_data)[20]
    # raw <-  CACend_raw
    # seq <- exp_am_data
    loc <- raw$loc[grep(paste0('^',name,'$'), raw$name_of_cell)]
    start_data <- paste0(unlist(strsplit(loc,'\\d+'))[1], 4)
    range_write(ss=ss_id, seq[,name, drop=F], sheet = locus, range=start_data, col_names = F, reformat = F)
    
    
  }, raw, locus_obj$peaks_nobis, ss_id, locus)
  
  Sys.sleep(10)
  
  
  #Статистики для nobis
  
  cat(paste0('Writing SD, SD after correction, supported nucleotides and Direction to locus: ', locus, '\n'))
  
  stat_df <- data.frame('nobisSD' = locus_obj$nobisSD, 'nobisSD_corr' = locus_obj$nobisSD_corr, 'nobisSup' = locus_obj$nobisSup, 'nobisDir' = locus_obj$nobisDir)
  
  l_ply(colnames(stat_df), function(name, raw, seq, ss_id, locus) {
    
    # name <- colnames(exp_am_data)[20]
    # raw <-  CACend_raw
    # seq <- exp_am_data
    loc <- raw$loc[grep(paste0('^',name,'$'), raw$name_of_cell)]
    start_data <- paste0(unlist(strsplit(loc,'\\d+'))[1], 4)
    range_write(ss=ss_id, seq[,name, drop=F], sheet = locus, range=start_data, col_names = F, reformat = F)
    
    
  }, raw, stat_df, ss_id, locus)
  
  
  Sys.sleep(10)
  
  cat(paste0('Writing nobis after correction to locus: ', locus, '\n'))
  
  #nobis после коррекции
  
  l_ply(colnames(locus_obj$peaks_nobis_corr), function(name, raw, seq, ss_id, locus) {
    
    # name <- colnames(exp_am_data)[20]
    # raw <-  CACend_raw
    # seq <- exp_am_data
    loc <- raw$loc[grep(paste0('^',name,'$'), raw$name_of_cell)]
    start_data <- paste0(unlist(strsplit(loc,'\\d+'))[1], 4)
    range_write(ss=ss_id, seq[,name, drop=F], sheet = locus, range=start_data, col_names = F, reformat = F)
    
    
  }, raw, locus_obj$peaks_nobis_corr, ss_id, locus)
  
  
  Sys.sleep(10)
  
  #blood после коррекции
  
  cat(paste0('Writing bloodN after correction to locus: ', locus, '\n'))
  
  l_ply(colnames(locus_obj$peaks_blood_corrected), function(name, raw, seq, ss_id, locus) {
    
    # name <- colnames(exp_am_data)[20]
    # raw <-  CACend_raw
    # seq <- exp_am_data
    loc <- raw$loc[grep(paste0('^',name,'$'), raw$name_of_cell)]
    start_data <- paste0(unlist(strsplit(loc,'\\d+'))[1], 4)
    range_write(ss=ss_id, seq[,name, drop=F], sheet = locus, range=start_data, col_names = F, reformat = F)
    
  }, raw, locus_obj$peaks_blood_corrected, ss_id, locus)
  
  
  Sys.sleep(10)
  
  df_meth_blood <- cbind(locus_obj$seq_df, locus_obj$methBlood)
  meth_blood <- df_meth_blood[,grep('bloodN', colnames(df_meth_blood))]
  df_mean_sd_blood <- data.frame('bloodNMean' = apply(meth_blood,1,mean, na.rm=TRUE), 'bloodNSD' = apply(meth_blood,1,sd,na.rm=TRUE))
  df_mean_sd_blood$bloodMean <- ifelse(df_meth_blood$CpG==1, df_mean_sd_blood$bloodMean, NA)
  df_mean_sd_blood$bloodSD <- ifelse(df_meth_blood$CpG==1, df_mean_sd_blood$bloodSD, NA)
  
  l_ply(colnames(df_mean_sd_blood), function(name, raw, seq, ss_id, locus) {
    
    # name <- colnames(exp_am_data)[20]
    # raw <-  CACend_raw
    # seq <- exp_am_data
    loc <- raw$loc[grep(paste0('^',name,'$'), raw$name_of_cell)]
    start_data <- paste0(unlist(strsplit(loc,'\\d+'))[1], 4)
    range_write(ss=ss_id, seq[,name, drop=F], sheet = locus, range=start_data, col_names = F, reformat = F)
    
  }, raw, df_mean_sd_blood, ss_id, locus)
  
}






