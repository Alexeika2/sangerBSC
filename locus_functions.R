# Connecting helpers functions

source('helpers.R')



get_batch_locus <- function(ss_id, loci_table) {
  
  path_proj <- '/media/HEAP-EPI/Claes-Jensen/'
  
  loci_gs <- load_sheet(ss_id = ss_id, sheet = loci_table)
  
  loci_gs <- as.data.frame(loci_gs)
  
  loci_gs <- loci_gs[loci_gs$plot==1,] # Только те у кого есть pathologic blood
  
  llply(loci_gs$Locus, function(locus, loci_df, path_proj) {
    
    
    chr <-  loci_df$chr[loci_df$Locus==locus]
    start_g <- as.numeric(loci_df$start[loci_df$Locus==locus])
    end_g <- as.numeric(loci_df$end[loci_df$Locus==locus])
    
    
    if(loci_df$bs_type[loci_df$Locus == locus] == 'OT') {
      
      refG <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, chr, 
                                 start_g,
                                 end_g)
      refBS <- Biostrings::chartr('C','T', refG)
      refBS <- as.character(refBS)
      refG <- as.character(refG)
      G_position <- paste0(chr,':',seq(start_g, end_g))
      coef <- c('A'=1,'G'=1,'T'=1)
      
    }
    else if (loci_df$bs_type[loci_df$Locus == locus] == 'CTOT') {
      
      refG <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, chr, 
                                 start_g,
                                 end_g)
      refBS <- Biostrings::chartr('C','T', refG)
      refBS <- as.character(Biostrings::reverseComplement(refBS))
      refG <- as.character(Biostrings::reverseComplement(refG))
      G_position <- paste0(chr,':',seq(end_g, start_g))
      coef <- c('A'=1,'C'=1,'T'=1)
      
    }
    else {stop('No such type in bs_type')}
    
    path_loc <- file.path(path_proj, locus)
    
    get_locus(locus, path_loc, refG, refBS, coef, chr, start_g, end_g, G_position)  
    
    
  }, loci_gs, path_proj)
  
  
}


get_locus <- function(locus, path, refG,refBS, coef, chr, start, end, G_position) {
  
  #' @description Main function to read ab1 files, calculates corrected amplification values for nucleotides and methylation values. 
  #' @param locus String Description for the Locus, "TUBGstart" for example
  #' @param path String path to files with three folders: 'nobis', 'bloodN', 'bloodP'
  #' @param refG Genomic Sequence 
  #' @param refBS Bisulfite treated genomic sequence 
  #' @param coef Nucleotide coefficients for optimization purposes, For "OT" it should be c('A'=1,'G'=1,'T'=1), whereas for CTOT  is c('A'=1,'C'=1,'T'=1). Feel Free to experiment with numbers for coefficients.
  #' @param chr Chromosome, 'chr1' for example
  #' @param start start position for locus 
  #' @param end end position for locus
  #' @param G_position for OT its genomic positions from start to end, for CTOT its positions from end to start. Example: paste0(chr,':',seq(end, start)) for CTOT.
  
  
  # фиксируем имена папок
  subfolders <- c('nobis', 'bloodN', 'bloodP')
  
  # пути до файлов nobis и blood
  path_nobis <- with(
    list(nobis = list.files(file.path(path, subfolders[1]), full.names = T, pattern = '\\.ab1$')),
    nobis[order(as.numeric(str_match(nobis, '(\\d+)-')[,2]))]  
  )
  
  path_blood <- with(
    list(blood = list.files(file.path(path, subfolders[2]), full.names = T, pattern = '\\.ab1$') ),
    blood[order(as.numeric(str_match(blood, '(\\d+)-')[,2]))]  
  )
  
  path_bloodP <- with(
    list(bloodP = list.files(file.path(path, subfolders[3]), full.names = T, pattern = '\\.ab1$') ),
    bloodP[order(as.numeric(str_match(bloodP, '(\\d+)-')[,2]))]  
  )
  
  cat(paste0('Getting abif files for nobis to locus: ', locus, '\n'))
  
  seq_list_nobis <- seq_list(path_nobis, subfolders[1])
  
  cat(paste0('Getting abif files for bloodN to locus: ', locus, '\n'))
  
  seq_list_blood <- seq_list(path_blood, subfolders[2])
  
  cat(paste0('Getting abif files for bloodP to locus: ', locus, '\n'))
  
  seq_list_bloodP <- seq_list(path_bloodP, subfolders[3])
  
  
  cat(paste0('Getting amplification positions for nobis to locus: ', locus, '\n'))
  
  exp_peaks_nobis <- exp_peaks_am(seq_list_nobis, refBS)
  
  cat(paste0('Getting amplification positions for bloodN to locus: ', locus, '\n'))
  
  exp_peaks_blood <- exp_peaks_am(seq_list_blood, refBS)
  
  cat(paste0('Getting amplification positions for bloodP to locus: ', locus, '\n'))
  
  exp_peaks_bloodP <- exp_peaks_am(seq_list_bloodP, refBS)
  
  
  cat(paste0('Collecting SDs, Supports and Direction for nobis to locus: ', locus, '\n'))
  
  nobisSD <- apply(exp_peaks_nobis, 1, sd, na.rm = TRUE)
  nobisSup <- ifelse(nobisSD < 50, 1, 0)
  
  cpg_mask <- mask_cpg(refG, refBS)
  cpg <- ifelse(cpg_mask==1,0,1)
  
  
  # чтобы CG не участовали в оптимизации
  nobisSup <- pmin(nobisSup, cpg_mask)
  
  nobisDir <- apply(exp_peaks_nobis, 1, median, na.rm=TRUE)
  nobisDir <- ifelse(nobisSup==1, nobisDir, NA)
  
  
  
  
  nobis_corrected <- optCorr(peaks = exp_peaks_nobis, nobisDir = nobisDir, refBS = refBS, coef = coef)
  blood_corrected <- optCorr(peaks = exp_peaks_blood, nobisDir = nobisDir, refBS = refBS, coef = coef)
  bloodP_corrected <- optCorr(peaks = exp_peaks_bloodP, nobisDir = nobisDir, refBS = refBS, coef = coef)
  
  
  colnames(nobis_corrected) <- paste('nobisD', seq(1, ncol(exp_peaks_nobis)), sep='')
  colnames(blood_corrected) <- colnames(exp_peaks_blood)
  colnames(bloodP_corrected) <- colnames(exp_peaks_bloodP)
  
  
  
  nobisDscale <- 100/apply(nobis_corrected, 1, median, na.rm = TRUE)
  
  meth_blood_corr <- 100 - (blood_corrected * t(nobisDscale))
  meth_bloodP_corr <- 100 - (bloodP_corrected * t(nobisDscale))
  

  
  nobisDirD <- apply(nobis_corrected, 1, median, na.rm=TRUE)
  
  
  
  cat(paste0('Getting corrected amplification positions for bloodN and nobis to locus: ', locus, '\n'))
  
  nobisSD_corr <- apply(nobis_corrected, 1, sd, na.rm=TRUE)
  bloodN_SD_corr <- apply(blood_corrected, 1, sd, na.rm=TRUE)
  
  seq_df <- data.frame('Position' = seq(1, length(unlist(strsplit(refG,'')))), 
                       'refG' = unlist(strsplit(refG,'')),
                       'refBS' = unlist(strsplit(refBS, '')),
                       'G_position' = G_position,
                       'CpG' = cpg)
  
  
  
  output <- list('Locus' = locus, 'refG' = refG, 'refBS'=refBS, 'nobis_seq_list' = seq_list_nobis, 'blood_seq_list' = seq_list_blood, 'peaks_nobis' = exp_peaks_nobis, 'peaks_blood' = exp_peaks_blood, 
                 'nobisSD' = nobisSD, 'nobisSup' = nobisSup, 'nobisDir' = nobisDir, 'peaks_nobis_corr' = nobis_corrected, 'peaks_blood_corrected' = blood_corrected, 
                 'nobisSD_corr' = nobisSD_corr, 'bloodN_SD_corr' = bloodN_SD_corr, 'nobisDirD' = nobisDirD, 'nobisDscale' = nobisDscale, 'seq_df' = seq_df, 'methBlood' = meth_blood_corr, 'methBloodP' = meth_bloodP_corr)
  
  
  class(output) <- 'Locus'
  
  cat(paste0('DONE for locus: ', locus, '\n'))
  
  
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






