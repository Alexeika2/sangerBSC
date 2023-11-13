#' @import googlesheets4
#' @import plyr
#' @import ape
#' @import reshape2
#' @import phangorn
#' @import stringr

require(googlesheets4)
require(plyr)
require(ape)
require(reshape2)
require(phangorn)
require(stringi)
require(stringr)

#bicondocutor packages
require(DECIPHER)
require(Biostrings)
require(sangerseqR)
require(sangeranalyseR)
require(BSgenome.Hsapiens.UCSC.hg19)



#devtools::load_all('/media/HEAP-EPI/R/CRAN/rUtils')


plot_raw_chromo <- function(abif,...) {
  nuc_mat <- make_nuc_mat(abif)
  matplot(nuc_mat, type='l', col=c('green','blue','black','red'), lty='solid', ylab='Fluorescence', xlab='Data points',...)
  abline(v=abif@data$PLOC.2)
  legend('topleft',legend=c('A','C','G','T'), col=c('green','blue','black','red'), fill = c('green','blue','black','red'))
}


get_seq_list <- function(fn_list, type){
  
  #fn <- blood_end_fn_list[1]
  res_list <- llply(fn_list, function(fn){ 
    ab1 <- read.abif(fn)
    idx_seq <- seq(1, length(ab1@data$P1AM.1))
    seq_ab1  <- unlist(strsplit(ab1@data$PBAS.2, ''))
    ab1@data$PBAS.2 <- paste(seq_ab1[idx_seq], collapse='')
    seq_corr_ab1 <- unlist(strsplit(ab1@data$PBAS.2, '')) %in% names(IUPAC_CODE_MAP)
    if (!all(seq_corr_ab1)) {
      stop('Some ab1 file has non supported character')
    }
    ab1
  })
  # names(res_list) <- paste0(type,str_match(basename(fn_list),'(\\d+)-')[,2])
  names(res_list) <- llply(res_list, function(ab1){ ab1@data$SMPL.1 })
  res_list
}


cost_funct <- function(coef, peaks, nobisDir, refBS) {
  # peaks <- peaks[,nobis]
  # nobisDir <- nobis_sd$nobisDir
  # seq <- ref_seq
  # coef <- c('A'=1,'C'=1,'T'=1)
  
  nobisCoef <- peakCorr(coef, peaks, refBS)
  nobisCoef_length <- sqrt(sum(nobisCoef^2, na.rm = TRUE))
  nobisDir_length <- sqrt(sum(nobisDir^2, na.rm = TRUE))
  
  # calculate cos of angle netween nobisCoef and nobisDir
  cosCoefDir <- sum(nobisCoef*nobisDir, na.rm = TRUE)/(nobisCoef_length*nobisDir_length)
  
  (1-cosCoefDir)^2 + (nobisCoef_length - nobisDir_length)^2/(nobisCoef_length*nobisDir_length)
}



peakCorr <- function(coef, peaks, refBS) {
  
  unlist(sapply(
    1:nchar(refBS),
    function(pos){
      # pos <- 1
      base <- as.character(substr(refBS, start = pos, stop = pos))
      peaks[pos]*coef[[base]]
    }))
  
  
}


optCorr <- function(peaks, nometh_Dir, refBS, coef) {

  # peaks <- exp_peaks_nobis

  peaks_corrected <- do.call(cbind.data.frame,
                             llply(colnames(peaks), function(nobis, peaks, nobisDir, refBS, coef) {
                               
                               # nobis <- "nobis1"
                               
                               coef_corr <- optim(coef, cost_funct, gr='CG', peaks = peaks[,nobis], nobisDir = nobisDir, refBS = refBS)$par
                               nobis_corr <- peakCorr(coef_corr, peaks[,nobis], refBS)
                               
                             }, peaks, nometh_Dir, refBS, coef))
  names(peaks_corrected) <- names(peaks)
  peaks_corrected
}


# если ОТ - маскируем С
# если СТОТ - маскируем G
# Возможно стоит брать маску пошире, только не очень ясно в какую сторону от CG сайта

mask_cpg <- function(refG, refBS) {
  
  cpg_content <- Biostrings::matchPattern('CG',DNAString(refG))
  
  cpg_mask <- rep(1, nchar(refG))
  for (cpg_i in 1:length(cpg_content@ranges)) {
    cpg <- cpg_content@ranges[cpg_i]
    for (pos_i in start(cpg):end(cpg)){
      if (substr(refG, start=pos_i, stop = pos_i) != substr(refBS, start=pos_i, stop=pos_i)) {
        cpg_mask[pos_i] <- 0
      }
    }
  }
  
  cpg_mask  
}

#' @description AM - AMplitude of basecall peaks (from ABIF file format specification)
exp_peaks_am <- function(seq_list, refBS){
  # matrix_align: last column - consensus, the one before last - reference.
  
  # seq_list <- seq_list_blood
  
  matrix_align <- {
    
    vec <- c()
    for (seqi in seq_list) {
      vec <- c(vec, seqi@data$PBAS.2)
    }
    
    vec <- c(vec, refBS)
    vec <- DNAStringSet(vec)
    names(vec) <- c(names(seq_list), 'ref_seq')
    res_align <- sangeranalyseR::merge.reads(vec)
    matrix_align <- t(as.matrix(res_align$alignment))
  }
  
  
  exp_col <- 1:(ncol(matrix_align) - 2)
  exp_pos <- rep(0, length(exp_col))
  exp_data <- data.frame(matrix(vector(), 0, 2, dimnames=list(c(), c('nobis1', 'nobis2'))))
  
  ref_col <- ncol(matrix_align)-1
  ref_pos <- 0
  
  for (i in 1:nrow(matrix_align)) {
    # i <- 16
    exp_base <- as.character(matrix_align[i, exp_col])
    ref_base <- as.character(matrix_align[i, ref_col])
    
    exp_pos <- exp_pos + ifelse(exp_base != '-', 1, 0)
    
    if( ref_base != '-'){
      ref_pos <- ref_pos + 1
      
      exp_am_row <- llply(
        exp_col,
        # exp_i <- 1
        function(exp_i, exp_base, ref_base, seq_list, exp_pos){
          if (exp_base[exp_i] != '-') {
            seq <- seq_list[[exp_i]]
            seq_pos <- exp_pos[exp_i]
            
            # pri - primary, sec - secondary
            seq_pri_base <- substr(seq@data$PBAS.2, seq_pos, seq_pos)
            seq_pri_base_iupac <- as.character(unlist(strsplit(IUPAC_CODE_MAP[seq_pri_base],'')))
            seq_pri_am <- seq@data$P1AM.1[seq_pos]
            seq_sec_base <- substr(seq@data$P2BA.1, seq_pos, seq_pos)
            seq_sec_am <- seq@data$P2AM.1[seq_pos]
            seq_pri_base <- seq_pri_base_iupac[!seq_pri_base_iupac %in% seq_sec_base]
            ifelse (seq_pri_base == ref_base, seq_pri_am, 
                    ifelse(seq_sec_base == ref_base, seq_sec_am, NA))
          }
          else NA
        }, 
        exp_base, ref_base, seq_list, exp_pos)
      
      exp_data <- rbind(exp_data, exp_am_row)
    }
  }
  names(exp_data) <- names(seq_list)
  rownames(exp_data) <- 1:nrow(exp_data)
  exp_data
}





plot_cpg_boxplot_old <- function(locus_obj, path) {
  
  
  # locus_obj <- batch_locus[[1]]
  # path <- '/media/HEAP-EPI/Claes-Jensen/L59mid/'
  
  locus <- locus_obj$Locus
  figures <- 'figures'
  dir.create(file.path(path, figures), showWarnings = FALSE)
  path_fig <- file.path(path, figures)
  
  
  scale_nobis_corr <- locus_obj$peaks_nobis_corr * t(locus_obj$nobisDscale)
  scale_blood_corr <- locus_obj$peaks_blood_corrected * t(locus_obj$nobisDscale)
  
  
  
  df_locus_sup <- cbind(locus_obj$seq_df, locus_obj$nobisSup, scale_nobis_corr, scale_blood_corr)
  df_sup <- subset(df_locus_sup, df_locus_sup$`locus_obj$nobisSup` == 1)
  df_sup$G_position <- paste0(df_sup$G_position,'_|',df_sup$refG)
  gather_cols <- colnames(df_sup)[7:ncol(df_sup)]
  keycol <- 'samples'
  valuecol <- 'RFU'
  df_sup_wide <- tidyr::gather_(df_sup, keycol,valuecol,gather_cols)
  df_sup_wide$groups <- gsub('[0-9]+','', df_sup_wide$samples)
  
  plot_all_supp <- ggplot(df_sup_wide, aes(x=G_position, y=RFU, fill=groups))+
    geom_boxplot()+
    theme(axis.text.x = element_text(size=15, angle = 55, hjust=1),
          axis.text.y = element_text(size=15), 
          legend.text = element_text(size=15),
          legend.title = element_text(size=15), 
          axis.title.x = element_blank(), 
          legend.position = 'top', 
          plot.title = element_text(size=17),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0('Locus: ',locus,'\n', 'All Sups'))+expand_limits(y=0)
  all_sup_plot <- file.path(path_fig, paste0(locus,'_all_sups.png'))
  ggsave(filename = all_sup_plot, plot=plot_all_supp, dpi=300, width = 35, height = 7)
  
  
  
  
  df_locus <- cbind(locus_obj$seq_df, scale_nobis_corr, scale_blood_corr)
  df_cpg <- subset(df_locus, df_locus$CpG == 1)
  df_cpg$G_position <- paste0(df_cpg$G_position,'_|',df_cpg$refG)
  gather_cols <- colnames(df_cpg)[6:ncol(df_cpg)]
  keycol <- 'samples'
  valuecol <- 'RFU'
  df_cpg_wide <- tidyr::gather_(df_cpg, keycol,valuecol,gather_cols)
  df_cpg_wide$groups <- gsub('[0-9]+','', df_cpg_wide$samples)
  df_cpg_wide$groups_genomics <- ifelse(df_cpg_wide$refG=='C', 
                                        paste0(df_cpg_wide$G_position,'_',dplyr::lead(df_cpg_wide$G_position,1)), 
                                        ifelse(df_cpg_wide$refG=='G', paste0(dplyr::lag(df_cpg_wide$G_position,1),'_',df_cpg_wide$G_position), 
                                               NA))
  
  plot_all <- ggplot(df_cpg_wide, aes(x=G_position, y=RFU, fill=groups))+
    geom_boxplot()+
    theme(axis.text.x = element_text(size=15, angle = 55, hjust=1),
          axis.text.y = element_text(size=15), 
          legend.text = element_text(size=15),
          legend.title = element_text(size=15), 
          axis.title.x = element_blank(), 
          legend.position = 'top', 
          plot.title = element_text(size=17),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0('Locus: ',locus,'\n', 'All CpGs'))+expand_limits(y=0)
  
  all_plot <- file.path(path_fig, paste0(locus,'_all_cpgs.png'))
  ggsave(filename = all_plot, plot=plot_all, dpi=300)
  
  
  df_cpg_split <- split(df_cpg_wide, df_cpg_wide$groups_genomics)  
  llply(df_cpg_split, function(df_cpgs, locus, path_fig) {
    
    
    # df_cpgs <- df_cpg_split[[2]]
    # locus <- locus_obj$Locus
    
    name_of_png <- paste0(locus,'_', unique(df_cpgs$groups_genomics),'.png')
    name_of_png <- gsub(':','_', name_of_png, fixed = TRUE)
    name_of_png <- gsub('|', '', name_of_png, fixed=TRUE)
    path_png <- file.path(path_fig, name_of_png)
    # png(filename=path_png)
    gplot <- ggplot(df_cpgs, aes(x=G_position, y=RFU, fill=groups))+
      geom_boxplot()+
      theme(axis.text.x = element_text(size=15),
            axis.text.y = element_text(size=15), 
            legend.text = element_text(size=15),
            legend.title = element_text(size=15), 
            axis.title.x = element_blank(), 
            legend.position = 'top', 
            plot.title = element_text(size=17),
            panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))+
      ggtitle(paste0('Locus: ',locus,'\n', 'CpGs: ', unique(df_cpgs$groups_genomics)))
    
    #setwd(path_fig)
    ggsave(filename = path_png, plot=gplot, dpi=300)
    
    #dev.off()
    
    
  }, locus, path_fig)
  
  #lapply(names(figs_list), function(x) {ggsave(filename = file.path(path_fig, paste0(x,'.png')), plot=figs_list[[x]], dpi=300)})
  
  
  
}

# Умеет рисовать вискеры с помощью кастомной функции

#' @export plot_cpg_boxplot
plot_cpg_boxplot <- function(locus_obj, path, ci=0.95) {
  
  if(class(locus_obj) != 'Locus'){
    stop("Only 'Locus' class allowed.")
  }
  
  # locus_obj <- batch_locuses[[6]]
  # path <- '/media/HEAP-EPI/Claes-Jensen/L59mid/'
  
  #path_loc <- file.path(path, locus)
  locus <- locus_obj$Locus
  path_loc <- file.path(path, locus)
  figures <- 'figures'
  dir.create(file.path(path, figures), showWarnings = FALSE)
  path_fig <- file.path(path_loc, figures)

  scale_nometh_corr <- 100-(locus_obj$groups_corrected[[1]] * t(locus_obj$nometh_Dscale))
  #scale_blood_corr <- locus_obj$peaks_blood_corrected * t(locus_obj$nobisDscale)
  
  cbind(scale_nometh_corr, locus_obj$groups_meth[[1]])
  
  # df_locus <- cbind(locus_obj$seq_df, scale_nometh_corr, locus_obj$methBlood, locus_obj$methBloodP)
  df_locus <- cbind(locus_obj$seq_df, locus_obj$groups_meth)
  df_cpg <- subset(df_locus, df_locus$CpG == 1)
  df_cpg$G_position <- paste0(df_cpg$G_position,'_|',df_cpg$refG)
  gather_cols <- colnames(df_cpg)[(ncol(locus_obj$seq_df)+1):ncol(df_cpg)]
  keycol <- 'samples'
  valuecol <- 'Methylation'
  df_cpg_wide <- tidyr::gather_(df_cpg, keycol, valuecol, gather_cols)
  df_cpg_wide$groups <- unlist(llply(str_split(df_cpg_wide$samples, '\\.'), '[[', 1))
  # gsub('[0-9]+','', df_cpg_wide$samples) 
  #df_cpg_wide$groups <- factor(df_cpg_wide$groups, levels = c('nobisD','bloodN','bloodP'))
  df_cpg_wide$groups_genomics <- ifelse(df_cpg_wide$refG=='C', 
                                        paste0(df_cpg_wide$G_position,'_',dplyr::lead(df_cpg_wide$G_position,1)), 
                                        ifelse(df_cpg_wide$refG=='G', paste0(dplyr::lag(df_cpg_wide$G_position,1),'_',df_cpg_wide$G_position), 
                                               NA))
  plot_all <- ggplot(df_cpg_wide, aes(x=G_position, y=Methylation, fill=groups, color=groups))+
    stat_summary(fun.data=custom_box, geom='boxplot', position = 'dodge2', fun.args = ci) + 
    stat_summary(fun = custom_outlier, geom="point", position = position_dodge(width=1.2)) + 
    stat_summary(fun.data = custom_box, geom = "errorbar", position='dodge2', fun.args = ci) + 
    scale_fill_manual(values = c('gray','green','red')) + 
    scale_color_manual(values = c('gray','green','red')) +
    theme(axis.text.x = element_text(size=15, angle = 55, hjust=1),
          axis.text.y = element_text(size=15), 
          legend.text = element_text(size=15),
          legend.title = element_text(size=15), 
          axis.title.x = element_blank(), 
          legend.position = 'top', 
          plot.title = element_text(size=17),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0('Locus: ',locus,'\n', 'All CpGs'))+expand_limits(y=-30)
  
  all_plot <- file.path(path_fig, paste0(locus,'_all_cpgs.png'))
  ggsave(filename = all_plot, plot=plot_all, dpi=300, width = 10, height = 6, units = 'in')
  
  # --
  # Собираем статистику для Вики Будадоржиевой
  # df_aggr <- aggregate(Methylation~G_position+groups, data=df_cpg_wide, FUN=customCI, ci=ci)
  # df_aggr$Locus <- locus
  # matr_aggr <- as.data.frame(df_aggr$Methylation)
  # df_aggr$Methylation <- NULL
  # df_aggr <- cbind(df_aggr, matr_aggr)
  # df_aggr <- df_aggr[,c(3,1,2,7,4,5,6)]
  # df_aggr_bloodN <- subset(df_aggr, df_aggr$groups == 'bloodN')
  # df_aggr_bloodP <- subset(df_aggr, df_aggr$groups == 'bloodP')
  # df_aggr_bloodP$groups <- NULL
  # df_aggr_bloodN$groups <- NULL
  # df_aggr_bloodP$Locus <- NULL
  # df_aggr_bloodP$lower <- NULL
  # df_aggr_bloodP$upper <- NULL
  # #df_aggr_bloodP$G_position <- NULL
  # df_aggr_bloodP$SD <- NULL
  # colnames(df_aggr_bloodP)[2] <- 'bloodP_Mean'
  # 
  # df_total <- merge(df_aggr_bloodN, df_aggr_bloodP, by='G_position')
  # df_total$significance <- ifelse(df_total$bloodP_Mean > df_total$lower & df_total$bloodP_Mean < df_total$upper, 0, 1)
  # 
  # df_total
  # --
  
  
  # df_cpg_split <- split(df_cpg_wide, df_cpg_wide$groups_genomics)  
  # llply(df_cpg_split, function(df_cpgs, locus, path_fig) {
  #   
  #   
  #   # df_cpgs <- df_cpg_split[[2]]
  #   # locus <- locus_obj$Locus
  #   
  #   name_of_png <- paste0(locus,'_', unique(df_cpgs$groups_genomics),'.png')
  #   name_of_png <- gsub(':','_', name_of_png, fixed = TRUE)
  #   name_of_png <- gsub('|', '', name_of_png, fixed=TRUE)
  #   path_png <- file.path(path_fig, name_of_png)
  #   # png(filename=path_png)
  #   gplot <- ggplot(df_cpgs, aes(x=G_position, y=RFU, fill=groups))+
  #     geom_boxplot()+
  #     theme(axis.text.x = element_text(size=15),
  #           axis.text.y = element_text(size=15), 
  #           legend.text = element_text(size=15),
  #           legend.title = element_text(size=15), 
  #           axis.title.x = element_blank(), 
  #           legend.position = 'top', 
  #           plot.title = element_text(size=17),
  #           panel.background = element_blank(), 
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank(), 
  #           axis.line = element_line(colour = "black"))+
  #     ggtitle(paste0('Locus: ',locus,'\n', 'CpGs: ', unique(df_cpgs$groups_genomics)))
  #   
  #   #setwd(path_fig)
  #   ggsave(filename = path_png, plot=gplot, dpi=300)
  #   
  #   #dev.off()
  #   
  #   
  # }, locus, path_fig)
  
  #lapply(names(figs_list), function(x) {ggsave(filename = file.path(path_fig, paste0(x,'.png')), plot=figs_list[[x]], dpi=300)})
  
}

custom_box <- function(x,ci) {
  r <- quantile(x, probs = c(0.25, 0.5, 0.75))
  yminmax <- Rmisc::CI(x, ci=ci)
  r <- c(yminmax[['lower']],r,yminmax[['upper']])
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

custom_outlier <- function(x) {
  subset(x, x < custom_box(x)[1] | custom_box(x)[5] < x)
}

customCI <- function (x, ci = 0.95) {
  a <- mean(x)
  s <- sd(x)
  n <- length(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(c(upper = a + error, mean = a, SD=s,lower = a - error))
}

