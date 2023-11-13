source('/media/HEAP-EPI/R/CRAN/rUtils/R/load_sheet.R')
devtools::load_all('/media/HEAP-EPI/R/CRAN/sangerBSC/')

# Здесь пихаешь ссылку на таблицу с инфой по точкам
ss_id <- ''

# Здесь читаешь лист с этой инфой и указываешь какой лист читать (loci)
locus_df <- as.data.frame(googlesheets4::read_sheet(ss=ss_id, sheet='loci'))

# Здесь выставляешь путь к папке проекта
proj_path <- ''

# здесь нужно из него data.frame сделать
locus_df <- as.data.frame(locus_df)

# здесь мутишь из него объект вида GRanges 

loci <- GRanges(locus_df$chr, IRanges(start=locus_df$start, end=locus_df$end, names=locus_df$Locus, bs_type=locus_df$bs_type))


# Получаешь результаты
# loci - granges с локусами
# group_list  - группы (должны совпадать с папками)
# proj_path - путь к папке с проектом
result_bis <- get_batch_loci(loci, 
                              proj_path, group_list = c('nobis','bloodN','bloodP'), genome = BSgenome.Hsapiens.UCSC.hg19)


# Чтобы были результаты по bloodU
# result_bisU <- get_batch_loci(loci, 
#                              proj_path, group_list = c('nobis','bloodN','bloodU','bloodP'), genome = BSgenome.Hsapiens.UCSC.hg19)

# присваиваем имена локусов в списке
names(result_bis) <- locus_df$Locus

#names(result_bisU) <- locus_df$Locus



# Смотреть можно на них с помощью View(test_df[[1]]), но вообще скорее этот объект нужен для дальнейшей загрузки в google sheets и удобного просмотра уже в них
test_df <- make_df_from_locus_list(result_bis)

#test_dfU <- make_df_from_locus_list_U(result_bisU)

names(test_df) <- locus_df$Locus

#names(test_dfU) <- locus_df$Locus

sotos_bis3_df <- test_df$`sotos-bis3`
sotos_bis643_df <- test_df$`sotos-bis643`
sotos_bis11_df <- test_df$`sotos-bis11`

#sotos_bis3_df <- test_dfU$`sotos-bis3`
#sotos_bis643_df <- test_dfU$`sotos-bis643`
#sotos_bis11_df <- test_dfU$`sotos-bis11`



#sotos_bis <- cbind(result_bis[[2]][["seq_df"]],result_bis[[2]][["groups_meth"]][["bloodN"]], result_bis[[2]][["groups_meth"]][["nobis"]], result_bis[[2]][["groups_meth"]][["bloodP"]])

sotos_bis3_plot <- subset(sotos_bis3_df, sotos_bis3_df$CpG==1)
melted_bis3_plot <- reshape2::melt(sotos_bis3_plot, measure.vars = colnames(sotos_bis3_plot)[6:37])
melted_bis3_plot$cat <- c(rep('bloodN', 10), rep('nometh', 40), rep('bloodP', 110))
melted_bis3_plot$G_position <- as.character(melted_bis3_plot$G_position)

theme_set(theme_classic())
ggplot2::ggplot(data=melted_bis3_plot, aes(x=G_position, y=value, fill=cat))+ggplot2::geom_boxplot() + labs(x='Координаты CpG', y='Метилирование (%)') +
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = 'top', axis.text = element_text(size=15), axis.title = element_text(size=16),
        legend.text = element_text(size=11), legend.title = element_text(size=14)) + scale_fill_discrete(name="Образцы",
                                                                                                         labels=c("Нормальный\nобразец (кровь)", "Патологический\nобразец (кровь)", "Нулевое метилирование")) 
ggsave('/path/bis3_boxplot.png', dpi=300, width=7.5, height = 5, units='in')



sotos_bis643_plot <- subset(sotos_bis643_df, sotos_bis643_df$CpG==1)
melted_bis643_plot <- reshape2::melt(sotos_bis643_plot, measure.vars = colnames(sotos_bis643_plot)[6:38])
melted_bis643_plot$cat <- c(rep('bloodN', 21), rep('nometh', 56), rep('bloodP', 154))
melted_bis643_plot$G_position <- as.character(melted_bis643_plot$G_position)

theme_set(theme_classic())
ggplot2::ggplot(data=melted_bis643_plot, aes(x=G_position, y=value, fill=cat))+ggplot2::geom_boxplot() + labs(x='Координаты CpG',y='Метилирование (%)') +
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = 'top', axis.text = element_text(size=15) ,axis.title = element_text(size=16),
        legend.text = element_text(size=11), legend.title = element_text(size=14)) + scale_fill_discrete(name="Образцы",
                                                                                                        labels=c("Нормальный\nобразец (кровь)", "Патологический\nобразец (кровь)", "Нулевое метилирование"))
ggsave('/path/bis643_boxplot.png', dpi=300, width=7.5, height = 5, units='in')


sotos_bis11_plot <- subset(sotos_bis11_df, sotos_bis11_df$CpG==1)
melted_bis11_plot <- reshape2::melt(sotos_bis11_plot, measure.vars = colnames(sotos_bis11_plot)[6:31])
melted_bis11_plot$cat <- c(rep('bloodN', 14), rep('nometh', 56), rep('bloodP', 112))
melted_bis11_plot$G_position <- as.character(melted_bis11_plot$G_position)

theme_set(theme_classic())
ggplot2::ggplot(data=melted_bis11_plot, aes(x=G_position, y=value, fill=cat))+ggplot2::geom_boxplot() + labs(x='Координаты CpG',y='Метилирование (%)') +
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = 'top', axis.text = element_text(size=15) ,axis.title = element_text(size=16),
        legend.text = element_text(size=11), legend.title = element_text(size=14)) + scale_fill_discrete(name="Образцы",
                                                                                                         labels=c("Нормальный\nобразец (кровь)", "Патологический\nобразец (кровь)", "Нулевое метилирование"))

ggsave('/path/bis11_boxplot.png', dpi=300, width=7.7, height = 5, units='in')



# Вспомогательные чтобы делать лист с датафреймами
make_df_from_locus_list <- function(locus_list) {
  
  res <- llply(1:length(locus_list), function(idx, locus_list) {
    
    frame <- cbind(locus_list[[idx]][["seq_df"]],
                   locus_list[[idx]][["groups_meth"]][["bloodN"]], 
                   locus_list[[idx]][["groups_meth"]][["nobis"]],
                   locus_list[[idx]][["groups_meth"]][["bloodP"]])
    
  }, locus_list)
  
  res
  
}

make_df_from_locus_list_U <- function(locus_list) {
  
  res <- llply(1:length(locus_list), function(idx, locus_list) {
    
    frame <- cbind(locus_list[[idx]][["seq_df"]],
                   locus_list[[idx]][["groups_meth"]][["bloodN"]], 
                   locus_list[[idx]][["groups_meth"]][["nobis"]],
                   locus_list[[idx]][["groups_meth"]][["bloodU"]],
                   locus_list[[idx]][["groups_meth"]][["bloodP"]])
    
  }, locus_list)
  
  res
  
}