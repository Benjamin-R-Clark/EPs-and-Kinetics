library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(readr)
library(lmtest)
library(ggfortify)
library(car)


#making a named vector where the names are the left anchor and values are right
left_anchor <- read.delim("data/left_anchor.bed", header = F)[,1]
right_anchor <- read.delim("data/right_anchor.bed", header = F)[,1]


#removing self linkages
anchors <- cbind(left_anchor, right_anchor)
anchors <- anchors[-which(anchors[,1] == anchors[,2]),]

bed_convert <- function(item){
  gsub("[:-]", "\t", item)
}

left_formatted_contacts <- sapply(anchors[,1], bed_convert, USE.NAMES = F)
right_formatted_contacts <- sapply(anchors[,2], bed_convert, USE.NAMES = F)


write(left_formatted_contacts, file = "data/left_anchor_formatted.bed")
write(right_formatted_contacts, file = "data/right_anchor_formatted.bed")

#bedtools intersect each with a promoter file


#promoter hits

left_promoters <- read.delim("data/left_promoters_filtered.bed", header = F)
right_promoters <-  read.delim("data/right_promoters_filtered.bed", header = F)

#turning coordinates into a string I can use as a key
toKey <- function(row){
  return(paste0(row[1], ":", row[2], "-", row[3]))
}

l_coors <- left_promoters[,1:3]
l_promoter_keys <- l_coors %>% dplyr::rowwise() %>% dplyr::mutate(key = toKey(c(V1,V2,V3)))
left_promoters <- cbind(left_promoters, key = l_promoter_keys$key)

r_coors <- right_promoters[,1:3]
r_promoter_keys <- r_coors %>% dplyr::rowwise() %>% dplyr::mutate(key = toKey(c(V1,V2,V3)))
right_promoters <- cbind(right_promoters, key = r_promoter_keys$key)
#here we try to identify promoter-promoter contacts.
#if the left key is a promoter, is its linked right key also a promoter?

promoter_promoter <- which(anchors[as.vector(l_promoter_keys$key)] %in% r_promoter_keys$key)
promoter_promoter <- anchors[promoter_promoter]
#84 possible promoter-promoter contacts

#which of these are actually enhancers in promoters?
right_enhancers <-  read.delim("data/right_enhancers.bed", header = F)
left_enhancers <- read.delim("data/left_enhancers.bed", header = F)

l_enh_keys <- left_enhancers[,1:3] %>% rowwise() %>% mutate(key = toKey(c(V1,V2,V3)))
r_enh_keys <- right_enhancers[,1:3] %>% rowwise() %>% mutate(key = toKey(c(V1,V2,V3)))



l_enhancer_in_promoter <-  l_enh_keys$key[which(l_enh_keys$key %in% l_promoter_keys$key)]
#One left enhancer is in a promoter

r_enhancer_in_promoter <-  r_enh_keys$key[which(r_enh_keys$key %in% r_promoter_keys$key)]
#one right enhancer is in a promoter


#Now lets check if they're in the promoter-promoter list

l_enhancer_in_promoter %in% names(promoter_promoter)
#FALSE
r_enhancer_in_promoter %in% promoter_promoter
#FALSE

#All promoter-promoter contacts seem to be genuine, lets remove them from our promoter list and then perform joins on the tables
filtered_l_promoters <- left_promoters[which(!left_promoters$key %in% names(promoter_promoter)),]
filtered_r_promoters <- right_promoters[which(!right_promoters$key %in% names(promoter_promoter)),]


clipname <- function(fullname){
  return(substr(fullname, 1, nchar(fullname) - 2))
}


#I'm noticing we have two promoter annotations which would inflate the contact frequency.
#I'm going to try to keep best annotation set (the one with the most interactions per promoter)
all_linked_promoters <- rbind(filtered_l_promoters,filtered_r_promoters)
all_linked_promoters["gene"] <- sapply(all_linked_promoters$V7, FUN = clipname)

genes <- all_linked_promoters$gene
unique_linked_promoters <- all_linked_promoters
for(gene_name in unique(genes)){
  gene_group <- all_linked_promoters %>% filter(gene == gene_name)
  if(length(unique(gene_group$V7)) > 1){
    not_gene_anno <- names(which.min(table(gene_group$V7)))
    unique_linked_promoters <- unique_linked_promoters %>% filter(V7 != not_gene_anno)
  }
}
unique_linked_genes <- unique_linked_promoters$gene
genereps <- table(unique_linked_genes)



#annotating to ens genes
ens <- read.delim("data/MGIBatchReport_annos.txt")
ens_unique <- ens %>% distinct(Input, .keep_all = T)

annotated_links <- unique_linked_promoters %>% inner_join(y = ens_unique, by = join_by(gene == Input))
unique_links <- annotated_links %>% distinct(Ensembl.ID, .keep_all = T)


ordered_gene_reps <-  genereps[order(match(names(genereps), unique_links$gene, nomatch = NA))]
ordered_gene_reps <-  ordered_gene_reps[which(names(ordered_gene_reps) %in% unique_links$gene)]
unique_links <- cbind(unique_links, ordered_gene_reps)
unique_links["anchor"] <- TRUE




#loading in kinetics
kinetics <- read_excel("data/journal.pcbi.1008772.s006.xlsx", sheet = "C57")
colnames(kinetics) <- kinetics[1,]
kinetics <- kinetics[-1,]

joined_df <- kinetics %>% dplyr::left_join(y = unique_links, by = join_by(gene_id == Ensembl.ID))
joined_df["anchor"] <- !is.na(joined_df["anchor"])



joined_df$k_on <- as.numeric(joined_df$k_on)
joined_df$k_off <- as.numeric(joined_df$k_off)
joined_df$k_syn <- as.numeric(joined_df$k_syn)

joined_df["bs"] <- joined_df$k_syn / joined_df$k_off
joined_df$Freq[which(is.na(joined_df$Freq))] <- 0

###########plotting anchor vs non-anchor##################
#library(ggridges)
ggplot(joined_df) + geom_density(mapping = aes(x = k_on, color = anchor, fill = anchor), alpha = 0.3) +
  scale_x_log10() + theme_bw()

bf_fig <- ggplot(joined_df, aes(x=factor(anchor), k_on, fill = anchor)) +
  #geom_violin() +
  geom_boxplot(show.legend = F, linewidth = 1.3, fatten = T) +
  theme_linedraw(base_size = 30) +
  theme(axis.text.x=element_text(size=rel(0.95))) +
  ylab( "Burst Frequency") +
  scale_y_log10() +
  xlab("") +
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(name = "Enhancer-Promoter Contact", type = c("#F0E442", "#56B4E9")) +
  #scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  stat_compare_means(method = "t.test", label.x.npc = 0.28, size = 9)

bs_fig <- ggplot(joined_df, aes(x=factor(anchor), bs, fill = anchor)) +
  #geom_violin() +
  geom_boxplot(linewidth = 1.3, fatten = T) +
  theme_linedraw(base_size = 30) +
  theme(axis.text.x=element_text(size=rel(0.95))) +
  ylab( "Burst Size") +
  scale_y_log10() +
  xlab("") +
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(name = "Enhancer-Promoter Contact", type = c("#F0E442", "#56B4E9")) +
  #scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  stat_compare_means(method = "t.test", label.x.npc = 0.28, size = 9)

ggarrange(bf_fig, bs_fig, common.legend = T)

ggplot(joined_df, aes(x=factor(anchor), k_on, fill = anchor)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_minimal(base_size = 21) +
  ylab( "Burst Frequency") +
  scale_y_log10() +
  xlab("") +
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(name = "Enhancer-Promoter \n Contact") +
  stat_compare_means(method = "t.test")


ggplot(joined_df, aes(x=factor(anchor), bs, fill = anchor)) +
  geom_violin() +
  stat_compare_means(method = "t.test") +
  geom_boxplot(width=0.1) +
  theme_minimal(base_size = 12) +
  ylab( "Burst Size") +
  scale_y_log10() +
  xlab("") +
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(name = "Enhancer-Promoter \n Contact")


ggplot(joined_df) + geom_density(mapping = aes(x = bs, color = anchor, fill = anchor), alpha = 0.3) +
  scale_x_log10() + theme_bw()

#t.test and measuring effect size
test_kon <- t.test(x = log10(joined_df$k_on[which(joined_df$anchor)]), y = log10(joined_df$k_on[which(!joined_df$anchor)]))

test_bs <-t.test(x = log10(joined_df$bs[which(joined_df$anchor)]), y = log10(joined_df$bs[which(!joined_df$anchor)]))

wilcox_effsize(joined_df, bs ~ anchor)
wilcox_effsize(joined_df, k_on ~ anchor)

## absolute effect size in hours= abs((1/10^-0.004669) - (1/10^-0.2052)) * 4.8
# 2.847241

## absolute effect size in number of transcripts per burst
#  abs((10^0.67764) - (10^0.612454))
# 0.6634751

#testing to see if multiple anchors do any to bf

ggplot(joined_df, aes(factor(Freq), k_on)) + geom_violin() + geom_jitter(size = 0.4) + scale_y_log10()

freq_means <- joined_df %>% dplyr::group_by(Freq) %>% dplyr::summarise(mean = mean(k_on),
                                                         median = median(k_on))


freq_means <- cbind(freq_means,  "n" = dplyr::count(joined_df, Freq)$n)




 ggplot(joined_df, aes(factor(Freq), k_on)) +
  geom_violin(scale = "width", linewidth = 1) +
  geom_dotplot(binaxis = "y", binwidth = 0.015,  mapping = aes(alpha = factor(Freq)), color = "blue") +
  scale_alpha_manual(values = c(0, seq(from = 0.2, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  geom_smooth(data = freq_means, mapping= aes(Freq, mean), color = "red") +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Log10(Burst Frequency)") +
  xlab("Number of Promoter Anchors")

bs_means <- joined_df %>% dplyr::group_by(Freq) %>% dplyr::summarise(mean = mean(bs),
                                                       median = median(bs))
bs_means <- cbind(bs_means, "n" = dplyr::count(joined_df, Freq)$n)

ggplot(joined_df, aes(factor(Freq), bs)) +
  geom_violin(linewidth = 1) +
  geom_dotplot(binaxis = "y", binwidth = 0.02,  mapping = aes(alpha = factor(Freq)), color = "blue") +
  scale_alpha_manual(values = c(0, seq(from = 0.15, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  geom_smooth(data = bs_means, mapping= aes(Freq, mean), color = "red") +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Log10(Burst Size)") +
  xlab("Number of Promoter Anchors")

## plots with anchors above 8 grouped
joined_df <- joined_df %>% mutate(string.Freq = as.character(Freq))

joined_df[which(joined_df$Freq > 8),]$string.Freq <- "8<"

freq_means <- joined_df %>% dplyr::group_by(string.Freq) %>% dplyr::summarise(mean = mean(k_on),
                                                         median = median(k_on),
                                                          sd = sd(k_on))


freq_means <- cbind(freq_means,  "n" = dplyr::count(joined_df, string.Freq)$n)
freq_means$string.Freq <- seq(0,9)



ggplot(joined_df, aes(factor(string.Freq), k_on)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.7) +
  geom_dotplot(binaxis = "y", binwidth = 0.015,
                mapping = aes(alpha = factor(string.Freq)), color = "blue", dotsize = 1.5, stackdir = "center") +
  scale_alpha_manual(values = c(0, seq(from = 0.2, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  geom_smooth(data = joined_df, mapping= aes(Freq, k_on), color = "red") +
  coord_cartesian(xlim = c(0,11)) +
  theme_minimal(base_size = 27) +
  theme(legend.position = "none", axis.text.x=element_text(size=rel(0.95))) +
  ylab("Burst Frequency") +
  xlab("Number of Enhancer Contacts")

ggplot(joined_df, aes(factor(string.Freq), bs)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.7) +
  geom_dotplot(binaxis = "y", binwidth = 0.015,
               mapping = aes(alpha = factor(string.Freq)), color = "blue", dotsize = 1.5, stackdir = "center") +
  scale_alpha_manual(values = c(0, seq(from = 0.2, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  geom_smooth(data = joined_df, mapping= aes(Freq, bs), color = "red") +
  coord_cartesian(xlim = c(0,11)) +
  theme_minimal(base_size = 27) +
  theme(legend.position = "none", axis.text.x=element_text(size=rel(0.95))) +
  ylab("Burst Size") +
  xlab("Number of Enhancer Contacts")



library(lmtest)
library(ggfortify)
library(aomisc)

kon_lm <- lm(log10(k_on) ~ Freq, joined_df)
bs_lm <- lm(log10(bs) ~ Freq, joined_df)





#Adding expression values
c57_ex <- read_excel("data/journal.pcbi.1008772.s007.xlsx")
c57_ex <- c57_ex %>% dplyr::select(1,14)
colnames(c57_ex) <- c("ens", "c57_exprs")
c57_ex <- c57_ex[-1,]
c57_ex <- c57_ex %>% dplyr::mutate(across(c57_exprs, as.numeric))
joined_df_exprs <- c57_ex %>% dplyr::left_join(y = joined_df, by = join_by(ens == "gene_id"))


#kinetics_exprs_ep <- joined_df_exprs %>% select()
#Adding loess expression curve to figures

ggplot(joined_df_exprs, aes(factor(string.Freq), k_on)) +
  geom_violin(scale = "width", linewidth = 1) +
  geom_dotplot(binaxis = "y", binwidth = 0.015,  mapping = aes(alpha = factor(string.Freq)), color = "blue") +
  scale_alpha_manual(values = c(0, seq(from = 0.2, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  geom_smooth(data = joined_df_exprs, mapping= aes(Freq, k_on), color = "red") +
  geom_smooth(data = joined_df_exprs, mapping = aes(Freq, c57_exprs), color = "orange") +
  coord_cartesian(xlim = c(0,11)) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none") +
  ylab("Burst Frequency") +
  xlab("Number of Anchored Enhancers")

ggplot(joined_df_exprs, aes(factor(string.Freq), bs)) +
  geom_boxplot(outliers = F) +
  geom_dotplot(binaxis = "y", binwidth = 0.015,  mapping = aes(alpha = factor(string.Freq)), color = "blue") +
  scale_alpha_manual(values = c(0, seq(from = 0.2, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  coord_cartesian(xlim = c(0,11)) +
  geom_smooth(data = joined_df_exprs, mapping= aes(Freq, bs), color = "red") +
  geom_smooth(data = joined_df_exprs, mapping = aes(Freq, c57_exprs), color = "orange") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none") +
  ylab("Burst Size") +
  xlab("Number of Anchored Enhancers")

exprs_f <- ggplot(joined_df_exprs, aes(factor(string.Freq), c57_exprs)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  geom_dotplot(binaxis = "y", binwidth = 0.012,
                mapping = aes(alpha = factor(string.Freq)), color = "blue",
                stackdir = "center", dotsize = 1.7) +
  scale_alpha_manual(values = c(0, seq(from = 0.1, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  coord_cartesian(xlim = c(0,11)) +
  geom_smooth(data = joined_df_exprs, mapping = aes(Freq, c57_exprs), color = "orange") +
  theme_gray(base_size = 28) +
  theme(legend.position = "none", axis.text.x=element_text(size=rel(0.95))) +
  ylab("Mean Expression (# UMI)") +
  xlab(NULL)

kon_mul <- ggplot(joined_df, aes(factor(string.Freq), k_on)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  geom_dotplot(binaxis = "y", binwidth = 0.012, dotsize = 1.5,
               mapping = aes(alpha = factor(string.Freq)), color = "blue", stackdir = "center") +
  scale_alpha_manual(values = c(0, seq(from = 0.1, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  geom_smooth(data = joined_df, mapping= aes(Freq, k_on), color = "red") +
  coord_cartesian(xlim = c(0,11)) +
  theme_gray(base_size = 28) +
  theme(legend.position = "none", axis.text.x=element_text(size=rel(0.95))) +
  ylab("Burst Frequency (1/min)") +
  xlab(NULL)

bs_mul <- ggplot(joined_df, aes(factor(string.Freq), bs)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  geom_dotplot(binaxis = "y", binwidth = 0.012,
               mapping = aes(alpha = factor(string.Freq)),
               color = "blue", dotsize = 1.5, stackdir = "center") +
  scale_alpha_manual(values = c(0, seq(from = 0.1, to = 1 , by = 0.9/20))) +
  scale_y_log10() +
  coord_cartesian(xlim = c(0,11)) +
  geom_smooth(data = joined_df, mapping= aes(Freq, bs), color = "red") +
  theme_gray(base_size = 28) +
  theme(legend.position = "none", axis.text.x=element_text(size=rel(0.95))) +
  ylab("Burst Size (# mRNA)") +
  xlab(NULL)


f <- gridExtra::arrangeGrob(kon_mul,bs_mul,exprs_f,  ncol = 3, nrow = 1 ,
                            bottom=grid::textGrob(label= "Number of Enhancer Contacts",
                                                  gp = grid::gpar(fontsize=20)))

grid::grid.newpage()
grid::grid.draw(f)

#####bootstrapping linear model for bs#########
binned_kinetics <- joined_df_exprs %>% arrange(Freq) %>% mutate(freq_bin = findInterval(Freq, c(1,2,3,4,5,6,7,8)),
                                                                c57_exprs = unlist(c57_exprs))
kon_lm_poly <- lm(formula = log10(k_on) ~ poly(freq_bin, 2), data = binned_kinetics)
exp_poly <- lm(formula = log10(binned_kinetics$c57_exprs) ~ poly(freq_bin, 2), data = binned_kinetics)
exp <- lm(formula = log10(binned_kinetics$c57_exprs) ~ freq_bin, data = binned_kinetics)
#bining data


binned_lm <- lm(log10(bs) ~ freq_bin, data = binned_kinetics)
binned_kinetics %>% dplyr::count(freq_bin)
bin_stats <- binned_kinetics %>% dplyr::group_by(factor(freq_bin)) %>%
  dplyr::summarise(range = str(range(Freq)))

get_coef <- function(data, i){
  lm <-  lm(log10(bs) ~ freq_bin, data = data[i,])
  return(c(lm$coefficients, summary(lm)$r.squared))
}

get_coef_kon <- function(data, i){
  lm <-  lm(log10(k_on) ~ poly(freq_bin,2), data = data[i,])
  return(c(lm$coefficients, summary(lm)$r.squared))
}


bs_ordinary_boot <- boot::boot(data = binned_kinetics, statistic = get_coef,
                               R = 5000)
bf_ordinary_boot <- boot::boot(data = binned_kinetics, statistic = get_coef_kon, R = 5000)


quantile(bs_ordinary_boot$t[,2], probs = c(0.025, 0.975), type = 6)
quantile(bf_ordinary_boot$t[,2], probs = c(0.025, 0.975), type = 6)
###### random effects model############


ggplot(binned_kinetics, aes(y = c57_exprs, x = bs)) +
  geom_point() +
  xlab("Exprs") +
  ylab("Kon") +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", formula = 'y ~ x', se=F,fullrange = T) +
  facet_wrap(~Freq)+
  scale_y_log10() +
  scale_x_log10()

lme4::lmer(log10(c57_exprs) ~ log10(k_on) + (k_on | Freq), data = binned_kinetics)


#########################Exploring promoter linkage distance################################

hichip <- read.csv("data/H3K27ac_hichip.csv")
joined_df_linked <- joined_df[which(joined_df$anchor),]
joined_keys <- joined_df_linked$key
joined_keys_left_indexes <- which(joined_keys %in% hichip$ATAC.summit..left.anchor.)
joined_keys_right_indexes <- which(joined_keys %in% hichip$ATAC.summit..right.anchor.)


filtered_hichip <- hichip %>% dplyr::filter(ATAC.summit..left.anchor. %in% joined_df_linked$key | ATAC.summit..right.anchor. %in% joined_df_linked$key)


get.coor <- function(string){
#there must be a prettier way of doing this
  return(as.numeric(strsplit(strsplit(string, ":")[[1]][2], "-")[[1]][1]))

}

stripped_filtered_hichip <- base::as.data.frame(apply(filtered_hichip, FUN = get.coor, c(1,2)))
chrom_dist <- abs(stripped_filtered_hichip[,1] - stripped_filtered_hichip[,2])

filtered_hichip["distance"] <- chrom_dist

joined_df_linked_dist_l <-  joined_df_linked %>% inner_join(filtered_hichip, by = join_by("key" == "ATAC.summit..left.anchor."))
joined_df_linked_dist_r <-  joined_df_linked %>% inner_join(filtered_hichip, by = join_by("key" == "ATAC.summit..right.anchor."))
#joined_df_linked_dist <-  joined_df_linked_dist_l %>% full_join(filtered_hichip, by = join_by("key" == "ATAC.summit..right.anchor."))

ones_df <- joined_df_linked_dist_l %>% filter(Freq == 1) %>%
  dplyr::mutate("ATAC.summit..left.anchor." = key)
ones_df_r <- joined_df_linked_dist_r %>% filter(Freq == 1) %>%
  dplyr::mutate("ATAC.summit..right.anchor." = key)

ones_df <- ones_df %>% union(ones_df_r)
ones_df$distance
#ones_df <- ones_df %>% inner_join(ones_df_r, by = join_by("gene_id"))

#ones_df <- ones_df %>% dplyr::full_join(ones_df_r, by = "gene_id")

#tdistance <- ones_df$distance
#tdistance[which(is.na(ones_df$distance.x))] <- ones_df$distance.y[which(is.na(ones_df$distance.x))]
#ones_df["distance"] <- tdistance


ones_df <- ones_df %>% dplyr::rename("promoter_start" = "V10","promoter_end" = "V11")

ones_left <- sapply(ones_df$ATAC.summit..left.anchor., bed_convert, USE.NAMES = F)
#im too tired to think of another way of doing this
write(ones_left, file = "data/left_ones_anchor.bed")
ones_left <- read.delim("data/left_ones_anchor.bed", header = F)
ones_left$gene <- ones_df$unique_linked_genes
write_delim(ones_left, delim = "\t", file = "data/left_ones_anchor.bed",col_names = F)

ones_right <- sapply(ones_df$ATAC.summit..right.anchor., bed_convert, USE.NAMES = F)
write(ones_right, file = "data/right_ones_anchor.bed")
ones_right <- read.delim("data/right_ones_anchor.bed", header = F)
ones_right$gene <- ones_df$unique_linked_genes
write_delim(ones_right, delim = "\t", file = "data/right_ones_anchor.bed", col_names = F)

left_ctcf <- read.delim("data/left_ones_ctcf.bed", header = F)
right_ctcf <- read.delim("data/right_ones_ctcf.bed", header = F)

double_ctcf <- left_ctcf %>% dplyr::full_join(right_ctcf, by=join_by(V4))
no_ctcf <- which(!(ones_left$gene %in% double_ctcf$V4))
no_ctcf <- ones_left[no_ctcf,]

double_ctcf$double <- FALSE
double_ctcf$double[which(complete.cases(double_ctcf))] <- T
double_ctcf$single <- FALSE
double_ctcf$single[which(!complete.cases(double_ctcf))] <- T

ones_df$ctcf <- 0
ones_df$ctcf[which(ones_df$Symbol %in% double_ctcf$V4[which(double_ctcf$double)])] <- 2
ones_df$ctcf[which(ones_df$Symbol %in% double_ctcf$V4[which(double_ctcf$single)])] <- 1

scale <- c("2" = "red", "1" = "blue","0" = "black" )
kon_d <- ggplot(ones_df) + geom_point(mapping = aes(x = log10(distance), y = k_on, color = factor(ctcf)), size = 2.5) +
 theme_classic(base_size = 18) + xlab(NULL) + ylab("Burst Frequency") + scale_y_log10()

bs_d <- ggplot(ones_df) + geom_point(size = 2.5, mapping = aes(x = log10(distance), y = bs, color = factor(ctcf))) +
  theme_classic(base_size = 18) + xlab(NULL) + ylab("Burst Size") + scale_y_log10()



f2 <- gridExtra::arrangeGrob(kon_d,bs_d,  ncol = 2, nrow = 1 ,
                            bottom=grid::textGrob(label= "log Distance (bp)",
                                                  gp = grid::gpar(fontsize=16)))

grid::grid.newpage()
grid::grid.draw(f2)

#df_out <- ones_df %>% dplyr::select(-c("V1", "V2", "V3","V4","V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "key",
#                                  "gene", "Input.Type","MGI.Gene.Marker.ID", "Name", "Feature.Type" )) %>%
#                       dplyr::rename("HiChIP_Intersect" = "key")

#write_delim(df_out, file = "data/linked_kinetics.bed", delim = "\t")


ones_exprs <- joined_df_exprs %>% dplyr::right_join(ones_df, by = join_by("ens" == "gene_id")) %>% na.omit()

mid3<-mean(log10(ones_exprs$c57_exprs))

ggplot(ones_exprs) + geom_point(size = 2.5, mapping = aes(x = log10(distance), y = k_on.x, color = log10(c57_exprs))) +
  theme_classic(base_size = 18) +
  xlab("log distance (bp)") +
  ylab("Burst Frequency") +
  scale_y_log10() +
scale_color_gradient2(midpoint=mid3, low="blue", mid="grey",
                      high="red", space ="Lab" )





#############GO analysis############
library(clusterProfiler)
library(org.Mm.eg.db)
library(huxtable)

group_8 <- unlist(joined_df %>% dplyr::filter(Freq >= 8) %>% dplyr::select(gene_id))
ego2 <- enrichGO(gene         = group_8,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2_table <- as.data.frame(ego2) %>% dplyr::select(-c("geneID", "Count", "pvalue"))
ego2_table <- ego2_table[1:20,]
dotplot(ego2, font.size = 12)
#huxtable::quick_latex(ego2_table)


######Enhancer By-Pass Analysis#############
coor_split <- function(item){
  strsplit("[:-]", x =  item, perl = T)
}

ones_linked_enhancers <- ones_df %>% dplyr::select(key, Symbol, gene_id)
enh_coors <- ones_linked_enhancers %>% rowwise() %>% dplyr::mutate(chr = coor_split(key)[[1]][1],
                                                                   start = coor_split(key)[[1]][2],
                                                                   end = coor_split(key)[[1]][3])
enh_coors <- enh_coors %>% dplyr::relocate(chr, start, end, Symbol, gene_id) %>% dplyr::select(-c(key))

write_delim(enh_coors, file = "data/singles_enhancer_coors.bed", delim = "\t", col_names = F)

#sorting bed
#bedtools sort -i singles_enhancer_coors.bed > sorted_singles_enh_coors.bed

#getting closest promoter
#bedtools closest -a sorted_singles_enh_coors.bed -b promoters_300.bed -io -t "first" > singles_enh_closest.bed

#reading file back

enh_with_closest <- read.delim("data/singles_enh_closest.bed", sep = "\t", header = F)
enh_with_closest <- enh_with_closest %>% dplyr::rowwise() %>% dplyr::mutate("Close_Promoter" = clipname(V9) )%>%
  dplyr::rename("Linked_Promoter" = "V4")


enh_with_closest$bypass <- FALSE
enh_with_closest$bypass[which(!enh_with_closest$Close_Promoter == enh_with_closest$Linked_Promoter)] <- T

close_promoter_id <- enh_with_closest %>% dplyr::filter(bypass) %>% dplyr::select("Close_Promoter")
write(close_promoter_id$Close_Promoter, "bypassed_id.txt")

bypass_annos <- read.delim("bypassed_annos.txt")

enh_with_closest_anno <- enh_with_closest %>% left_join(bypass_annos, by = join_by("Close_Promoter" == "Input")) %>%
  dplyr::rename("bypass_id" = "Ensembl.ID") %>% dplyr::mutate(bypass_id = ifelse(is.na(bypass_id), V5, bypass_id))

enh_with_closest_anno_kin <- enh_with_closest_anno %>% left_join(kinetics, by = join_by(bypass_id == gene_id)) %>%
  dplyr::mutate(k_on = as.numeric(k_on), k_syn = as.numeric(k_syn), k_off = as.numeric(k_off)) %>%
  distinct(V5, .keep_all = T) %>%
  dplyr::mutate(no_exprs = ifelse(is.na(k_on), T, F)) #%>%
  #dplyr::mutate(k_on = tidyr::replace_na(k_on, 0.006))


table(enh_with_closest_anno_kin$no_exprs, enh_with_closest_anno_kin$bypass, dnn = c("no.exprs", "bypass"))

ggplot(enh_with_closest_anno_kin, aes(y = k_on, x = factor(bypass))) +
  geom_boxplot() +
  scale_y_log10() +
  stat_compare_means()

ggplot(enh_with_closest_anno_kin, aes(x = bypass, fill = no_exprs)) +
  geom_bar(position = "fill")
