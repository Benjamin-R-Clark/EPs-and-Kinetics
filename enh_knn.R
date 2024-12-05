library(Biostrings)
library(bedr)
library(dplyr)
library(org.Mm.eg.db)
library(readxl)
library(FNN)
library(M3C)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(randomForest)

generateKmers <- function(length) {
  nucleotides <- c("A", "C", "T", "G")
  sequences <- expand.grid(rep(list(nucleotides), length))
  result <- apply(sequences, 1, paste, collapse = "")
  return(result)
}

fasta <- readDNAStringSet("data/enh.fasta")
name <- names(fasta)
seq <- paste(fasta)
df.fasta <- data.frame(name, seq)

threemers <- generateKmers(3)


getVec <- function (seq, kmers){
  outtable <- c()
  for(k in kmers){
    match <- gregexpr(k, df.fasta$seq[1], ignore.case = T)
    outtable[k] <- length(match[[1]]) / (nchar(seq) - nchar(k) + 1)
  }
  return(outtable)
}

freqs <- t(sapply(df.fasta$seq, getVec, kmers = threemers))
row.names(freqs) <- df.fasta$name
freqs <- data.frame(freqs)
en_codes <- str_extract(rownames(freqs), "e[0-9]*")
freqs["en_code"] <- en_codes
unique.freqs.ids <- freqs %>% distinct(en_code)
unique.freqs <- freqs[rownames(unique.freqs.ids),]



### load in enhancer promoter contacts
enh_prom <- read.delim("data/enh_promoter_contacts.bed", header = F)

library(stringr)
ennames <- df.fasta$name

en_codes <- str_extract(ennames, "e[0-9]*")

df.fasta["enhancer"] <- en_codes
unique.df.fasta <- dplyr::distinct(df.fasta)

enh_prom <- enh_prom %>% dplyr::select(V12, V17)

linked.df.fasta <- df.fasta %>% left_join(enh_prom, by = join_by("enhancer" == V12))

## annotating genes with ensemble Ids
annos <- AnnotationDbi::select(org.Mm.eg.db,
                               keys = linked.df.fasta$V17,
                               columns = c("ENSEMBL", "GENENAME"),
                               keytype = c("SYMBOL"))

annos <- distinct(annos)

linked.df.fasta <- linked.df.fasta %>% left_join(y = annos, by = join_by(V17 == SYMBOL))
unique.links <- linked.df.fasta %>% distinct() %>% na.omit()


#loading in kinetics
kinetics <- read_excel("data/journal.pcbi.1008772.s006.xlsx", sheet = "C57")
colnames(kinetics) <- kinetics[1,]
kinetics <- kinetics[-1,]

kinetics <- kinetics %>% dplyr::select(gene_id, k_on, k_off, k_syn)

#join tables
unique.links.kinetics <- unique.links %>% left_join(kinetics, by = join_by("ENSEMBL" == "gene_id")) %>% na.omit()
unique.links.kinetics.freqs <- unique.links.kinetics %>%
  right_join(unique.freqs, by = join_by("enhancer" == "en_code")) %>%
  mutate(id = row_number()) %>% na.omit() %>% mutate("k_on" = as.numeric(k_on), "k_off" = as.numeric(k_off), "k_syn" = as.numeric(k_syn))

three_train <- unique.links.kinetics.freqs %>% sample_frac(.70)
three_test <- anti_join(unique.links.kinetics.freqs, three_train, by = 'id')

#  knn3 <- knn.reg(train = three_train[,c(11:(ncol(three_train) - 1))],
#                 y = as.numeric(three_train$k_on),
#                 test = three_test[,c(11:(ncol(three_test) - 1))],
#                 k=10)
#
# mse.knn3.datapred <- knn.reg(train = three_train[,c(11:(ncol(three_train) - 1))],
#                              y = three_train$k_on,
#                              test = three_train[,c(11:(ncol(three_train) - 1))],
#                              k = 10)
#
# mse.knn3 <- mean((mse.knn3.datapred$pred - three_train$k_on) ^ 2)
# r2.knn3 <- 1 - mse.knn3/(var(three_train$k_on)* nrow(three_train) - 1 / nrow(three_train))
tsne.table.all <- as.data.frame((distinct(unique.links.kinetics.freqs[, c(7,11:ncol(unique.links.kinetics.freqs) - 1)])))


#features <- featurefilter(tsne.table)
#tplot <- tsne(tsne.table.u[2:nrow(tsne.table.u),], perplex = 5, labels = as.vector(scale(as.numeric(tsne.table.u[1,]))), scale = 2, controlscale = TRUE)
#tplot
tsne.table <- unique(tsne.table.all[,1:ncol(tsne.table.all)])
Z_scaled_kon <- as.vector(scale(tsne.table[,1]))
tsne <- Rtsne(tsne.table, perplexity = 17, check_duplicates = T )
tsne.data <- as.data.frame(tsne$Y)
ggplot(tsne.data) + geom_point(mapping = aes(x = V2, y = V1, color = Z_scaled_kon))+ theme_light()


#Random Forest Regression

three.fit <-  randomForest(k_on ~ ., data = unique.links.kinetics.freqs[,c(7,10:(ncol(three_train) - 1))],
                           importance = T, ntree = 500)
three.fit

#very poor result with threemers

fourmers <- generateKmers(4)



fourfreqs <- t(sapply(df.fasta$seq, getVec, kmers = fourmers))
row.names(fourfreqs) <- df.fasta$name
fourfreqs <- data.frame(fourfreqs)
en_codes <- str_extract(rownames(fourfreqs), "e[0-9]*")
fourfreqs["en_code"] <- en_codes
unique.fourfreqs.ids <- freqs %>% distinct(en_code)
unique.fourfreqs <- fourfreqs[rownames(unique.fourfreqs.ids),]



unique.links.kinetics.fourfreqs <- unique.links.kinetics %>%
  right_join(unique.fourfreqs, by = join_by("enhancer" == "en_code")) %>%
  mutate(id = row_number()) %>% na.omit() %>% mutate("k_on" = as.numeric(k_on), "k_off" = as.numeric(k_off), "k_syn" = as.numeric(k_syn))



four.fit <-  randomForest(k_on ~ ., data = unique.links.kinetics.fourfreqs[,c(7,10:(ncol(unique.links.kinetics.fourfreqs) - 1))],
                           importance = T, ntree = 500, mtry = 200)
four.fit