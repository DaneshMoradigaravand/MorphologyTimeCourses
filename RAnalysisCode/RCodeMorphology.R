#compute difference in time series example section manuscript
library(tidyverse)
library(sqldf)
library(Hmisc) 
library(pheatmap)
library(tidyverse)

growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints_summarised.csv") %>% 
        group_by(label, timepoint) %>% 
        summarise(morphology_average=median(morphology.score.fixed.circles)) %>%
        pivot_wider(names_from =timepoint , values_from = morphology_average) %>% 
  filter(`44h` >0) %>%
  arrange(`44h`)
growth$rank<-seq(dim(growth)[1])

diff_frame_tmp <- growth[-1,4] - growth[-nrow(growth),4]
growth$diff_frame_39h<- c(NA, diff_frame_tmp$`39h` )

diff_frame_tmp <- growth[-1,5] - growth[-nrow(growth),5]
growth$diff_frame_44h<- c(NA, diff_frame_tmp$`44h` )

growth_shortened<-growth[which(growth$diff_frame_44h==0),]
nth_highest<-7
rank_holder<- which( growth_shortened$diff_frame_39h==sort(growth_shortened$diff_frame_39h, TRUE)[nth_highest])
rbind(growth_shortened[rank_holder,],growth[growth_shortened$rank[rank_holder]-1,])

#Two distributions 
growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints_summarised.csv") %>% 
  group_by(label, timepoint) %>% 
  summarise(morphology_average=median(morphology.score.fixed.circles)) %>%
  filter(timepoint=="44h")

biofilm_genes<-c("PA14_09480", "PA14_39910","PA14_09470","PA14_42760","PA14_13250", "PA14_13240", "PA14_13230", "PA14_24900", "PA14_13850", "PA14_51410", "PA14_51430", "PA14_24480", "PA14_24490", "PA14_24510", "PA14_24530", "PA14_24550", "PA14_24560",
                 "PA14_16500",
                 "PA14_64050",
                 "PA14_12810",
                 "PA14_56790",
                 "PA14_50060",
                 "PA14_02110",
                 "PA14_42220",
                 "PA14_66320",
                 "PA14_21190",
                 "PA14_07730",
                 "PA14_31290",
                 "PA14_45940",
                 "PA14_52260",
                 "PA14_54430",
                 "PA14_09150",
                 "PA14_61200","PA14_12530")

genes<-c("phzA", "phzE","phzB","aroC", "moaE", "moaD", "moaC", "moaB", "moaA", "pqsC", "pqsA", "pelA","pelB", "pelD","pelE", "pelF","pelG",
         "wspR",
         "gcbA",
         "rocR",
         "bifA",
         "roeA",
         "siaD",
         "mucR",
         "dipA",
         "nbdA",
         "rsmA",
         "lecA",
         "lasI",
         "gacS",
         "algU",
         "katA",
         "cdrA","top")

gene_id<-paste0(biofilm_genes,"_",genes)

shortened<-growth[match(biofilm_genes, growth$label),] %>% drop_na()

growth %>%
  ggplot( aes(x=growth$morphology_average)) +
 # geom_density(fill="blue") +
  geom_histogram(fill="blue")+
  ggplot2::annotate("text", x = shortened$morphology_average, y=1500, label = gene_id[match(shortened$label,biofilm_genes )] , angle = 90)+
  geom_vline(xintercept = shortened$morphology_average, linetype="dotted", size = 0.3)+
  geom_vline(xintercept = c(quantile(growth$morphology_average, 0.95), quantile(growth$morphology_average, 0.05)), linetype="solid", col="red", size = 0.3)+
  #geom_hline(yintercept = 1, linetype="solid", col="black", size = 0.3)+
  theme_bw()+ 
  xlim(range(0,300))+
  ylim(range(0,2000))+
  ylab("Frequency")+
  xlab("Morphology at 44h")+
  theme(axis.text.x = element_text( size=14, angle=45,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))


bonferroni_correction<-function(x, number_of_test){
  tmp<-length(which(x>growth$morphology_average))/length(growth$morphology_average)
  tmp<-ifelse(tmp>0.5, (1-tmp), tmp)
  ifelse(tmp<0.05/number_of_test, "sig", "nonsig")
}
bonferroni_correction(shortened$morphology_average[2],1)
shortened$sig<-sapply(shortened$morphology_average, function(x) bonferroni_correction(x,1))
shortened$full_label<-gene_id[match(shortened$label, biofilm_genes)]


library(edgeR)
library(gmodels)
library(tidyverse)
growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv") %>% 
  select(time_point, gene, value) %>% 
  group_by(gene, time_point) %>% 
  summarise(color_average=median(value)) %>%
  filter(time_point=="44h")


shortened<-growth[match(biofilm_genes, growth$gene),] %>% drop_na()

growth %>%
  ggplot( aes(x=growth$color_average)) +
  geom_histogram(fill="blue") +
  ggplot2::annotate("text", x = shortened$color_average, y=3000, label = gene_id[match(shortened$gene,biofilm_genes )] , angle = 90)+
  geom_vline(xintercept = shortened$color_average, linetype="dotted", size = 0.3)+
  geom_vline(xintercept = c(quantile(growth$color_average, 0.95), quantile(growth$color_average, 0.05)), linetype="solid", col="red", size = 0.3)+
  #geom_hline(yintercept = 1, linetype="solid", col="black", size = 0.3)+
  theme_bw()+ 
  xlim(range(0,9000000))+
  ylim(range(0,4000))+
  ylab("Frequency")+
  xlab("Morphology at 44h")+
  theme(axis.text.x = element_text( size=14, angle=45,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))

####################################
#Build the placeholder
library(RSQL) 
library(RSQLite)
library(tidyverse)
library(DBI)

growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints.csv") 
con <- dbConnect(drv = RSQLite::SQLite(), dbname = ":memory:")
dbWriteTable(conn = con, name = "growth", value = growth)

tmp<-dbGetQuery(con,"SELECT [PA14.ID], timepoint, AVG([morphology.score.fixed.circles]) AS average_morphology, VARIANCE([morphology.score.fixed.circles]) AS variance_morphology, VARIANCE([colony.size]) AS variance_colony_size ,  AVG([colony.size]) AS average_colony_size, COUNT(*) AS replicas
           FROM growth
           GROUP BY [PA14.ID],timepoint 
           ")
dbWriteTable(conn = con, name = "tmp", value = tmp, overwrite=T)

tmp<-dbGetQuery(con,"SELECT *, variance_morphology/replicas AS average_variance_morphology_per_replica,  variance_colony_size/replicas AS average_variance_colony_size_per_replica 
           FROM tmp
           ORDER BY average_variance_morphology_per_replica  DESC
           ")
dbWriteTable(conn = con, name = "tmp", value = tmp, overwrite=T)

df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv") 
dbWriteTable(conn = con, name = "df_colour", value = df_colour)
tmp1<-dbGetQuery(con,"SELECT gene, time_point, AVG([value]) AS average_colony_colour, VARIANCE([value]) AS variance_colony_colour, COUNT(*) AS replicas
           FROM df_colour
           GROUP BY gene,time_point 
           ")
dbWriteTable(conn = con, name = "tmp1", value = tmp1, overwrite=T)
tmp1<-dbGetQuery(con,"SELECT *, variance_colony_colour/replicas AS average_variance_colony_colour_per_replica
           FROM tmp1
           ORDER BY average_variance_colony_colour_per_replica  DESC
           ")
dbWriteTable(conn = con, name = "tmp1", value = tmp1, overwrite=T)

dbExecute(con,"DROP TABLE temp_table")
dbExecute(con,"CREATE TABLE temp_table AS
               SELECT *, [PA14.ID] As gene  
               FROM tmp
          ")
dbExecute(con,"ALTER TABLE temp_table
           DROP COLUMN [PA14.ID]")
dbExecute(con,"SELECT * FROM temp_table
            INNER JOIN tmp1 
            ON temp_table.gene= tmp1.gene
           ")
tmp3<-dbGetQuery(con,"SELECT * FROM temp_table")
#write_csv(tmp3, "/Users/moradigd/Documents/chemicalgenomic/summary_table_time_points.csv" )



####################################

load("/Users/moradigd/Documents/Plasmid/BioInfProject/Data/abundetc.RData", ex <- new.env())
ls() #returns a list of all the objects you just loaded (and anything else in your environment)


#Correlation
library(tidyverse)
library(corrplot)
inp<-read.csv("/Users/moradigd/Documents/Plasmid/BioInfProject/Data/aptable.csv") %>% select(colnames(.)[3:24])
#inp[!is.na(inp)]<-1
#inp[is.na(inp)]<-0

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  for(i in seq(n)){
    for(j in seq(n)){
      tmp <- cor.test(as.numeric(as.character(gsub(" ","",mat[, j]))) , as.numeric(as.character(gsub(" ","",mat[, i]))) ,use="pairwise.complete.obs",method ="spearman")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  diag(p.mat) <- 0
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

cor_df<-cor(inp, use = "pairwise.complete.obs",method ="spearman")

ordered_tags<-c(
  "DK062017_1A_PKJK", "DK062017_1B_PKJK", "DK062017_2_PKJK" ,  
  "DK062017_4_PKJK" ,
  "DK062017_6_PKJK" ,  "DK062017_7_PKJK"  ,
  "DK022018_1A_PKJK", "DK022018_1B_PKJK", "DK022018_2_PKJK",
  "DK022018_4_PKJK",
  "UK082017_1_PKJK" ,  "UK082017_2_PKJK",   
  "UK082017_4_PKJK" ,  "UK082017_6_PKJK", "UK082017_7_PKJK" ,
  "DK062017_6_RP" ,   "UK082017_6_RP",
  "DK062017_6_PB"   ,     "UK082017_6_PB" 
)

cor_df_ordered<-cor_df[match(ordered_tags, colnames(cor_df)),match(ordered_tags, colnames(cor_df))] 
inp_ordered<-inp[, match(ordered_tags, colnames(inp))]



summerised_ordered_tags<-c("1A", "1B", "2" ,  "4" ,"6" ,  "7"  ,"1A", "1B", "2","4",
                           "1" ,  "2",   "4" ,  "6", "7" ,"6" ,   "6", "6"   ,     "6" 
)
colnames(cor_df_ordered)<-summerised_ordered_tags
row.names(cor_df_ordered)<-summerised_ordered_tags
#[1] "DK062017_1A_PKJK"  "DK062017_1B_PKJK"  "DK062017_2_PKJK"   "DK062017_4_PKJK"  
#[5] "DK062017_4.2_PKJK" "DK062017_6_PKJK"   "DK062017_7_PKJK"   "DK062017_6_RP"    
#[9] "DK062017_6_PB"     "UK082017_1_PKJK"   "UK082017_2_PKJK"   "UK082017_2.2_PKJK"
#[13] "UK082017_4_PKJK"   "UK082017_6_PKJK"   "UK082017_6_PB"     "UK082017_6_RP"    
#[17] "UK082017_7_PKJK"   "UK082017_7.2_PKJK" "DK022018_1A_PKJK"  "DK022018_1B_PKJK" 
#[21] "DK022018_2_PKJK"   "DK022018_4_PKJK"  

corrplot(cor_df_ordered, method = "circle", type="upper",p.mat = cor.mtest(inp_ordered), sig.level = 0.01,tl.col="black", tl.cex = 1)
corrplot(cor_df_ordered, method = "circle", type="upper",p.mat = cor.mtest(inp_ordered), sig.level = 0.01,tl.col="black", tl.cex = 1.5)

mean(cor_df_ordered[cor_df_ordered!=1])
range(cor_df_ordered[cor_df_ordered<1])
max(cor_df_ordered[cor_df_ordered!=1])
?cor.test

mean(cor_df_ordered[lower.tri(cor_df_ordered, diag = FALSE)])
range(cor_df_ordered[lower.tri(cor_df_ordered, diag = FALSE)])
#PCA
library(ggbiplot)
plasmidpca <- prcomp(cor(cor_df_ordered), center = TRUE,scale. = TRUE)
ggbiplot(plasmidpca,labels = as.character(sapply(colnames(inp_ordered), function(x) str_split(x,"_")[[1]][3])), labels.size = 2.5, groups = as.character(sapply(colnames(inp_ordered), function(x) str_split(x,"_")[[1]][1])),  var.axes=FALSE )+ 
  geom_point(aes(colour=as.character(sapply(colnames(inp_ordered), function(x) str_split(x,"_")[[1]][1]))), size = 10,pch=22) +
  theme_bw() +
  theme(axis.text.x = element_text( size=11, angle=0,hjust = 1),
        axis.text.y = element_text( size=11, hjust = 1),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
        #,legend.position="none"
  )

exp<-as.character(sapply(colnames(inp_ordered), function(x) str_split(x,"_")[[1]][1]))
ggbiplot(plasmidpca, labels.size = 2.5, groups =  as.character(sapply(exp, function(x) substr(x,1,2)))   ,  var.axes=FALSE )+ 
  geom_point(aes(colour=as.character(sapply(exp, function(x) substr(x,1,2)))), size = 10,pch=22) +
  theme_bw() +
  theme(axis.text.x = element_text( size=14, angle=0,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
        #,legend.position="none"
  )  

exp<-as.character(sapply(colnames(inp_ordered), function(x) str_split(x,"_")[[1]][1]))
ggbiplot(plasmidpca, labels.size = 2.5, groups =  as.character(sapply(exp, function(x) substr(x,3,8)))   ,  var.axes=FALSE )+ 
  geom_point(aes(colour=as.character(sapply(exp, function(x) substr(x,3,8)))), size = 10,pch=22) +
  theme_bw() +
  theme(axis.text.x = element_text( size=14, angle=0,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
        #,legend.position="none"
  ) 

exp<-as.character(sapply(colnames(inp_ordered), function(x) str_split(x,"_")[[1]][3]))
ggbiplot(plasmidpca, labels.size = 2.5, groups =  exp  ,  var.axes=FALSE )+ 
  geom_point(aes(colour=exp), size = 10,pch=22) +
  theme_bw() +
  theme(axis.text.x = element_text( size=14, angle=0,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
        #,legend.position="none"
  ) 

cor.test(inp$DK062017_4_PKJK, inp$DK062017_4.2_PKJK, method= "spearman")


#unceratinty 
inp<-read.csv("/Users/moradigd/Documents/Plasmid/BioInfProject/Data/aptable.csv") %>% select(colnames(.)[3:24])
#inp[!is.na(inp)]<-1
inp[is.na(inp)]<-0
data.frame(table(inp$DK062017_1A_PKJK, inp$DK062017_1B_PKJK))
datatable(table(inp$DK062017_1A_PKJK, inp$DK062017_1B_PKJK))
datatable(table(inp$UK082017_7_PKJK, inp$UK082017_7.2_PKJK))


#Permissiveness
inp<-read.csv("/Users/moradigd/Documents/Plasmid/BioInfProject/Data/aptable.csv") %>% select(colnames(.)[3:24])
inp[is.na(inp)]<-0
ggplot(inp, aes(x=log(DK062017_1A_PKJK), y=log(DK062017_1B_PKJK)))+
  geom_point()+
  ylab("")+
  geom_smooth(method=lm, se=TRUE)+
  theme_bw()

ggplot(inp, aes(x=log(DK062017_4_PKJK), y=log(DK062017_4.2_PKJK)))+
  geom_point()+
  ylab("")+
  geom_smooth(method=lm, se=TRUE)+
  theme_bw()

ggplot(inp, aes(x=log(UK082017_2_PKJK), y=log(UK082017_2.2_PKJK)))+
  geom_point()+
  ylab("")+
  geom_smooth(method=lm, se=TRUE)+
  theme_bw()





#Predictioion
library(tidyverse)
inp<-read_csv("/Users/moradigd/Documents/Plasmid/BioInfProject/Data/aptable.csv") %>% select(colnames(.)[1:25]) %>% select(-colnames(.)[1])
#inp[2:dim(inp)[2]]<-ifelse(is.na(inp[2:dim(inp)[2]]), 0, 1)
inp[2:dim(inp)[2]]<-ifelse(is.na(inp[2:dim(inp)[2]]), 0, inp[2:dim(inp)[2]])
dim(inp)

length_seq<-as.numeric(as.character(sapply(inp$seq, nchar)))   
inp<-inp[which(length_seq>380),]

out_seq<-c()
for(i in seq(length(inp$seq))){out_seq<-c(out_seq, paste0(">seq",i), inp$seq[i])}
#writeLines(out_seq, "/Users/moradigd/Documents/Plasmid/BioInfProject/inp_seq")

#inp_sum<-apply(inp[,which(grepl("_PKJK",colnames(inp) ) & grepl("_7",colnames(inp) ))], 1 , sum) 
inp_sum<-apply(inp[,which(grepl("_PKJK",colnames(inp) ))], 1 , sum) 
inp_sum_PKJK<-ifelse(inp_sum==0, 0, 1)
table(inp_sum_PKJK)

inp_sum<-apply(inp[,which(grepl("_RP",colnames(inp) ))], 1 , sum) 
inp_sum_RP<-ifelse(inp_sum==0, 0, 1)
table(inp_sum_RP)

inp_sum<-apply(inp[,which(grepl("_PB",colnames(inp) ))], 1 , sum) 
inp_sum_PB<-ifelse(inp_sum==0, 0, 1)
table(inp_sum_PB)


vcf<-read_tsv("/Users/moradigd/Documents/Plasmid/BioInfProject/vcf_out.vcf") %>% select(colnames(.)[10:dim(.)[2]])
vcf_t<-data.frame(t(vcf))
dim(vcf_t)
vcf_t$label<-factor(inp_sum_PKJK)

library(randomForest)
require(caTools)
sample = sample.split(vcf_t$label, SplitRatio = .5)
train = subset(vcf_t, sample == TRUE)
test  = subset(vcf_t, sample == FALSE)
dim(train)
dim(test)
rf <- randomForest(
  label~ .,
  data=train,
  mtry=9,
  ntree=1500
)
pred = predict(rf, newdata=test[-483])
table(test[,483], pred)
rf

#Grid search 
control <- trainControl(method="repeatedcv", number=2, repeats=2, search="grid")
seed <- 7
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(8,9,10))
metric <- "Accuracy"
rf_gridsearch <- train(label~., data=train, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)



#Multiple tuning
customRF <- list(type = "Classification",
                 library = "randomForest",
                 loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"),
                                  class = rep("numeric", 2),
                                  label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
  randomForest(x, y,
               mtry = param$mtry,
               ntree=param$ntree)
}




#Predict label
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

#Predict prob
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes


control <- trainControl(method="repeatedcv", 
                        number=2, 
                        repeats=1,
                        allowParallel = TRUE)

tunegrid <- expand.grid(.mtry=c(12,15,18),.ntree=c(1000))

set.seed(123)

custom <- train(label~., data=train, 
                method=customRF, 
                metric=metric, 
                tuneGrid=tunegrid, 
                trControl=control)

summary(custom)
plot(custom)


# 1 ROC curve, mock vs non mock
library(pROC)
library(randomForest)
require(caTools)

rf <- randomForest(
  label~ .,
  data=train,
  mtry=12,
  ntree=1000
)
pred = predict(rf, newdata=test[-483])
table(test[,483], pred)
rf

predictions <- as.data.frame(predict(rf, test, type = "prob"))
predictions$predict<-names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- as.character(test$label)

roc.one <- roc(ifelse(predictions$observed=="1", "1", "0"), as.numeric(predictions$`1`))
roc.zero <- roc(ifelse(predictions$observed=="0", "0", "1"), as.numeric(predictions$`0`))
plot(roc.zero, col = "red")
lines(roc.one, col = "blue")
auc(roc.one)

#vcf<-read_tsv("/Users/moradigd/Documents/Plasmid/BioInfProject/vcf_out.vcf") %>% select(colnames(.)[10:dim(.)[2]])
#vcf_t<-data.frame(t(vcf))
#dim(vcf_t)
#vcf_t$label<-factor(inp_sum_PB)

#sample = sample.split(vcf_t$label, SplitRatio = .5)
#train = subset(vcf_t, sample == TRUE)
#test  = subset(vcf_t, sample == FALSE)
library(caret) 
randomized<-c()
for(j in seq(50)){
  print(j)
  train_sampled<-train
  train_sampled$label<-sample(train_sampled$label)
  
  dim(train)
  dim(test)
  rf <- randomForest(
    label~ .,
    data=train_sampled
  )
  pred = predict(rf, newdata=test[-483])
  res<-confusionMatrix(table(test[,483], pred))
  randomized<-c(randomized, res$byClass[11])
  print(res$byClass[11])
}

rf <- randomForest(
  label~ .,
  data=train
)
pred = predict(rf, newdata=test[-483])
res<-confusionMatrix(table(test[,483], pred))

hist(randomized, breaks = 100, col = "black", xlab = "", main="", xlim = range(0.4, 0.8))
abline(v = res$byClass[11], col="red", lwd=3, lty=2,las=2)
box()



n = nrow(vcf_t)
train.index = sample(n,floor(0.75*n))
train.data = as.matrix(vcf_t[train.index,])
train.label = as.factor(inp_sum_PB[train.index])
test.data = as.matrix(vcf_t[-train.index,])
test.label = as.factor(inp_sum_PB[-train.index])

library(xgboost)
xgb.train = xgb.DMatrix(data=train.data,label=train.label)
xgb.test = xgb.DMatrix(data=test.data,label=test.label)

num_class = length(levels(as.factor(inp_sum_PB)))
params = list(
  booster="gbtree",
  eta=0.001,
  max_depth=5,
  gamma=3,
  subsample=0.75,
  colsample_bytree=1,
  objective="multi:softprob",
  eval_metric="mlogloss",
  num_class=num_class
)


xgb.fit=xgb.train(
  params=params,
  data=xgb.train,
  nrounds=100,
  nthreads=1,
  early_stopping_rounds=10,
  watchlist=list(val1=xgb.train,val2=xgb.test),
  verbose=0
)

# xboost
grid_default <- expand.grid(
  nrounds = 1500,
  max_depth = 3,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.8
)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_base <- caret::train(
  x = as.matrix(train[,-483]),
  y = as.matrix(train[,483]),
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE
)
xgb_base
pred = predict(xgb_base, newdata=test[-483])
table(test[,483], pred)
res<-confusionMatrix(table(test[,483], pred))
res
predictions <- as.data.frame(predict(xgb_base, test, type = "prob"))
predictions$predict<-names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- as.character(test$label)

roc.one <- roc(ifelse(predictions$observed=="1", "1", "0"), as.numeric(predictions$`1`))
roc.zero <- roc(ifelse(predictions$observed=="0", "0", "1"), as.numeric(predictions$`0`))
lines(roc.zero, col = "red")
lines(roc.one, col = "blue")
auc(roc.one)


pred = predict(rf, newdata=test[-483])
table(test[,483], pred)


# note to start nrounds from 200, as smaller learning rates result in errors so
# big with lower starting points that they'll mess the scales
nrounds <- 1000
tune_grid <- expand.grid(
  nrounds = seq(from = 200, to = nrounds, by = 50),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(2, 3, 4, 5, 6),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

tune_control <- caret::trainControl(
  method = "cv", # cross-validation
  number = 3, # with n folds 
  #index = createFolds(tr_treated$Id_clean), # fix the folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_tune <- caret::train(
  x = as.matrix(train[,-483]),
  y = as.matrix(train[,483]),
  trControl = tune_control,
  tuneGrid = tune_grid,
  method = "xgbTree",
  verbose = TRUE
)

# helper function for the plots
tuneplot <- function(x, probs = .90) {
  ggplot(x) +
    coord_cartesian(ylim = c(quantile(x$results$RMSE, probs = probs), min(x$results$RMSE))) +
    theme_bw()
}

tuneplot(xgb_tune)


#IToL 
library(rncl)
library(phangorn)
library(ape)

dna<-read.dna("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/seq_PKJK_aligned.fasta",format= "fasta")
distdna<-dist.dna(dna,model ="N",pairwise.deletion = TRUE,as.matrix = TRUE)
#rnd<-round(distdna*117889 )
write.csv(distdna, file="/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/PKJK_dist-Core.csv",quote=FALSE, row.names = FALSE)
tre<-nj(distdna)
tre<-midpoint(tre)
#write.nexus(tre, file="/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/seq_PKJK_aligned.fasta.tre")
write.tree(tre, file="/Users/moradigd/Documents/Plasmid/seq_PKJK_aligned.fasta.tre")



library(tidyverse)
library(tidymodels)
library(AmesHousing)
library(gt)
# function copied from here:
# https://github.com/rstudio/gt/issues/613#issuecomment-772072490 
# (simpler solution should be implemented in future versions of {gt})
fmt_if_number <- function(..., digits = 2) {
  input <- c(...)
  fmt <- paste0("%.", digits, "f")
  if (is.numeric(input))   return(sprintf(fmt, input))
  return(input)
}

ames <- make_ames() %>% 
  mutate(Years_Old = Year_Sold - Year_Built,
         Years_Old = ifelse(Years_Old < 0, 0, Years_Old))
set.seed(4595)
data_split <- initial_split(ames, strata = "Sale_Price", prop = 0.75)
ames_train <- training(data_split)
ames_holdout  <- testing(data_split) 
dim(ames_train)
dim(ames_holdout)



#RF models require comparably less pre-processing to linear models
rf_recipe <- 
  recipe(
    Sale_Price ~ Lot_Area + Neighborhood  + Years_Old + Gr_Liv_Area + Overall_Qual + Total_Bsmt_SF + Garage_Area, 
    data = ames_train
  ) %>%
  step_log(Sale_Price, base = 10) %>%
  step_other(Neighborhood, Overall_Qual, threshold = 50) %>% 
  step_novel(Neighborhood, Overall_Qual) %>% 
  step_dummy(Neighborhood, Overall_Qual) 


rf_mod <- rand_forest() %>%
  set_engine("ranger", importance = "impurity", seed = 63233, quantreg = TRUE) %>%
  set_mode("regression")
set.seed(63233)
rf_wf <- workflows::workflow() %>% 
  add_model(rf_mod) %>% 
  add_recipe(rf_recipe) %>% 
  fit(ames_train)



preds_bind <- function(data_fit, lower = 0.05, upper = 0.95){
  predict(
    rf_wf$fit$fit$fit, 
    workflows::pull_workflow_prepped_recipe(rf_wf) %>% bake(data_fit),
    type = "quantiles",
    quantiles = c(lower, upper, 0.50)
  ) %>% 
    with(predictions) %>% 
    as_tibble() %>% 
    set_names(paste0(".pred", c("_lower", "_upper",  ""))) %>% 
    mutate(across(contains(".pred"), ~10^.x)) %>% 
    bind_cols(data_fit) %>% 
    select(contains(".pred"), Sale_Price, Lot_Area, Neighborhood, Years_Old, Gr_Liv_Area, Overall_Qual, Total_Bsmt_SF, Garage_Area)
}
rf_preds_test <- preds_bind(ames_holdout)

set.seed(1234)
rf_preds_test %>% 
  mutate(pred_interval = ggplot2::cut_number(Sale_Price, 10)) %>% 
  group_by(pred_interval) %>% 
  sample_n(2) %>% 
  ggplot(aes(x = .pred))+
  geom_point(aes(y = .pred, color = "prediction interval"))+
  geom_errorbar(aes(ymin = .pred_lower, ymax = .pred_upper, color = "prediction interval"))+
  geom_point(aes(y = Sale_Price, color = "actuals"))+
  scale_x_log10(labels = scales::dollar)+
  scale_y_log10(labels = scales::dollar)+
  labs(title = "90% Prediction intervals on a holdout dataset",
       subtitle = "Random Forest Model",
       y = "Sale_Price prediction intervals and actuals")+
  theme_bw()+
  coord_fixed()


#Correlation test 
library(tidyverse)
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 
#write_csv(df, "/Users/moradigd/Documents/Plasmid/input_agg.csv")

df<-df[!is.na(df$PB),]
output<-c()
for(i in seq(dim(df)[1])){output<-c(output, paste0(">seq",i), df$seq[i])}
writeLines(output,"/Users/moradigd/Documents/seq_PB.fasta")
#380 
mySequences <- readAAStringSet("/Users/moradigd/Documents/seq_PB.fasta")
mySequences
myFirstAlignment <- msa(mySequences)

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}
alignment2Fasta(myFirstAlignment, '/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/seq_PB_aligned.fasta')

library(rhierbaps)

fasta.file.name <- system.file("extdata", "seqs.fa", package = "rhierbaps")
snp.matrix <- load_fasta("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/seq_PB_aligned.fasta")
hb.results <- hierBAPS(snp.matrix, max.depth = 50, n.pops = 50, quiet = FALSE)
hb.results
results_baps<-hb.results$partition.df

results_baps$`level 49`==results_baps$`level 50`
which()
all(results_baps$`level 5`==results_baps$`level 6`)
#HierBAPS based prediction 
baps<-results_baps[,1:7]
for(i in seq(2, dim(baps)[2])){
  baps[,i]<-paste0("level_",i,"_",baps[,i] )
}
write_csv(baps, "/Users/moradigd/Documents/Plasmid/PB_BAPS.csv")

#importace analysis for BAPS 
results<-c()
for(i in seq(0,9)){
  file_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/PB_baps_importance_",i,"_PB")) %>% arrange(desc(importance))
  file_tmp$label<-paste0("scale_",i)
  file_tmp$rank<-seq(dim(file_tmp)[1])
  results<-rbind(results, file_tmp)
}

importance_average<-results %>% 
  group_by(tag) %>%
  summarise(importance_average=mean(importance), rank_average=mean(rank)) 

importance_average<-results %>% 
  group_by(tag) %>%
  count(tag) %>% inner_join(importance_average, by=c("tag"="tag")) %>%
  arrange(desc(n))

#Families 
taxa<-read_csv("/Users/moradigd/Documents/Plasmid/TaxaFamily.csv")
table(taxa$Genus)

sequences<-readLines("/Users/moradigd/Documents/seq_PB.fasta")
baps<-read_csv("/Users/moradigd/Documents/Plasmid/PB_BAPS.csv")


column_id<-sapply(importance_average$tag, function(x) as.numeric(strsplit(x,"_")[[1]][2])) 
names(column_id)

for(k in seq(1)){
  baps_tmp<- baps[, c(1,column_id[k])]
  seqs_tmp<-sequences[match(paste0(">seq",which(baps_tmp[,2]==names(column_id)[k])),sequences)+1]
  print(taxa$Genus[match(seqs_tmp, taxa$seq)])
}

#Results from prediction 
library(gmodels)
PB_results<-read_csv("/Users/moradigd/Documents/Plasmid/PB_baps_results_screening_agg_scale.csv")
colnames(PB_results)<-c("index","correlation","pvalue","plasmid","model","run")

RP_results<-read_csv("/Users/moradigd/Documents/Plasmid/RP_baps_results_screening_agg_scale.csv")
colnames(RP_results)<-c("index","correlation","pvalue","plasmid","model","run")

PKJK_results<-read_csv("/Users/moradigd/Documents/Plasmid/PKJK_baps_results_screening_agg_scale.csv")
colnames(PKJK_results)<-c("index","correlation","pvalue","plasmid","model","run")

results<-rbind(PB_results,RP_results, PKJK_results) %>% 
  drop_na() %>%
  group_by(plasmid, model) %>%
  summarise(mean_correlation=mean(correlation) , lowCI = ci(correlation)[2],
            hiCI = ci(correlation)[3], sd_correlation= ci(correlation)[4]) %>%
  mutate(lowCI = replace(lowCI, which(lowCI<0), 0))

#
read_csv("")







cor.test(corr$DK062017_4_PKJK,corr$`DK062017_4-2_PKJK`)

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) length(which(!is.na(x))))

tmp<-corr[, which(grepl("_PKJK",colnames(corr)))]


hist(as.matrix(tmp[788,]))

dim(corr[, which(grepl("_PKJK",colnames(corr)))])

#alpha beta diveristy 
#speices frequency 

#prediction results

#interpretable machine learning what speices are overrepresentated? 

#how much population strcture matters? Done 

#Clusters Done

fasta.file.name <- system.file("extdata", "seqs.fa", package = "rhierbaps")
snp.matrix <- load_fasta(fasta.file.name)


#Sample diversity 
library(vegan)
data(BCI)
H <- diversity(BCI)
simp <- diversity(BCI, "simpson")
invsimp <- diversity(BCI, "inv")


#sample augmentation 
#Include ancestral values



#Interval prediction 

library(tidyverse)
library(tidymodels)
library(AmesHousing)
library(gt)
# function copied from here:
# https://github.com/rstudio/gt/issues/613#issuecomment-772072490 
# (simpler solution should be implemented in future versions of {gt})
fmt_if_number <- function(..., digits = 2) {
  input <- c(...)
  fmt <- paste0("%.", digits, "f")
  if (is.numeric(input))   return(sprintf(fmt, input))
  return(input)
}


ames <- make_ames() %>% 
  mutate(Years_Old = Year_Sold - Year_Built,
         Years_Old = ifelse(Years_Old < 0, 0, Years_Old))
set.seed(4595)

ames <- read.csv("/Users/moradigd/Documents/Plasmid/interval.csv")
#ames <-ames[1000:1026] 

ames <- apply(ames, 2, function(x) scale(x))
#ames  <- ames %>% mutate_at(., ~(scale(.) %>% as.vector))
ames <- data.frame(ames)
#ames_1 <- scale(ames)

data_split <- initial_split(ames, strata = "label", prop = 0.75)
ames_train <- training(data_split)
ames_holdout  <- testing(data_split) 


#RF models require comparably less pre-processing to linear models
rf_recipe <- 
  recipe(
    label ~ . , 
    data = ames_train
  )  
#%>%
# step_log(label, base = 2) 
# %>%
#  step_other(Neighborhood, Overall_Qual, threshold = 50) %>% 
#  step_novel(Neighborhood, Overall_Qual) %>% 
#  step_dummy(Neighborhood, Overall_Qual) 


rf_mod <- rand_forest() %>%
  set_engine("ranger", importance = "impurity", seed = 63233, quantreg = TRUE) %>%
  set_mode("regression")
set.seed(63233)
rf_wf <- workflows::workflow() %>% 
  add_model(rf_mod) %>% 
  add_recipe(rf_recipe) %>% 
  fit(ames_train)


preds_bind <- function(data_fit, lower = 0.01, upper = 0.99){
  predict(
    rf_wf$fit$fit$fit, 
    workflows::pull_workflow_prepped_recipe(rf_wf) %>% bake(data_fit),
    type = "quantiles",
    quantiles = c(lower, upper, 0.50)
  ) %>% 
    with(predictions) %>% 
    as_tibble() %>% 
    set_names(paste0(".pred", c("_lower", "_upper",  ""))) 
  #%>% 
  #mutate(across(contains(".pred"), ~10^.x)) %>% 
  #bind_cols(data_fit) %>% 
  #select(contains(".pred"), label)
}
rf_preds_test <- preds_bind(ames_holdout)
rf_preds_test$label<-ames_holdout$label
rf_preds_test
hist(rf_preds_test$.pred_lower)

table(apply(rf_preds_test, 1, function(x) ifelse(x[1]<x[4] & x[4]<x[2],"I" ,"O")))


set.seed(1234)
rf_preds_test %>% 
  mutate(pred_interval = ggplot2::cut_number(label, 2)) %>% 
  group_by(pred_interval) %>% 
  sample_n(70) %>% 
  ggplot(aes(x = .pred))+
  geom_point(aes(y = .pred, color = "prediction interval"))+
  geom_errorbar(aes(ymin = .pred_lower, ymax = .pred_upper, color = "prediction interval"))+
  geom_point(aes(y = label, color = "actuals"))+
  #scale_x_log10(labels = scales::dollar)+
  #scale_y_log10(labels = scales::dollar)+
  labs(title = "",
       subtitle = "Random Forest Model",
       y = "Sale_Price prediction intervals and actuals")+
  theme_bw()

#Phytools 
packageVersion("phytools")
library(phytools)
anole.tree<-read.tree("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/tree.newick.txt") ## or
plotTree(anole.tree,type="fan",ftype="i")

identify(anole.tree)

library(ape)
library(tidytree)
set.seed(2017)
tree <- rtree(4)
tree

x <- as_tibble(tree)
x
as.phylo(x)

## read tree from file
anole.tree<-read.tree("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/tree.newick_PKJK.txt") ## or
#anole.tree<-read.tree("http://www.phytools.org/eqg2015/data/anole.tre")
## plot tree
plotTree(anole.tree,type="fan",ftype="i")

df_tmp<-data.frame(df$RP)
row.names(df_tmp)<- paste0("seq", seq(1, dim(df)[1]))       
svl<-as.matrix(df_tmp)[,1]

fit<-fastAnc(anole.tree,svl,vars=TRUE,CI=TRUE)
fit

##Real analysis 
library(tidyverse)
library(ape)
library(tidytree)
library(phangorn)

anole.tree<-read.tree("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/tree.newick.txt") 

#Correlation test 
library(tidyverse)
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")
PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 
write_csv(df, "/Users/moradigd/Documents/Plasmid/input_agg.csv")


#Scoary Features 
df_RP<-df[!is.na(df$RP),]
df_RP$Name<-paste0("seq", seq(1, dim(df_RP)[1])) 
df_RP$RP_GWAS<- ifelse(df_RP$RP>=median(df_RP$RP, na.rm = TRUE),1,0)
write_csv(df_RP[,c('Name','RP_GWAS')], "/Users/moradigd/Documents/Plasmid/input_RP_GWAS.csv")


#Scoary Features 
df_PKJK<-df[!is.na(df$PKJK),]
df_PKJK$Name<-paste0("seq", seq(1, dim(df_PKJK)[1])) 
df_PKJK$PKJK_GWAS<- ifelse(df_PKJK$PKJK>=median(df_PKJK$PKJK, na.rm = TRUE),1,0)
write_csv(df_PKJK[,c('Name','PKJK_GWAS')], "/Users/moradigd/Documents/Plasmid/input_PKJK_GWAS.csv")


#Scoary Features 
df_PB<-df[!is.na(df$PB),]
df_PB$Name<-paste0("seq", seq(1, dim(df_PB)[1])) 
df_PB$PB_GWAS<- ifelse(df_PB$PB>=median(df_PB$PB, na.rm = TRUE),1,0)
write_csv(df_PB[,c('Name','PB_GWAS')], "/Users/moradigd/Documents/Plasmid/input_PB_GWAS.csv")





df<-df[!is.na(df$RP),]
df_tmp<-data.frame(df$RP)
row.names(df_tmp)<- paste0("seq", seq(1, dim(df)[1]))   
svl<-as.matrix(df_tmp)[,1]

anole.tree<-midpoint(anole.tree)
tmp<-as_tibble(anole.tree)
estimate<-ace(svl,anole.tree, type="continuous", method="pic")
tmp_ace<-tmp[ which(grepl("seq",tmp$label)) ,]
length(estimate$ace)


#estimate_mat<-data.frame(list(labels=tmp_ace$label, estimates= as.numeric(as.character(estimate$ace)) ))
estimate_mat<-data.frame(list(labels=tmp_ace$label, estimates= svl[match(tmp_ace$label   ,names(svl))] ))



dna<-read.dna("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/sample_out.fasta",format= "fasta")
distdna<-dist.dna(dna,model ="N",pairwise.deletion = TRUE,as.matrix = TRUE)
distdna_ace<- distdna[ which(grepl("seq",row.names(distdna))),which(grepl("seq",row.names(distdna)))]


dist_m<-matrix(0, nrow =  dim(distdna_ace)[1], ncol = dim(distdna_ace)[1])
for(i in seq(dim(distdna_ace)[1])){
  for(j in seq(dim(distdna_ace)[1])){
    dist_m[i,j]<-estimate_mat$estimates[i]/estimate_mat$estimates[j]
  }
}
colnames(dist_m)<-colnames(distdna_ace)
row.names(dist_m)<-row.names(distdna_ace)

intersect<-intersect(colnames(dist_m), row.names(distdna_ace))
dist_m_1<-dist_m[match(intersect,colnames(dist_m)),match(intersect,colnames(dist_m))]
distdna_1<-distdna_ace[match(intersect,colnames(distdna_ace)),match(intersect,colnames(distdna_ace))]

distdna_1[2,]
distdna_1[3,]
distdna_1[4,]

cor.test(dist_m_1[,5],distdna_1[,5])


plot(log(dist_m_1[,5]),distdna_1[,5])






library(ape)
library(phangorn)

anole.tree<-read.tree("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/tree.newick.txt") 
svl<-read.csv("/Users/moradigd/Documents/Plasmid/svl.csv",row.names=1) ## or
## change this into a vector
svl<-as.matrix(svl)[,1]
svl
df_tmp<-data.frame(df$PKJK)
row.names(df_tmp)<- paste0("seq", seq(1, dim(df)[1]))       
svl<-as.matrix(df_tmp)[,1]

tmp<-ace(svl, anole.tree, type="continuous", method="pic")

tmp$ace
getDescendants(anole.tree,node=3089)

anole.tree$tip.label

fit<-fastAnc(tree_tmp,svl,vars=TRUE,CI=TRUE)
termite_tree2<-multi2di(anole.tree,random=TRUE) 
tmp<-ace(svl, tree_tmp, type="continuous", method="pic")

tmp$ace
tree_tmp$node.label

estimate<-data.frame(list(labels=c(names(svl), tree_tmp$node.label), 
                          estimates= as.numeric(as.character(c(svl,tmp$ace ))) ))



library(ape)
library(phytools)
M<-anole.tree$Nnode
N<-length(anole.tree$tip.label)
x<-fastBM(anole.tree)
anc<-vector()
for(i in 1:M+N){
  print(i)
  anc[i-N]<-ace(x,multi2di(root(anole.tree,node=i)), method="pic")$ace[1]
  names(anc)[i-N]<-i
}

hd<-c()
for(i in 1:M+N){
  hd<-c(hd, i)
}


length(tmp$ace)
which(tmp$ace==0)
which(svl==0)

plot(tree_tmp, show.tip.label=FALSE)
tiplabels( cex=(svl-min(svl))/(max(svl)-min(svl)))
nodelabels( cex=(svl-min(svl))/(max(svl)-min(svl)))


#dna_m<-read.dna("seq.marginal.txt")

dna<-read.dna("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/sample_out.fasta",format= "fasta")
distdna<-dist.dna(dna,model ="N",pairwise.deletion = TRUE,as.matrix = TRUE)

dist_m<-matrix(0, nrow =  dim(estimate)[1], ncol = dim(estimate)[1])
for(i in seq(dim(estimate)[1])){
  for(j in seq(dim(estimate)[1])){
    dist_m[i,j]<-estimate$estimates[i]/estimate$estimates[j]
  }
}
colnames(dist_m)<-estimate$labels
row.names(dist_m)<-estimate$labels

intersect<-intersect(colnames(dist_m), row.names(distdna))
dist_m_1<-dist_m[match(intersect,colnames(dist_m)),match(intersect,colnames(dist_m))]
distdna_1<-distdna[match(intersect,colnames(distdna)),match(intersect,colnames(distdna))]

cor.test(dist_m_1[,4:5],distdna_1[,4:5])

plot(log10(dist_m_1[,1]), distdna_1[,1])

?read.dna
fastAnc_w(anole.tree,svl)

fastAnc_w<-function(tree,x){
  if(!is.binary.tree(tree)) btree<-multi2di(tree)
  else btree<-tree
  M<-btree$Nnode
  N<-length(btree$tip)
  anc<-vector()
  for(i in 1:M+N){
    anc[i-N]<-ace(x,multi2di(root(btree,node=i)),
                  method="pic")$ace[1]
    names(anc)[i-N]<-i
  }
  if(!is.binary.tree(tree)){
    ancNames<-matchNodes(tree,btree)
    anc<-anc[as.character(ancNames[,2])]
    names(anc)<-ancNames[,1]
  }
  return(anc)
}

fit
is.binary(anole.tree)
anole.tree$tip.label

#fit<-fastAnc(anole.tree,svl,vars=TRUE,CI=TRUE)
#fit
#HierBA
#anole.tree$node.label
library(phangorn)
tree_tmp<-midpoint(anole.tree)

anole.tree$node.label[1]
obj<-contMap(midpoint(anole.tree),svl,plot=FALSE)
plot(obj,type="fan",legend=0.7*max(nodeHeights(anole.tree)),
     fsize=c(0.7,0.9))

library(geiger)
nodelabel.phylo(anole.tree)
nodelabels(anole.tree)
library(tidytree)
nodeid(midpoint(anole.tree), names(x))
nodelab(anole.tree, names(svl))
?nodeid
require(devtools) 
as.treedata(anole.tree)
tree <- rtree(4)
x <- as_tibble(tree)

as.phylo(anole.tree)
install_version("tidytree", version = "0.3.4", repos = "http://cran.us.r-project.org")

tmp<-ace(svl, anole.tree, type="continuous", method="pic")

getSisters(anole.tree, 2593, mode=c("number","label"))

nd <- data.frame(node = names(fit$ace), trait = fit$ace)

require(remotes)
install_version("tidytree", version = "0.3.4")

nodelabels(anole.tree)

my_tree=anole.tree
node_labels_in_edge <- my_tree$tip.label[my_tree$edge[,1]-Ntip(my_tree)]
tips_nodes <- my_tree$edge[,2]

select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
}

edge_table <- data.frame(
  "parent" = my_tree$edge[,1],
  "par.name" = sapply(my_tree$edge[,1], select.tip.or.node, tree = my_tree),
  "child" = my_tree$edge[,2],
  "chi.name" = sapply(my_tree$edge[,2], select.tip.or.node, tree = my_tree)
)

library(DECIPHER)

library(msa)
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet("/Users/moradigd/Documents/Plasmid/seq.fasta")
mySequences
myFirstAlignment <- msa(mySequences)

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}
alignment2Fasta(myFirstAlignment, '/Users/moradigd/Documents/Plasmid/seq_1.fasta')

write.msa(myFirstAlignment, file="/Users/moradigd/Documents/Plasmid/seq_1.fasta",format= "FASTA")
library(phangorn)

primates <- read.phyDat("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/seq_PKJK_aligned.fasta", format =  "fasta")
tree <- pratchet(primates, trace=0) %>% acctran(primates)
tree_midpoint <- midpoint(tree)
plotTree(tree_midpoint,type="fan",ftype="i")


parsimony(tree, primates)

anole.tree<-read.tree("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/tree.newick_PKJK.txt") ## or

anc.acctran <- ancestral.pars(tree_midpoint, primates, "ACCTRAN")
anc.mpr <- ancestral.pars(tree, primates, "MPR")

t(subset(anc.acctran, getRoot(tree), 1:20)[[1]])

library(seqLogo)
seqLogo( t(subset(anc.acctran, getRoot(tree), 1:20)[[1]]), ic.scale=FALSE)

svl<-read.csv("/Users/moradigd/Documents/Plasmid/input_agg.csv",row.names=1)
svl<-svl[,1]
svl<-svl[complete.cases(svl)]
## or
svl<-as.matrix(svl)[,1]
names(svl)<-paste0("seq", seq(1, length(svl)))  
#df_tmp<-data.frame(df$PKJK)
#row.names(df_tmp)<- paste0("seq", seq(1, dim(df)[1]))       
#svl<-as.matrix(df_tmp)[,1]
write.tree(tree_midpoint,"/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/tree.newick_PKJK_m.txt")
anole.tree_1<-read.tree("/Users/moradigd/Documents/Plasmid/tree_1.newick") ## or
anole.tree_1$node.label<-paste0("N",seq(length(anole.tree_1$node.label)))

fit<-fastAnc(tree_midpoint ,svl)
fit
seq(anole.tree_1$node.label)

tree_midpoint$tip.label
tree_midpoint$node.label<-paste0("N",seq(length(tree_midpoint$node.label)))

tmp<-ace(svl, tree_midpoint, type="continuous", method="pic")
tail(names(tmp$ace))

#chemical genomcis 
library(tidyverse)
chem<-read_csv("/Users/moradigd/Documents/chemicalgenomic/Klebsiella.csv") %>% 
  filter(Confirmation=="Co") %>%
  filter(!grepl("intergenic", Notes)) 


dim(chem)
length(chem$`KP Locus`)
length( unique(chem$`KP Locus`))
hist(chem$From, breaks = 1000)

insertaion<-as.numeric(sapply(chem$`Position relative to locus`, function(x) strsplit(x[1],"\\(")[[1]][1] )) 
length_gene<-apply(chem, 1, function(x) as.numeric(as.character(x[11]))-as.numeric(as.character(x[10]))+1)

hist(insertaion/length_gene, breaks = 500)
gene_length<-c()
for(j in seq(0,1,0.02)){
  gene_length<-c(gene_length, length(unique(chem$`KP Locus`[which((insertaion/length_gene)>j)])))
}
barplot(gene_length, names.arg = seq(0,1,0.02), las=2, col="black")
box()

df<-data.frame(list(gene_length=gene_length, cutoff=seq(0,1,0.02)))
library(highcharter) 
hc <- df %>%
  hchart('column', hcaes(x = cutoff, y = gene_length))
hc 

#TB chemical genomics

#chemical genomcis 
library(tidyverse)
chem<-read_csv("/Users/moradigd/Documents/TB/TB.csv") %>% 
  filter(`Tn% after start codon`!=".") 
#%>%filter(!grepl("intergenic", Notes)) 
chem$`name Danish`<-chem$`locus Danish`
chem$`Tn% after start codon`<-as.numeric(as.character(chem$`Tn% after start codon`))
chem$gene_affected<-chem$`Tn% after start codon`/100

dim(chem)
length(chem$`KP Locus`)
length( unique(chem$`KP Locus`))
hist(chem$gene_affected, breaks = 1000)

insertaion<-as.numeric(sapply(chem$`Position relative to locus`, function(x) strsplit(x[1],"\\(")[[1]][1] )) 
length_gene<-apply(chem, 1, function(x) as.numeric(as.character(x[11]))-as.numeric(as.character(x[10]))+1)

hist(insertaion/length_gene, breaks = 500)
gene_length<-c()
for(j in seq(0,1,0.02)){
  gene_length<-c(gene_length, length(unique(chem$`name Danish`[which((chem$gene_affected)>j)])))
}
barplot(gene_length, names.arg = seq(0,1,0.02), las=2, col="black")
box()

df<-data.frame(list(gene_length=gene_length, cutoff=seq(0,1,0.02)))
library(highcharter) 
hc <- df %>%
  hchart('column', hcaes(x = cutoff, y = gene_length))
hc 



cutoffs<-seq(0,0.98,0.02)
results<-matrix(0, nrow = length(cutoffs), ncol = 5)
results_counts<-matrix(0, nrow = length(cutoffs), ncol = 1)

for(k in seq(1, length(cutoffs))){
  print(cutoffs[k])
  gene<-chem$`name Danish`[which((chem$gene_affected)>cutoffs[k])] 
  freq_table<-table(sort(table(gene), decreasing = T))
  for(j in seq(1, length(freq_table))){
    if(as.numeric(names(freq_table[j]))>=5){
      results[k, 5]=results[k, 5]+freq_table[j]
    }
    else{
      results[k, as.numeric(names(freq_table[j]))]=freq_table[j]
    }
  }
  
  results_counts[k,1]<-length(gene)
}
results<-data.frame(results)
row.names(results)<-cutoffs
colnames(results)<-c("1 mut","2 mut","3 mut","4 mut","≥5 mut")
results<-results%>%
  pivot_longer(
    everything(),
    values_to = "value"
  )
results$cutoff<-as.character(as.vector(sapply(cutoffs, function(x) rep(x,5) )))
hchart(results,'column', hcaes(x = 'cutoff', y = 'value', group = 'name'),stacking = "normal")%>%
  hc_colors(c("#FF00FF", "#00FFFF", "#FF0000","#00FF00","#0000FF"))%>%
  hc_yAxis(title = list(text = "Count unique loci included"))

results_counts<-data.frame(results_counts)
results_counts$cutoff<-as.character(cutoffs)
#row.names(results_counts)<-cutoffs
hchart(results_counts,'column', hcaes(x='cutoff', y = 'results_counts')) %>%
  hc_colors(c("#000000")) %>%
  hc_yAxis(title = list(text = "Count mutants included"))






cutoffs<-seq(0,0.98,0.02)
results<-matrix(0, nrow = length(cutoffs), ncol = 2)
results_counts<-matrix(0, nrow = length(cutoffs), ncol = 1)

for(k in seq(1, length(cutoffs))){
  print(cutoffs[k])
  gene<-chem$`name Danish`[which((chem$gene_affected)>cutoffs[k])] 
  freq_table<-table(sort(table(gene), decreasing = T))
  for(j in seq(1, length(freq_table))){
    if(as.numeric(names(freq_table[j]))>=2){
      results[k, 2]=results[k, 2]+freq_table[j]
    }
    else{
      results[k, as.numeric(names(freq_table[j]))]=freq_table[j]
    }
  }
  
  results_counts[k,1]<-results[k, 1]+2*results[k, 2]
}
results<-data.frame(results)
row.names(results)<-cutoffs
colnames(results)<-c("1 mut","≥2 mut")
results<-results%>%
  pivot_longer(
    everything(),
    values_to = "value"
  )
results$cutoff<-as.character(as.vector(sapply(cutoffs, function(x) rep(x,2) )))
hchart(results,'column', hcaes(x = 'cutoff', y = 'value', group = 'name'),stacking = "normal")%>%
  hc_colors(c("#FF00FF", "#00FFFF")) %>%
  hc_yAxis(title = list(text = "Count unique loci included"))

results_counts<-data.frame(results_counts)
results_counts$cutoff<-as.character(cutoffs)
#row.names(results_counts)<-cutoffs
hchart(results_counts,'column', hcaes(x='cutoff', y = 'results_counts')) %>%
  hc_colors(c("#000000")) %>%
  hc_yAxis(title = list(text = "Count mutants included single+two mutants"))



#Iterative random forest 

#Today and tomorrow 
#Simulation the day after tomorrow 
#Next week write-up plasmid

#iRF
library(tidyverse)
input<-read_csv("/Users/moradigd/Documents/Plasmid/test.csv") %>% 
  select(-X1)
sum_columns<- as.numeric(as.character(apply(input, 2, function(x) sum(x)/dim(input)[1])))  
input[,which(sum_columns>0.4 & sum_columns<0.6)]

#Correlation test 
library(tidyverse)
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 

PKJK<-PKJK[complete.cases(PKJK)]
hist(PKJK, breaks = 100)




library(iRF)
n <- 500
p <- 250
X <- matrix(rnorm(n * p), nrow=n)
hist(X)

Y <- (X[,1] > 0 & X[,2] > 0 & X[,3] > 0)
Y <- as.factor(as.numeric(Y))

y <- as.numeric(Y) + rnorm(n, sd=0.25)


train.id <- 1:(n / 2)
test.id <- setdiff(1:n, train.id)


sel.prob <- rep(1/p, p)

# iteratively grow RF, use Gini importance of features as weights
rf <- list()
for (iter in 1:4){
  rf[[iter]] <- randomForest(x=X[train.id,], y=Y[train.id], 
                             xtest=X[test.id,], ytest=Y[test.id], 
                             mtry.select.prob=sel.prob)
  
  # update selection probabilities for next iteration
  sel.prob <- rf[[iter]]$importance
}

library(AUC)
plot(0:1, 0:1, type='l', lty = 2, xlab = 'FPR', ylab = 'TPR', main='ROC Curve')
for (iter in 1:4){
  # performance on test set
  cat(paste('iter = ', iter, ':: '))
  roc.info <- roc(rf[[iter]]$test$votes[,2], Y[test.id])
  lines(roc.info$fpr, roc.info$tpr, type='l', col=iter, lwd=2)
  cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep=''))
} 
legend('bottomright', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')


sel.prob <- rep(1/p, p)

# iteratively grow RF, use Gini importance of features as weights
rf <- list()
for (iter in 1:4){
  rf[[iter]] <- randomForest(x=X[train.id,], y=Y[train.id], 
                             xtest=X[test.id,], ytest=Y[test.id], 
                             mtry.select.prob=sel.prob, 
                             keep.forest=TRUE)
  
  # update selection probabilities for next iteration
  sel.prob <- rf[[iter]]$importance / sum(rf[[iter]]$importance)
}
y <- as.numeric(Y) + rnorm(n, sd=0.25)

fit <- iRF(x=X[train.id,], 
           y=y[train.id], 
           xtest=X[test.id,], 
           ytest=y[train.id], 
           n.iter=5, 
           interactions.return=3,
           select.iter = TRUE,
           n.bootstrap=10
)

library(dplyr)
head(fit$importance)
fit$rf.list


library(dplyr)
fit$interaction
head(fit$importance)


class1.nodes <- rforest$tree.info$prediction - 1 == 1
wt <- rforest$tree.info$size.node[class1.nodes]
RIT(rforest$node.feature[class1.nodes,], weights=wt,
    depth=5, branch=2, n_trees=100)


# set up cores for parallel implementation
library(doMC)
registerDoMC()
n.cores <- 4
options(cores=n.cores)
n.tree.per.core <- 30

rfpar <- foreach(n.tree=rep(n.tree.per.core, n.cores), 
                 .combine=combine, .multicombine=TRUE) %dopar% {
                   randomForest(x=X[train.id,], y=Y[train.id], ntree=n.tree)
                 }



sel.prob <- rep(1/p, p)

# iteratively grow RF, use Gini importance of features as weights
rf <- list()
for (iter in 1:4){
  rf[[iter]] <- randomForest(x=X[train.id,], y=Y[train.id], 
                             xtest=X[test.id,], ytest=Y[test.id], 
                             mtry.select.prob=sel.prob, 
                             keep.forest=TRUE)
  
  # update selection probabilities for next iteration
  sel.prob <- rf[[iter]]$importance / sum(rf[[iter]]$importance)
}
rforest <- readForest(rfobj=rf[[3]], x=X[train.id,])
head(rforest$tree.info, n=10)


class1.nodes <- rforest$tree.info$prediction - 1 == 1
wt <- rforest$tree.info$size.node[class1.nodes]

RIT(rforest$node.feature[class1.nodes,], weights=wt,
    depth=5, branch=2, n_trees=100)



#K-mer Simulation 
#for i in 0.01 0.05 0.01 0.25 0.5 0.25 1; do for j in 100 500 1000 2000; do ./SimBac -N ${j} -T ${i} -B 700 -o genomes_${i}_${j}.fasta ; done ; done
#for i in 0.01 0.05 0.01 0.25 0.5 0.25 1; do for j in 100 500 1000 2000; do snp-sites genomes_${i}_${j}.fasta > variant_genomes_${i}_${j}.vcf; done ; done
#for i in 0.0001 0.001 0.01 0.1; do for j in 100 500 1000 2000; do ./SimBac -N ${j} -T ${i} -B 700 -o genomes_${i}_${j}.fasta; done ; done


#Figure 1 
library(tidyverse)
inp<-read_csv("/Users/moradigd/Documents/Plasmid/results_screening_agg_scale.csv") %>%
  group_by(col_2, col_3, col_4, col_6) %>% 
  drop_na()%>% 
  summarise(mean = mean(col_0), sd = sd(col_0))


inp$col_0[inp$col_3=="RP" & inp$col_4=="lasso" & inp$col_6=="Kmers"]
inp$col_0[inp$col_3=="RP" & inp$col_4=="lasso" & inp$col_6=="BAPS"]




#GWAS
inp<-read_csv("/Users/moradigd/Documents/Plasmid/PKJK_GWAS_10_08_2021_0855.results.csv")%>% 
  filter(Best_pairwise_comp_p < 0.05 & Worst_pairwise_comp_p<0.05 & Bonferroni_p<0.05) 


inp$Bonferroni_p[which(sapply(inp$Gene, function(x) nchar(x))>15)]
inp$Odds_ratio[which(sapply(inp$Gene, function(x) nchar(x))>15)]
tmp_gene<-inp$Gene[which(sapply(inp$Gene, function(x) nchar(x))>10)]

inp$Bonferroni_p[which(sapply(inp$Gene, function(x) nchar(x))>1),]

genes<-read_csv("/Users/moradigd/Documents/Plasmid/pan_genome_agg_PB.csv")

tmp_mat<-as.data.frame(genes[ which(genes$Gene %in% tmp_gene),-1])
for(i in seq(dim(tmp_mat)[1])){
  print(sum(as.numeric(tmp_mat[i,])))
  for(j in seq(dim(tmp_mat)[1])){
    print(sum(as.numeric(tmp_mat[i,])))
    print(cor(as.numeric(tmp_mat[i,]),as.numeric(tmp_mat[j,])))
  }
}



genes<-read_csv("/Users/moradigd/Documents/Plasmid/RP_GWAS_09_08_2021_2251.results.csv")
#group_by(tags_res) %>% 
#  summarise(pca= pca,test_mean = mean(test_mae),train_mean = mean(train_mae),pvalue= t.test(test_sp,mu=0)$p.value , sd_test = sd(train_mae),sd_test = sd(train_mae),mean_sp_test=mean(test_sp),mean_sp_train=mean(train_sp), sd_sp_test=sd(test_sp),sd_sp_train=sd(train_sp))

#Interval prediction cleans

#Interval prediction 

library(tidyverse)
library(tidymodels)
library(gt)
# function copied from here:
# https://github.com/rstudio/gt/issues/613#issuecomment-772072490 
# (simpler solution should be implemented in future versions of {gt})
fmt_if_number <- function(..., digits = 2) {
  input <- c(...)
  fmt <- paste0("%.", digits, "f")
  if (is.numeric(input))   return(sprintf(fmt, input))
  return(input)
}


predictor<-read_csv("/Users/moradigd/Documents/Plasmid/pan_genome_agg_PB.csv")

gene_length<-sapply(predictor$Gene, function(x) nchar(x))   

predictor<-predictor[which(gene_length==9),]
predictor<-predictor[,15:dim(predictor)[2]]
predictor<-data.frame(t(predictor))

tmp_pca<-prcomp(predictor, center = TRUE,scale. = TRUE)
predictor<-data.frame(tmp_pca$x)

#Correlation test 
library(tidyverse)
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

#PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
#RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
#PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
#df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 
#write_csv(df, "/Users/moradigd/Documents/Plasmid/input_agg.csv")
outcome<-read_csv("/Users/moradigd/Documents/Plasmid/input_agg.csv")

outcome<-outcome[!is.na(outcome$PB),]

predictor$label<-outcome$PB
predictor$label<-scale(predictor$label)
predictor<- data.frame(predictor)

data_split <- initial_split(predictor, strata = "label", prop = 0.8)
ames_train <- training(data_split)
ames_holdout  <- testing(data_split) 


#RF models require comparably less pre-processing to linear models
options(expressions = 5e5)
rf_recipe <- 
  recipe(
    label ~ . , 
    data = ames_train
  )  

preds_bind <- function(data_fit, lower = 0.01, upper = 0.99){
  predict(
    rf_wf$fit$fit$fit, 
    workflows::pull_workflow_prepped_recipe(rf_wf) %>% bake(data_fit),
    type = "quantiles",
    quantiles = c(lower, upper, 0.50)
  ) %>% 
    with(predictions) %>% 
    as_tibble() %>% 
    set_names(paste0(".pred", c("_lower", "_upper",  ""))) 
}

mtry_vals<-c(20)
trees_vals<-c(100)

for(k in seq(length(mtry_vals))){
  
  for(m in seq(length(trees_vals))){
    print(mtry_vals[k])
    print(trees_vals[m])
    
    rf_mod <- rand_forest(mtry = mtry_vals[k], trees = trees_vals[m]) %>%
      set_engine("ranger", importance = "impurity", seed = 63233, quantreg = TRUE) %>%
      set_mode("regression")
    set.seed(63233)
    rf_wf <- workflows::workflow() %>% 
      add_model(rf_mod) %>% 
      add_recipe(rf_recipe) %>% 
      fit(ames_train)
    
    lower_q<-seq(1,25)*0.01
    upper_q<-1-seq(1,25)*0.01
    report_output_test<-c()
    report_output_train<-c()
    
    for(i in seq(1,length(lower_q))){
      rf_preds_test <- preds_bind(ames_holdout, lower_q[i], upper_q[i])
      rf_preds_test$label<-ames_holdout$label
      report_final<-table(apply(rf_preds_test, 1, function(x) ifelse(x[1]<x[4] & x[4]<x[2],"I" ,"O")))
      report_output_test<-c(report_output_test,report_final[1]/sum(report_final))
      
      rf_preds_train <- preds_bind(ames_train, lower_q[i], upper_q[i])
      rf_preds_train$label<-ames_train$label
      report_final<-table(apply(rf_preds_train, 1, function(x) ifelse(x[1]<x[4] & x[4]<x[2],"I" ,"O")))
      report_output_train<-c(report_output_train,report_final[1]/sum(report_final))
    }
    print(as.numeric(as.character(report_output_train[1])))
    print("")
  }
}

rf_mod <- rand_forest() %>%
  set_engine("ranger", importance = "impurity", seed = 63233, quantreg = TRUE) %>%
  set_mode("regression")
set.seed(63233)
rf_wf <- workflows::workflow() %>% 
  add_model(rf_mod) %>% 
  add_recipe(rf_recipe) %>% 
  fit(ames_train)

lower_q<-seq(1,25)*0.01
upper_q<-1-seq(1,25)*0.01
report_output_test<-c()
report_output_train<-c()

for(i in seq(1,length(lower_q))){
  print(i)
  rf_preds_test <- preds_bind(ames_holdout, lower_q[i], upper_q[i])
  rf_preds_test$label<-ames_holdout$label
  report_final<-table(apply(rf_preds_test, 1, function(x) ifelse(x[1]<x[4] & x[4]<x[2],"I" ,"O")))
  report_output_test<-c(report_output_test,report_final[1]/sum(report_final))
  
  rf_preds_train <- preds_bind(ames_train, lower_q[i], upper_q[i])
  rf_preds_train$label<-ames_train$label
  report_final<-table(apply(rf_preds_train, 1, function(x) ifelse(x[1]<x[4] & x[4]<x[2],"I" ,"O")))
  report_output_train<-c(report_output_train,report_final[1]/sum(report_final))
}
print(report_output_test)
print(report_output_train)

list(value=c(report_output_test, report_output_train), 
     dataset=c(rep("test", length(report_output_test)), rep("train", length(report_output_train))),
     interval=c(paste0(lower_q,"-",upper_q), paste0(lower_q,"-",upper_q))
) %>% 
  data.frame(.) %>%
  ggplot(aes(fill=dataset, x=interval,y=value))+
  geom_bar( stat = "identity", position = "dodge")+
  theme_bw() +
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )+
  ylab("True detection rate")+
  xlab("Interval")+ 
  ylim(range(0,1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("#998899", "#E69F00"))

set.seed(1234)
rf_preds_test <- preds_bind(ames_holdout, 0.05, 0.95)
rf_preds_test$label<-ames_holdout$label
report_final<-table(apply(rf_preds_test, 1, function(x) ifelse(x[1]<x[4] & x[4]<x[2],"I" ,"O")))

rf_preds_test %>% 
  mutate(pred_interval = ggplot2::cut_number(label, 2)) %>% 
  group_by(pred_interval) %>% 
  sample_n(30) %>% 
  ggplot(aes(x = .pred))+
  geom_point(aes(y = .pred, color = "prediction interval"))+
  geom_errorbar(aes(ymin = .pred_lower, ymax = .pred_upper, color = "prediction interval"))+
  geom_point(aes(y = label, color = "actuals"))+
  #scale_x_log10(labels = scales::dollar)+
  #scale_y_log10(labels = scales::dollar)+
  labs(title = "",
       subtitle = "Random Forest Model",
       y = "Permissiveness prediction interval",
       x="Permissiveness prediction")+
  theme_bw()+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))

install.packages(vcfR)
library(vcfR)
vcf = read.vcfR("/Users/moradigd/Desktop/Y", verbose=TRUE, nrows = 10)
vcf

#Preparation of samples 
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")
PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))


PKJK_RP_vals_train<-PKJK[which(!is.na(PKJK) & is.na(RP))]
PKJK_RP_seq_train<-corr$seq[which(!is.na(PKJK) & is.na(RP))]

PKJK_RP_vals_test<-RP[which(is.na(PKJK) & !is.na(RP))]
PKJK_RP_seq_test<-corr$seq[which(is.na(PKJK) & !is.na(RP))]

PKJK_RP_train<-data.frame(list("seq"=PKJK_RP_seq_train, "value"=PKJK_RP_vals_train, "label"=rep("train", length(PKJK_RP_vals_train))))
PKJK_RP_test<-data.frame(list("seq"=PKJK_RP_seq_test, "value"=PKJK_RP_vals_test, "label"=rep("test", length(PKJK_RP_vals_test))))

write_csv(rbind(PKJK_RP_train,PKJK_RP_test), "/Users/moradigd/Documents/Plasmid/input_PKJK_RP_data.csv")


PKJK_PB_vals_train<-PKJK[which(!is.na(PKJK) & is.na(PB))]
PKJK_PB_seq_train<-corr$seq[which(!is.na(PKJK) & is.na(PB))]

PKJK_PB_vals_test<-PB[which(is.na(PKJK) & !is.na(PB))]
PKJK_PB_seq_test<-corr$seq[which(is.na(PKJK) & !is.na(PB))]

PKJK_PB_train<-data.frame(list("seq"=PKJK_PB_seq_train, "value"=PKJK_PB_vals_train, "label"=rep("train", length(PKJK_PB_vals_train))))
PKJK_PB_test<-data.frame(list("seq"=PKJK_PB_seq_test, "value"=PKJK_PB_vals_test, "label"=rep("test", length(PKJK_PB_vals_test))))

write_csv(rbind(PKJK_PB_train,PKJK_PB_test), "/Users/moradigd/Documents/Plasmid/input_PKJK_PB_data.csv")


RP_PB_vals_train<-RP[which(!is.na(RP) & is.na(PB))]
RP_PB_seq_train<-corr$seq[which(!is.na(RP) & is.na(PB))]

RP_PB_vals_test<-PB[which(is.na(RP) & !is.na(PB))]
RP_PB_seq_test<-corr$seq[which(is.na(RP) & !is.na(PB))]

RP_PB_train<-data.frame(list("seq"=RP_PB_seq_train, "value"=RP_PB_vals_train, "label"=rep("train", length(RP_PB_vals_train))))
RP_PB_test<-data.frame(list("seq"=RP_PB_seq_test, "value"=RP_PB_vals_test, "label"=rep("test", length(RP_PB_vals_test))))

write_csv(rbind(RP_PB_train,RP_PB_test), "/Users/moradigd/Documents/Plasmid/input_RP_PB_data.csv")


#Preparation of samples 
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")
PKJK_DK<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP_DK<-apply(corr[, which(grepl("_RP",colnames(corr)) &  grepl("DK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB_DK<-apply(corr[, which(grepl("_PB",colnames(corr)) &  grepl("DK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))

PKJK_UK<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("UK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP_UK<-apply(corr[, which(grepl("_RP",colnames(corr)) &  grepl("UK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB_UK<-apply(corr[, which(grepl("_PB",colnames(corr)) &  grepl("UK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))

PKJK_vals_train<-PKJK_DK[which(!is.na(PKJK_DK) & is.na(PKJK_UK))]
PKJK_seq_train<-corr$seq[which(!is.na(PKJK_DK) & is.na(PKJK_UK))]

PKJK_vals_test<-PKJK_UK[which(is.na(PKJK_DK) & !is.na(PKJK_UK))]
PKJK_seq_test<-corr$seq[which(is.na(PKJK_DK) & !is.na(PKJK_UK))]

PKJK_train<-data.frame(list("seq"=PKJK_seq_train, "value"=PKJK_vals_train, "label"=rep("train", length(PKJK_vals_train))))
PKJK_test<-data.frame(list("seq"=PKJK_seq_test, "value"=PKJK_vals_test, "label"=rep("test", length(PKJK_vals_test))))

write_csv(rbind(PKJK_train,PKJK_test), "/Users/moradigd/Documents/Plasmid/input_PKJK_country_data.csv")


RP_vals_train<-RP_DK[which(!is.na(RP_DK) & is.na(RP_UK))]
RP_seq_train<-corr$seq[which(!is.na(RP_DK) & is.na(RP_UK))]

RP_vals_test<-RP_UK[which(is.na(RP_DK) & !is.na(RP_UK))]
RP_seq_test<-corr$seq[which(is.na(RP_DK) & !is.na(RP_UK))]

RP_train<-data.frame(list("seq"=RP_seq_train, "value"=RP_vals_train, "label"=rep("train", length(RP_vals_train))))
RP_test<-data.frame(list("seq"=RP_seq_test, "value"=RP_vals_test, "label"=rep("test", length(RP_vals_test))))

write_csv(rbind(RP_train,RP_test), "/Users/moradigd/Documents/Plasmid/input_RP_country_data.csv")


PB_vals_train<-PB_DK[which(!is.na(PB_DK) & is.na(PB_UK))]
PB_seq_train<-corr$seq[which(!is.na(PB_DK) & is.na(PB_UK))]

PB_vals_test<-PB_UK[which(is.na(PB_DK) & !is.na(PB_UK))]
PB_seq_test<-corr$seq[which(is.na(PB_DK) & !is.na(PB_UK))]

PB_train<-data.frame(list("seq"=PB_seq_train, "value"=PB_vals_train, "label"=rep("train", length(PB_vals_train))))
PB_test<-data.frame(list("seq"=PB_seq_test, "value"=PB_vals_test, "label"=rep("test", length(PB_vals_test))))

write_csv(rbind(PB_train,PB_test), "/Users/moradigd/Documents/Plasmid/input_PB_country_data.csv")






df<-data.frame(list("seq"=corr$seq,"PKJK_DK"=PKJK_DK, "RP_DK"=RP_DK , "PB_DK"=PB_DK , "PKJK_UK"=PKJK_UK ,"RP_UK"=RP_UK,"PB_UK"=PB_UK)) 
write_csv(df, "/Users/moradigd/Documents/Plasmid/input_agg_country_plasmid.csv")

#Preparation of samples 
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")
PKJK_DK_2017<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("2017",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_2018<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("2018",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))

PKJK_DK_2017_vals<-PKJK_DK_2017[which(is.na(PKJK_DK_2018) & !is.na(PKJK_DK_2017))]
PKJK_DK_2017_seq<-corr$seq[which(is.na(PKJK_DK_2018) & !is.na(PKJK_DK_2017))]

PKJK_DK_2018_vals<-PKJK_DK_2018[which(!is.na(PKJK_DK_2018) & is.na(PKJK_DK_2017))]
PKJK_DK_2018_seq<-corr$seq[which(!is.na(PKJK_DK_2018) & is.na(PKJK_DK_2017))]

tmp_train<-data.frame(list("seq"=PKJK_DK_2017_seq, "value"=PKJK_DK_2017_vals, "label"=rep("train", length(PKJK_DK_2017_vals))))
tmp_test<-data.frame(list("seq"=PKJK_DK_2018_seq, "value"=PKJK_DK_2018_vals, "label"=rep("test", length(PKJK_DK_2018_vals))))

write_csv(rbind(tmp_train,tmp_test), "/Users/moradigd/Documents/Plasmid/input_agg_country_year.csv")

#sites
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")
PKJK_DK_site_1A_test<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_1A",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_1B_test<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_1B",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_2_test<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_2",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_4_test<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_4",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
#PKJK_DK_site_5<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_5",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_6_test<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_6",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_7_test<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_7",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))

PKJK_DK_site_1A_train<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  !grepl("_1A",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_1B_train<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  !grepl("_1B",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_2_train<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  !grepl("_2",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_4_train<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  !grepl("_4",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
#PKJK_DK_site_5<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  grepl("_5",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_6_train<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  !grepl("_6",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK_DK_site_7_train<-apply(corr[, which(grepl("_PKJK",colnames(corr)) &  grepl("DK",colnames(corr))  &  !grepl("_7",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))

PKJK_DK_site_1A_train_vals<-PKJK_DK_site_1A_train[which(is.na(PKJK_DK_site_1A_test) & !is.na(PKJK_DK_site_1A_train))]
PKJK_DK_site_1A_train_seq<-corr$seq[which(is.na(PKJK_DK_site_1A_test) & !is.na(PKJK_DK_site_1A_train))]
tmp_train<-data.frame(list("seq"=PKJK_DK_site_1A_train_seq, "value"=PKJK_DK_site_1A_train_vals, "label"=rep("train", length(PKJK_DK_site_1A_train_vals))))

PKJK_DK_site_1A_test_vals<-PKJK_DK_site_1A_test[which(!is.na(PKJK_DK_site_1A_test) & is.na(PKJK_DK_site_1A_train))]
PKJK_DK_site_1A_test_seq<-corr$seq[which(!is.na(PKJK_DK_site_1A_test) & is.na(PKJK_DK_site_1A_train))]
tmp_test<-data.frame(list("seq"=PKJK_DK_site_1A_test_seq, "value"=PKJK_DK_site_1A_test_vals, "label"=rep("test", length(PKJK_DK_site_1A_test_vals))))

write_csv(rbind(tmp_train,tmp_test), "/Users/moradigd/Documents/Plasmid/input_agg_sites_PKJK_DK_site_1A.csv")


PKJK_DK_site_7_train_vals<-PKJK_DK_site_7_train[which(is.na(PKJK_DK_site_7_test) & !is.na(PKJK_DK_site_7_train))]
PKJK_DK_site_7_train_seq<-corr$seq[which(is.na(PKJK_DK_site_7_test) & !is.na(PKJK_DK_site_7_train))]
tmp_train<-data.frame(list("seq"=PKJK_DK_site_7_train_seq, "value"=PKJK_DK_site_7_train_vals, "label"=rep("train", length(PKJK_DK_site_7_train_vals))))

PKJK_DK_site_7_test_vals<-PKJK_DK_site_7_test[which(!is.na(PKJK_DK_site_7_test) & is.na(PKJK_DK_site_7_train))]
PKJK_DK_site_7_test_seq<-corr$seq[which(!is.na(PKJK_DK_site_7_test) & is.na(PKJK_DK_site_7_train))]
tmp_test<-data.frame(list("seq"=PKJK_DK_site_7_test_seq, "value"=PKJK_DK_site_7_test_vals, "label"=rep("test", length(PKJK_DK_site_7_test_vals))))

write_csv(rbind(tmp_train,tmp_test), "/Users/moradigd/Documents/Plasmid/input_agg_sites_PKJK_DK_site_7.csv")


PKJK_DK_site_6_train_vals<-PKJK_DK_site_6_train[which(is.na(PKJK_DK_site_6_test) & !is.na(PKJK_DK_site_6_train))]
PKJK_DK_site_6_train_seq<-corr$seq[which(is.na(PKJK_DK_site_6_test) & !is.na(PKJK_DK_site_6_train))]
tmp_train<-data.frame(list("seq"=PKJK_DK_site_6_train_seq, "value"=PKJK_DK_site_6_train_vals, "label"=rep("train", length(PKJK_DK_site_6_train_vals))))

PKJK_DK_site_6_test_vals<-PKJK_DK_site_6_test[which(!is.na(PKJK_DK_site_6_test) & is.na(PKJK_DK_site_6_train))]
PKJK_DK_site_6_test_seq<-corr$seq[which(!is.na(PKJK_DK_site_6_test) & is.na(PKJK_DK_site_6_train))]
tmp_test<-data.frame(list("seq"=PKJK_DK_site_6_test_seq, "value"=PKJK_DK_site_6_test_vals, "label"=rep("test", length(PKJK_DK_site_6_test_vals))))

write_csv(rbind(tmp_train,tmp_test), "/Users/moradigd/Documents/Plasmid/input_agg_sites_PKJK_DK_site_6.csv")


PKJK_DK_site_2_train_vals<-PKJK_DK_site_2_train[which(is.na(PKJK_DK_site_2_test) & !is.na(PKJK_DK_site_2_train))]
PKJK_DK_site_2_train_seq<-corr$seq[which(is.na(PKJK_DK_site_2_test) & !is.na(PKJK_DK_site_2_train))]
tmp_train<-data.frame(list("seq"=PKJK_DK_site_2_train_seq, "value"=PKJK_DK_site_2_train_vals, "label"=rep("train", length(PKJK_DK_site_2_train_vals))))

PKJK_DK_site_2_test_vals<-PKJK_DK_site_2_test[which(!is.na(PKJK_DK_site_2_test) & is.na(PKJK_DK_site_2_train))]
PKJK_DK_site_2_test_seq<-corr$seq[which(!is.na(PKJK_DK_site_2_test) & is.na(PKJK_DK_site_2_train))]
tmp_test<-data.frame(list("seq"=PKJK_DK_site_2_test_seq, "value"=PKJK_DK_site_2_test_vals, "label"=rep("test", length(PKJK_DK_site_2_test_vals))))

write_csv(rbind(tmp_train,tmp_test), "/Users/moradigd/Documents/Plasmid/input_agg_sites_PKJK_DK_site_2.csv")



PKJK_DK_site_4_train_vals<-PKJK_DK_site_4_train[which(is.na(PKJK_DK_site_4_test) & !is.na(PKJK_DK_site_4_train))]
PKJK_DK_site_4_train_seq<-corr$seq[which(is.na(PKJK_DK_site_4_test) & !is.na(PKJK_DK_site_4_train))]
tmp_train<-data.frame(list("seq"=PKJK_DK_site_4_train_seq, "value"=PKJK_DK_site_4_train_vals, "label"=rep("train", length(PKJK_DK_site_4_train_vals))))

PKJK_DK_site_4_test_vals<-PKJK_DK_site_4_test[which(!is.na(PKJK_DK_site_4_test) & is.na(PKJK_DK_site_4_train))]
PKJK_DK_site_4_test_seq<-corr$seq[which(!is.na(PKJK_DK_site_4_test) & is.na(PKJK_DK_site_4_train))]
tmp_test<-data.frame(list("seq"=PKJK_DK_site_4_test_seq, "value"=PKJK_DK_site_4_test_vals, "label"=rep("test", length(PKJK_DK_site_4_test_vals))))

write_csv(rbind(tmp_train,tmp_test), "/Users/moradigd/Documents/Plasmid/input_agg_sites_PKJK_DK_site_4.csv")


PKJK_DK_site_1B_train_vals<-PKJK_DK_site_1B_train[which(is.na(PKJK_DK_site_1B_test) & !is.na(PKJK_DK_site_1B_train))]
PKJK_DK_site_1B_train_seq<-corr$seq[which(is.na(PKJK_DK_site_1B_test) & !is.na(PKJK_DK_site_1B_train))]
tmp_train<-data.frame(list("seq"=PKJK_DK_site_1B_train_seq, "value"=PKJK_DK_site_1B_train_vals, "label"=rep("train", length(PKJK_DK_site_1B_train_vals))))

PKJK_DK_site_1B_test_vals<-PKJK_DK_site_1B_test[which(!is.na(PKJK_DK_site_1B_test) & is.na(PKJK_DK_site_1B_train))]
PKJK_DK_site_1B_test_seq<-corr$seq[which(!is.na(PKJK_DK_site_1B_test) & is.na(PKJK_DK_site_1B_train))]
tmp_test<-data.frame(list("seq"=PKJK_DK_site_1B_test_seq, "value"=PKJK_DK_site_1B_test_vals, "label"=rep("test", length(PKJK_DK_site_1B_test_vals))))

write_csv(rbind(tmp_train,tmp_test), "/Users/moradigd/Documents/Plasmid/input_agg_sites_PKJK_DK_site_1B.csv")




#test train site 


write_csv(df, "/Users/moradigd/Documents/Plasmid/input_agg_sites_DK_test.csv")

#Error prediction 

corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

corr_short<-corr[which(is.na(corr[,6]) |  is.na(corr[,7])),] 
corr_short<-corr_short[,-c(6,7)]
PKJK<-apply(corr_short[, which(grepl("_PKJK",colnames(corr_short)) )], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr_short$seq,"PKJK"=PKJK))
df<-df[complete.cases(df),]
df$label<-"train"

df_test_1<-corr[which(!is.na(corr[,6]) & !is.na(corr[,7])),c(2,6)]
df_test_1$label<-"test"
colnames(df_test_1)<-colnames(df)
df_tot_1<-rbind(df,df_test_1)
write_csv(df_tot_1,"/Users/moradigd/Documents/Plasmid/PKJK_replicate_1_Rep1.csv")

df_test_2<-corr[which(!is.na(corr[,6]) & !is.na(corr[,7])),c(2,7)]
df_test_2$label<-"test"
colnames(df_test_2)<-colnames(df)
df_tot_2<-rbind(df,df_test_2)
write_csv(df_tot_2,"/Users/moradigd/Documents/Plasmid/PKJK_replicate_1_Rep2.csv")

corr_short<-corr[which(is.na(corr[,13]) |  is.na(corr[,14])),] 
corr_short<-corr_short[,-c(13,14)]
PKJK<-apply(corr_short[, which(grepl("_PKJK",colnames(corr_short)) )], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr_short$seq,"PKJK"=PKJK))
df<-df[complete.cases(df),]
df$label<-"train"

df_test_1<-corr[which(!is.na(corr[,13]) & !is.na(corr[,14])   ),c(2,13)]
df_test_1$label<-"test"
colnames(df_test_1)<-colnames(df)
df_tot_1<-rbind(df,df_test_1)
write_csv(df_tot_1,"/Users/moradigd/Documents/Plasmid/PKJK_replicate_2_Rep1.csv")

df_test_2<-corr[which(!is.na(corr[,13])  & !is.na(corr[,14])),c(2,14)]
df_test_2$label<-"test"
colnames(df_test_2)<-colnames(df)
df_tot_2<-rbind(df,df_test_2)
write_csv(df_tot_2,"/Users/moradigd/Documents/Plasmid/PKJK_replicate_2_Rep2.csv")


#Figure 1
library(gmodels)
library(tidyverse)

iteration<-c(5,8)
plasmids<-c("RP","PKJK","PB")
#New_results_screening_agg_scale_8_PB.csv
output<-c()
for(i in seq(length(iteration))){
  for(j in seq(length(plasmids))){
    print(plasmids[j])
    tmp_file<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/New_results_screening_agg_scale_",iteration[i],"_", plasmids[j],".csv"))
    output<-rbind(output, tmp_file)
  }
}

colnames(output)<-c("index","correlation","pvalue","kmer","plasmid","model","run")

output<-output %>%
  drop_na() %>%
  group_by(  kmer, plasmid ,model) %>%
  summarise(mean_correlation=mean(correlation) , lowCI = ci(correlation)[2],
            hiCI = ci(correlation)[3], sd_correlation= ci(correlation)[4]) %>%
  mutate(lowCI = replace(lowCI, which(lowCI<0), 0))


output%>%
  ggplot(aes(fill=model, x=plasmid,y=mean_correlation))+
  geom_bar( stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin= output$lowCI, ymax=output$hiCI), width=.2,
                position=position_dodge(.9)) +
  facet_grid(~ kmer)+
  theme_bw()+
  ylab("Spearman's ρ")+
  xlab("Concentration")+
  scale_fill_manual(values=c("#998899", "#E69F00","#FF0000"))

#Figure 1 new 
#Figure 1
library(gmodels)
library(tidyverse)

output<-read_csv("/Users/moradigd/Documents/Plasmid/new_results_screening_agg_scale_tot.csv")
colnames(output)<-c("index","correlation","pvalue","kmer","plasmid","model","run","dataset")
output$plasmid[output$plasmid=="PB"]<-"pB10"
output$plasmid[output$plasmid=="PKJK"]<-"pKJK5"
output$plasmid[output$plasmid=="RP"]<-"RP4"

output<-output %>%
  drop_na() %>%
  dplyr::group_by(  kmer, plasmid ,model,dataset) %>%
  dplyr::summarise(mean_correlation=mean(correlation) , lowCI = ci(correlation)[2],
                   hiCI = ci(correlation)[3], sd_correlation= ci(correlation)[4]) %>%
  mutate(lowCI = replace(lowCI, which(lowCI<0), 0))


output%>%
  filter(dataset=="test") %>%
  ggplot(aes(fill=model, x=plasmid,y=mean_correlation))+
  geom_bar( stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin= lowCI, ymax=hiCI), width=.2,
                position=position_dodge(.9)) +
  facet_grid(~ kmer)+
  theme_bw()+
  ylab("Spearman's ρ")+
  xlab("Plasmid")+
  ylim(range(0,1))+
  scale_fill_manual(values=c("#998899", "#E69F00","#FF0000"))+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))




#+
#  facet_grid(~ Drug.y)+
#  theme_bw() +
#  theme(axis.text.x = element_text( size=18, angle=45,hjust = 1),
#        axis.text.y = element_text( size=18, hjust = 1),
#        axis.title.x = element_text(color="black", size=20, face="bold"),
#        axis.title.y = element_text(color="black", size=20, face="bold"),
#        strip.text.x = element_text(
#          size = 15, color = "black", face = "bold"
#        )
#  )+
#  ylab("MAE")+
#  xlab("Concentration")+ 
#  ylim(range(0,1.5))+
#ylim(range(0,0.4))+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#  scale_fill_manual(values=c("#998899", "#E69F00"))


#Simulated data
library(tidyverse)

corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")
PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PKJK<-PKJK[complete.cases(PKJK)]

files<-dir("/Users/moradigd/Documents/Plasmid/SimBac-master", pattern = "^table_g5_")
files<-files[which(grepl("_2000.fasta.csv",files ))]

for(j in 1:length(files)){
  print(files[j])
  simulated<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/SimBac-master/",files[j]))
  simulated_median<-apply(simulated, 2, function(x) length(which(x==0))/as.numeric(strsplit(gsub(".fasta.csv","",files[j]),"_")[[1]][5]))
  
  alphas<-c(-100,100,1000000)
  gene_ids<-simulated[,which(simulated_median> 0.05 & simulated_median<0.95)]
  
  df_dist<-c()
  for(alpha in alphas){
    for(id in 1:min(dim(gene_ids)[2],10)){
      print(min(dim(gene_ids)[2],10))
      zeros<-which(gene_ids[,id] == 0)
      ones<-which(gene_ids[,id] != 0)
      absence<-sample(PKJK[PKJK<(median(PKJK, na.rm=T)-(range(PKJK)[2]-range(PKJK)[1])/alpha)],length(zeros),replace = T)
      presence<-sample(PKJK[PKJK>(median(PKJK, na.rm=T)-(range(PKJK)[2]-range(PKJK)[1])/alpha)],length(ones),replace = T)
      output_df<-data.frame(list(ids=c(zeros,ones ), values=c(absence, presence))) 
      write_csv(output_df, paste0("/Users/moradigd/Documents/Plasmid/XX_new_simulation_input_",alpha,"_",id,"_",files[j],".csv" ))
      
      df_dist<-rbind(df_dist, data.frame(list( label=c(rep("ab", length(absence)), rep("pr", length(presence))), value=c(absence, presence), alpha=rep( paste0("alpha_",alpha) , length(c(absence, presence))))))
    }
  }
  
  zeros<-which(gene_ids[,id] == 0)
  ones<-which(gene_ids[,id] != 0)
  absence<-sample(PKJK[PKJK<(median(PKJK, na.rm=T))],length(zeros),replace = T)
  presence<-sample(PKJK[PKJK>(median(PKJK, na.rm=T))],length(ones),replace = T)
  df_dist<-rbind(df_dist, data.frame(list( label=c(rep("ab", length(absence)), rep("pr", length(presence))), value=c(absence, presence), alpha=rep( paste0("alpha_or") , length(c(absence, presence))))))
  
  
  
}

ggplot(df_dist, aes(x=alpha, y=value, fill=label)) +
  #geom_violin()+
  geom_boxplot( color="grey") +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  ylim(range(0,50))

#plotting 
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

input<-c()
#for(i in 1:length(alphas)  ){
for(i in 1:3){
  input_tmp<-read_csv( paste0("/Users/moradigd/Documents/Plasmid/XX_new_simulation_input_",alphas[i], "_2_table_g5_genomes_1_2000.fasta.csv.csv"))
  input_tmp$tag<- paste0("tag_", i) 
  input<-rbind(input,input_tmp)
}
median(input$values[input$tag=="Original"])
median(input$values[input$tag=="tag_1"])
median(input$values[input$tag=="tag_2"])
median(input$values[input$tag=="tag_3"])


input_org<-data.frame(list(ids=seq(length(PKJK)), values=PKJK, tag="Original"))
input<-rbind(input, input_org)
hist(input$values)

ggplot(input, aes(x=tag, y=values, fill=tag)) +
  #geom_violin()+
  geom_boxplot( color="grey") +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  ylim(range(0,50))








#Chemical genomics clustring results 
library(solitude)
library(tidyverse)
filenames <- list.files("/Users/moradigd/Documents/chemicalgenomic/errors", pattern="*.jpg", full.names=TRUE)
filenames<-gsub("/Users/moradigd/Documents/chemicalgenomic/errors/image_output_","",filenames )
filenames<-gsub(".jpg","",filenames )

tmp<-read_csv("/Users/moradigd/Documents/chemicalgenomic/similarity_morphology.csv")
tmp<-tmp[,-1]

tmp<-tmp[which(!colnames(tmp)  %in% filenames), which(!colnames(tmp)  %in% filenames)]
tmp<-as.matrix(tmp)

tmp<-log(tmp)
head(tmp[1:10,])
tmp[lower.tri(tmp)]<-tmp[upper.tri(tmp)]
tmp<-1/tmp
diag(tmp) <- 0
head(tmp[,1:10])
write_csv(data.frame(tmp), "/Users/moradigd/Documents/chemicalgenomic/similarity_morphology_new.csv" )

hist(tmp[lower.tri(tmp)], breaks = 100)
range(tmp[lower.tri(tmp)])

labels<-read_csv("/Users/moradigd/Documents/chemicalgenomic/anomaly_labels.csv")

library(edgeR)
keg<-kegga(labels$names[labels$anomaly==-1], species.KEGG="pau")
go<-goana(labels$names[labels$anomaly==-1], species.KEGG="pau")
topKEGG(keg)
topGO(go)

keg<-kegga(filenames, species.KEGG="pau")
go<-goana(filenames, species.KEGG="pau")
topKEGG(keg)
topGO(go)


#contamination range screening 
tmp_0<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/anomaly_morphology_labels",0,".csv"))
cont_mat<-matrix(0, nrow = dim(tmp_0)[1], ncol = 9)
cont_mat[which(tmp_0$anomaly==-1),1]<-1

for(i in 2:9){
  tmp_1<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/anomaly_morphology_labels",i-1,".csv"),show_col_types = FALSE)
  tmp_2<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/anomaly_morphology_labels",i,".csv"),show_col_types = FALSE)
  cont_mat[which((tmp_2$anomaly == -1) & (tmp_1$anomaly == 1)),(i)]<-1
  print(table(tmp_2$anomaly, tmp_1$anomaly))
}





apply(cont_mat, 2, sum)
sum(cont_mat[which(tmp_0$anomaly==-1),1])

library(gmodels)
abnormality<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/anomaly_colony_size_labels",0,".csv"))
dist_matrix<-read_csv("/Users/moradigd/Documents/chemicalgenomic/similarity_colony_size_matrix.csv")
abnormal_tags<-abnormality$names[abnormality$anomaly==-1]

system(paste0("mkdir /Users/moradigd/Documents/chemicalgenomic/abnormal_group_",0))
for(i in 1:length(abnormal_tags)){
  system(paste0("cp /Users/moradigd/Documents/chemicalgenomic/morphology/image_output_",abnormal_tags[i],".jpg /Users/moradigd/Documents/chemicalgenomic/abnormal_group_",0))
}

df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
for(i in 1:9){
  system(paste0("mkdir /Users/moradigd/Documents/chemicalgenomic/size_group_",i))
}
for(i in 1:length(df_colour$isolate)){
  system(paste0("cp /Users/moradigd/Documents/chemicalgenomic/morphology/image_output_",df_colour$isolate[i], ".jpg /Users/moradigd/Documents/chemicalgenomic/size_group_",df_colour$cluster[i]))
}
#dist_matrix_red<-dist_matrix[match(abnormal_tags,   colnames(dist_matrix) ),match(abnormal_tags,   colnames(dist_matrix) )]

#hist(as.matrix(data.frame(dist_matrix)), breaks = 200 )
#kmeans.re <- kmeans(as.matrix(data.frame(dist_matrix_red)), centers = 2, nstart = 20,iter.max = 100)
#matrix_dist_red<- as.matrix(data.frame(dist_matrix_red))
#matrix_tmp_1<-matrix_dist_red[which(kmeans.re$cluster==1),which(kmeans.re$cluster==1)]
#matrix_tmp_2<-matrix_dist_red[which(kmeans.re$cluster==2),which(kmeans.re$cluster==2)]
#hist(matrix_tmp_1[lower.tri(matrix_tmp_1, diag = FALSE)], breaks = 200)
#hist(matrix_tmp_2[lower.tri(matrix_tmp_2, diag = FALSE)], breaks = 200)
#mean(matrix_tmp_1[lower.tri(matrix_tmp_1, diag = FALSE)])
#mean(matrix_tmp_2[lower.tri(matrix_tmp_2, diag = FALSE)])

growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints_summarised.csv")
growth$cluster<- as.numeric(abnormality$anomaly[match(growth$label, abnormality$names )])+2
growth<-growth[complete.cases(growth),]
colnames(growth)
growth_short<-growth %>% group_by(timepoint,cluster) %>% 
  summarise(importance_average=mean(colony.size), importance_low=ci(colony.size)[2], importance_high=ci(colony.size)[3]) 

growth_short<-growth %>% group_by(timepoint,cluster) %>% 
  summarise(importance_average=mean(colony.size), importance_low=min(colony.size), importance_high=max(colony.size)) 


growth_short$timepointsequence<-array(sapply(seq(1,6), function(x) rep(x,2)))
growth_short$cluster<-as.character(growth_short$cluster)

ggplot(growth_short, aes(x =  timepointsequence, y = importance_average, color = cluster)) +
  geom_ribbon(aes(ymin=importance_low,ymax=importance_high),fill="grey80")+
  geom_line(linetype = 9,lwd = 1.1)+
  theme_bw()+
  ylab("Morphology Score")+
  xlab("Time Pionts")+ 
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))







concat_df<-c()
for(i in 1:9){
  print(i)
  keg<-kegga(tmp_1$names[which(cont_mat[,i]==1)], species.KEGG="pau")
  print(head(tmp_1$names[which(cont_mat[,i]==1)]))
  top_tmp<-topKEGG(keg) %>% filter(P.DE<0.05) %>% mutate(contamination=paste("contamination",i))
  concat_df<-rbind(concat_df, top_tmp)
  print(top_tmp)
}

for(j in 1:9){
  print(j)
  abnormal_tags<-tmp_1$names[which(cont_mat[,j]==1)]
  system(paste0("mkdir /Users/moradigd/Documents/chemicalgenomic/abnormal_morphology_group_",j))
  for(i in 1:length(abnormal_tags)){
    system(paste0("cp /Users/moradigd/Documents/chemicalgenomic/morphology/image_output_",abnormal_tags[i],".jpg /Users/moradigd/Documents/chemicalgenomic/abnormal_morphology_group_",j))
  }
}

#copy images 
df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")

for(j in 1:8){
  print(j)
  tags<-  df_colour$isolate[which(df_colour$cluster==j)]
  system(paste0("mkdir /Users/moradigd/Documents/chemicalgenomic/morphology_group_",j))
  for(i in 1:length(tags)){
    system(paste0("cp /Users/moradigd/Documents/chemicalgenomic/morphology/image_output_",tags[i],".jpg /Users/moradigd/Documents/chemicalgenomic/morphology_group_",j))
  }
}

#dictribution of clusters aand roups 

df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")

shared<-intersect(intersect(df_morphology$isolate,df_size$isolate), df_colour$isolate )
df_morphology<-df_morphology[match(shared, df_morphology$isolate), ]
df_size<-df_size[match(shared, df_size$isolate), ]
df_colour<-df_colour[match(shared, df_colour$isolate), ]

table(df_size$cluster[which(df_morphology$cluster==6)],  df_colour$cluster[which(df_morphology$cluster==6)] )
table(df_size$cluster[which(df_morphology$cluster==6)])/sum(table(df_size$cluster[which(df_morphology$cluster==6)]))
table(df_colour$cluster[which(df_morphology$cluster==6)])
table(df_colour$cluster[which(df_morphology$cluster==6)])/sum(table(df_colour$cluster[which(df_morphology$cluster==6)]))


table(df_size$cluster[which(df_morphology$cluster==8)],  df_colour$cluster[which(df_morphology$cluster==8)] )
table(df_size$cluster[which(df_morphology$cluster==8)])/sum(table(df_size$cluster[which(df_morphology$cluster==8)]))
table(df_colour$cluster[which(df_morphology$cluster==8)])/sum(table(df_colour$cluster[which(df_morphology$cluster==8)]))

table(df_size$cluster)/sum(table(df_size$cluster))
table(df_colour$cluster)/sum(table(df_colour$cluster))
par(mar=c(15, 5, 2, 4))
barplot(concat_df$N, col = cols[sapply(concat_df$contamination, function(x) as.numeric(strsplit(x, "")[[1]][15])) ], las=2, names.arg = concat_df$Pathway, cex.names = 0.7 , ylab="Frequency")
box()


for(i in 1:10){
  tmp_1<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/anomaly_labels",i-1,".csv"),show_col_types = FALSE)
  cont_mat[which((tmp_2$anomaly == -1) & (tmp_1$anomaly == 1)),(i)]<-1
}


#Chemical genomics 
library(usethis) 
usethis::edit_r_environ()


library(tidyverse)
library(mclust)
library(randomcoloR)
library(solitude)

filenames <- list.files("/Users/moradigd/Documents/chemicalgenomic/errors", pattern="*.jpg", full.names=TRUE)
filenames<-gsub("/Users/moradigd/Documents/chemicalgenomic/errors/image_output_","",filenames )
filenames<-gsub(".jpg","",filenames )

distanec_matrix<-read_csv( "/Users/moradigd/Documents/chemicalgenomic/morphology_distance_matrix.csv")
distanec_matrix<-distanec_matrix[which(!colnames(distanec_matrix)  %in% filenames), which(!colnames(distanec_matrix)  %in% filenames)]
distanec_matrix<-as.matrix(distanec_matrix)
distanec_matrix<-distanec_matrix[which(apply(distanec_matrix, 2, var) != 0) , which(apply(distanec_matrix, 2, var) != 0)]
#write_csv(data.frame(distanec_matrix), "/Users/moradigd/Documents/chemicalgenomic/morphology_distance_mat_shortened.csv")
#distanec_matrix<-read_csv( "/Users/moradigd/Documents/chemicalgenomic/morphology_distance_mat_shortened.csv")


BIC <- mclustBIC(distanec_matrix,verbose = interactive(),G=1:9)
plot(BIC, las=1)
summary(BIC)
BIC[2,2]

for(i in 2:dim(BIC)[1]){
  for(j in 1:dim(BIC)[2]){
    print((-BIC[i,j]+BIC[i-1,j])/BIC[i-1,j])
  }
}
#?mclustBIC

model <- Mclust(distanec_matrix, x = BIC)
#summary(model , parameters = TRUE)[[9]]

df<-data.frame(list(isolate=colnames(distanec_matrix), cluster=model$classification)) 
#write_csv(df,"/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_size_new_5.csv" )

require(Rtsne)
clustering<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")

shared<-intersect( colnames(distanec_matrix), clustering$isolate    )
clustering<-clustering[match(shared, clustering$isolate),]
distanec_matrix<-distanec_matrix[match(shared,colnames(distanec_matrix)),match(shared,colnames(distanec_matrix))]

tsne <- Rtsne(distanec_matrix, check_duplicates = FALSE, pca = T, perplexity=100, theta=0.5, dims=2)

cols <- rainbow(9)
par(mar=c(4, 4, 2,2))
plot(tsne$Y,col= alpha(cols[clustering$cluster], 0.55) , pch = 16,alpha = 0.5, xlab="t-SNE Dim(1)", ylab="t-SNE Dim(2)")

#clsuter distances 
library(clv)
kid<-cls.scatt.data(distanec_matrix, clustering$cluster)

#box plot 
inp1<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_colour_annotation.csv")
inp2<-read_csv("/Users/moradigd/Documents/chemicalgenomic/property_colony_colour_tot.csv")

boxplot( inp2$abnormality[match(inp1$isolate, inp2$Gene)] ~ inp1$cluster ) 
#
#size
#mod1$cluster[mod1$cluster==9]<-8
#mod1$cluster[mod1$cluster==6]<-5
#mod1$cluster[mod1$cluster==4]<-2
# display the results of t-SNE
?Rtsne
#res.km <- kmeans(tmp_red, 10, nstart = 50,iter.max = 15)
# K-means clusters showing the group of each individuals
#res.km$cluster

#text(tsne$Y,  col=cols[mod1$cluster])

#silouhte 
require(cluster)
library(tidyverse)
clustering<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_1.csv")

filenames <- list.files("/Users/moradigd/Documents/chemicalgenomic/errors", pattern="*.jpg", full.names=TRUE)
filenames<-gsub("/Users/moradigd/Documents/chemicalgenomic/errors/image_output_","",filenames )
filenames<-gsub(".jpg","",filenames )

distance_matrix<-read_csv( "/Users/moradigd/Documents/chemicalgenomic/morphology_distance_matrix.csv")


#for colour
#distance_matrix<-distance_matrix[,-1]


#for(i in 1:dim(distance_matrix)[1]){
#  print(i)
#  for(j in i:dim(distance_matrix)[1]){
#    distance_matrix[j,i]<-distance_matrix[i,j]
#  }
#}
#write_csv(distance_matrix,  "/Users/moradigd/Documents/chemicalgenomic/colour_distance_matrix.csv")
distance_matrix<-distance_matrix[which(!colnames(distance_matrix)  %in% filenames), which(!colnames(distance_matrix)  %in% filenames)]
distance_matrix<-as.matrix(distance_matrix)
distance_matrix<-distance_matrix[which(apply(distance_matrix, 2, var) != 0) , which(apply(distance_matrix, 2, var) != 0)]

shared<-intersect(colnames(distance_matrix),clustering$isolate )
clustering<-clustering[match(shared, clustering$isolate),]
distance_matrix<-distance_matrix[match(shared, colnames(distance_matrix)),match(shared, colnames(distance_matrix))]
D <- daisy(distance_matrix)
sil_results<- silhouette(clustering$cluster, distance_matrix)
plot(silhouette(clustering$cluster, D), col=1:9,border=NA)
summary(sil_results)

cluster_SIL <- t(rbind(summary(sil_results)[["clus.sizes"]], summary(sil_results)[["clus.avg.widths"]]))
colnames(cluster_SIL) <- c("No. of Obs", "Avg. Silh. Width")
#write_csv(data.frame(cluster_SIL), "/Users/moradigd/Documents/chemicalgenomic/morphology/sum_sil_colony_colour.csv")

cluster_full<-as.data.frame(cbind(gene=shared,cluster=sil_results[,1],neighbor=sil_results[,2],sil_width=sil_results[,3]))
#write_csv(cluster_full, "/Users/moradigd/Documents/chemicalgenomic/morphology/full_sil_colony_colour.csv")

#CIvalid
library(clValid)
intern <- clValid(distance_matrix, 2:9, clMethods=c("hierarchical","kmeans","pam"),validation="internal")

summary(intern)


require(Rtsne)
tsne <- Rtsne(distance_matrix, check_duplicates = FALSE, pca = T, perplexity=100, theta=0.5, dims=2)

cols <- rainbow(9)
par(mar=c(4, 4, 2,2))
plot(tsne$Y,col= alpha(cols[clustering$cluster], 0.55) , pch = 16,alpha = 0.5, xlab="t-SNE Dim(1)", ylab="t-SNE Dim(2)")


#tSNE abnormality 
tsne_cluster<-read_csv("/Users/moradigd/Documents/chemicalgenomic/anomaly_colony_size_labels9.csv")
tsne_cluster<-tsne_cluster[match(colnames(tmp_red) ,tsne_cluster$names ),]
tsne_cluster$anomaly<-ifelse(tsne_cluster$anomaly==1, 2, 1)

cols <- c("red","blue")
par(mar=c(4, 4, 2,2))
plot(tsne$Y,col= alpha(cols[tsne_cluster$anomaly], 0.55) , pch = 16,alpha = 0.5, xlab="t-SNE Dim(1)", ylab="t-SNE Dim(2)")


#COG Analysis 
library(edgeR)
library(gmodels)
library(tidyverse)

df_circularity<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv")
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_1.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")

#annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/annotation.csv")
#df_colour$annotation<-annotation$`"Active" Gene Description`[match(df_colour$isolate, annotation$`"Active" Gene Locus`)]
#write_csv(df_colour, "/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_colour_annotation.csv" )

df<-df_size
output_KEGG<-c()
output_go<-c()
for(i in 1:9){
  print(i)
  label_tmp<-df$isolate[df$cluster==i]
  keg<-kegga(label_tmp, species.KEGG="pau")
  go<-goana(label_tmp, species.KEGG="pau")
  top_kegg<-topKEGG(keg)
  top_kegg$cluster<-paste0("cluster_", i)
  print(top_kegg[top_kegg$P.DE<0.8,])
  
  output_KEGG<-rbind(output_KEGG,top_kegg[top_kegg$P.DE<0.05,] )
  
  top_go<-topGO(go)
  top_go$cluster<-paste0("cluster_", i)
  output_go<-c()
  print(top_go[top_go$P.DE<0.05,])
  output_go<-rbind(output_go,top_kegg[top_go$P.DE<0.05,] )
}
table(df$cluster)

cols <- rainbow(9)
par(mar=c(15, 5, 2, 4))
barplot( -log10(output_KEGG$P.DE)   , col = cols[match( output_KEGG$cluster, unique(output_KEGG$cluster))], las=2, names.arg = output_KEGG$Pathway, cex.names = 0.7 , ylab="pvalue")
box()

#Runs 1 and 2
df_morphology_1<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_4.csv")
df_morphology_2<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_3.csv")

reordered_id<-apply(table(df_morphology_2$cluster,df_morphology_1$cluster), 2, function(x) which.max(x)   )
df_morphology_reordered<-sapply(df_morphology_1$cluster, function(x) reordered_id[x])

table(df_morphology_2$cluster,df_morphology_1$cluster)

table(df_morphology_2$cluster,df_morphology_reordered)

table(df_morphology_1$cluster[match(df_morphology_2$isolate[df_morphology_2$cluster==7],df_morphology_1$isolate)])

#Robust morphoology 
df_1<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_1.csv")
df_2<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_2.csv")
df_3<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_3.csv")
df_4<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_4.csv")
df_5<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_new_5.csv")

df_tot<-cbind(df_1$cluster, df_2$cluster, df_3$cluster, df_4$cluster, df_5$cluster  )


#GO terms 
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")
annotaion<-read_csv("/Users/moradigd/Documents/chemicalgenomic/annotation.csv")
df_morphology$annotation<-annotaion$`"Active" Gene Description`[match(df_morphology$isolate, annotaion$`"Active" Gene Locus`)]
df_morphology<-df_morphology[df_morphology$cluster==1,]
#df_morphology<-df_morphology[which(grepl("Argi", df_morphology$annotation)),]

for(i in 1:dim(df_morphology)[1]){
  
  if(file.exists(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/KEGG_",df_morphology$isolate[i],".csv"))){
    df<- read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/KEGG_",df_morphology$isolate[i],".csv"),show_col_types = FALSE) %>% filter(P.DE<0.05)
    if( length(df$P.DE)>0){
      print(df_morphology$isolate[i])
      print(head(df))
    }
  }
}

freq<-sort(table(paste0(df_circularity$cluster,"_",df_morphology$cluster, "_", df_size$cluster)))
length(freq)

table(df_circularity$cluster)
table(df_morphology$cluster)
table(df_size$cluster)

growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints_summarised.csv")
growth$cluster<-df$cluster[match(growth$label, df$isolate )]
growth<-growth[complete.cases(growth),]
colnames(growth)
growth_short<-growth %>% group_by(timepoint,cluster) %>% 
  summarise(importance_average=mean(colony.circularity), importance_low=ci(colony.circularity)[2], importance_high=ci(colony.circularity)[3]) 
growth_short$timepointsequence<-array(sapply(seq(1,6), function(x) rep(x,9)))
growth_short$cluster<-as.character(growth_short$cluster)

ggplot(growth_short, aes(x =  timepointsequence, y = importance_average, color = cluster)) +
  geom_ribbon(aes(ymin=importance_low,ymax=importance_high),fill="grey80")+
  geom_line(linetype = 9,lwd = 1.1)+
  theme_bw()+
  ylab("Morphology Score")+
  xlab("Time Pionts")+ 
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))

#colour 
library(edgeR)
library(gmodels)
library(tidyverse)
df<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_colour_annotation.csv")
growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv")
growth$cluster<-df$cluster[match(growth$gene, df$isolate )]
growth<-growth[complete.cases(growth),]
colnames(growth)
growth_short<-growth %>% group_by(time_point,cluster) %>% 
  summarise(importance_average=mean(value), importance_low=ci(value)[2], importance_high=ci(value)[3]) 
growth_short$timepointsequence<-array(sapply(seq(1,6), function(x) rep(x,9)))
growth_short$cluster<-as.character(growth_short$cluster)

ggplot(growth_short, aes(x =  timepointsequence, y = importance_average, color = cluster)) +
  geom_ribbon(aes(ymin=importance_low,ymax=importance_high),fill="grey80")+
  geom_line(linetype = 9,lwd = 1.1)+
  theme_bw()+
  ylab("Colony Colour Score")+
  xlab("Time Pionts")+ 
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))


#Add protein interaction data 
library(tidyverse)
library(edgeR)
df_circularity<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv")
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")

shared<-intersect(intersect(df_circularity$isolate, df_morphology$isolate), df_size$isolate)
df_circularity<-df_circularity[match(shared,df_circularity$isolate),]
df_morphology<-df_morphology[match(shared,df_morphology$isolate),]
df_size<-df_size[match(shared,df_size$isolate),]
df_colour<-df_size[match(shared,df_size$isolate),]

annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/annotation.csv")
genes_orthologous<-annotation$`PAO1 Orthologs`[match(shared, annotation$`"Active" Gene Locus`)]
genes_function<-annotation$`"Active" Gene Description`[match(shared, annotation$`"Active" Gene Locus`)]


df_circularity$orthologous<-genes_orthologous
df_morphology$orthologous<-genes_orthologous
df_size$orthologous<-genes_orthologous
df_colour$orthologous<-genes_orthologous

df_circularity$description<-genes_function
df_morphology$description<-genes_function
df_size$description<-genes_function
df_colour$description<-genes_function

write_csv(df_circularity, "/Users/moradigd/Documents/chemicalgenomic/df_circularity_cluster_annotation.csv")
write_csv(df_morphology, "/Users/moradigd/Documents/chemicalgenomic/df_morphology_cluster_annotation.csv")
write_csv(df_size, "/Users/moradigd/Documents/chemicalgenomic/df_size_cluster_annotation.csv")
write_csv(df_colour, "/Users/moradigd/Documents/chemicalgenomic/df_colour_cluster_annotation.csv")


#Prediction 
library(tidyverse)
library(edgeR)
df_circularity<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv")
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")

shared<-intersect(intersect(df_circularity$isolate, df_morphology$isolate), df_size$isolate)
df_circularity<-df_circularity[match(shared,df_circularity$isolate),]
df_morphology<-df_morphology[match(shared,df_morphology$isolate),]
df_size<-df_size[match(shared,df_size$isolate),]

pattern<-paste0(df_circularity$cluster,"_",df_morphology$cluster, "_", df_size$cluster)

freq<-sort(table(pattern), decreasing = T)
length(freq)

freq[which(freq>1)]
barplot(freq[which(freq>1)], las=2, cex.names = 0.5, ylab = "Freqeucny", xlab="Growth groups", col="black")
box()

output_KEGG<-c()
for(j in 1:length(names(freq))){
  print(j)
  print(names(freq)[j])
  tmp_pattern<-df_circularity$isolate[which(pattern==names(freq)[j])]
  keg<-kegga(tmp_pattern, species.KEGG="pau")
  top_kegg<-topKEGG(keg)
  
  top_kegg$tag<-names(freq)[j]
  print(top_kegg[top_kegg$P.DE<0.05,])
  output_KEGG<-rbind(output_KEGG,top_kegg[top_kegg$P.DE<0.05,] )
}
#write_csv(output_KEGG,"/Users/moradigd/Documents/chemicalgenomic/profile_KEGG_all_genes.csv" )
PA14_34210


annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/annotation.csv")
genes_annotation<-annotation$`"Active" Gene Description`[match(shared, annotation$`"Active" Gene Locus`)]
output_total<-c()
for(i in 1:length(shared)){
  tmp<-output_KEGG_red[which(output_KEGG_red$tag %in% pattern[i]),]
  if(dim(tmp)[1]!=0){
    tmp$gene_annotation<-genes_annotation[i]
    tmp$locus<-shared[i]
    output_total<-rbind(output_total, tmp)
  }
}
#write_csv(output_KEGG,"/Users/moradigd/Documents/chemicalgenomic/profile_KEGG_all_genes_train.csv" )
#write_csv(output_total,"/Users/moradigd/Documents/chemicalgenomic/profile_KEGG_all_genes_all.csv" )

#funxtional analysis 
test_data_id<-sample(seq(length(pattern)),300)
train_data_id<-which(!seq(length(pattern)) %in% test_data_id)

df_circularity_red<-df_circularity[train_data_id,]
df_size_red<-df_size[train_data_id,]
df_morphology_red<-df_morphology[train_data_id,]

pattern_red<-paste0(df_circularity_red$cluster,"_",df_morphology_red$cluster, "_", df_size_red$cluster)

freq_red<-sort(table(pattern_red), decreasing = T)
length(freq)

output_KEGG_red<-c()
for(j in 1:length(names(freq_red))){
  print(j)
  print(names(freq_red)[j])
  tmp_pattern<-df_circularity_red$isolate[which(pattern_red==names(freq_red)[j])]
  keg<-kegga(tmp_pattern, species.KEGG="pau")
  top_kegg<-topKEGG(keg)
  
  top_kegg$tag<-names(freq)[j]
  print(top_kegg[top_kegg$P.DE<0.05,])
  output_KEGG_red<-rbind(output_KEGG_red,top_kegg[top_kegg$P.DE<0.05,] )
}
#write_csv(output_KEGG,"/Users/moradigd/Documents/chemicalgenomic/profile_KEGG_all_genes_train.csv" )

annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/annotation.csv")
genes_annotation<-annotation$`"Active" Gene Description`[match(shared[test_data_id], annotation$`"Active" Gene Locus`)]
output_test<-c()
for(i in 1:300){
  tmp<-output_KEGG_red[which(output_KEGG_red$tag %in% pattern[test_data_id][i]),]
  if(dim(tmp)[1]!=0){
    tmp$gene_annotation<-genes_annotation[i]
    tmp$locus<-shared[test_data_id][i]
    output_test<-rbind(output_test, tmp)
  }
}

output_test<-output_test[complete.cases(output_test),]
write_csv(output_test,"/Users/moradigd/Documents/chemicalgenomic/profile_KEGG_all_genes_test.csv" )

#protein inetraction
library(tidyverse)
library(edgeR)
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")

shared<-intersect(intersect(df_colour$isolate, df_morphology$isolate), df_size$isolate)
df_morphology<-df_morphology[match(shared,df_morphology$isolate),]
df_size<-df_size[match(shared,df_size$isolate),]
df_colour<-df_colour[match(shared,df_colour$isolate),]

unique(df_morphology$cluster)
interaction<-readLines("/Users/moradigd/Documents/chemicalgenomic/PA14_links_score.csv")

genes_id<-df_morphology$isolate
for(i in 1:length(genes_id)){
  print(i)
  tmp<-which(grepl(genes_id[i], interaction, fixed = TRUE))
  write_csv(data.frame(tmp), paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/",genes_id[i],".csv" ))
}


for(k in 1:9){
  print(k)
  genes_id<-df_colour$isolate[df_colour$cluster==k]
  genes_inter<-c()
  for(i in genes_id){
    print(i)
    print(match(i,genes_id )/length(genes_id))
    gene_tmp<-i
    tmp1<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/",i,".csv" ),show_col_types = FALSE)
    
    inter="N"
    for(j in genes_id){
      if(i != j){
        tmp2<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/",j,".csv" ),show_col_types = FALSE)
        if(length(intersect(tmp1$tmp, tmp2$tmp))>0){
          inter="Y"
          break
        }
      }
    }
    genes_inter<-c(genes_inter, inter)
  }
  print(table(genes_inter))
  df_out<-data.frame(list(gene= genes_id, interaction=genes_inter ))
  write_csv(df_out, paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/colour_interaction_",k,".csv" ))
}

df_output<-c()
pvalue<-c()

for(i in 1:9){
  df_tmp<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/colour_interaction_",i,".csv" ), show_col_types = FALSE) %>% count(interaction)
  df_tmp$cluster<-paste0("cluster_",i)
  df_output<-rbind(df_output, df_tmp)
  pvalue<-c(pvalue, -log10(chisq.test(data.frame("cnt"=c(82,118), "sam"=df_tmp$n, row.names = c("N", "Y")))$p.value) )
}
barplot(pvalue, col = "black", ylab="-log10(pvalue)", xlab="cluster")
abline(h=c(1.3,2), col=c("blue", "red"), lty=c(2,2), lwd=c(3, 3))
box()


#baseline distribution 
print(k)
genes_id<-sample(df_colour$isolate,200) 
genes_inter<-c()
for(i in genes_id){
  print(i)
  print(match(i,genes_id )/length(genes_id))
  gene_tmp<-i
  tmp1<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/",i,".csv" ),show_col_types = FALSE)
  
  inter="N"
  for(j in genes_id){
    if(i != j){
      tmp2<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/",j,".csv" ),show_col_types = FALSE)
      if(length(intersect(tmp1$tmp, tmp2$tmp))>0){
        inter="Y"
        break
      }
    }
  }
  genes_inter<-c(genes_inter, inter)
}
print(table(genes_inter))
df_out<-data.frame(list(gene= genes_id, interaction=genes_inter ))
#write_csv(df_out, paste0("/Users/moradigd/Documents/chemicalgenomic/inetraction/colour_interaction_",k,".csv" ))


#1.439024




ggplot(df_output, aes(fill=interaction, y=n, x=cluster)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_bw()+
  ylab("Frequency")+
  xlab("Clusters")+ 
  theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
        axis.text.y = element_text( size=11, hjust = 1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))

#Operon analysis 
library(tidyverse)
library(edgeR)
df_circularity<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv")
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")

shared<-intersect(intersect(df_circularity$isolate, df_morphology$isolate), df_size$isolate)
df_circularity<-df_circularity[match(shared,df_circularity$isolate),]
df_morphology<-df_morphology[match(shared,df_morphology$isolate),]
df_size<-df_size[match(shared,df_size$isolate),]

operon_names<-read_tsv("/Users/moradigd/Documents/chemicalgenomic/operon_PA14.txt")

operon_names$cluster1<-df_morphology$cluster[match( operon_names$SysName1,  df_morphology$isolate)]
operon_names$cluster2<-df_morphology$cluster[match( operon_names$SysName2,  df_morphology$isolate)]
operon_names<-select(operon_names, cluster1, cluster2, bOp, pOp, COGSim)
operon_names<-operon_names[complete.cases(operon_names),]
operon_names_red<-operon_names[operon_names$cluster1==1 & operon_names$cluster2==1, ]

table(operon_names_red$bOp)/sum(table(operon_names_red$bOp))
table(operon_names$bOp)/sum(table(operon_names$bOp))

table(operon_names_red$COGSim)/sum(table(operon_names_red$COGSim))
table(operon_names$COGSim)/sum(table(operon_names$COGSim))

#Expression analysis 
library(tidyverse)
library(edgeR)
df_circularity<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv")
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")

shared<-intersect(intersect(df_circularity$isolate, df_morphology$isolate), df_size$isolate)
df_circularity<-df_circularity[match(shared,df_circularity$isolate),]
df_morphology<-df_morphology[match(shared,df_morphology$isolate),]
df_size<-df_size[match(shared,df_size$isolate),]

proteome<-read_csv("/Users/moradigd/Documents/chemicalgenomic/proteome.csv")
df_morphology
proteome$cluster<-df_morphology$cluster[match(proteome$PA14_ID,df_morphology$isolate)]
proteome<-proteome[complete.cases(proteome),]
boxplot( proteome$`ratio bf vs pl` ~ proteome$cluster )

anomaly<-read_csv("/Users/moradigd/Documents/chemicalgenomic/anomaly_colony_size_labels0.csv")
proteome$cluster<-anomaly$anomaly[match(proteome$PA14_ID,anomaly$names)]
boxplot( proteome$`ratio bf vs pl` ~ proteome$cluster )



#Features associations

library(inflection)
library(nlstools)
library(tidyverse)
sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

start_val_sigmoid <- function(x, y) {
  fit <- lm(log(y[which.max(x)] - y + 1e-6) ~ x)
  list(
    a = y[which.max(x)],
    b = unname(-coef(fit)[2]),
    c = unname(-coef(fit)[1] / coef(fit)[2]))
}


#input_red<-input[which(!input$PA14.ID %in% filenames),]
#write_csv(input_red, "/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv")
length(which( is.na(input$morphology.score.fixed.circles) ))
length(which( is.na(input$colony.size) ))

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv")
unique_loci<-unique(input$PA14.ID)

properties_tot<-c()
length_unique<-c()
for(ko in 1:20){
  properties<-c()
  #  for(j in seq(length(unique_loci))){
  for(j in seq(length(unique_loci))){
    set.seed(j) 
    print(j)
    print(unique_loci[j])
    timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
    x_timepoints<-match(input$timepoint[input$PA14.ID==unique_loci[j]],timepoints )
    x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
    y_growth<-input$morphology.score.fixed.circles[input$PA14.ID==unique_loci[j]]
    print(length(y_growth))
    plot(x_timepoints,y_growth)
    check="T"
    tryCatch({
      fitmodel <- nls(y_growth ~ a / ( 1 + exp(-b * (x_timepoints - c))),start=start_val_sigmoid(x_timepoints, y_growth), algorithm = "port")
      params=coef(fitmodel)
      check="F"
      
      x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
      y_growth_predicted <- sigmoid(params,x_timepoint_range)
      
      inflection_point_slope<-y_growth_predicted[ which(abs(diff(y_growth_predicted))==max(abs(diff(y_growth_predicted))) )]
      
      inflection_point=bese(x_timepoint_range,y_growth_predicted,0)
      
      data <- data.frame(
        x = x_timepoints,
        y = y_growth
      )
      
      model = nls(y ~ SSlogis(x, a, b, c), data = data)
      
      rate_increase<-(-median(y_growth[1:4])+median(y_growth[(length(y_growth)-3):length(y_growth)]))/5
      properties<-rbind(properties,c(inflection_point_slope,params,as.vector(confint2(fitmodel)),inflection_point$iplast,  summary(model)$coefficients[1],   summary(model)$coefficients[2] ,    summary(model)$coefficients[3] , model$m$deviance() , rate_increase ,unique_loci[j]))
      
    }, error = function(e) {
      
      
      an.error.occured <<- TRUE},
    finally = {
      if(check=="T"){
        rate_increase<-(-median(y_growth[1:4])+median(y_growth[(length(y_growth)-3):length(y_growth)]))/5
        properties<-rbind(properties,c(rep(NA, 15) , rate_increase ,unique_loci[j]))
      }
    }
    )
  }
  properties<-data.frame(properties)
  properties$round<-ko
  unique_loci<-as.character(properties$V17[which(is.na(properties$V15))]) 
  properties_tot<-rbind(properties_tot, properties)
  print(length(unique_loci))
  length_unique<-c(length_unique,length(unique_loci))
  
}

properties_tot_modified<-c()
unique_genes<-unique(properties_tot$V17)
for(i in 1:length(unique_genes)){
  tmp<-as.character(unique_genes[i])
  if(length(which(!is.na( properties_tot$V15[properties_tot$V17==tmp]   )))>0){
    properties_tot_modified<-rbind(properties_tot_modified, properties_tot[properties_tot$V17==tmp,][which(!is.na( properties_tot$V15[properties_tot$V17==tmp]   ))[1],])
  }
  else{
    properties_tot_modified<-rbind(properties_tot_modified,properties_tot[properties_tot$V17==tmp,][1,])
  }
}

colnames(properties_tot_modified)<-c("inflection_point_slope", "a","b","c","a_min","b_min","c_min","a_max","b_max","c_max","inflection","a_curve", "b_curve", "c_curve", "R_squared" ,"growth rate",   "locus id")
for(j in 1:16){
  properties_tot_modified[,j]<-as.numeric(as.character(properties_tot_modified[,j]))
}

hist(properties$inflection_point_slope)

write_csv(properties_tot_modified, "/Users/moradigd/Documents/chemicalgenomic/new_new_Growth_properties_morphology.csv")

#color
library(tidyverse)
input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv")
unique_loci<-unique(input$gene)

#input_red<-input[which(!input$PA14.ID %in% filenames),]
#write_csv(input_red, "/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv")
#input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv")
#unique_loci<-unique(input$PA14.ID)

properties_tot<-c()
length_unique<-c()
for(ko in 1:20){
  properties<-c()
  #  for(j in seq(length(unique_loci))){
  for(j in seq(length(unique_loci))){
    set.seed(j) 
    print(j)
    print(unique_loci[j])
    timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
    x_timepoints<-match(input$time_point[input$gene==unique_loci[j]],timepoints )
    x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
    y_growth<-input$value[input$gene==unique_loci[j]]
    plot(x_timepoints,y_growth)
    check="T"
    tryCatch({
      
      fitmodel <- nls(y_growth ~ a / ( 1 + exp(-b * (x_timepoints - c))),start=start_val_sigmoid(x_timepoints, y_growth), algorithm = "port")
      params=coef(fitmodel)
      check="F"
      
      x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
      y_growth_predicted <- sigmoid(params,x_timepoint_range)
      
      inflection_point_slope<-y_growth_predicted[ which(abs(diff(y_growth_predicted))==max(abs(diff(y_growth_predicted))) )]
      
      inflection_point=bese(x_timepoint_range,y_growth_predicted,0)
      
      data <- data.frame(
        x = x_timepoints,
        y = y_growth
      )
      
      model = nls(y ~ SSlogis(x, a, b, c), data = data)
      #rate_increase<-(-median(y_growth[1:4])+median(y_growth[21:24]))/5
      rate_increase<-(-median(y_growth[1:4])+median(y_growth[(length(y_growth)-3):length(y_growth)]))/5
      properties<-rbind(properties,c(inflection_point_slope,params,as.vector(confint2(fitmodel)),inflection_point$iplast,  summary(model)$coefficients[1],   summary(model)$coefficients[2] ,    summary(model)$coefficients[3] , model$m$deviance() , rate_increase ,unique_loci[j]))
      
    }, error = function(e) {
      
      
      an.error.occured <<- TRUE},
    finally = {
      if(check=="T"){
        rate_increase<-(-median(y_growth[1:4])+median(y_growth[(length(y_growth)-3):length(y_growth)]))/5
        
        # rate_increase<-(-median(y_growth[1:4])+median(y_growth[21:24]))/5
        properties<-rbind(properties,c(rep(NA, 15) , rate_increase ,unique_loci[j]))
      }
    }
    )
  }
  properties<-data.frame(properties)
  properties$round<-ko
  unique_loci<-as.character(properties$V17[which(is.na(properties$V15))]) 
  properties_tot<-rbind(properties_tot, properties)
  print(length(unique_loci))
  length_unique<-c(length_unique,length(unique_loci))
}




properties_tot_modified<-c()
unique_genes<-unique(properties_tot$V17)
for(i in 1:length(unique_genes)){
  tmp<-as.character(unique_genes[i])
  if(length(which(!is.na( properties_tot$V15[properties_tot$V17==tmp]   )))>0){
    properties_tot_modified<-rbind(properties_tot_modified, properties_tot[properties_tot$V17==tmp,][which(!is.na( properties_tot$V15[properties_tot$V17==tmp]   ))[1],])
  }
  else{
    properties_tot_modified<-rbind(properties_tot_modified,properties_tot[properties_tot$V17==tmp,][1,])
  }
}
properties_tot_modified<-data.frame(properties_tot_modified)
colnames(properties_tot_modified)<-c("inflection_point_slope", "a","b","c","a_min","b_min","c_min","a_max","b_max","c_max","inflection","a_curve", "b_curve", "c_curve", "R_squared" ,"growth rate",   "locus id")
for(j in 1:16){
  properties_tot_modified[,j]<-as.numeric(as.character(properties_tot_modified[,j]))
}

write_csv(properties_tot_modified, "/Users/moradigd/Documents/chemicalgenomic/new_new_Growth_properties_colony_colour.csv")


#Cicularity 
properties<-c()
for(j in seq(length(unique_loci))){
  set.seed(j) 
  print(j)
  print(unique_loci[j])
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$timepoint[input$PA14.ID==unique_loci[j]],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$colony.circularity[input$PA14.ID==unique_loci[j]]
  plot(x_timepoints,y_growth)
  tryCatch({
    #fitmodel <- nls(y_growth ~ a / ( 1 + exp(-b * (x_timepoints - c))),start=start_val_sigmoid(x_timepoints, y_growth), algorithm = "port")
    #params=coef(fitmodel)
    
    #x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
    #y_growth_predicted <- sigmoid(params,x_timepoint_range)
    
    #inflection_point_slope<-y_growth_predicted[ which(abs(diff(y_growth_predicted))==max(abs(diff(y_growth_predicted))) )]
    
    #inflection_point=bese(x_timepoint_range,y_growth_predicted,0)
    
    #data <- data.frame(
    #  x = x_timepoints,
    #  y = y_growth
    #)
    
    #model = nls(y ~ SSlogis(x, a, b, c), data = data)
    
    overal_rate_increase<-(-median(y_growth[1:4])+median(y_growth[21:24]))/5
    first_rate_increase<-(-median(y_growth[1:4])+median(y_growth[5:8]))/5
    first_rate_decline<-(-median(y_growth[9:12])+median(y_growth[5:8]))/5
    second_rate_increase<-(-median(y_growth[9:12])+median(y_growth[13:16]))/5
    variance_total<-(var(y_growth[9:12])+var(y_growth[1:4])+var(y_growth[5:8])+var(y_growth[13:16])+var(y_growth[17:20])+var(y_growth[21:24]))
    
    properties<-rbind(properties,c(overal_rate_increase, first_rate_increase, first_rate_decline, second_rate_increase, variance_total ,unique_loci[j]))
    
  }, error = function(e) {an.error.occured <<- TRUE}
  )
}
properties<-data.frame(properties)
colnames(properties)<-c("overal_increase", "first_rate_increase","first_rate_decline","second_rate_increase","variance_total", "locus id")
for(j in 1:5){
  properties[,j]<-as.numeric(as.character(properties[,j]))
}

hist(properties$inflection_point_slope)

#write_csv(properties, "/Users/moradigd/Documents/chemicalgenomic/Growth_properties_circularity.csv")

#Properties and features 
properties<-read_csv("/Users/moradigd/Documents/chemicalgenomic/Growth_properties_colony_size.csv")


plotting<-function(j){
  library(tidyverse)
  input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv",show_col_types = FALSE)
  unique_loci<-unique(input$PA14.ID)
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$timepoint[input$PA14.ID==unique_loci[j]],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$morphology.score.fixed.circles[input$PA14.ID==unique_loci[j]]
  
  fitmodel <- nls(y_growth ~ a / ( 1 + exp(-b * (x_timepoints - c))),start=start_val_sigmoid(x_timepoints, y_growth), algorithm = "port")
  params=coef(fitmodel)
  
  x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
  y_growth_predicted <- sigmoid(params,x_timepoint_range)
  
  
  plot(x_timepoints,y_growth)
  df<-data.frame(list(x=x_timepoints, y=y_growth))
  df_predict<-data.frame(list(x=x_timepoint_range, y=y_growth_predicted ))
  
  output_plot<-ggplot(df,aes(x,y))+
    geom_point(col="red")+
    geom_line(data=df_predict,size=2)+ 
    ggtitle(paste0("Gene ", unique_loci[j]))+
    theme_bw()+
    #ylim(range(0,30000))+
    # ylab("Colony Size")+
    ylab("Morphology")+
    xlab("Time Pionts")+ 
    theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
          axis.text.y = element_text( size=11, hjust = 1),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          strip.text.x = element_text(
            size = 13, color = "black", face = "bold"
          ))+
    scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
  
  return(output_plot)
}



plotting_color<-function(j){
  library(tidyverse)
  input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv",show_col_types = FALSE)
  unique_loci<-unique(input$gene)
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$time_point[input$gene==unique_loci[j]],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$value[input$gene==unique_loci[j]]
  
  fitmodel <- nls(y_growth ~ a / ( 1 + exp(-b * (x_timepoints - c))),start=start_val_sigmoid(x_timepoints, y_growth), alg = "plinear")
  params=coef(fitmodel)
  
  x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
  y_growth_predicted <- sigmoid(params,x_timepoint_range)
  
  
  plot(x_timepoints,y_growth)
  df<-data.frame(list(x=x_timepoints, y=y_growth))
  df_predict<-data.frame(list(x=x_timepoint_range, y=y_growth_predicted ))
  
  output_plot<-ggplot(df,aes(x,y))+
    geom_point(col="red")+
    geom_line(data=df_predict,size=2)+ 
    ggtitle(paste0("Gene ", unique_loci[j]))+
    theme_bw()+
    ylim(range(0,30000))+
    # ylab("Colony Size")+
    ylab("Morphology")+
    xlab("Time Pionts")+ 
    theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
          axis.text.y = element_text( size=11, hjust = 1),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          strip.text.x = element_text(
            size = 13, color = "black", face = "bold"
          ))+
    scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
  
  return(output_plot)
}

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv",show_col_types = FALSE)
unique_loci<-unique(input$PA14.ID)
plotting(which(unique_loci=="PA14_42220"))


library(corrplot)
library(tidyverse)
corrplot( cor(properties[,1:14], use="pairwise.complete.obs"), method = "circle", type="upper")

annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/annotation.csv")
properties<-read_csv("/Users/moradigd/Documents/chemicalgenomic/new_Growth_properties_morphology.csv")
properties$annotation<-annotation$`"Active" Gene Description`[match(properties$`locus id`, annotation$`"Active" Gene Locus`)]
#properties$R_squared<-as.numeric(as.character(properties$R_squared))
head(properties[order(-properties$a),])


filenames <- list.files("/Users/moradigd/Documents/chemicalgenomic/errors", pattern="*.jpg", full.names=TRUE)
filenames<-gsub("/Users/moradigd/Documents/chemicalgenomic/errors/image_output_","",filenames )
filenames<-gsub(".jpg","",filenames )

#input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology_new.csv",show_col_types = FALSE)
#four_colony<-names(table(input$gene))[which(table(input$gene)==24)]

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv",show_col_types = FALSE)
four_colony<-names(table(input$PA14.ID))[which(table(input$PA14.ID)==24)]

properties<-properties[which(! properties$`locus id` %in% filenames),]

properties<-properties[order(-properties$inflection_point_slope),]
properties<-properties[!is.na(properties$inflection_point_slope),]
properties<-properties[which(properties$`locus id` %in%  four_colony),]
head(properties$`locus id`)
tail(properties$`locus id`,50)




j<-50
library(inflection)
library(nlstools)
sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv",show_col_types = FALSE)
properties<-read.csv("/Users/moradigd/Documents/chemicalgenomic/new_Growth_properties_morphology.csv")
input_gene<-unique_loci[j]
plot_gene_from_file<-function(input_gene,input,properties){
  unique_loci<-unique(input$PA14.ID)
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$timepoint[input$PA14.ID==input_gene],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$morphology.score.fixed.circles[input$PA14.ID==input_gene]
  
  params<-properties[properties$locus.id==input_gene,c(2,3,4)]
  as.numeric(params)
  
  x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
  y_growth_predicted <- sigmoid(as.numeric(params),x_timepoint_range)
  
  plot(x_timepoints,y_growth)
  df<-data.frame(list(x=x_timepoints, y=y_growth))
  df_predict<-data.frame(list(x=x_timepoint_range, y=y_growth_predicted ))
  
  output<-ggplot(df,aes(x,y))+
    geom_point(col="red")+
    geom_line(data=df_predict,size=2)+ 
    ggtitle(paste0("Gene ", unique_loci[j]))+
    theme_bw()+
    ylim(range(0,300))+
    ylab("Morphology")+
    xlab("Time Pionts")+ 
    theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
          axis.text.y = element_text( size=11, hjust = 1),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          strip.text.x = element_text(
            size = 13, color = "black", face = "bold"
          ))+
    scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
  
  return(output)
}

plot_gene_from_file(input_gene,input,  properties )

#Abnormality 
abnormality<-read_csv("/Users/moradigd/Documents/chemicalgenomic/abnormality_colour.csv")
clustering<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")
clustering$abnormality<-abnormality$percentage[match( clustering$isolate, abnormality$tag)]
boxplot(clustering$abnormality~ clustering$cluster, xlab="Clusters", ylab = "Abnormality")

#top bottom properties  
library(tidyverse)
library(edgeR)
properties<-read_csv("/Users/moradigd/Documents/chemicalgenomic/new_new_Growth_properties_colony_colour.csv")

abnormality<-read_csv("/Users/moradigd/Documents/chemicalgenomic/abnormality_colour.csv")
properties$abnormality<-abnormality$percentage[match(properties$`locus id`, abnormality$tag)]

properties<-properties[,c(1,2,3,4,15, 16,19,17)]
colnames(properties)<-c("Slope at Midpoint", "Projected Max", "Projected slope", "Midpoint", "Stochasticity", "Average increase", "Abnormality","locus id")

keg_tot<-c()
for(i in 1:7){
  print(i)
  
  properties_tmp<-properties[which(!is.na(properties[,i])),]
  locus_tags<-properties_tmp$`locus id`[order(properties_tmp[,i])][1:floor(nrow(properties)*0.05)]
  write_csv(  data.frame(list(tags=locus_tags)), paste0("/Users/moradigd/Documents/chemicalgenomic/Bottom_color_", colnames(properties)[i], ".csv" ))
  keg<-kegga(locus_tags, species.KEGG="pau")
  keg_tmp<-topKEGG(keg) %>% filter(P.DE<0.05)
  if(nrow(keg_tmp)>0){
    keg_tmp$tag<-"bot"
    keg_tmp$name<-colnames(properties_tmp)[i]
    keg_tot<-rbind(keg_tot, keg_tmp)
  }
  
  properties_tmp<-properties[which(!is.na(properties[,i])),]
  locus_tags<-properties_tmp$`locus id`[order(-properties_tmp[,i])][1:floor(nrow(properties)*0.05)]
  write_csv(  data.frame(list(tags=locus_tags)), paste0("/Users/moradigd/Documents/chemicalgenomic/Top_color_", colnames(properties)[i], ".csv" ))
  keg<-kegga(locus_tags, species.KEGG="pau")
  keg_tmp<-topKEGG(keg) %>% filter(P.DE<0.05)
  if(nrow(keg_tmp)>0){
    keg_tmp$tag<-"top"
    keg_tmp$name<-colnames(properties_tmp)[i]
    keg_tot<-rbind(keg_tot, keg_tmp)
  }
}

#b rate of incrrease 
#a projected max
#c inflectiion point 
#rate of increase 


keg_tot %>% filter(tag=="top" & P.DE<0.01) %>%
  ggplot(aes(x=name, y=-log10(P.DE), fill=Pathway))+
  geom_bar(stat = 'identity',position=position_dodge(), colour="black")+
  theme_bw()+
  ylim(range(0,10))+
  ylab("PValue")+
  xlab("Characteristics")+ 
  theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
        axis.text.y = element_text( size=11, hjust = 1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        strip.text.x = element_text(
          size = 13, color = "black", face = "bold"
        ))


keg_tot %>% filter(tag=="top" & P.DE<0.01) %>%
  group_by(name) %>% 
  mutate(position = rank(log10(P.DE))) %>%
  ggplot(aes(x=name, y=-log10(P.DE), fill= reorder(Pathway, -log10(P.DE)), group=position ))+
  geom_bar(stat = 'identity',position="dodge", colour="black")+
  theme_bw()+
  ylim(range(0,6))+
  ylab("-log10(P-value)")+
  xlab("Characteristics")+ 
  guides(fill=guide_legend(title="Pathway"))+
  theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
        axis.text.y = element_text( size=11, hjust = 1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        strip.text.x = element_text(
          size = 13, color = "black", face = "bold"
        ))


length(which(is.na(properties$a)))

locus_tags<-properties$`locus id`[which(properties$signficance=="L")]
keg<-kegga(locus_tags, species.KEGG="pau")
go<-goana(locus_tags, species.KEGG="pau")
topKEGG(keg)
topGO(go)
?goana
go<-goana(locus_tags)

#extract kegg and goana top bottom
tags_biofilm<-c()
attribute<-c()
gene_biofilm<-c()
for(i in 1:7){
  print(i)
  bottom<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/Bottom_size_", colnames(properties)[i], ".csv" ))
  for(j in 1:length(bottom$tags)){
    print(j)
    kegg<-topKEGG(kegga(bottom$tags[j], species.KEGG="pau"))
    ids_tmp<-which(grepl("Biofilm",kegg$Pathway[kegg$P.DE<1] ))
    if (length(ids_tmp)>0){
      tags_biofilm<-c(tags_biofilm, "bottom")
      attribute<-c(attribute, colnames(properties)[i])
      gene_biofilm<-c(gene_biofilm, bottom$tags[j] )
    }
  }
  
  top<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/Top_colour_", colnames(properties)[i], ".csv" ))
  for(j in 1:length(bottom$tags)){
    print(j)
    kegg<-topKEGG(kegga(top$tags[j], species.KEGG="pau"))
    ids_tmp<-which(grepl("Biofilm",kegg$Pathway[kegg$P.DE<1] ))
    if (length(ids_tmp)>0){
      tags_biofilm<-c(tags_biofilm, "top")
      attribute<-c(attribute, colnames(properties)[i])
      gene_biofilm<-c(gene_biofilm, bottom$tags[j] )
    }
  }
}

biofilm_genes<-data.frame(list(tags_biofilm=tags_biofilm,attribute=attribute,gene_biofilm=gene_biofilm  ))
biofilm_genes<-data.frame(list(tags_biofilm=tags_biofilm,attribute=attribute,gene_biofilm=gene_biofilm,pvalue=pvalue  ))
biofilm_genes$annotation<-annotation$`"Active" Gene Description`[match(biofilm_genes$gene_biofilm, annotation$`"Active" Gene Locus`)]
biofilm_genes$gene_names<-annotation$`"Active" Gene Name`[match(biofilm_genes$gene_biofilm, annotation$`"Active" Gene Locus`)]
write_csv(biofilm_genes, "/Users/moradigd/Documents/chemicalgenomic/biofilm_genes_morphology.csv")

#extract kegg and goana top bottom from files 
tags_biofilm<-c()
attribute<-c()
gene_biofilm<-c()
pvalue<-c()
for(i in 1:7){
  print(i)
  
  bottom<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/Bottom_size_", colnames(properties)[i], ".csv" ),show_col_types = FALSE)
  for(j in 1:length(bottom$tags)){
    print(j)
    if(file.exists(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/KEGG_", bottom$tags[j], ".csv" ))){
      
      kegg<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/KEGG_", bottom$tags[j], ".csv" ),show_col_types = FALSE)
      #topKEGG(kegga(bottom$tags[j], species.KEGG="pau"))
      ids_tmp<-which(grepl("Porphyrin",kegg$Pathway[kegg$P.DE<1] ))
      if (length(ids_tmp)>0){
        tags_biofilm<-c(tags_biofilm, "bottom")
        attribute<-c(attribute, colnames(properties)[i])
        gene_biofilm<-c(gene_biofilm, bottom$tags[j] )
        pvalue<-c(pvalue, kegg$P.DE[kegg$P.DE<1][ids_tmp])
      }
    }
  }
  
  top<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/Top_size_", colnames(properties)[i], ".csv" ),show_col_types = FALSE)
  for(j in 1:length(top$tags)){
    print(j)
    if(file.exists(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/KEGG_", top$tags[j], ".csv" ))){
      
      kegg<-read_csv(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/KEGG_", top$tags[j], ".csv" ),show_col_types = FALSE)
      #kegg<-topKEGG(kegga(top$tags[j], species.KEGG="pau"))
      ids_tmp<-which(grepl("Porphyrin",kegg$Pathway[kegg$P.DE<1] ))
      if (length(ids_tmp)>0){
        tags_biofilm<-c(tags_biofilm, "top")
        attribute<-c(attribute, colnames(properties)[i])
        gene_biofilm<-c(gene_biofilm, top$tags[j] )
        pvalue<-c(pvalue, kegg$P.DE[kegg$P.DE<1][ids_tmp])
      }
    }
  }
}

annotation<-read_csv( "/Users/moradigd/Documents/chemicalgenomic/annotation.csv")
biofilm_genes<-data.frame(list(tags_biofilm=tags_biofilm,attribute=attribute,gene_biofilm=gene_biofilm  ))
biofilm_genes$annotation<-annotation$`"Active" Gene Description`[match(biofilm_genes$gene_biofilm, annotation$`"Active" Gene Locus`)]
biofilm_genes$gene_names<-annotation$`"Active" Gene Name`[match(biofilm_genes$gene_biofilm, annotation$`"Active" Gene Locus`)]
write_csv(biofilm_genes, "/Users/moradigd/Documents/chemicalgenomic/Porphyrin_genes_size.csv")

tmp<-read_csv("/Users/moradigd/Documents/chemicalgenomic/biofilm_genes_morphology.csv")
unique(tmp$gene_biofilm)

#extract kegg and goana all gene 
library(tidyverse)
library(edgeR)
properties<-read_csv("/Users/moradigd/Documents/chemicalgenomic/new_new_Growth_properties_morphology.csv")
for(i in 1:length(properties$`locus id`)){
  print(i/length(properties$`locus id`))
  write_csv(topKEGG(kegga(properties$`locus id`[i], species.KEGG="pau")), paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/KEGG_",properties$`locus id`[i],".csv" )) 
  write_csv(topGO(goana(properties$`locus id`[i], species.KEGG="pau")), paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/GO_",properties$`locus id`[i],".csv" )) 
}



#isolation forests 
#import packages
library(ggplot2)
library(solitude)

library(usethis) 
usethis::edit_r_environ()


library(tidyverse)
library(mclust)
library(randomcoloR)
library(solitude)

filenames <- list.files("/Users/moradigd/Documents/chemicalgenomic/errors", pattern="*.jpg", full.names=TRUE)
filenames<-gsub("/Users/moradigd/Documents/chemicalgenomic/errors/image_output_","",filenames )
filenames<-gsub(".jpg","",filenames )

distanec_matrix<-read_csv( "/Users/moradigd/Documents/chemicalgenomic/aggregated_morphology.csv")
distanec_matrix<-distanec_matrix[which(!colnames(distanec_matrix)  %in% filenames), which(!colnames(distanec_matrix)  %in% filenames)]
distanec_matrix<-as.matrix(distanec_matrix)
distanec_matrix<-distanec_matrix[which(apply(distanec_matrix, 2, var) != 0) , which(apply(distanec_matrix, 2, var) != 0)]

iforest<- isolationForest$new()
iforest$fit(distanec_matrix)

prediction<- iforest$predict(distanec_matrix)

require(Rtsne)
tsne <- Rtsne(distanec_matrix, check_duplicates = FALSE, pca = T, perplexity=100, theta=0.5, dims=2)

cols <- rainbow(2)
par(mar=c(4, 4, 2,2))

percentage<-(prediction$anomaly_score-min(prediction$anomaly_score))/(max(prediction$anomaly_score)-min(prediction$anomaly_score))
outlier <- as.factor(ifelse(percentage >=0.5, "outlier", "normal"))
plot(tsne$Y,col= alpha(cols[outlier], 0.55) , pch = 16,alpha = 0.35, xlab="t-SNE Dim(1)", ylab="t-SNE Dim(2)")
write_csv(data.frame(list(tag=colnames(distanec_matrix),percentage=percentage )), "/Users/moradigd/Documents/chemicalgenomic/abnormality_morphology.csv")







#create and plot sample data

#create data
n = 1000
Var1 = c(rnorm(n, 0, 0.5), rnorm(n*0.1, -2, 1))
Var2 = c(rnorm(n, 0, 0.5), rnorm(n*0.1,  2, 1))
outliers = c(rep(0, n), rep(1, (0.1*n))) + 3
data = data.frame(Var1, Var2)

#plot data
ggplot(data, aes(x = Var1, y = Var2)) + 
  geom_point(shape = 1, alpha = 0.5) +
  labs(x = "x", y = "y") +
  labs(alpha = "", colour="Legend")

#create isolation forest using isolationForest function from solitude package with default parameters

n = 1000
Var1 = c(rnorm(n, 0, 0.5), rnorm(n*0.1, -2, 1))
Var2 = c(rnorm(n, 0, 0.5), rnorm(n*0.1,  2, 1))
outliers = c(rep(0, n), rep(1, (0.1*n))) + 3
data = data.frame(Var1, Var2)

iforest<- isolationForest$new()

iforest$fit(data)
#print(iforest$scores)
#predict outliers within dataset
data$pred <- iforest$predict(data)
data$outlier <- as.factor(ifelse(data$pred$anomaly_score >=0.80, "outlier", "normal"))

range(data$pred$anomaly_score)
#plot data again with outliers identified
ggplot(data, aes(x = Var1, y = Var2, color = outlier)) + 
  geom_point(shape = 1, alpha = 0.5) +
  labs(x = "x", y = "y") +
  labs(alpha = "", colour="Legend")


rep<-c()
for(i in properties$`locus id`){
  rep<-c(rep, length(which(input$PA14.ID %in% i)))
}
properties<-properties[which(rep==24),]
table(rep)

#biofilm genes
library(Hmisc) 
library(pheatmap)
library(tidyverse)

#biofilm_genes<-c("PA14_16500", "PA14_12810", "PA14_56790","PA14_16500","PA14_64050",  "PA14_16500")

biofilm_genes<-c("PA14_09480", "PA14_39910","PA14_09470","PA14_42760","PA14_13250", "PA14_13240", "PA14_13230", "PA14_24900", "PA14_13850", "PA14_51410", "PA14_51430", "PA14_24480", "PA14_24490", "PA14_24510", "PA14_24530", "PA14_24550", "PA14_24560",
                 "PA14_16500",
                 "PA14_64050",
                 "PA14_12810",
                 "PA14_56790",
                 "PA14_50060",
                 "PA14_02110",
                 "PA14_42220",
                 "PA14_66320",
                 "PA14_21190",
                 "PA14_07730",
                 "PA14_31290",
                 "PA14_45940",
                 "PA14_52260",
                 "PA14_54430",
                 "PA14_09150",
                 "PA14_61200")

genes<-c("phzA", "phzE","phzB","aroC", "moaE", "moaD", "moaC", "moaB", "moaA", "pqsC", "pqsA", "pelA","pelB", "pelD","pelE", "pelF","pelG",
         "wspR",
         "gcbA",
         "rocR",
         "bifA",
         "roeA",
         "siaD",
         "mucR",
         "dipA",
         "nbdA",
         "rsmA",
         "lecA",
         "lasI",
         "gacS",
         "algU",
         "katA",
         "cdrA")

#df<-data.frame(list(gene=genes, llocud_is=biofilm_genes))
#write_csv(df,"/Users/moradigd/Documents/chemicalgenomic/biofilm_genes.csv" )

gene_id<-paste0(biofilm_genes,"_",genes)

properties<-read_csv("/Users/moradigd/Documents/chemicalgenomic/new_new_Growth_properties_colony_size.csv")
properties<-properties[,c(1,2,3,4,15, 16,17)]
colnames(properties)<-c("Slope at Midpoint", "Projected Max", "Projected slope", "Midpoint", "Stochasticity", "Average increase", "locus id")

abnormality<-read_csv("/Users/moradigd/Documents/chemicalgenomic/abnormality_morphology.csv")
properties$abnormality<-abnormality$percentage[match(properties$`locus id`, abnormality$tag)]
#"/Users/moradigd/Documents/chemicalgenomic/abnormality_morphology.csv"

properties_quant<-data.frame(apply(properties[,c(1:6,8)],2 , function(x) as.numeric(cut2(x, g=20))))
row.names(properties_quant) <-properties$`locus id`
properties_quant_tot<-(properties_quant-1)/20
properties_quant_tot$Gene <-properties$`locus id`
#write.csv(properties_quant_tot, "/Users/moradigd/Documents/chemicalgenomic/property_colony_colour_tot.csv")



properties_quant<-properties_quant[match(biofilm_genes, row.names(properties_quant)), ] #%>% filter(!which(row.names(.)=="NA"))
properties_quant_short<-properties_quant[which(!is.na(match(biofilm_genes, row.names(properties_quant))) ),]
row.names(properties_quant_short)<-gene_id[which(!is.na(match(biofilm_genes, row.names(properties_quant))) )]


pheatmap((properties_quant_short-1)/20, display_numbers = T, cluster_rows = FALSE)
row.names(properties_quant_short)
?pheatmap
#b rate of incrrease 
#a projected max
#c inflectiion point 
#rate of increase 

dim(properties_quant_short)


#distribution of variance 
library(tidyverse)
input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv",show_col_types = FALSE)
frequency<-data.frame(table(input$PA14.ID)/6) 

properties<-read_csv("/Users/moradigd/Documents/chemicalgenomic/new_Growth_properties_morphology.csv")
abnormality<-read_csv("/Users/moradigd/Documents/chemicalgenomic/abnormality_morphology.csv")
properties$abnormality<-abnormality$percentage[match(properties$`locus id`, abnormality$tag)]
properties<-properties[!is.na(properties$R_squared),] %>% select(`locus id`,R_squared )
properties<-inner_join(properties,frequency, by=c(`locus id`="Var1") )
properties$R_squared<-properties$R_squared/properties$Freq

shortened<-properties[match(biofilm_genes,properties$`locus id`),] %>% drop_na()
shortened$tag<-gene_id[match(shortened$`locus id`,biofilm_genes)] 

properties %>%
  ggplot( aes(x=log(R_squared))) +
  geom_histogram(fill="blue") +
  ggplot2::annotate("text", x = log(shortened$R_squared), y=500, label = shortened$tag, angle = 90)+
  geom_vline(xintercept = log(shortened$R_squared), linetype="dotted", size = 0.3)+
  geom_vline(xintercept = c(log(quantile(properties$R_squared, 0.95)), log(quantile(properties$R_squared, 0.05))), linetype="solid", col="red", size = 0.3)+
  #geom_hline(yintercept = 1, linetype="solid", col="black", size = 0.3)+
  theme_bw()+ 
  ylim(range(0,1500))+
  ylab("Frequency")+
  xlab("R squared")+
  theme(axis.text.x = element_text( size=14, angle=45,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+ 
  scale_x_continuous(labels = scales::comma)+ scale_x_continuous(labels = scales::comma)



#distribution 
library(tidyverse)
properties<-read_csv("/Users/moradigd/Documents/chemicalgenomic/new_Growth_properties_morphology.csv")
abnormality<-read_csv("/Users/moradigd/Documents/chemicalgenomic/abnormality_morphology.csv")
properties$abnormality<-abnormality$percentage[match(properties$`locus id`, abnormality$tag)]

#properties$abnormality<-scale(properties$abnormality)
properties<-properties[!is.na(properties$abnormality),]
properties$signficance<-ifelse(properties$a< sort(properties$a)[floor(length(properties$a)*0.05)], "L",      ifelse(properties$a> sort(properties$a)[floor(length(properties$a)*0.95)] , "H", "M" ))    

shortened<-properties[match(biofilm_genes,properties$`locus id`),]
shortened$tag<-gene_id[match(shortened$`locus id`,biofilm_genes)]
shortened<-shortened[!is.na(shortened$abnormality),]
shortened$`growth rate`


properties %>%
  ggplot( aes(x=abnormality)) +
  geom_density(fill="blue") +
  ggplot2::annotate("text", x = shortened$abnormality, y=4.5, label = shortened$tag, angle = 90)+
  geom_vline(xintercept = shortened$abnormality, linetype="dotted", size = 0.3)+
  geom_vline(xintercept = c(quantile(properties$abnormality, 0.95), quantile(properties$abnormality, 0.05)), linetype="solid", col="red", size = 0.3)+
  #geom_hline(yintercept = 1, linetype="solid", col="black", size = 0.3)+
  theme_bw()+ 
  ylim(range(0,max(properties$abnormality)+5))+
  ylab("Frequency")+
  xlab("Average increase")+
  theme(axis.text.x = element_text( size=14, angle=45,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+ 
  scale_x_continuous(labels = scales::comma)+ scale_x_continuous(labels = scales::comma)
shortened$tag

#abnormality 

properties %>%
  ggplot( aes(x=abnormality)) +
  geom_density(fill="blue") +
  ggplot2::annotate("text", x = shortened$abnormality, y=4.5, label = shortened$tag, angle = 90)+
  geom_vline(xintercept = shortened$abnormality, linetype="dotted", size = 0.3)+
  geom_vline(xintercept = c(quantile(properties$abnormality, 0.95)), linetype="solid", col="red", size = 0.3)+
  #geom_hline(yintercept = 1, linetype="solid", col="black", size = 0.3)+
  theme_bw()+ 
  ylim(range(0,max(properties$abnormality)+5))+
  ylab("Frequency")+
  xlab("Abnormality")+
  theme(axis.text.x = element_text( size=14, angle=45,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+ 
  scale_x_continuous(labels = scales::comma)+ scale_x_continuous(labels = scales::comma)

shortened$abnormality[order(shortened$abnormality)]
shortened$tag[order(shortened$abnormality)]

properties %>%
  ggplot( aes(x=`growth rate`)) +
  geom_histogram(data=subset(properties,signficance=="H"), fill="red",binwidth = 5) +
  geom_histogram(data=subset(properties,signficance=="M"), fill="blue",binwidth = 5) +
  geom_histogram(data=subset(properties,signficance=="L"), fill="red",binwidth = 5)+ 
  ggplot2::annotate("text", x = shortened$`growth rate`, y=550, label = shortened$tag, angle = 90)+
  geom_vline(xintercept = shortened$`growth rate`, linetype="dotted", size = 0.3)+
  theme_bw()+ 
  ylim(range(0,700))+
  ylab("Frequency")+
  xlab("Maximum Antipitated Growth")+
  theme(axis.text.x = element_text( size=14, angle=45,hjust = 1),
        axis.text.y = element_text( size=14, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+ 
  scale_x_continuous(labels = scales::comma)+ scale_x_continuous(labels = scales::comma)

library(topGO)
library(edgeR)
locus_tags<-properties$`locus id`[which(properties$signficance=="L")]
keg<-kegga(locus_tags, species.KEGG="pau")
go<-goana(locus_tags, species.KEGG="pau")
topKEGG(keg)
topGO(go)
topGO(goana("PA14_34210", species.KEGG="pau"))
topKEGG(kegga("PA14_34210", species.KEGG="pau"))

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv",show_col_types = FALSE)
unique_loci<-unique(input$PA14.ID)
plotting(which(unique_loci=="PA14_72560"))
properties[which(properties$`locus id`=="PA14_56530"),]




tail(properties[order(-properties$`growth rate`),])
tail(properties)
properties
library(edgeR)
chunk_size<-50
for(i in seq(floor(dim(properties)[1]/chunk_size),dim(properties)[1],floor(dim(properties)[1]/chunk_size))){
  isolates<-properties$`locus id`[(i-floor(dim(properties)[1]/chunk_size)+1):i]
  keg<-kegga(isolates, species.KEGG="pau")
  top_kegg<-topKEGG(keg)
  top_kegg$id<-paste0("chunk_",i/floor(dim(properties)[1]/chunk_size) )
  print(top_kegg[top_kegg$P.DE<0.05,])
  print(head(isolates))
}

#Abnormaly prediction 

#Biofilm formation 
#gcbA 
#NicD

#Clsuters unsupervised learning
library(tidyverse)
library(edgeR)
library(circlize)

df_circularity<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv")
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")

shared<-intersect(intersect(df_circularity$isolate, df_morphology$isolate), df_size$isolate)
df_morphology<-df_morphology[match(shared,df_morphology$isolate),]
df_size<-df_size[match(shared,df_size$isolate),]
df_colour<-df_colour[match(shared,df_size$isolate),]

#write_csv(df_morphology,"/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology_shortened.csv" )
#write_csv(df_size,"/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv" )
#write_csv(df_circularity,"/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv" )



library(tidyverse)
library(mclust)
library(randomcoloR)
library(solitude)
library(ClassDiscovery)

#Worked
df_morphology_clusters<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints_summarised.csv")
growth$cluster<-df_morphology_clusters$cluster[match(growth$label, df_morphology_clusters$isolate )]
growth<-growth[complete.cases(growth),]
colnames(growth)
growth_short<-growth %>% group_by(timepoint,cluster) %>% 
  summarise(average=mean(colony.size)) 
growth_short$timepointsequence<-array(sapply(seq(1,6), function(x) rep(x,9)))
growth_short$cluster<-as.character(growth_short$cluster)


#For circularity 
#df_morphology_clusters$cluster[which(df_morphology_clusters$cluster==5)]<-6

distance_groups<-matrix(0, nrow = length(unique(df_morphology_clusters$cluster)), ncol = length(unique(df_morphology_clusters$cluster)))

for(i in 1:length(unique(df_morphology_clusters$cluster))){
  for(j in 1:length(unique(df_morphology_clusters$cluster))){
    tmp1<-growth_short$average[growth_short$cluster==unique(df_morphology_clusters$cluster)[i]]
    tmp2<-growth_short$average[growth_short$cluster==unique(df_morphology_clusters$cluster)[j]]
    distance_groups[i,j]<-mean(tmp1)/mean(tmp2)
  }
}
colnames(distance_groups)<-unique(df_morphology_clusters$cluster)
row.names(distance_groups)<-unique(df_morphology_clusters$cluster)
hclust_avg <- hclust(dist(distance_groups))
plot(hclust_avg, hang = -1, cex = 1.5, las=2, xlab="", )
#morphology
barplot(table(growth$cluster)[c(2,1,5,4,3,7,8,6)], col="black", las=2,ylim=range(0,8000), ylab="Frequency")
box()
#circularity
barplot(table(df_morphology_clusters$cluster)[c(3,4,2,7,1,5,6,8)], col="black", las=2,ylim=range(0,1000), ylab="Frequency")
box()
#size
barplot(table(growth$cluster)[c(2,4,7,3,1,8,9,6,5)], col="black", las=2,ylim=range(0,6000), ylab="Frequency")
box()


#colony coulour 
df_morphology_clusters<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_colour_annotation.csv")
growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv")
growth$cluster<-df_morphology_clusters$cluster[match(growth$gene, df_morphology_clusters$isolate )]
growth<-growth[complete.cases(growth),]
colnames(growth)
growth_short<-growth %>% group_by(time_point,cluster) %>% 
  summarise(average=mean(value)) 
growth_short$timepointsequence<-array(sapply(seq(1,6), function(x) rep(x,9)))
growth_short$cluster<-as.character(growth_short$cluster)

unique(df_morphology_clusters$cluster)
#For circularity 
#df_morphology_clusters$cluster[which(df_morphology_clusters$cluster==5)]<-6

distance_groups<-matrix(0, nrow = length(unique(df_morphology_clusters$cluster)), ncol = length(unique(df_morphology_clusters$cluster)))

for(i in 1:length(unique(df_morphology_clusters$cluster))){
  for(j in 1:length(unique(df_morphology_clusters$cluster))){
    tmp1<-growth_short$average[growth_short$cluster==unique(df_morphology_clusters$cluster)[i]]
    tmp2<-growth_short$average[growth_short$cluster==unique(df_morphology_clusters$cluster)[j]]
    distance_groups[i,j]<-mean(tmp1)/mean(tmp2)
  }
}
colnames(distance_groups)<-unique(df_morphology_clusters$cluster)
row.names(distance_groups)<-unique(df_morphology_clusters$cluster)
hclust_avg <- hclust(dist(distance_groups))
plot(hclust_avg, hang = -1, cex = 1.5, las=2, xlab="", )
#morphology
barplot(table(df_morphology_clusters$cluster)[c(7,1,3,8,5,6,9,2,4)], col="black", las=2,ylim=range(0,1000), ylab="Frequency")
box()
#circularity
barplot(table(df_morphology_clusters$cluster)[c(3,4,2,7,1,5,6,8)], col="black", las=2,ylim=range(0,1000), ylab="Frequency")
box()
#size
barplot(table(growth$cluster)[c(2,4,7,3,1,8,9,6,5)], col="black", las=2,ylim=range(0,6000), ylab="Frequency")
box()





for(j in 1:8){
  annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/annotation.csv",show_col_types = FALSE)
  ids<-annotation[match( df_morphology_clusters$isolate[df_morphology_clusters$cluster==j] ,    annotation$`"Active" Gene Locus`),]
  ids<-ids[, c(23, 20, 14, 19)]
  ids$`PAO1 Orthologs`<- gsub('"', '', ids$`PAO1 Orthologs`)
  colnames(ids)<-c("PAO1_loci", "Description","PA14_loci","Gene_name")
  write_csv(ids, paste0("/Users/moradigd/Documents/chemicalgenomic/string_circularity_", j, ".csv"))
}

#concordance between clusters
df1<-data.frame(table(df_colour$cluster, df_morphology$cluster)) 
df1$Var1<-paste0("C_",df1$Var1 )
df1$Var2<-paste0("M_",df1$Var2 )

grid.col = c(C_1 = "brown", C_2 = "brown1",C_3 = "brown2", C_4 = "brown3",C_5 = "brown4",C_6 = "darkred",C_7 = "firebrick", C_8 = "firebrick2",C_9 = "firebrick4",
             M_1 = "chartreuse", M_2 = "chartreuse1",M_3 = "chartreuse2", M_4 = "chartreuse3",M_5 = "chartreuse4",M_6 = "darkolivegreen",M_7 = "darkolivegreen1", M_8 = "darkolivegreen2",M_9 = "darkolivegreen3"
)

chordDiagram(df1, grid.col = grid.col)
chordDiagram(df1)



df1<-data.frame(table(df_colour$cluster, df_size$cluster)) 
df1$Var1<-paste0("C_",df1$Var1 )
df1$Var2<-paste0("S_",df1$Var2 )

grid.col = c(C_1 = "brown", C_2 = "brown1",C_3 = "brown2", C_4 = "brown3",C_5 = "brown4",C_6 = "darkred",C_7 = "firebrick", C_8 = "firebrick2",C_9 = "firebrick4",
             S_1 = "darkorange", S_2 = "darkorange1",S_3 = "darkorange2", S_4 = "darkorange3",S_5 = "darkorange4",S_6 = "chocolate",S_7 = "chocolate1", S_8 = "chocolate2",S_9 = "chocolate3"
)

chordDiagram(df1, grid.col = grid.col)
chordDiagram(df1)


df1<-data.frame(table(df_size$cluster, df_morphology$cluster)) 
df1$Var1<-paste0("S_",df1$Var1 )
df1$Var2<-paste0("M_",df1$Var2 )

grid.col = c(S_1 = "darkorange", S_2 = "darkorange1",S_3 = "darkorange2", S_4 = "darkorange3",S_5 = "darkorange4",S_6 = "chocolate",S_7 = "chocolate1", S_8 = "chocolate2",S_9 = "chocolate3",
             M_1 = "chartreuse", M_2 = "chartreuse1",M_3 = "chartreuse2", M_4 = "chartreuse3",M_5 = "chartreuse4",M_6 = "darkolivegreen",M_7 = "darkolivegreen1", M_8 = "darkolivegreen2",M_9 = "darkolivegreen3"
)

chordDiagram(df1, grid.col = grid.col)
chordDiagram(df1)





set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)
mat
col_mat = rand_color(length(mat), transparency = 0.5)


mat2 = matrix(sample(100, 35), nrow = 5)
rownames(mat2) = letters[1:5]
colnames(mat2) = letters[1:7]
mat2


chordDiagram(mat2,  directional = 1, row.col = 1:5)


annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated_colony_size.csv")
df<-data.frame(list(cluster=mod1$classification, go=annotation$go[match(colnames(tmp_red), annotation$label)]  ))
df$go<-ifelse(is.na(df$go),"A",df$go)
freq_df<-table(df$cluster, df$go)
freq_df<-freq_df[,which(sapply(colnames(freq_df), function(x) nchar(x)) ==1)]

freq_df_rel<-freq_df/apply(freq_df, 1, sum)
#freq_df_rel<-freq_df
palette <- distinctColorPalette(dim(freq_df_rel)[2])
barplot(t(freq_df_rel), col=palette, las=2, type='n')
barplot(t(freq_df_rel), col="black", las=2, type='n')
legend("topright", colnames(freq_df), fill = palette, bty = "n")
box()

annotation$go
annotation$go<-ifelse(is.na(annotation$go),"A",annotation$go)
ann_table<-table(annotation$go)/sum(table(annotation$go))
ann_table<-ann_table[match(colnames(freq_df_rel),names(ann_table))]

for(i in 1:dim(freq_df_rel)[1]){
  freq_df_rel[i,][which(freq_df_rel[i,]/ann_table<1)]<-0
}
barplot(t(freq_df_rel), col=palette, las=2, type='n')




?Rtsne

f <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}
m<-f(tmp)
is.na(m) <- sapply(m, is.infinite)

fit <- prcomp(na.omit(m), center = TRUE, scale. = T)

df<-data.frame(summary(fit)$importance)
head(df[,1:3])
#stdData1 <- cbind(fit$x, as.factor(res.km$cluster))

wss <- sapply(1:10, function(k){kmeans(scale(na.omit(m)), k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:10, wss,
     type="b", pch = 19, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares", las=2)

write_csv(data.frame(wss), "/Users/moradigd/Documents/chemicalgenomic/screeplot_kmeans_colony_size.csv") 
plot(wss)

res.km <- kmeans(scale(na.omit(m)), 4, nstart = 50,iter.max = 15)
# K-means clusters showing the group of each individuals
res.km$cluster
fviz_pca_ind(fit, geom="point",habillage=res.km$cluster)

gene_list<-read_csv("/Users/moradigd/Documents/chemicalgenomic/similarity_colony_size.csv")
length(res.km$cluster)


gene_cluster<-data.frame(cbind(gene_list$X1,res.km$cluster))
colnames(gene_cluster)<-c("gene","cluster")
write_csv(gene_cluster, "/Users/moradigd/Documents/chemicalgenomic/cluster_similarity_colony_size.csv")

growth<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints_summarised.csv")
growth$cluster<-gene_cluster$cluster[match(growth$label, gene_cluster$gene )]

growth_short<-growth %>% group_by(timepoint,cluster) %>% 
  summarise(importance_average=mean(colony.size)) 
growth_short$timepointsequence<-array(sapply(seq(1,6), function(x) rep(x,4)))

ggplot(growth_short, aes(x = timepointsequence, y = importance_average, fill = cluster)) +
  geom_line() +
  geom_point(size = 4, shape = 21)+
  theme_bw()+
  ylab("Morphology Score")+
  xlab("Time Pionts")

annotation<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/PA14_computational_v2_updated_annotations.csv")
Goterm<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/ListCOGGOTerms.csv")
annotation$`"Active" Gene Locus`

cogs<-annotation$COG[match(growth$label,annotation$`Active Gene Locus` )]
growth$go<-Goterm$func[match(cogs, Goterm$`# COG`)]
table(growth$go, growth$cluster)/apply(table(growth$go, growth$cluster), 2, sum)

apply(table(growth$go, growth$cluster),2,function(x){x/sum(x)})

apply(table(growth$go, growth$cluster), 2, sum)

tmp<-cbind(annotation[match(growth$label,annotation$`Active Gene Locus`),], growth)
write_csv(tmp,"/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated_colony_size.csv")



plot(tmp$colony.size[tmp$`Active Gene Locus`=="PA14_70850"],
     type="b", pch = 19, 
     xlab="Number of clusters K",
     ylab="Colony Circularity", las=2)
, ylim=range(0.5,0.9)


plot(growth_short$importance_average,type="b", pch = 19,ylim=range(0,200), las=2, ylab="Morphology Score", xlab="Times")

plot( growth_short$timepoint ~ growth_short$importance_average , type="b" , bty="l" , xlab="value of a" , ylab="value of b" , col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17 )
lines(c ~a , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=19 , type="b" )


growth_short<-growth[which(growth$label %in% gene_cluster$gene[gene_cluster$cluster=="4"]),]



boxplot( growth_short$morphology.score.fixed.circles ~ growth_short$timepoint , ylim=range(0,1500), las=2 , ylab="Morphology Score", xlab="Times")


plot(1:10, wss/1000000,las=2,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
box()

library(factoextra)
fviz_cluster(res.km , data = scale(na.omit(m)),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

for(i in seq(dim(tmp[,-1])[1])){
  print(i)
  for(j in seq(dim(tmp[,-1])[1])){
    tmp[,-1][j,i]<-tmp[,-1][i,j]
  }
}
is.na(tmp[,-1]) <- sapply(tmp[,-1], is.infinite)

hist(tmp[upper.tri(tmp, diag = F)])

fit <- prcomp(na.omit(m), cor = T)

princomp(~., data=tmp[,-1], na.action=na.omit)

# do PCA
pca(mydata,colvec=c('gold'),printres=TRUE)

#Fit sigmoid 
library(nlstools)
sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

start_val_sigmoid <- function(x, y) {
  fit <- lm(log(y[which.max(x)] - y + 1e-6) ~ x)
  list(
    a = y[which.max(x)],
    b = unname(-coef(fit)[2]),
    c = unname(-coef(fit)[1] / coef(fit)[2]))
}

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints.csv")
unique_loci<-unique(input$PA14.ID)

properties<-c()
for(j in seq(length(unique_loci))){
  set.seed(j) 
  print(j)
  print(unique_loci[j])
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$timepoint[input$PA14.ID==unique_loci[j]],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$colony.size[input$PA14.ID==unique_loci[j]]
  plot(x_timepoints,y_growth)
  tryCatch({
    fitmodel <- nls(y_growth ~ a / ( 1 + exp(-b * (x_timepoints - c))),start=start_val_sigmoid(x_timepoints, y_growth), algorithm = "port")
    params=coef(fitmodel)
    
    x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
    y_growth_predicted <- sigmoid(params,x1)
    
    inflection_point_slope<-y_growth_predicted[ which(abs(diff(y_growth_predicted))==max(abs(diff(y_growth_predicted))) )]
    
    inflection_point=bese(x_timepoint_range,y_growth_predicted,0)
    properties<-rbind(properties,c(inflection_point_slope,params,as.vector(confint2(fitmodel)),inflection_point$iplast,unique_loci[j]))
    
  }, error = function(e) {an.error.occured <<- TRUE}
  )
}
properties<-data.frame(properties)
colnames(properties)<-c("inflection_point_slope", "a","b","c","a_min","b_min","c_min","a_max","b_max","c_max","inflection","locus id")
write_csv(properties, "/Users/moradigd/Documents/chemicalgenomic/Growth_properties_colony_size.csv")

#Error & confidence interval 
library(gmodels)
library(tidyverse)
library(corrplot)

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints.csv")
unique_loci<-unique(input$PA14.ID)

#lowCI = ci(variable_of_interest)[2],
#hiCI = ci(variable_of_interest)[3], 
#sd = ci (variable_of_interest)[4]

input_sum<- input %>%
  group_by(PA14.ID, timepoint) %>%
  summarise(mean_cl=mean(morphology.score.fixed.circles), min_cl=min(morphology.score.fixed.circles), max_cl=max(morphology.score.fixed.circles), lowCI = ci(morphology.score.fixed.circles, confidence=0.99)[2],hiCI = ci(morphology.score.fixed.circles, confidence=0.99)[3])

#"13h", "19h", "39h", "44h", "69h", "89h"

input_sum_morphology_tmp_13h<-input_sum[input_sum$timepoint=="13h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_morphology_tmp_13h)[8]<-"ranking_morphology_13h"
input_sum_morphology_tmp_19h<-input_sum[input_sum$timepoint=="19h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_morphology_tmp_19h)[8]<-"ranking_morphology_19h"
input_sum_morphology_tmp_39h<-input_sum[input_sum$timepoint=="39h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_morphology_tmp_39h)[8]<-"ranking_morphology_39h"
input_sum_morphology_tmp_44h<-input_sum[input_sum$timepoint=="44h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_morphology_tmp_44h)[8]<-"ranking_morphology_44h"
input_sum_morphology_tmp_69h<-input_sum[input_sum$timepoint=="69h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_morphology_tmp_69h)[8]<-"ranking_morphology_69h"
input_sum_morphology_tmp_89h<-input_sum[input_sum$timepoint=="89h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_morphology_tmp_89h)[8]<-"ranking_morphology_89h"

correlation_summary<-list(input_sum_morphology_tmp_13h, input_sum_morphology_tmp_19h, input_sum_morphology_tmp_39h, input_sum_morphology_tmp_44h,input_sum_morphology_tmp_69h, input_sum_morphology_tmp_89h ) %>% reduce(inner_join, by = "PA14.ID") %>%
  select(ranking_morphology_13h,ranking_morphology_19h,ranking_morphology_39h, ranking_morphology_44h,ranking_morphology_69h, ranking_morphology_89h ) 
correlation_summary<-correlation_summary[,-1]
cor(correlation_summary) %>% corrplot(.)

correlation_summary<-list(input_sum_morphology_tmp_13h, input_sum_morphology_tmp_19h, input_sum_morphology_tmp_39h, input_sum_morphology_tmp_44h,input_sum_morphology_tmp_69h, input_sum_morphology_tmp_89h ) %>% reduce(inner_join, by = "PA14.ID") %>%
  select(ranking_morphology_13h,ranking_morphology_19h,ranking_morphology_39h, ranking_morphology_44h,ranking_morphology_69h, ranking_morphology_89h ) 

slope<-c()
for(j in seq(dim(correlation_summary)[1])){
  tmp<-as.numeric(as.character(coef(lm(unlist(correlation_summary[j,-1],use.names = FALSE) ~ seq(1,6)))[2]))
  slope<-c(slope, tmp)
  print(j)
}
correlation_summary$slope<-slope


data.frame(Time_point=seq(1,6), Ranking= t(correlation_summary[which(correlation_summary$PA14.ID=="PA14_35690"),-1])[-7]) %>%
  ggplot( aes(x = Time_point, y = Ranking) ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15)+
  theme_bw()

genes_functional_group<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated_colony_size.csv")

g<-correlation_summary %>% arrange(desc(slope)) %>% select(PA14.ID, slope) %>% inner_join(genes_functional_group, by=c("PA14.ID"="Active Gene Locus")) %>%
  select('Active Gene Name', 'Active Gene Description','go') %>%
  distinct() 
grid.table(g[1:10,])
dev.off()

#abnormal 
abnormal<-read_csv("/Users/moradigd/Documents/chemicalgenomic/aggegated_morphology_mean_anomaly.csv")
abnormal <- arrange(abnormal, desc(anomaly_score))
boxplot(  abnormal$anomaly_score ~ abnormal$anomaly)
head(abnormal[abnormal$anomaly_score> 0.15 & abnormal$anomaly_score< 0.7,],100)
#0.6875 mean
#0.6485 sd
hist(abnormal$anomaly_score, breaks = 500)

#extract colonies functions 
library(tidyverse)
library(magick)
image_producer<-function(locus_tag){
  genes<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/PA14_computational_v2_updated_annotations.csv")
  row_id<-genes$`Row 384`[which(genes$`"Active" Gene Locus`==locus_tag)]
  column_id<-genes$`Column 384`[which(genes$`"Active" Gene Locus`==locus_tag)]
  plate_id<-genes$`384 plate`[which(genes$`"Active" Gene Locus`==locus_tag)][2]
  plate_id_range<-seq((plate_id-1)*4+1, plate_id*4)
  time_points<-c("13h","19h","39h","44h","69h","89h")
  
  image_all_timepoints<-c()
  for(k in seq(length(time_points))){
    #print(time_points[k])
    image_cropped_annotated_set<-c()
    for(j in seq(length(plate_id_range))){
      #print(plate_id_range[j])
      tag<-ifelse(plate_id_range[j]<10,paste0("00",plate_id_range[j]),paste0("0",plate_id_range[j]) )
      plate <- image_read(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/Morphologyscreen/TimeImages/", time_points[k],"/TB-1p_",tag,".JPG.grid.jpg"), density =2000)
      plate_scale<- image_scale(plate, "x2000") 
      image_scale(plate, "x2000") 
      
      column<-24
      row<-16
      
      lambda_column<-3011/column
      lambda_row<-2000/row
      crop_info<-paste0(lambda_column,"x",lambda_row , "+", (column_id-1)*lambda_column,"+", (row_id-1)*lambda_row)
      
      image_crop(plate_scale,crop_info)
      image_scale(image_crop(plate_scale, crop_info), "x150")
      image_cropped<-image_crop(plate_scale, crop_info)
      image_cropped_annotated<-image_annotate(image_cropped, time_points[k], size = 15, color = "red", boxcolor = "pink", degrees = 0, location = "+45")
      image_cropped_annotated_set<-c(image_cropped_annotated_set, image_cropped_annotated)
    }
    appended_image_row<-image_append(c(image_cropped_annotated_set[[1]], image_cropped_annotated_set[[2]], image_cropped_annotated_set[[3]], image_cropped_annotated_set[[4]]))
    image_write(appended_image_row, path = paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/tmp_", time_points[k],".jpg"), format = "jpeg")
    image_tmp <- image_read(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/tmp_", time_points[k],".jpg"), density =2000)
    image_all_timepoints<-c(image_all_timepoints,image_tmp )
  }
  image_output<-image_append(c(image_all_timepoints[[1]],image_all_timepoints[[2]],image_all_timepoints[[3]],image_all_timepoints[[4]],image_all_timepoints[[5]],image_all_timepoints[[6]]), stack = T)
  image_write(image_output, path = paste0("/Users/moradigd/Documents/chemicalgenomic/new_2_image_output_", locus_tag,".jpg"), format = "jpeg")
}

image_producer("PA14_42220")

genes_unique<-unique(genes$`"Active" Gene Locus`)[-1]

#Until 70%
for(i in seq(1,length(genes_unique))){
  print(i/length(genes_unique))
  print(genes_unique[i])
  image_producer(genes_unique[i])
}


#Extract data 
library(tidyverse)

genes<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/PA14_computational_v2_updated_annotations.csv")
unique_genes<-unique(genes$`"Active" Gene Locus`)

value_extractor<-function(locus_tag){
  genes<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/PA14_computational_v2_updated_annotations.csv",show_col_types = FALSE)
  row_id<-genes$`Row 384`[which(genes$`"Active" Gene Locus`==locus_tag)]
  column_id<-genes$`Column 384`[which(genes$`"Active" Gene Locus`==locus_tag)]
  plate_id<-genes$`384 plate`[which(genes$`"Active" Gene Locus`==locus_tag)][1]
  plate_id_range<-seq((plate_id-1)*4+1, plate_id*4)
  time_points<-c("13h","19h","39h","44h","69h","89h")
  
  data_all_timepoints<-c()
  for(k in seq(length(time_points))){
    for(j in seq(length(plate_id_range))){
      tag<-ifelse(plate_id_range[j]<10,paste0("00",plate_id_range[j]),paste0("0",plate_id_range[j]) )
      plate <- read_tsv(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/Morphologyscreen/TimeImages/", time_points[k],"/TB-1p_",tag,".JPG.iris"), skip = 6,show_col_types = FALSE)
      for(m in 1:length(row_id)){
        tmp_value<-filter(plate, row==row_id[m]) %>% filter(column==column_id[m]) %>% select(`biofilm color intensity`) %>% pull(1)
        data_all_timepoints<-rbind( data_all_timepoints,  c(time_points[k],locus_tag, row_id[m],  column_id[m],tmp_value ))
      }
    }
  }
  data_all_timepoints<-data.frame(data_all_timepoints)
  colnames(data_all_timepoints)<-c("time_point","gene","column","row","value")
  return(data_all_timepoints)
}


final_outputs<-c()
for(j in 2:length(unique_genes)){
  print(j/length(unique_genes))
  final_outputs<-rbind(final_outputs, value_extractor(unique_genes[j]))
}
unique(final_outputs$value)
#write_csv(final_outputs, "/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv")




resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")




#final_outputs<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_size_new.csv")
input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/PAmorphology.csv",show_col_types = FALSE) %>% select(PA14.ID, colony.size, timepoint) %>% group_by(PA14.ID, timepoint) %>% summarise(data=median(colony.size))
final_outputs <- final_outputs %>% group_by(time_point , gene) %>% summarise(new_data=median(as.numeric(  as.character(value))))

plot(final_outputs$new_data, input$data[  match(  paste0(final_outputs$gene, final_outputs$time_point) , paste0(input$PA14.ID, input$timepoint) ) ]  )
cor(final_outputs$new_data, input$data[  match(  paste0(final_outputs$gene, final_outputs$time_point) , paste0(input$PA14.ID, input$timepoint) ) ] , use="pairwise.complete.obs" )
final_outputs$value[final_outputs$gene=="PA14_37530"]

#Extract colonies
library(magick)
library(tidyverse)

locus_tag<-"PA14_03370"
genes<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/PA14_computational_v2_updated_annotations.csv")
row_id<-genes$`Row 384`[which(genes$`"Active" Gene Locus`==locus_tag)]
column_id<-genes$`Column 384`[which(genes$`"Active" Gene Locus`==locus_tag)]
plate_id<-genes$`384 plate`[which(genes$`"Active" Gene Locus`==locus_tag)][1]
plate_id_range<-seq((plate_id-1)*4+1, plate_id*4)
time_points<-c("13h","19h","39h","44h","69h","89h")

image_all_timepoints<-c()
for(k in seq(length(time_points))){
  print(time_points[k])
  image_cropped_annotated_set<-c()
  for(j in seq(length(plate_id_range))){
    print(plate_id_range[j])
    tag<-ifelse(plate_id_range[j]<10,paste0("00",plate_id_range[j]),paste0("0",plate_id_range[j]) )
    plate <- image_read(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/Morphologyscreen/TimeImages/", time_points[k],"/TB-1p_",tag,".JPG.grid.jpg"), density =2000)
    plate_scale<- image_scale(plate, "x2000") 
    image_scale(plate, "x2000") 
    
    column<-24
    row<-16
    
    lambda_column<-3011/column
    lambda_row<-2000/row
    crop_info<-paste0(lambda_column,"x",lambda_row , "+", (column_id-1)*lambda_column,"+", (row_id-1)*lambda_row)
    
    image_crop(plate_scale,crop_info)
    image_scale(image_crop(plate_scale, crop_info), "x150")
    image_cropped<-image_crop(plate_scale, crop_info)
    image_cropped_annotated<-image_annotate(image_cropped, time_points[k], size = 15, color = "red", boxcolor = "pink", degrees = 0, location = "+45")
    image_cropped_annotated_set<-c(image_cropped_annotated_set, image_cropped_annotated)
  }
  appended_image_row<-image_append(c(image_cropped_annotated_set[[1]], image_cropped_annotated_set[[2]], image_cropped_annotated_set[[3]], image_cropped_annotated_set[[4]]))
  image_write(appended_image_row, path = paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/tmp_", time_points[k],".jpg"), format = "jpeg")
  image_tmp <- image_read(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/tmp_", time_points[k],".jpg"), density =2000)
  image_all_timepoints<-c(image_all_timepoints,image_tmp )
}
image_output<-image_append(c(image_all_timepoints[[1]],image_all_timepoints[[2]],image_all_timepoints[[3]],image_all_timepoints[[4]],image_all_timepoints[[5]],image_all_timepoints[[6]]), stack = T)
image_write(image_output, path = paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/output_", locus_tag,".jpg"), format = "jpeg")

#morphology 

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints.csv")
input<-input[input$PA14.ID!="EMPTY",]
input<-input %>% 
  group_by(PA14.ID,timepoint) %>% 
  summarise(colonysize_mean=mean(colony.size),circularity_mean=mean(colony.circularity), morphology_mean=mean(morphology.score.fixed.circles) )%>%
  group_by(PA14.ID) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = c(colonysize_mean, circularity_mean, morphology_mean), values_fill=0) %>%
  summarise(across(everything(), list(sum))) %>%
  select(-row_1)%>%
  rename_all(~ str_remove(., "h_1")) %>%
  write_csv("/Users/moradigd/Documents/chemicalgenomic/aggegated_morphology_mean.csv")


input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints.csv")
input<-input[input$PA14.ID!="EMPTY",]
input<-input %>% 
  group_by(PA14.ID,timepoint) %>% 
  summarise(colonysize_mean=sd(colony.size),circularity_mean=sd(colony.circularity), morphology_mean=sd(morphology.score.fixed.circles) )%>%
  group_by(PA14.ID) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = c(colonysize_mean, circularity_mean, morphology_mean), values_fill=0) %>%
  summarise(across(everything(), list(sum))) %>%
  select(-row_1)%>%
  rename_all(~ str_remove(., "h_1")) %>%
  write_csv("/Users/moradigd/Documents/chemicalgenomic/aggegated_morphology_sd.csv")


input_tmp<- input[input$PA14.ID==locus_tag & input$timepoint==time_points[1],]
output<-c()
output<-c(output, input_tmp$colony.size, input_tmp$colony.circularity, input_tmp$morphology.score.fixed.circles)
paste0("seq_",seq(1,4),"_",colnames(input_tmp)[7])
paste0("seq_",seq(1,4),"_", seq(1,4))
head(input)







frink <- image_read(paste0("/Users/moradigd/Documents/chemicalgenomic/morphology/tmp_", time_point,".jpg"), density =2000)
print(frink )

im2<-image_append(c(im1,im1,im1,im1,im1))

im3<-image_annotate(im1, "13h", size = 30, color = "red", boxcolor = "pink", degrees = 0, location = "+17")
image_append(c(im3,im3,im1,im1,im1))

image_write(image_append(c(im3,im3,im1,im1,im1)), path = "/Users/moradigd/Documents/chemicalgenomic/morphology/tiger.jpg", format = "jpeg")

frink <- image_read("/Users/moradigd/Documents/chemicalgenomic/morphology/tiger.jpg", density =2000)
print(frink )
image_append(c(frink,frink), stack = T)



g<-correlation_summary %>% arrange(slope) %>% select(PA14.ID, slope) %>% inner_join(genes_functional_group, by=c("PA14.ID"="Active Gene Locus")) %>%
  select('Active Gene Name', 'Active Gene Description','go') %>%
  distinct() 
grid.table(g[1:10,])
dev.off()



#correlaation/Figure
library(tidyverse)
library(gmodels)
library(corrplot)

input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints.csv")
input<-input[input$PA14.ID!="EMPTY",]

input_sum<- input %>%
  group_by(PA14.ID, timepoint) %>%
  summarise(mean_cl=mean(colony.size), min_cl=min(colony.size), max_cl=max(colony.size), lowCI = ci(colony.size, confidence=0.99)[2],hiCI = ci(colony.size, confidence=0.99)[3])

#"13h", "19h", "39h", "44h", "69h", "89h"

input_sum_colonysize_tmp_13h<-input_sum[input_sum$timepoint=="13h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_colonysize_tmp_13h)[8]<-"ranking_colonysize_13h"
input_sum_colonysize_tmp_19h<-input_sum[input_sum$timepoint=="19h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_colonysize_tmp_19h)[8]<-"ranking_colonysize_19h"
input_sum_colonysize_tmp_39h<-input_sum[input_sum$timepoint=="39h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_colonysize_tmp_39h)[8]<-"ranking_colonysize_39h"
input_sum_colonysize_tmp_44h<-input_sum[input_sum$timepoint=="44h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_colonysize_tmp_44h)[8]<-"ranking_colonysize_44h"
input_sum_colonysize_tmp_69h<-input_sum[input_sum$timepoint=="69h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_colonysize_tmp_69h)[8]<-"ranking_colonysize_69h"
input_sum_colonysize_tmp_89h<-input_sum[input_sum$timepoint=="89h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_colonysize_tmp_89h)[8]<-"ranking_colonysize_89h"


correlation_summary<-list(input_sum_colonysize_tmp_13h, input_sum_colonysize_tmp_19h, input_sum_colonysize_tmp_39h, input_sum_colonysize_tmp_44h,input_sum_colonysize_tmp_69h, input_sum_colonysize_tmp_89h ) %>% reduce(inner_join, by = "PA14.ID") %>%
  select(ranking_colonysize_13h,ranking_colonysize_19h,ranking_colonysize_39h, ranking_colonysize_44h,ranking_colonysize_69h, ranking_colonysize_89h ) 
correlation_summary<-correlation_summary[,-1]
cor(correlation_summary) %>% corrplot(.)

input_sum<- input %>%
  group_by(PA14.ID, timepoint) %>%
  summarise(mean_cl=mean(colony.circularity), min_cl=min(colony.circularity), max_cl=max(colony.circularity), lowCI = ci(colony.circularity, confidence=0.99)[2],hiCI = ci(colony.circularity, confidence=0.99)[3])

#"13h", "19h", "39h", "44h", "69h", "89h"
input_sum_circularity_tmp_13h<-input_sum[input_sum$timepoint=="13h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_circularity_tmp_13h)[8]<-"ranking_circularity_13h"
input_sum_circularity_tmp_19h<-input_sum[input_sum$timepoint=="19h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_circularity_tmp_19h)[8]<-"ranking_circularity_19h"
input_sum_circularity_tmp_39h<-input_sum[input_sum$timepoint=="39h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_circularity_tmp_39h)[8]<-"ranking_circularity_39h"
input_sum_circularity_tmp_44h<-input_sum[input_sum$timepoint=="44h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_circularity_tmp_44h)[8]<-"ranking_circularity_44h"
input_sum_circularity_tmp_69h<-input_sum[input_sum$timepoint=="69h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_circularity_tmp_69h)[8]<-"ranking_circularity_69h"
input_sum_circularity_tmp_89h<-input_sum[input_sum$timepoint=="89h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$timepoint=="69h",])[1]))
colnames(input_sum_circularity_tmp_89h)[8]<-"ranking_circularity_89h"


correlation_summary<-list(input_sum_circularity_tmp_13h, input_sum_circularity_tmp_19h, input_sum_circularity_tmp_39h, input_sum_circularity_tmp_44h,input_sum_circularity_tmp_69h, input_sum_circularity_tmp_89h ) %>% reduce(inner_join, by = "PA14.ID") %>%
  select(ranking_circularity_13h,ranking_circularity_19h,ranking_circularity_39h, ranking_circularity_44h,ranking_circularity_69h, ranking_circularity_89h ) 
correlation_summary<-correlation_summary[,-1]
cor(correlation_summary) %>% corrplot(.)

#list(input_sum_morphology_tmp_13h, input_sum_morphology_tmp_19h, input_sum_morphology_tmp_39h, input_sum_morphology_tmp_44h,input_sum_morphology_tmp_69h, input_sum_morphology_tmp_89h ) %>% reduce(inner_join, by = "PA14.ID") %>%
#  select(ranking_morphology_13h,ranking_morphology_19h,ranking_morphology_39h, ranking_morphology_44h,ranking_morphology_69h, ranking_morphology_89h ) 

correlation_summary<-list(input_sum_morphology_tmp_13h, input_sum_morphology_tmp_19h, input_sum_morphology_tmp_39h, input_sum_morphology_tmp_44h,input_sum_morphology_tmp_69h, input_sum_morphology_tmp_89h, input_sum_colonysize_tmp_13h, input_sum_colonysize_tmp_19h, input_sum_colonysize_tmp_39h, input_sum_colonysize_tmp_44h,input_sum_colonysize_tmp_69h, input_sum_colonysize_tmp_89h, input_sum_circularity_tmp_13h, input_sum_circularity_tmp_19h, input_sum_circularity_tmp_39h, input_sum_circularity_tmp_44h,input_sum_circularity_tmp_69h, input_sum_circularity_tmp_89h ) %>% reduce(inner_join, by = "PA14.ID") %>%
  select(ranking_morphology_13h,ranking_morphology_19h,ranking_morphology_39h, ranking_morphology_44h,ranking_morphology_69h, ranking_morphology_89h, ranking_colonysize_13h,ranking_colonysize_19h,ranking_colonysize_39h, ranking_colonysize_44h,ranking_colonysize_69h, ranking_colonysize_89h, ranking_circularity_13h,ranking_circularity_19h,ranking_circularity_39h, ranking_circularity_44h,ranking_circularity_69h, ranking_circularity_89h ) 
correlation_summary<-correlation_summary[,-1]
cor(correlation_summary) %>% corrplot(.)

#color


input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colony_colour_new.csv")
input<-input[input$gene!="EMPTY",]

input_sum<- input %>%
  group_by(gene, time_point) %>%
  summarise(mean_cl=median(value), min_cl=min(value), max_cl=max(value), lowCI = ci(value, confidence=0.99)[2],hiCI = ci(value, confidence=0.99)[3])

#"13h", "19h", "39h", "44h", "69h", "89h"

input_sum_color_tmp_13h<-input_sum[input_sum$time_point=="13h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$time_point=="69h",])[1])) %>% rename("ranking_color_13h" = colnames(.)[8] )
input_sum_color_tmp_19h<-input_sum[input_sum$time_point=="19h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$time_point=="69h",])[1])) %>% rename("ranking_color_19h" = colnames(.)[8] )
input_sum_color_tmp_39h<-input_sum[input_sum$time_point=="39h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$time_point=="69h",])[1])) %>% rename("ranking_color_39h" = colnames(.)[8] )
input_sum_color_tmp_44h<-input_sum[input_sum$time_point=="44h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$time_point=="69h",])[1])) %>% rename("ranking_color_44h" = colnames(.)[8] )
input_sum_color_tmp_69h<-input_sum[input_sum$time_point=="69h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$time_point=="69h",])[1])) %>% rename("ranking_color_69h" = colnames(.)[8] )
input_sum_color_tmp_89h<-input_sum[input_sum$time_point=="89h",] %>% arrange(desc(mean_cl))  %>% add_column(seq(dim(input_sum[input_sum$time_point=="69h",])[1])) %>% rename("ranking_color_89h" = colnames(.)[8] )

input_sum_color_tmp_13h


correlation_summary<-list(input_sum_color_tmp_13h, input_sum_color_tmp_19h, input_sum_color_tmp_39h, input_sum_color_tmp_44h,input_sum_color_tmp_69h, input_sum_color_tmp_89h ) %>% reduce(inner_join, by = "gene") %>%
  select(ranking_color_13h,ranking_color_19h,ranking_color_39h, ranking_color_44h,ranking_color_69h, ranking_color_89h ) 
correlation_summary<-correlation_summary[,-1]
cor(correlation_summary, method="spearman") %>% corrplot(.)
?cor

#Autocorrelation color 

library(tseries)
gene_ids<-unique(input_sum$gene)
autocorr_df<-c()
for(gene_id in gene_ids){
  #print(match(gene_id, gene_ids))
  tmp<-input_sum$mean_cl[input_sum$gene==gene_id]
  autocorr_df<-rbind(autocorr_df, cbind(acf(tmp, pl=FALSE, lag=2)$acf[2:3], seq(2,3), rep(gene_id,2)))
}
autocorr_df<-data.frame(autocorr_df)
autocorr_df$X1<-as.numeric(as.character(autocorr_df$X1))
ggplot(autocorr_df, aes(x=X2, y=X1)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  #geom_jitter(height = 0, width = 0.4,alpha = 1/50, col="red")+
  ylim(range(-1,1))+
  theme_bw()+
  theme(axis.text.x = element_text( size=15, angle=90,hjust = 1),
        axis.text.y = element_text( size=15, hjust = 1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  ylab("Autocorrelation")+
  xlab("Lag")


acf(c(0.1, 0.5, 0.85, 0.99))
acf(c(1, 400, 3, 400))

#Autocorrelation 
input_sum<- input %>%
  group_by(PA14.ID, timepoint) %>%
  summarise(mean_cl=mean(colony.circularity), min_cl=min(morphology.score.fixed.circles), max_cl=max(colony.size), lowCI = ci(colony.size, confidence=0.99)[2],hiCI = ci(colony.size, confidence=0.99)[3])

library(tseries)
gene_ids<-unique(input_sum$PA14.ID)
autocorr_df<-c()
for(gene_id in gene_ids){
  print(match(gene_id, gene_ids))
  tmp<-input_sum$mean_cl[input_sum$PA14.ID==gene_id]
  autocorr_df<-rbind(autocorr_df, cbind(acf(tmp, pl=FALSE, lag=2)$acf[2:3], seq(2,3), rep(gene_id,2)))
}
autocorr_df<-data.frame(autocorr_df)
autocorr_df$X1<-as.numeric(as.character(autocorr_df$X1))
ggplot(autocorr_df, aes(x=X2, y=X1)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  #geom_jitter(height = 0, width = 0.4,alpha = 1/50, col="red")+
  ylim(range(-1,1))+
  theme_bw()+
  theme(axis.text.x = element_text( size=15, angle=90,hjust = 1),
        axis.text.y = element_text( size=15, hjust = 1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  ylab("Autocorrelation")+
  xlab("Lag")
#calculate autocorrelations










ranking<-c()
for(i in seq(dim(input_sum_tmp)[1])){
  print(i/dim(input_sum_tmp)[1])
  tmp<-length(which(input_sum_tmp$lowCI[i]<input_sum_tmp$hiCI[seq(dim(input_sum_tmp)[1])[-1]]))
  ranking<-c(ranking, tmp)
}

ranking<-1
while(i < (dim(input_sum_tmp)[1]-1)){
  j<-i+1
  while(input_sum_tmp$lowCI[i]<input_sum_tmp$hiCI[j]){
    input_sum_tmp$ranking[j]<-ranking
    j<-j+1
  }
  ranking<-ranking+1
  i<-i+1
}


ranking<-c()
for(i in seq(1,dim(input_sum_tmp)[1])){
  count<-0
  for(j in seq(1,dim(input_sum_tmp)[1])){
    if(input_sum_tmp$lowCI[i]<input_sum_tmp$hiCI[j]){
      count<-count+1
    }
  }
  ranking<-c(ranking,count)
}

predict(as.lm(fitmodel), interval = "confidence")


library(inflection)
cc=bese(x1,y2,0)
cc$iplast
plot(x1,y2,cex=0.3,pch=19,type="l", ylim=range(0,100))
points(x_timepoints,y_growth)
abline(v=cc$iplast,col='blue')

unique(input$timepoint[input$PA14.ID==unique_loci[2]])

x=c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6,
    6, 6)
x<-x+runif(n=length(x), min=1e-12, max=.01)


y = c(0,  0,  0,  0, 22, 11, 17, 12,  5, 23, 16, 14,  3, 64, 72, 75, 49,
      45, 51, 25, 93, 59, 22, 20)

sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

start_val_sigmoid <- function(x, y) {
  fit <- lm(log(y[which.max(x)] - y + 1e-6) ~ x)
  list(
    a = y[which.max(x)],
    b = unname(-coef(fit)[2]),
    c = unname(-coef(fit)[1] / coef(fit)[2]))
}
start_val_sigmoid(x, y)

plot(y)
fitmodel <- nls(y ~ a / ( 1 + exp(-b * (x - c))),start=start_val_sigmoid(x, y), algorithm = "port")

#nls(y~a/(1 + exp(-b * (x-c))))
# visualization code
# get the coefficients using the coef function
params=coef(fitmodel)



x1<-seq(0,6,0.1)
y2 <- sigmoid(params,x1)


f=expression( 48 / (1 + exp(-8.9* (x - 3.1))))
library(inflection)
cc=bese(x1,y2,0)
cc$iplast
plot(x1,y2,cex=0.3,pch=19,type="l")
abline(v=cc$iplast,col='blue')

y2[ which(abs(diff(y2))==max(abs(diff(y2))) )]

#variance between trends 
#variance between time series for one gene

#install.packages("forecast")
library(forecast)
#install.packages("fpp")
library(fpp)
y = c( 0, 12, 14, 25, 20,75)
plot(y)
data(ausbeer)
timeserie_beer = tail(head(ausbeer, 17*4+2),17*4-4)
plot(as.ts(timeserie_beer))

trend_beer = ma(y, order = 2, centre = T)
plot(as.ts(y))
lines(trend_beer)
plot(as.ts(y))

#predictive power 

library(caret)
library(tidyverse)
tmp<-read_csv("/Users/moradigd/Documents/chemicalgenomic/colour_distance_mat.csv")
tmp<-tmp[,-1]
tmp<-as.matrix(tmp)

f <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}
m<-f(tmp)
is.na(m) <- sapply(m, is.infinite)
hist(m)

write_csv(data.frame(m), "/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated_colony_colour_new.csv")

genes_functional_group<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated_colony_size.csv")
genes_functional_group_filtered<-genes_functional_group %>% 
  group_by(`Active Gene Locus`) %>%
  summarise(variance=var(colony.size, na.rm = T)) %>%
  filter(variance>0.1)

genes_functional_group$colony.size[which(genes_functional_group$`Active Gene Locus`=="PA14_55770")]
genes_functional_group_filtered[which(genes_functional_group_filtered$`Active Gene Locus`=="PA14_55770"),]

hist(genes_functional_group_filtered$variance, breaks = 500)

genes_functional_group<-genes_functional_group[! is.na( genes_functional_group$go),]

#filter out low variance 
#gene_list<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated_colony_size.csv")




gene_names<-read_csv("/Users/moradigd/Documents/chemicalgenomic/similarity_colony_size.csv")
shared_genes<-intersect(gene_names$X1, genes_functional_group$`Active Gene Locus`)
m_shortened<-m[match(shared_genes, gene_names$X1),match(shared_genes, gene_names$X1)]
genes_functional_group_shortened<-genes_functional_group[match(shared_genes,genes_functional_group$`Active Gene Locus`),]

labels<-matrix("F", nrow = dim(m_shortened)[1], ncol = dim(m_shortened)[1])
for(i in seq(dim(m_shortened)[1])){
  print(i/dim(m_shortened)[1])
  for(j in seq(dim(m_shortened)[1])){
    if(genes_functional_group_shortened$go[i]==genes_functional_group_shortened$go[j]){
      labels[i,j]<-"T"
    }
  }
}

#operon
gene_names<-read_csv("/Users/moradigd/Documents/chemicalgenomic/similarity_colony_size.csv")
operon_names<-read_tsv("/Users/moradigd/Documents/chemicalgenomic/operon_PA14.txt")
operon_genes<- unique(c(operon_names$SysName1, operon_names$SysName2))
shared_genes<-intersect(operon_genes, genes_functional_group$`Active Gene Locus`)
m_shortened<-m[match(shared_genes, gene_names$X1),match(shared_genes, gene_names$X1)]

labels<-matrix("F", nrow = dim(m_shortened)[1], ncol = dim(m_shortened)[1])
for(i in seq(dim(m_shortened)[1])){
  print(i/dim(m_shortened)[1])
  for(j in seq(dim(m_shortened)[1])){
    id<-which(operon_names$SysName1 %in% shared_genes[i] & operon_names$SysName2 %in% shared_genes[j])
    if(length(id)>0 && operon_names$bOp[id]==TRUE){
      labels[i,j]<-"T"
      labels[j,i]<-"T"
    }
    
    id<-which(operon_names$SysName1 %in% shared_genes[j] & operon_names$SysName2 %in% shared_genes[i])
    if(length(id)>0 && operon_names$bOp[id]==TRUE){
      labels[i,j]<-"T"
      labels[j,i]<-"T"
    }
  }
}
row.names(labels)<-shared_genes
colnames(labels)<-shared_genes
write_csv(data.frame(labels), "/Users/moradigd/Documents/chemicalgenomic/Operon_table.csv")
labels_df<-data.frame(labels)
which(labels_df=="T", arr.ind = T)
median(m_shortened[which(labels_df=="T", arr.ind = T)])
hist(m_shortened[which(labels_df=="F", arr.ind = T)])
range(m_shortened[which(labels_df=="F", arr.ind = T)],na.rm = T)
mean(m_shortened,na.rm = T )


m_shortened[which(m_shortened == sort(m_shortened,decreasing=T)[4], arr.ind = TRUE)]
shared_genes[which(m_shortened == sort(m_shortened,decreasing=T)[50], arr.ind = TRUE)]

sort(m_shortened,decreasing=T)[4]



table(as.vector(labels))
cutoff_range<-seq(0.01,0.95,0.05)*max(1000)
TPR<-c()
FPR<-c()
for(j in cutoff_range){
  print(j)
  m_shortened_cutoff<-ifelse(m_shortened>j,"T","F")
  print(table(as.vector(m_shortened_cutoff)))
  conf_tmp<-confusionMatrix(as.factor(as.vector(labels)), as.factor(as.vector(m_shortened_cutoff)), positive = "T")
  TPR<-c(TPR, conf_tmp[[4]][1])
  FPR<-c(FPR, 1-conf_tmp[[4]][2])
  print(conf_tmp[[4]])
}
plot(FPR,TPR, ylim=range(0,1))


groups<-sort(table(genes_functional_group_shortened$go), decreasing=TRUE)
go_term<-c()
similarity_score<-c()
for(i in seq(20)){
  tmp<-m_shortened[which(genes_functional_group_shortened$go %in% names(groups)[i]),which(genes_functional_group_shortened$go %in% names(groups)[i])]
  similarity_score<-c(similarity_score, as.vector(tmp)) 
  go_term<-c(go_term, rep(names(groups)[i], length(as.vector(tmp))))
}

table_sum<-table(genes_functional_group_shortened$go, genes_functional_group_shortened$cluster)


boxplot(similarity_score ~ go_term)

table_freq<- apply(table_sum[which( rownames(table_sum) %in% names(groups)[1:20] ),],2,function(x){x/sum(x)})
table_freq<-data.frame(table_freq )
table_freq$GOlabel<-row.names(table_freq)
pivot_longer(table_freq, 
             cols = c("X1","X2","X3","X4"),
             names_to = "Clusters",
             values_to = "values"
) %>% ggplot( aes(x=GOlabel, y=values, fill=Clusters)) + 
  geom_bar(stat = 'identity',position=position_dodge()) +
  coord_cartesian(ylim = c(0,0.25))+
  theme_light()+
  ylab("Relative Frequency")+
  theme_light()


#Summary morphology
genes_functional_group_shortened<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated.csv")


groups<-sort(table(genes_functional_group_shortened$go), decreasing=TRUE)
go_term<-c()
similarity_score<-c()
for(i in seq(20)){
  tmp<-genes_functional_group_shortened[which(genes_functional_group_shortened$go %in% names(groups)[i]),which(genes_functional_group_shortened$go %in% names(groups)[i])]
  go_term<-c(go_term, rep(names(groups)[i], length(as.vector(tmp))))
}

table_sum<-table(genes_functional_group_shortened$go, genes_functional_group_shortened$cluster)


boxplot(similarity_score ~ go_term)

table_freq<- apply(table_sum[which( rownames(table_sum) %in% names(groups)[1:20] ),],2,function(x){x/sum(x)})
table_freq<-data.frame(table_freq )
table_freq$GOlabel<-row.names(table_freq)
pivot_longer(table_freq, 
             cols = c("X1","X2","X3","X4","X5","X6"),
             names_to = "Clusters",
             values_to = "values"
) %>% ggplot( aes(x=GOlabel, y=values, fill=Clusters)) + 
  geom_bar(stat = 'identity',position=position_dodge()) +
  coord_cartesian(ylim = c(0,0.25))+
  theme_light()+
  ylab("Relative Frequency")+
  theme_light()

g<-genes_functional_group_shortened[genes_functional_group_shortened$cluster=="2",c(20,21,44)] %>% distinct()
library(gridExtra)
grid.table(g)

#colony_circularity
genes_functional_group_shortened<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology/aggregated_colony_circularity.csv.csv")

groups<-sort(table(genes_functional_group_shortened$go), decreasing=TRUE)
go_term<-c()
similarity_score<-c()
for(i in seq(20)){
  tmp<-genes_functional_group_shortened[which(genes_functional_group_shortened$go %in% names(groups)[i]),which(genes_functional_group_shortened$go %in% names(groups)[i])]
  go_term<-c(go_term, rep(names(groups)[i], length(as.vector(tmp))))
}

table_sum<-table(genes_functional_group_shortened$go, genes_functional_group_shortened$cluster)


boxplot(similarity_score ~ go_term)

table_freq<- apply(table_sum[which( rownames(table_sum) %in% names(groups)[1:20] ),],2,function(x){x/sum(x)})
table_freq<-data.frame(table_freq )
table_freq$GOlabel<-row.names(table_freq)
pivot_longer(table_freq, 
             cols = c("X1","X2","X3","X4"),
             names_to = "Clusters",
             values_to = "values"
) %>% ggplot( aes(x=GOlabel, y=values, fill=Clusters)) + 
  geom_bar(stat = 'identity',position=position_dodge()) +
  coord_cartesian(ylim = c(0,0.25))+
  theme_light()+
  ylab("Relative Frequency")+
  theme_light()

g<-genes_functional_group_shortened[genes_functional_group_shortened$cluster=="4",c(20,21,44)] %>% distinct()
library(gridExtra)
grid.table(g[25:33,])
dev.off()

#detect zero
input<-read_csv("/Users/moradigd/Documents/chemicalgenomic/all.PA.morphology.timepoints.csv")





###plasmid project 
#sample size 
library(tidyverse)
timepoints<-read_csv("/Users/moradigd/Documents/Plasmid/Subset_New_results_screening_agg_scale_5_RP.csv") %>%
  mutate(across('6', as.character))
colnames(timepoints)<-c("index","rho","pvalue","kmer","plasmid","model","counter","subsample","dataset")
ggplot(timepoints, aes(x=subsample, y=rho, fill=dataset)) + 
  geom_boxplot()+
  ylim(range(0,1))+
  theme_bw()+
  theme(axis.text.x = element_text( size=18, angle=90,hjust = 1),
        axis.text.y = element_text( size=18, hjust = 1),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  ylab("Spearman's ρ")+
  xlab("Subsample size")


#idfferent plasmid types
library(tidyverse)
library(tidymodels)
library(tidyverse)
library(gmodels)

plasmid_tmp_1<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_PB_RP_data.csv") %>%
  mutate(plasmid_test="RP",plasmid_train="PB" )

plasmid_tmp_2<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_PB_PKJK_data.csv") %>%
  mutate(plasmid_test="PKJK", plasmid_train="PB")

plasmid_tmp_3<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_PKJK_RP_data.csv") %>%
  mutate(plasmid_test="RP",plasmid_train="PKJK" )

plasmid_tmp_4<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_PKJK_PB_data.csv") %>%
  mutate(plasmid_test="PB", plasmid_train="PKJK")

plasmid_tmp_5<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_RP_PB_data.csv") %>%
  mutate(plasmid_test="PB",plasmid_train="RP" )

plasmid_tmp_6<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_RP_PKJK_data.csv") %>%
  mutate(plasmid_test="PKJK", plasmid_train="RP")

plasmid_tot<-rbind(plasmid_tmp_1, plasmid_tmp_2,plasmid_tmp_3, plasmid_tmp_4,plasmid_tmp_5, plasmid_tmp_6)
colnames(plasmid_tot)<-c("index","rho","pvalue","kmer","model","replicate","plasmid_test","plasmid_train")

plasmid_tot_sum<-group_by(plasmid_tot, plasmid_test,plasmid_train) %>%
  summarise(sd=sd(rho), mean=mean(rho))


ggplot(plasmid_tot_sum, aes(x = plasmid_train, y=mean, fill=plasmid_test))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = plasmid_tot_sum$mean-plasmid_tot_sum$sd, ymax =  plasmid_tot_sum$mean+plasmid_tot_sum$sd), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=23, angle=45,hjust = 1),
        axis.text.y = element_text( size=23, hjust = 1),
        axis.title.x = element_text(color="black", size=25, face="bold"),
        axis.title.y = element_text(color="black", size=25, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Training plasmid")+
  ylim(range(0,0.7))




output<-read_csv("/Users/moradigd/Documents/Plasmid/new_results_screening_agg_scale_tot.csv")
colnames(output)<-c("index","correlation","pvalue","kmer","plasmid","model","run","dataset")

output<-output %>%
  drop_na() %>%
  group_by(  kmer, plasmid ,model,dataset) %>%
  summarise(mean_correlation=mean(correlation) , lowCI = ci(correlation)[2],
            hiCI = ci(correlation)[3], sd_correlation= ci(correlation)[4]) %>%
  mutate(lowCI = replace(lowCI, which(lowCI<0), 0)) %>% 
  filter(dataset=="test") %>%
  group_by(plasmid ,model,dataset) %>%
  summarise(max_correlation=max(mean_correlation),lowCI =min(lowCI) , hiCI=max(hiCI), sd_correlation= max(sd_correlation))



plasmid_baps_sum$mean_correlation<-output$max_correlation-plasmid_baps_sum$mean_correlation
plasmid_baps_sum$lowCI<-plasmid_baps_sum$mean_correlation- sqrt(plasmid_baps_sum$sd_correlation^2+ output$sd_correlation^2)
plasmid_baps_sum$hiCI<-plasmid_baps_sum$mean_correlation+ sqrt(plasmid_baps_sum$sd_correlation^2+ output$sd_correlation^2 )
plasmid_baps_sum$lowCI<-ifelse(plasmid_baps_sum$lowCI<0,0,plasmid_baps_sum$lowCI)

max_ids<-c(3, 6, 9)
output<-output[max_ids,]

plasmid_tot_sum<-inner_join(plasmid_tot_sum,output, by=c("plasmid_train"="plasmid")) %>%
  mutate(mean_corr=max_correlation-mean) %>%
  mutate(lowCI_new=mean_corr-sqrt(sd^2+ sd_correlation^2), highCI_new=mean_corr+sqrt(sd^2+ sd_correlation^2)) %>%
  select(plasmid_test, plasmid_train, mean_corr,lowCI_new, highCI_new )


mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="PB"])
mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="RP"])
mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="PKJK"])

ggplot(plasmid_tot_sum, aes(x = plasmid_train, y=mean_corr, fill=plasmid_test))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin =lowCI_new , ymax = highCI_new), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Plasmid")+
  ylim(range(0,0.5))+
  ylab(expression("Spearman's ρ [Same] - Spearman's ρ [Diff.]"))+
  scale_fill_manual(values=c("#AA0099", "#E63F80","#FFBAA9"))






#Country


plasmid_country_PB<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_PB_country_data.csv") %>%
  mutate(plasmid="PB")
plasmid_country_RP<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_RP_country_data.csv") %>%
  mutate(plasmid="RP")
plasmid_country_PKJK<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_PKJK_country_data.csv") %>%
  mutate(plasmid="PKJK")

plasmid_country_tot<- rbind(plasmid_country_PB,plasmid_country_RP, plasmid_country_PKJK)
colnames(plasmid_country_tot)<-c("index","rho","pvalue","kmer","model","replicate","plasmid")
plasmid_country_tot_sum<-group_by(plasmid_country_tot, plasmid) %>%
  summarise(sd=sd(rho), mean=mean(rho))

ggplot(plasmid_country_tot_sum, aes(x = plasmid, y=mean))+
  geom_bar(stat = "identity", position = "dodge",fill=c("red","green","blue"))+
  geom_errorbar(aes(ymin = plasmid_country_tot_sum$mean-plasmid_country_tot_sum$sd, ymax =  plasmid_country_tot_sum$mean+plasmid_country_tot_sum$sd), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=23, angle=45,hjust = 1),
        axis.text.y = element_text( size=23, hjust = 1),
        axis.title.x = element_text(color="black", size=25, face="bold"),
        axis.title.y = element_text(color="black", size=25, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Plasmid")+
  ylim(range(0,0.7))




output<-read_csv("/Users/moradigd/Documents/Plasmid/new_results_screening_agg_scale_tot.csv")
colnames(output)<-c("index","correlation","pvalue","kmer","plasmid","model","run","dataset")

output<-output %>%
  drop_na() %>%
  group_by(  kmer, plasmid ,model,dataset) %>%
  summarise(mean_correlation=mean(correlation) , lowCI = ci(correlation)[2],
            hiCI = ci(correlation)[3], sd_correlation= ci(correlation)[4]) %>%
  mutate(lowCI = replace(lowCI, which(lowCI<0), 0)) %>% 
  filter(dataset=="test") %>%
  group_by(plasmid ,model,dataset) %>%
  summarise(max_correlation=max(mean_correlation),lowCI =min(lowCI) , hiCI=max(hiCI), sd_correlation= max(sd_correlation))


max_ids<-c(3, 6, 9)
output<-output[max_ids,]


plasmid_tot_sum<-inner_join(plasmid_country_tot_sum,output, by=c("plasmid"="plasmid")) %>%
  mutate(mean_corr=max_correlation-mean) %>%
  mutate(lowCI_new=mean_corr-sqrt(sd^2+ sd_correlation^2), highCI_new=mean_corr+sqrt(sd^2+ sd_correlation^2)) %>%
  select( plasmid, mean_corr,lowCI_new, highCI_new )


mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="PB"])
mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="RP"])
mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="PKJK"])

ggplot(plasmid_tot_sum, aes(x = plasmid, y=mean_corr, fill=plasmid))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin =lowCI_new , ymax = highCI_new), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Plasmid")+
  ylim(range(0,0.5))+
  ylab(expression("Spearman's ρ [Same] - Spearman's ρ [Diff.]"))+
  scale_fill_manual(values=c("#AA0099", "#E63F80","#FFBAA9"))



ggplot(plasmid_tot_sum, aes(x = plasmid, y=mean_corr, fill=plasmid))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin =lowCI_new , ymax = highCI_new), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Plasmid")+
  ylim(range(0,0.5))+
  ylab(expression(""))+
  scale_fill_manual(values=c("#AA0099", "#E63F80","#FFBAA9"))













plasmid_time<-read_csv("/Users/moradigd/Documents/Plasmid/results_input_agg_year.csv") %>%
  mutate(year="year")
colnames(plasmid_time)<-c("index","rho","pvalue","kmer","model","replicate","year")

plasmid_time_sum<-group_by(plasmid_time, year) %>%
  summarise(sd=sd(rho), mean=mean(rho))

plasmid_time_sum<-rbind(plasmid_time_sum,c("mixed", 0.04,0.38))

#sites 
tags<-c("1A","1B","2","4","6","7")
plasmid_sites_tot<-c()

for(i in seq(length(tags))){
  plasmid_sites<- read_csv(paste0("/Users/moradigd/Documents/Plasmid/results_input_agg_sites_PKJK_DK_site_", tags[i], ".csv")) %>%
    mutate(site=tags[i])
  plasmid_sites_tot<-rbind(plasmid_sites_tot, plasmid_sites)
}
colnames(plasmid_sites_tot)<-c("index","rho","pvalue","kmer","model","replicate","site")

plasmid_sites_tot_sum<-group_by(plasmid_sites_tot, site) %>%
  summarise(sd=sd(rho), mean=mean(rho))


ggplot(plasmid_sites_tot_sum, aes(x = site, y=mean))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = plasmid_sites_tot_sum$mean-plasmid_sites_tot_sum$sd, ymax =  plasmid_sites_tot_sum$mean+plasmid_sites_tot_sum$sd), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=23, angle=45,hjust = 1),
        axis.text.y = element_text( size=23, hjust = 1),
        axis.title.x = element_text(color="black", size=25, face="bold"),
        axis.title.y = element_text(color="black", size=25, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Plasmid")+
  ylim(range(0,0.7))


output<-output[2, ]

plasmid_sites_tot_sum$mean_corr<-output$max_correlation-plasmid_sites_tot_sum$mean
plasmid_sites_tot_sum$lowCI_new<-plasmid_sites_tot_sum$mean_corr-sqrt(plasmid_sites_tot_sum$sd^2+ output$sd_correlation^2)
plasmid_sites_tot_sum$highCI_new<-plasmid_sites_tot_sum$mean_corr+sqrt(plasmid_sites_tot_sum$sd^2+ output$sd_correlation^2)

ggplot(plasmid_sites_tot_sum, aes(x = site, y=mean_corr))+
  geom_bar(stat = "identity", position = "dodge",fill="#E63F80")+
  geom_errorbar(aes(ymin =lowCI_new , ymax = highCI_new), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Plasmid")+
  ylim(range(0,0.5))+
  ylab(expression("Spearman's ρ [Same] - Spearman's ρ [Diff.]"))


ggplot(plasmid_sites_tot_sum, aes(x = site, y=mean_corr))+
  geom_bar(stat = "identity", position = "dodge",fill="#E63F80")+
  geom_errorbar(aes(ymin =lowCI_new , ymax = highCI_new), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("WWTP sites")+
  ylim(range(0,0.5))+
  ylab(expression(""))



#BAPS
library(tidyverse)
library(gmodels)
library(tidyverse)

plasmid_baps_RP<-read_csv("/Users/moradigd/Documents/Plasmid/RP_baps_results_screening_agg_scale.csv") 
plasmid_baps_PB<-read_csv("/Users/moradigd/Documents/Plasmid/PB_baps_results_screening_agg_scale.csv") 
plasmid_baps_PKJK<-read_csv("/Users/moradigd/Documents/Plasmid/PKJK_baps_results_screening_agg_scale.csv") 
plasmid_baps<-rbind(plasmid_baps_RP, plasmid_baps_PB, plasmid_baps_PKJK)

colnames(plasmid_baps)<-c("index","rho","pvalue","plasmid","model","replicate")

plasmid_baps_sum<- plasmid_baps %>%
  dplyr::group_by(plasmid, model) %>%
  dplyr::summarise(mean_correlation=mean(rho, na.rm=T) , lowCI = ci(rho, na.rm=T)[2],
                   hiCI = ci(rho, na.rm=T)[3], sd_correlation= ci(rho, na.rm=T)[4]) %>%
  dplyr::mutate(lowCI = replace(lowCI, which(lowCI<0), 0))


output<-read_csv("/Users/moradigd/Documents/Plasmid/new_results_screening_agg_scale_tot.csv")
colnames(output)<-c("index","correlation","pvalue","kmer","plasmid","model","run","dataset")

output<-output %>%
  drop_na() %>%
  dplyr::group_by(  kmer, plasmid ,model,dataset) %>%
  dplyr::summarise(mean_correlation=mean(correlation) , lowCI = ci(correlation)[2],
                   hiCI = ci(correlation)[3], sd_correlation= ci(correlation)[4]) %>%
  dplyr::mutate(lowCI = replace(lowCI, which(lowCI<0), 0)) %>% 
  dplyr::filter(dataset=="test") %>%
  dplyr::group_by(plasmid ,model,dataset) %>%
  dplyr::summarise(max_correlation=max(mean_correlation),lowCI =min(lowCI) , hiCI=max(hiCI), sd_correlation= max(sd_correlation))


plasmid_baps_sum$mean_correlation<-output$max_correlation-plasmid_baps_sum$mean_correlation
plasmid_baps_sum$lowCI<-plasmid_baps_sum$mean_correlation- sqrt(plasmid_baps_sum$sd_correlation^2+ output$sd_correlation^2)
plasmid_baps_sum$hiCI<-plasmid_baps_sum$mean_correlation+ sqrt(plasmid_baps_sum$sd_correlation^2+ output$sd_correlation^2 )
plasmid_baps_sum$lowCI<-ifelse(plasmid_baps_sum$lowCI<0,0,plasmid_baps_sum$lowCI)

mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="PB"])
mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="RP"])
mean(plasmid_baps_sum$mean_correlation[plasmid_baps_sum$plasmid=="PKJK"])

plasmid_baps_sum$plasmid_new<-plasmid_baps_sum$plasmid
plasmid_baps_sum$plasmid_new[plasmid_baps_sum$plasmid_new=="PB"]<-"pB10"
plasmid_baps_sum$plasmid_new[plasmid_baps_sum$plasmid_new=="PKJK"]<-"pKJK5"
plasmid_baps_sum$plasmid_new[plasmid_baps_sum$plasmid_new=="RP"]<-"RP4"

ggplot(plasmid_baps_sum, aes(x = plasmid_new, y=mean_correlation, fill=model))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin =lowCI , ymax = hiCI), width=.2,position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text( size=18, angle=45,hjust = 1),
        axis.text.y = element_text( size=18, hjust = 1),
        axis.title.x = element_text(color="black", size=17, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Plasmid")+
  ylim(range(0,0.4))+
  ylab(expression("Spearman's ρ [kmer] - Spearman's ρ [BAPS]"))+
  scale_fill_manual(values=c("#998899", "#E69F00","#FF0000"))

importance<-c()
for(j in seq(9)){
  importance<-rbind(importance, read_csv(paste0("/Users/moradigd/Documents/Plasmid/RF_RP_baps_importance_",j,"_RP")) %>%
                      arrange(desc(importance))%>% 
                      mutate(rank = dense_rank(desc(importance))))
}
importance_RP<- importance %>%
  dplyr::group_by(tag) %>%
  dplyr::summarise(mean_importance=mean(importance), mean_rank=mean(rank), count_n=n()) %>%
  dplyr::arrange(mean_rank) %>%
  dplyr::filter(count_n>8)


#Correlation test 
library(tidyverse)
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 

df<-df[!is.na(df$RP),]
baps<-read_csv("/Users/moradigd/Documents/Plasmid/RP_BAPS.csv")

#count the number of occurrence 
sum_list<-c()
for(i in importance_RP$tag){
  sum_list<-c(sum_list, sum(apply(X=baps,2,FUN=function(x) length(which(x==i)))))
}
importance_RP$sum<-sum_list
importance_RP<-importance_RP[importance_RP$sum>5,]

ggplot(importance_RP, aes(x=mean_importance, y=mean_rank, size =sum/10 , col= count_n/5 )) +
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 24))+
  ylim(range(-5,60))+
  xlim(range(-0.01,0.25))+
  theme_bw()


mean(df[which(baps==importance_RP$tag[1], arr.ind=TRUE)[,1],3])
mean(df[which(baps!=importance_RP$tag[1], arr.ind=TRUE)[,1],3])
#compare significance 
significance<-c()
mean_diff<-c()

results<-c()
for(j in seq(length(importance_RP$tag))){
  print(j)
  with<-df[which(baps==importance_RP$tag[j], arr.ind=TRUE)[,1],3]
  without<-df[which(baps!=importance_RP$tag[j], arr.ind=TRUE)[,1],3]
  list_tmp<-cbind(c(with, without), rep(importance_RP$tag[j], length(c(with, without))), c(rep("with", length(with)), rep("without", length(without))))
  
  results<-rbind(results,list_tmp)
  significance<-c(significance, wilcox.test(with, without)$p.value)
  mean_diff<-c(mean_diff, mean(with)/mean(without))
}
importance_RP$mean_diff<-mean_diff
importance_RP$significance<-significance
importance_RP_filtered<- importance_RP[importance_RP$significance<0.05 & importance_RP$mean_diff>1 ,]
#importance_RP_filtered<- importance_RP


importance_RP[importance_RP$significance<0.01 ,]

similarity<-matrix(0,nrow = dim(importance_RP_filtered)[1], ncol = dim(importance_RP_filtered)[1] )
for(i in 1:dim(importance_RP_filtered)[1]){
  for(j in 1:dim(importance_RP_filtered)[1]){
    tmp1<-which(baps==importance_RP$tag[i], arr.ind=TRUE)[,1]
    tmp2<-which(baps==importance_RP$tag[j], arr.ind=TRUE)[,1]
    similarity[i,j]<-length(intersect(tmp1, tmp2))/length(union(tmp1, tmp2))
  }
}

library(corrplot)
corrplot( similarity, method = "circle", type="upper")




results<-data.frame(results)
results$X1<-as.numeric(as.character(results$X1))
results<-results[which(results$X2 %in% importance_RP_filtered$tag),]
results$X2<-paste0("", seq(length(unique(results$X2))))[match(results$X2, unique(results$X2))]
ggplot(results, aes(fill=X3, y=X1, x=X2)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent")+
  theme_bw()+
  theme(axis.text.x = element_text( size=18, angle=45,hjust = 1),
        axis.text.y = element_text( size=18, hjust = 1),
        axis.title.x = element_text(color="black", size=17, face="bold"),
        axis.title.y = element_text(color="black", size=17, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black"
        ))+
  xlab("BAPS Groups")+
  ylim(range(0,125))+
  ylab(expression("Permissiveness"))

barplot(importance_RP_filtered$sum, col="black", names.arg = importance_RP_filtered$tag, ylab = "Frequency", las=2)
box()

results_frequency<-c()
for(j in seq(dim(importance_RP_filtered)[1])){
  #for(j in seq(1,1)){
  print(importance_RP_filtered$tag[j])
  seq_RP<-df$seq[as.numeric(gsub("seq","",baps$Isolate[which(baps[,as.numeric(strsplit(importance_RP_filtered$tag[j],"_")[[1]][2])]%>% .[[1]] %in% importance_RP_filtered$tag[j],)]))]
  taxa<-read_csv("/Users/moradigd/Documents/Plasmid/taxa_seq.csv", col_types = cols())
  #print(taxa[match(seq_RP, taxa$seq),6] %>% .[[1]] %>% unique(.))
  #print(sort(table(taxa[match(seq_RP, taxa$seq),6] %>% .[[1]] )))
  #print(sort(table(taxa[match(seq_RP, taxa$seq),5] %>% .[[1]] )))
  results_frequency<-rbind(results_frequency, taxa[match(seq_RP, taxa$seq),c(1,5)] )
  #print(sort(table(taxa[match(seq_RP, taxa$seq),4] %>% .[[1]] )))
}



taxa_order<-read_csv("/Users/moradigd/Documents/Plasmid/taxa_seq.csv", col_types = cols()) 
taxa_order<-taxa_order[match(df$seq , taxa_order$seq), ] %>% pull(Order)
table_2<-sort(table(pull(distinct(results_frequency),2)),decreasing = T)/sum(table(pull(distinct(results_frequency),2)))

#RP
table_2<-table_2[-9]

#PKJK
table_2<-table_2[1:14]
table_2<-table_2[-13]

#fileteretd<-c("Betaproteobacteriales", "Alteromonadales", "Pseudomonadales" ,"Caulobacterales")
#table_2<-table_2[match(fileteretd, names(table_2) )]

table_1<-table(taxa_order)[match(names(table_2) , names(table(taxa_order)))]/sum(table(taxa_order))


table_1<-data.frame(table_1)
table_1$label<-"Base"

table_2<-data.frame(table_2)
table_2$label<-"Target"
colnames(table_2)[1]<-"taxa_order"

table_res<-rbind(table_1,table_2)

ggplot(table_res, aes(x=taxa_order, y=Freq, fill=label )) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text( size=12, angle=45,hjust = 1),
        axis.text.y = element_text( size=12, hjust = 1),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("Order")+
  ylim(range(0,0.7))+
  ylab(expression("Frequency"))+
  scale_fill_manual(values=c("#009FBB","#FF00FF"))

#table_res$plasmid<-"RP4"
#write_csv(table_res, "/Users/moradigd/Documents/Plasmid/Figure_2_RP4.csv")

pkjk<-read_csv("/Users/moradigd/Documents/Plasmid/Figure_2_pKJK5.csv")
rp<-read_csv("/Users/moradigd/Documents/Plasmid/Figure_2_RP4.csv")
pb<-read_csv("/Users/moradigd/Documents/Plasmid/Figure_2_pB10.csv")
total<-rbind(pkjk, rp,pb )

ggplot(total, aes(fill=label, x=taxa_order ,y=Freq))+
  geom_bar( stat = "identity", position = "dodge")+
  facet_grid(~ plasmid, , scales="free_x") +
  theme_bw() +
  xlab("Order")+
  ylab("Frequency")+
  ylim(range(0,1))+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )+
  scale_fill_manual(values=c("#00FF00","#FF00FF"))


#BAPS iTOL
labels<-rep("#FFFFFF",  dim(baps)[1] )
labels[match(results_frequency$seq ,df$seq)]<-"#0000FF"

lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","PKJK"), paste0("COLOR",",","#0000FF",",","#FFFFFF"))
lines<-c(lines, "DATA",paste0(paste0("seq", seq(dim(baps)[1])) ,",", labels ) )
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/New_BAPS_PB_output.txt"))


#BAPS iTOL
for(j in seq( length(importance_RP_filtered))) {
  print(j)
  ids<-which(baps==importance_RP_filtered$tag[j], arr.ind=TRUE)[,1]
  ids_ex<-which(! seq(dim(baps)[1]) %in% ids)
  baps_tmp<-baps[,c(1,which(baps==importance_RP_filtered$tag[j], arr.ind=TRUE)[,2][1])]
  baps_tmp[ids,2]<-"#0000FF"
  baps_tmp[ids_ex,2]<-"#FFFFFF"
  
  lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",",importance_RP$tag[j]), paste0("COLOR",",","#0000FF",",","#FFFFFF"))
  lines<-c(lines, "DATA",paste0(pull(baps_tmp[,1],1) ,",", pull(baps_tmp[,2],1) ) )
  writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/BAPS_PKJK_output_",importance_RP$tag[j],".txt"))
}

with<-df[which(baps==importance_RP$tag[j], arr.ind=TRUE)[,1],2]
without<-df[which(baps!=importance_RP$tag[j], arr.ind=TRUE)[,1],2]
list_tmp<-cbind(c(with, without), rep(importance_RP$tag[j], length(c(with, without))), c(rep("with", length(with)), rep("without", length(without))))





ggplot(importance_RP, aes(x=mean_importance, y=mean_rank, size = count_n/5, col=count_n/10)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 24))+
  ylim(range(-5,60))+
  xlim(range(-0.01,0.25))+
  theme_bw()


#kmer analysis
library(tidyverse)
iterator<-seq(1,9)
kmer_size<-c(5,6,7,8)
plasmids<-c("PB","RP","PKJK")

importance<-c()
for(m in 1:length(plasmids)){
  for(k in 1:length(iterator)){
    for(j in 1:length(kmer_size)){
      tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/randomforests_importance",kmer_size[j],"_",iterator[k],"_",plasmids[m])) %>% arrange(desc(importance))
      tmp$plasmid<-rep(plasmids[m], dim(tmp)[1])
      tmp$kmer<-rep(as.character(kmer_size[j]), dim(tmp)[1])
      tmp$iterator<-rep(iterator[k], dim(tmp)[1])
      tmp$rank<-seq(1, dim(tmp)[1])
      importance<-rbind(importance, tmp)
    }
  }
}
importance_kmer_agg<- importance %>%
  group_by(tag, plasmid, kmer) %>%
  summarise(mean_importance=mean(importance), mean_rank=mean(rank), count_n=n()) %>%
  arrange(mean_rank) %>%
  filter(count_n>8) %>%
  filter(plasmid=="RP")

inclusion<-c()
for(i in 1:length(importance_kmer_agg$tag)){
  print(i)
  tmp="U"
  for(j in 1:length(importance_kmer_agg$tag)){
    if(i!=j){
      if(grepl(importance_kmer_agg$tag[i], importance_kmer_agg$tag[j], fixed=TRUE)){
        tmp="N"
        break
      }
    }
  }
  inclusion<-c(inclusion, tmp)
}
importance_kmer_agg<-importance_kmer_agg[which(inclusion=="U"),]
#write_csv(importance_kmer_agg, '/Users/moradigd/Documents/Plasmid/PB_informative_kmer.csv')
length(which(importance_kmer_agg$mean_rank<140))/length(importance_kmer_agg$mean_rank)

importance_kmer_agg<-importance_kmer_agg[importance_kmer_agg$mean_rank<150,]



hist(importance_kmer_agg$mean_rank)
library(tidyverse)
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 

ids<-which(!is.na(df$RP))
df<-df[!is.na(df$RP),]


results<-c()
significance<-c()
mean_diff<-c()
included<-c()
count<-c()

for(j in 1:dim(importance_kmer_agg)[1]){
  print(j)
  kmer<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/kmer_",importance_kmer_agg$kmer[j],".csv"), col_types = cols()) 
  kmer_tmp<-kmer[ids,match(importance_kmer_agg$tag[j], colnames(kmer))]
  
  with<-df[which(kmer_tmp!=0),3]
  without<-df[which(kmer_tmp==0),3]
  list_tmp<-cbind(c(with, without), rep(importance_kmer_agg$tag[j], length(c(with, without))), c(rep("with", length(with)), rep("without", length(without))))
  
  results<-rbind(results,list_tmp)
  significance<-c(significance, wilcox.test(with, without)$p.value)
  mean_diff<-c(mean_diff, mean(with)/mean(without))
  count<-c(count, length(which(kmer_tmp!=0)))
}  

importance_kmer_agg$mean_diff<-mean_diff
importance_kmer_agg$significance<-significance
#write_csv(importance_kmer_agg, file="/Users/moradigd/Documents/Plasmid/importance_kmer_agg_RP.csv")

#kmer analysis
kmer<-read_csv("/Users/moradigd/Documents/Plasmid/importance_kmer_agg_PKJK.csv")
kmer<-kmer[kmer$significance<0.01,]
presence_pattern<-c()
for(j in 1:dim(kmer)[1]){
  print(j)
  kmer_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/kmer_",kmer$kmer[j],".csv"), col_types = cols()) 
  kmer_tmp<-kmer_tmp[ids,match(kmer$tag[j], colnames(kmer_tmp))]
  presence_pattern<-rbind(presence_pattern, pull(kmer_tmp,1))
}
#write_csv(data.frame(presence_pattern),"/Users/moradigd/Documents/Plasmid/RP_important_kmers.csv")

kmer<-read_csv("/Users/moradigd/Documents/Plasmid/importance_kmer_agg_PKJK.csv")
kmer<-kmer[kmer$significance<0.01,]
presence_pattern<-read.csv("/Users/moradigd/Documents/Plasmid/PKJK_important_kmers.csv")
matrix_similarity<-matrix(0, nrow=dim(presence_pattern)[1], ncol=dim(presence_pattern)[1])

for(j in 1:dim(matrix_similarity)[1]){
  print(sum(presence_pattern[j,]))
  for(i in 1:dim(matrix_similarity)[1]){
    matrix_similarity[i,j]<-cor(as.numeric(presence_pattern[i,]),as.numeric(presence_pattern[j,]))
    #matrix_similarity[i,j]<-length(intersect(which(presence_pattern[i,]==1),which(presence_pattern[j,]==1)))/length(union(which(presence_pattern[i,]==1),which(presence_pattern[j,]==1)))
  }
}
library(corrplot)
colnames(matrix_similarity)<-kmer$tag
row.names(matrix_similarity)<-kmer$tag
corrplot(matrix_similarity, method = "circle", type="upper",tl.col="black", tl.cex = 1)


hist(matrix_similarity[lower.tri(matrix_similarity)], col = "black", breaks = 200, las=2, xlab = "Pairwise Pearson Correlation (kmer presence patterns)", main="")
box()

corr_ind<-which(matrix_similarity>0.90, arr.ind = T)
holder<-c()
for(j in 1:dim(corr_ind)[1]){
  if(corr_ind[j,1]!=corr_ind[j,2]){
    holder<-c(holder, j)
  }
}
holder

#Boxplot
box_df<-c()
for(j in 1:length(kmer$tag)){
  print(j)
  kmer_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/kmer_",kmer$kmer[j],"_PB",".csv"), col_types = cols()) 
  kmer_tmp<-kmer_tmp[,match(kmer$tag[j], colnames(kmer_tmp))]
  box_df<-rbind(box_df, cbind(rep(kmer$tag[j], dim(kmer_tmp)[1]) ,pull(kmer_tmp,1), df$PKJK))
}
box_df<-data.frame(box_df)
box_df$X3<-as.numeric(as.character(box_df$X3))
box_df$X2

#write_csv(box_df,"/Users/moradigd/Documents/Plasmid/PB_boxplot_kmers.csv")
box_df<-read_csv("/Users/moradigd/Documents/Plasmid/PKJK_boxplot_kmers.csv")
box_df$X1 <- factor(box_df$X1)
#levels(box_df$X1)<-kmer$tag[order(kmer$mean_diff)]
box_df$X1 <- factor(box_df$X1 , levels=kmer$tag[order(kmer$mean_diff)])
#box_df$X1 <- factor(box_df$X1 , levels=kmer$tag)
box_df$X2<-as.character(box_df$X2)

kmer$tag[order(kmer$mean_diff)]
kmer[kmer$tag=="CCTTTGAG",]

ggplot(box_df, aes(fill=X2, y=X3, x=X1)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent", outlier.alpha = 0)+
  theme_bw()+
  theme(axis.text.x = element_text( size=11, angle=90,hjust = 1),
        axis.text.y = element_text( size=15, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=18, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("kmers")+
  ylim(range(0,75))+
  ylab(expression("Permissiveness"))



ggplot(box_df, aes(fill=X2, y=X3, x=X1)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent")+
  theme_bw()+
  theme(
    #axis.text.x = element_text( size=13, angle=45,hjust = 1),
    axis.text.y = element_text( size=13, hjust = 1),
    #axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank()
  )+
  ylim(range(0,75))+
  ylab(expression("Permissiveness"))



#taxa kmer 
library(corrplot)

taxa<-read_csv("/Users/moradigd/Documents/Plasmid/taxa_seq.csv", col_types = cols())
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 

df<-df[!is.na(df$PKJK),]

ids<-colnames(pan_genome)[which(pan_genome[match(significant$Gene, pan_genome$Gene)[2],]!=0)]
ids<-as.numeric(gsub("seq","",ids))

for(j in seq(dim(kmer)[1])){
  kmer_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/kmer_",kmer$kmer[j],"_PKJK",".csv"), col_types = cols()) 
  kmer_tmp_alleles<-kmer_tmp[,match(kmer$tag[j], colnames(kmer_tmp))]  %>% .[[1]]
  print(table(taxa$Order[match(df$seq[which(kmer_tmp_alleles!=0)],taxa$seq)]))
}

for(j in seq(dim(kmer)[1])){
  print(j)
  kmer_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/kmer_",kmer$kmer[j],"_PKJK",".csv"), col_types = cols()) 
  kmer_tmp_alleles<-kmer_tmp[,match(kmer$tag[j], colnames(kmer_tmp))]  %>% .[[1]]
  #itol_mat<-t(pan_genome[match(significant$Gene, pan_genome$Gene),15:dim(pan_genome)[2]])
  lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",",kmer$tag[j]), paste0("COLOR",",","#0000FF",",","#FFFFFF"))
  lines<-c(lines, "DATA",paste0(paste0("seq",seq(1, length(kmer_tmp_alleles))) ,",", ifelse(kmer_tmp_alleles==1,"#0000FF", "#FFFFFF")) )
  writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/output_kmer_",kmer$tag[j],".txt"))
}

#Taxa analysis 
library(tidyverse)
kmer<-read_csv("/Users/moradigd/Documents/Plasmid/importance_kmer_agg_RP.csv")
kmer<-kmer[kmer$significance<0.01,]
frequency_orders<-c()
significance<-c()

output_spieces<-c()
for(j in seq(length(kmer$tag))){
  print(kmer$tag[j])
  taxa<-read_csv("/Users/moradigd/Documents/Plasmid/TaxaFamily.csv", col_types = cols())
  
  kmer_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/kmer_",kmer$kmer[j],"_PB",".csv"), col_types = cols()) 
  kmer_tmp_alleles<-kmer_tmp[,match(kmer$tag[j], colnames(kmer_tmp))]  %>% .[[1]]
  
  ids<-which(kmer_tmp_alleles!=0)
  shared<-intersect(df$seq[ids],taxa$seq)
  
  hld_table_or<-sort(table(taxa[match(shared,taxa$seq),7]),decreasing = TRUE)
  filetered_frequent<-which(hld_table_or>0)
  hld_table<-hld_table_or/sum(hld_table_or)
  
  hld_table_tmp<- data.frame(hld_table)
  
  hld_table<-hld_table[filetered_frequent]
  
  shared<-intersect(df$seq,taxa$seq)
  hld_table_tot<-sort(table(taxa[match(shared,taxa$seq),7]),decreasing = TRUE)
  hld_table_tot<-hld_table_tot/sum(hld_table_tot)
  hld_table_freq<-hld_table/hld_table_tot[match( names(hld_table) ,names(hld_table_tot) )]
  significance<-c(significance, hld_table_freq)
  
  hld_table_tmp$label<-kmer$tag[j]
  if(dim(hld_table_tmp)[1]==1){
    output_spieces<-rbind(output_spieces,c(row.names(hld_table_tmp),as.numeric(hld_table_tmp[1]), as.character(hld_table_tmp[2]) ))
  }else{output_spieces<-rbind(output_spieces,hld_table_tmp)}
  
  frequency_orders<-c(frequency_orders,hld_table_freq)
}
output_spieces$sig<-significance
output_spieces$species<-"PB"
#write_csv(output_spieces,'/Users/moradigd/Documents/Plasmid/New_PB_kmer_Taxa.csv' )


#PKJK_plasmid<-read_csv('/Users/moradigd/Documents/Plasmid/PKJK_GWAS_Taxa.csv')
#PB_plasmid<-read_csv('/Users/moradigd/Documents/Plasmid/PB_GWAS_Taxa.csv')
RP_plasmid<-read_csv('/Users/moradigd/Documents/Plasmid/New_RP_kmer_Taxa.csv')

tot_plasmid<-rbind(RP_plasmid)
names(sort(table(tot_plasmid$Var1)))
names_freq<-names(sort(table(RP_plasmid$Var1), decreasing = T))[1:10]



tot_plasmid$new_lab<-names_freq[match(tot_plasmid$Var1, names_freq)]
tot_plasmid$new_lab<-ifelse( is.na(tot_plasmid$new_lab),"Others", tot_plasmid$new_lab   )
tot_plasmid<-tot_plasmid[tot_plasmid$sig>1,]
tot_plasmid$label <- factor(tot_plasmid$label , levels=kmer$tag[order(kmer$mean_diff)])
levels(tot_plasmid$label)<-kmer$tag
colours<-c("#D1C3DA","#8084D4","#C7E4D8","#D6D755","#6DE6D3","#DF5070","#D767C9","#91D794","#D78A47","#D89593","#72B5D3","#D99DDD","#D546E1","#7A46DE","#E2D4A0","#7CE45B","#0000FF","#000000")
orders<-c("Pseudomonadales","Flavobacteriales","Betaproteobacteriales","Aeromonadales","Rhizobiales","Xanthomonadales","Enterobacteriales","Alteromonadales","Sphingomonadales","Sphingobacteriales","Micrococcales","Cytophagales","Caulobacterales","Rhodobacterales","Clostridiales","Corynebacteriales","Others","unclassified")

tot_plasmid$label[tot_plasmid$label=="Enterobacteriales"]<-"Enterobacterales"

ggplot(tot_plasmid, aes(fill=new_lab , x=factor(label) ,y=Freq))+
  geom_bar( stat = "identity")+
  theme_bw() +
  xlab("kmers")+
  ylab("Frequency")+
  ylim(range(0,1))+
  theme(axis.text.x = element_text( size=13, angle=90,hjust = 1),
        axis.text.y = element_text( size=15, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=18),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )+
  scale_fill_manual(values=colours[match(sort(unique(tot_plasmid$new_lab)),orders)])


ggplot(tot_plasmid, aes(fill=new_lab , x=factor(label) ,y=Freq))+
  geom_bar( stat = "identity")+
  theme_bw() +
  xlab("kmers")+
  ylab("")+
  ylim(range(0,1))+
  theme(
    #axis.text.x = element_text( size=13, angle=45,hjust = 1),
    axis.text.y = element_text( size=13, hjust = 1),
    #axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank()
  )+
  ylab("Frequency")



library(randomcoloR)


taxa<-read_csv("/Users/moradigd/Documents/Plasmid/taxa_seq.csv", col_types = cols())
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 

df<-df[!is.na(df$PKJK),]


palette <- distinctColorPalette(length(kmer$tag))
itol<-matrix("#FFFFFF", nrow = dim(df)[1], ncol = 2)
itol[,1]<-paste0("seq",seq(1,dim(df)[1])) 
for(j in seq( length(kmer$tag))){
  print(j)
  if(kmer$mean_diff[j]>1){
    kmer_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/kmer_",kmer$kmer[j],"_PKJK",".csv"), col_types = cols()) 
    kmer_tmp<-kmer_tmp[,match(kmer$tag[j], colnames(kmer_tmp))]
    ids<-which(kmer_tmp[,1]!=0)
    itol[ids,2]<-palette[j]
  }
}
itol<-data.frame(itol)
palette<-c("#FFFFFF",palette)
lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","GAWS"), paste0(c("COLOR",paste0(",",palette)),collapse = "" ) )
lines<-c(lines, "DATA", paste0(itol[,1] ,",", itol[,2] ))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/Pos_Predictive_kmer_PKJK_output.txt"))










#Correlation test 
library(tidyverse)
taxa<-read_csv("/Users/moradigd/Documents/Plasmid/TaxaFamily.csv")

corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 

shared_isolates<-intersect(taxa$seq, df$seq)
shared_isolates_ids<-match(shared_isolates, df$seq)
df<-df[match(shared_isolates, df$seq),]
#write_csv(df, "/Users/moradigd/Documents/Plasmid/input_agg.csv")

#df<-df[!is.na(df$PB),]
PKJK_label<-ntile(df$PKJK,4)
PKJK_label<-ifelse(is.na(PKJK_label),5,PKJK_label)
RP_label<-ntile(df$RP,4)
RP_label<-ifelse(is.na(RP_label),5,RP_label)
PB_label<-ntile(df$PB,4)
PB_label<-ifelse(is.na(PB_label),5,PB_label)

colors<-c("#FF6A6A", "#EE6363", "#CD5555", "#8B3A3A","#FFFFFF")

lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","RP"), paste0("COLOR",",","#FF6A6A,", "#EE6363,", "#CD5555,", "#8B3A3A,","#FFFFFF"))
lines<-c(lines, "DATA",paste0(paste0("seq",seq(1, dim(df)[1])) ,",", colors[match(RP_label, seq(1,5))]) )
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/Total_RP_output_kmer.txt"))


library(rncl)
library(phangorn)
library(ape)


kmer<-read_csv("/Users/moradigd/Documents/Plasmid/Total_kmer_8_PKJK.csv")
kmer<-kmer[,-1]
distance_kmer<-matrix(0, nrow = dim(kmer)[1], ncol = dim(kmer)[1])
for(i in 1:dim(kmer)[1]){
  print(i)
  for(j in i:dim(kmer)[1]){
    vector_1<-kmer[i,] %>% pull(.)
    vector_2<-kmer[j,] %>% pull(.)
    distance_kmer[i,j]<-dist(rbind(vector_1, vector_2))
  }
}
colnames(distance_kmer)<-paste0("seq",seq(1,dim(kmer)[1])) 
row.names(distance_kmer)<-paste0("seq",seq(1,dim(kmer)[1])) 
#write_csv(data.frame(distance_kmer), "/Users/moradigd/Documents/Plasmid/Distance.csv")
#distance_kmer<-data.frame(distance_kmer)
distance_kmer<-read_csv("/Users/moradigd/Documents/Plasmid/Distance.csv")
distance_kmer<-distance_kmer[shared_isolates_ids,shared_isolates_ids]

colnames(distance_kmer)<-paste0("seq",seq(1,dim(distance_kmer)[1])) 
row.names(distance_kmer)<-paste0("seq",seq(1,dim(distance_kmer)[1])) 
distance_kmer<-as.matrix(distance_kmer)
distance_kmer[lower.tri(distance_kmer)]<-distance_kmer[upper.tri(distance_kmer)]
dist_distance_kmer<-dist(distance_kmer)

tre<-nj(dist_distance_kmer)
tre<-midpoint(tre)
write.tree(tre, file="/Users/moradigd/Documents/Plasmid/iTOL/kmerfive.tre")


#iTOL spieces 
library(tidyverse)
taxa<-read_csv("/Users/moradigd/Documents/Plasmid/TaxaFamily.csv")
isolates<-read_csv("/Users/moradigd/Documents/Plasmid/input_agg.csv")
taxa<-taxa[match(shared_isolates, taxa$seq), ]

names(sort(-table(taxa$Order)))[which(sort(-table(taxa$Order))/sum(sort(-table(taxa$Order)))>0.01)]
majot_orders<- c("Pseudomonadales"  ,     "Flavobacteriales"    ,  "Betaproteobacteriales", "Aeromonadales"      ,  
                 "Rhizobiales"       ,    "Xanthomonadales"     ,  "Enterobacteriales"    , "Alteromonadales"   ,   
                 "Sphingomonadales"   ,   "Sphingobacteriales" ,   "Micrococcales"        ,
                 "Cytophagales"        ,  "Caulobacterales"   ,    "Rhodobacterales"   ,    "Clostridiales"        ,
                 "Corynebacteriales" )

library(randomcoloR)
palette <- distinctColorPalette(length(majot_orders))

colours_assignment<-palette[match(taxa$Order,majot_orders)]
colours_assignment[is.na(colours_assignment)]<-"#FFFFFF"

lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","Order"), paste0("COLOR",paste0(paste0(",",c(palette,"#FFFFFF")), collapse = "")))
lines<-c(lines, "DATA",paste0(  paste0("seq",seq(1,dim(taxa)[1])) ,",",  colours_assignment))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/Total_Order_permissveness_output_p.txt"))


library(randomcoloR)
names(sort(-table(taxa$Class)))[which(sort(-table(taxa$Class))/sum(sort(-table(taxa$Class)))>0.01)]
major_classes<-c( "Gammaproteobacteria", "Bacteroidia"       ,  "Alphaproteobacteria" ,"Actinobacteria" ,     
                  "Bacilli"           ,  "Clostridia")
palette <- distinctColorPalette(length(major_classes))

colours_assignment<-palette[match(taxa$Class,major_classes)]
colours_assignment[is.na(colours_assignment)]<-"#FFFFFF"

lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","Class"), paste0("COLOR",paste0(paste0(",",c(palette,"#FFFFFF")), collapse = "")))
lines<-c(lines, "DATA",paste0(paste0("seq",seq(1,dim(taxa)[1])) ,",",  colours_assignment))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/Total_Class_permissveness_output_p.txt"))



library(randomcoloR)
names(sort(-table(taxa$Phylum )))[which(sort(-table(taxa$Phylum))/sum(sort(-table(taxa$Phylum)))>0.01)]
major_classes<-c( "Proteobacteria", "Bacteroidetes" , "Actinobacteria" ,"Firmicutes")
palette <- distinctColorPalette(length(major_classes))

colours_assignment<-palette[match(taxa$Phylum,major_classes)]
colours_assignment[is.na(colours_assignment)]<-"#FFFFFF"

lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","Phylum"), paste0("COLOR",paste0(paste0(",",c(palette,"#FFFFFF")), collapse = "")))
lines<-c(lines, "DATA",paste0(paste0("seq",seq(1,dim(taxa)[1])) ,",",  colours_assignment))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/Total_Phylum_permissveness_output_p.txt"))





#GWAS 
library(tidyverse)
significant<-read_csv("/Users/moradigd/Documents/Plasmid/PKJK_GWAS.csv") %>% filter(Best_pairwise_comp_p <0.05 & Worst_pairwise_comp_p<0.05 & Bonferroni_p<0.05)
pan_genome<-read_csv("/Users/moradigd/Documents/Plasmid/pan_genome_agg_PKJK.csv")
significant$Odds_ratio

taxa<-read_csv("/Users/moradigd/Documents/Plasmid/taxa_seq.csv", col_types = cols())
corr<-read_csv("/Users/moradigd/Documents/Plasmid/input.csv")

PKJK<-apply(corr[, which(grepl("_PKJK",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
RP<-apply(corr[, which(grepl("_RP",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
PB<-apply(corr[, which(grepl("_PB",colnames(corr)))], 1, function(x) mean(x,na.rm=TRUE))
df<-data.frame(list("seq"=corr$seq,"PKJK"=PKJK, "RP"=RP, "PB"=PB)) 

df<-df[!is.na(df$PKJK),]

ids<-colnames(pan_genome)[which(pan_genome[match(significant$Gene, pan_genome$Gene)[2],]!=0)]
ids<-as.numeric(gsub("seq","",ids))

#Correlation between significant features
library(corrplot)


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  for(i in seq(n)){
    for(j in seq(n)){
      tmp <- cor.test(as.numeric(as.character(gsub(" ","",mat[, j]))) , as.numeric(as.character(gsub(" ","",mat[, i]))) ,use="pairwise.complete.obs",method ="spearman")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  diag(p.mat) <- 0
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

cor_mat<-cor(t(pan_genome[match(significant$Gene, pan_genome$Gene),15:dim(pan_genome)[2]]))
corrplot(cor_mat, method = "circle", type="upper",p.mat = cor.mtest(cor_mat), sig.level = 0.01,tl.col="black", tl.cex = 1)


#itol file maker 
#DATASET_COLORSTRIP
#SEPARATOR TAB
#DATASET_LABEL	r_CIP.0.004
#COLOR	#FFBFBF	#FF7F7F	#FF3F3F	#FF0000	#FFFFFF
#DATA

kmers<-significant$Gene
i<-3
itol_mat<-t(pan_genome[match(significant$Gene, pan_genome$Gene),15:dim(pan_genome)[2]])
lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",",kmers[i]), paste0("COLOR",",","#0000FF",",","#FFFFFF"))
lines<-c(lines, "DATA",paste0(row.names(itol_mat),",", ifelse(itol_mat[,i]==1,"#0000FF", "#FFFFFF")) )
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/output_",kmers[i],".txt"))



itol_mat<-t(pan_genome[match(significant$Gene, pan_genome$Gene),15:dim(pan_genome)[2]])
lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","r_CIP.0.004"), paste0("COLOR",",","#FF0000",",","#FFFFFF"))
lines<-c(lines, "DATA",paste0(row.names(itol_mat),",",  ifelse(df$PKJK> median(df$PKJK),"#FF0000", "#FFFFFF")))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/PB_permissveness_output.txt"))

taxa<-read_csv("/Users/moradigd/Documents/Plasmid/TaxaFamily.csv")

names(sort(-table(taxa$Order)))[which(sort(-table(taxa$Order))/sum(sort(-table(taxa$Order)))>0.01)]
majot_orders<- c("Pseudomonadales"  ,     "Flavobacteriales"    ,  "Betaproteobacteriales", "Aeromonadales"      ,  
                 "Rhizobiales"       ,    "Xanthomonadales"     ,  "Enterobacteriales"    , "Alteromonadales"   ,   
                 "Sphingomonadales"   ,   "Sphingobacteriales" ,   "Micrococcales"        ,
                 "Cytophagales"        ,  "Caulobacterales"   ,    "Rhodobacterales"   ,    "Clostridiales"        ,
                 "Corynebacteriales" )

library(randomcoloR)
#palette <- distinctColorPalette(length(majot_orders))
palette <-c("#DBB1D1", "#CEE25C" ,"#7AC5DD", "#E1B95D", "#7763DA", "#D85A86", "#77E1C9", "#838F7A",
            "#CFE1A9" ,"#7C93D5", "#7DDE8A", "#DF8766" ,"#D583D8", "#C942D8", "#80EC4D", "#D6DCD8",
            "#FFFFFF"
)
#[1] "#DBB1D1" "#CEE25C" "#7AC5DD" "#E1B95D" "#7763DA" "#D85A86" "#77E1C9" "#838F7A"
#[9] "#CFE1A9" "#7C93D5" "#7DDE8A" "#DF8766" "#D583D8" "#C942D8" "#80EC4D" "#D6DCD8"

colours_assignment<-palette[match(taxa$Order[match(df$seq, taxa$seq)],majot_orders)]
colours_assignment[is.na(colours_assignment)]<-"#FFFFFF"

itol_mat<-t(pan_genome[match(significant$Gene, pan_genome$Gene),15:dim(pan_genome)[2]])
lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","Order"), paste0("COLOR",paste0(paste0(",",c(palette,"#FFFFFF")), collapse = "")))
lines<-c(lines, "DATA",paste0(row.names(itol_mat),",",  colours_assignment))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/RP_Order_permissveness_output_p.txt"))

library(randomcoloR)
names(sort(-table(taxa$Class)))[which(sort(-table(taxa$Class))/sum(sort(-table(taxa$Class)))>0.01)]
major_classes<-c( "Gammaproteobacteria", "Bacteroidia"       ,  "Alphaproteobacteria" ,"Actinobacteria" ,     
                  "Bacilli"           ,  "Clostridia")
palette <- distinctColorPalette(length(major_classes))

colours_assignment<-palette[match(taxa$Class[match(df$seq, taxa$seq)],major_classes)]
colours_assignment[is.na(colours_assignment)]<-"#FFFFFF"

itol_mat<-t(pan_genome[match(significant$Gene, pan_genome$Gene),15:dim(pan_genome)[2]])
lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","Order"), paste0("COLOR",paste0(paste0(",",c(palette,"#FFFFFF")), collapse = "")))
lines<-c(lines, "DATA",paste0(row.names(itol_mat),",",  colours_assignment))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/PB_Class_permissveness_output_p.txt"))




#Box plot 
results<-c()
significance<-c()
mean_diff<-c()

for(j in seq(length(significant$Gene))){
  ids<-colnames(pan_genome)[which(pan_genome[match(significant$Gene, pan_genome$Gene)[j],]!=0)]
  ids<-as.numeric(gsub("seq","",ids))[-1]
  with<-df[which(seq(dim(df)[1]) %in% ids),2]
  without<-df[which(!seq(dim(df)[1]) %in% ids),2]
  list_tmp<-cbind(c(with, without), rep(significant$Gene[j], length(c(with, without))), c(rep("with", length(with)), rep("without", length(without))))
  results<-rbind(results,list_tmp)
  significance<-c(significance, wilcox.test(with, without)$p.value)
  mean_diff<-c(mean_diff, mean(with)/mean(without))
}
significant$mean_diff<-mean_diff
significant$significance<-significance
significant_filtered<- significant[significant$significance<0.05 & significant$mean_diff>1 ,]


holder_uniqe<-c()
for(i in unique(significant_filtered$Gene)){
  for(j in unique(significant_filtered$Gene)){
    if(i !=j & grepl( i, j, fixed = TRUE)){
      holder_uniqe<-c(holder_uniqe, i)
      break
    }
  }
}
significant_filtered<-significant_filtered[which(! significant_filtered$Gene %in% holder_uniqe),]

significant<-significant[match(significant_filtered$Gene,significant$Gene),]

significant_filtered_gene<- significant$Gene[match(unique(significant$significance),significant$significance )]
significant<-significant[match(significant_filtered_gene,significant$Gene),]


significant_filtered<-significant[match(significant_filtered_gene,significant$Gene),]


#corrplot(, method = "circle", type="upper")
frequency_orders<-c()
significance<-c()


output_spieces<-c()
for(j in seq(length(significant$Gene))){
  print(significant$Gene[j])
  taxa<-read_csv("/Users/moradigd/Documents/Plasmid/TaxaFamily.csv", col_types = cols())
  shared<-intersect(df$seq[ids],taxa$seq)
  ids<-colnames(pan_genome)[which(pan_genome[match(significant$Gene[j], pan_genome$Gene),]!=0)]
  ids<-as.numeric(gsub("seq","",ids))[-1]
  hld_table_or<-sort(table(taxa[match(shared,taxa$seq),7]),decreasing = TRUE)
  filetered_frequent<-which(hld_table_or>0)
  hld_table<-hld_table_or/sum(hld_table_or)
  
  hld_table_tmp<- data.frame(hld_table)
  
  hld_table<-hld_table[filetered_frequent]
  
  shared<-intersect(df$seq,taxa$seq)
  hld_table_tot<-sort(table(taxa[match(shared,taxa$seq),7]),decreasing = TRUE)
  hld_table_tot<-hld_table_tot/sum(hld_table_tot)
  hld_table_freq<-hld_table/hld_table_tot[match( names(hld_table) ,names(hld_table_tot) )]
  significance<-c(significance, hld_table_freq)
  #hld_table_freq<-hld_table_freq[which(hld_table_freq>1)]
  #print(hld_table_or[match( names(hld_table_freq), names(hld_table_or))])
  #print(hld_table_freq)
  
  hld_table_tmp$label<-significant$Gene[j]
  if(dim(hld_table_tmp)[1]==1){
    output_spieces<-rbind(output_spieces,c(row.names(hld_table_tmp),as.numeric(hld_table_tmp[1]), as.character(hld_table_tmp[2]) ))
  }else{output_spieces<-rbind(output_spieces,hld_table_tmp)}
  
  
  frequency_orders<-c(frequency_orders,hld_table_freq)
  
  
  
  
  print("---------------")
}
output_spieces$sig<-significance
output_spieces$species<-"PKJK"
write_csv(output_spieces,'/Users/moradigd/Documents/Plasmid/PKJK_GWAS_Taxa.csv' )
#PB
#[1] "Flavobacteriales"   "Bacillales"         "Cytophagales"       "Sphingobacteriales"
#[5] "Chitinophagales"    "Xanthomonadales" 



#PB

significant_filtered$Number_pos_present_in
significant_filtered$Number_neg_present_in
significant_filtered$Number_pos_not_present_in
significant_filtered$Number_neg_not_present_in

results<-data.frame(results)
results$X1<-as.numeric(as.character(results$X1))
results<-results[which(results$X2 %in% significant$Gene),]
write_csv(results, "/Users/moradigd/Documents/Plasmid/PKJK_boxplot.csv")
ggplot(results, aes(fill=X3, y=X1, x=X2)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent")+
  theme_bw()+
  theme(axis.text.x = element_text( size=13, angle=90,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))+
  xlab("kmer")+
  ylim(range(0,50))+
  ylab(expression("Permissiveness"))


library(randomcoloR)
palette <- distinctColorPalette(length(significant$Gene))
palette<-c("#FFFFFF",palette)
itol<-matrix("#FFFFFF", nrow = dim(df)[1], ncol = 2)
itol[,1]<-paste0("seq",seq(1,dim(df)[1])) 
for(j in seq( length(significant_filtered$Gene))) {
  print(j)
  ids<-colnames(pan_genome)[which(pan_genome[match(significant$Gene, pan_genome$Gene)[j],]!=0)]
  ids<-as.numeric(gsub("seq","",ids))[-1]
  ids_ex<-which(! seq(dim(df)[1]) %in% ids)
  
  #baps_tmp<-baps[,c(1,which(baps==significant_filtered$tag[j], arr.ind=TRUE)[,2][1])]
  
  itol[ids,2]<-palette[j]
}
itol<-data.frame(itol)
lines<-c("DATASET_COLORSTRIP", "SEPARATOR COMMA",paste0("DATASET_LABEL",",","GAWS"), paste0(c("COLOR",paste0(",",palette)),collapse = "" ) )
lines<-c(lines, "DATA", paste0(itol[,1] ,",", itol[,2] ))
writeLines(lines,paste0("/Users/moradigd/Documents/Plasmid/iTOL/GAWS_kmer_RP_output.txt"))



#plot frequency and boxplot 
PKJK_plasmid<-read_csv('/Users/moradigd/Documents/Plasmid/PKJK_GWAS_Taxa.csv')
PB_plasmid<-read_csv('/Users/moradigd/Documents/Plasmid/PB_GWAS_Taxa.csv')
RP_plasmid<-read_csv('/Users/moradigd/Documents/Plasmid/RP_GWAS_Taxa.csv')

tot_plasmid<-rbind(PKJK_plasmid,PB_plasmid,RP_plasmid)
names(sort(table(tot_plasmid$Var1)))
names_freq<-c("Corynebacteriales"    , "Micrococcales"        ,
              "Sphingomonadales"   ,   "Cellvibrionales"     ,  "Aeromonadales" ,       
              "Betaproteobacteriales" ,"Enterobacteriales"  ,   "Xanthomonadales",      
              "Pseudomonadales") 

tot_plasmid$species[tot_plasmid$species=="PB"]<-"pB10"
tot_plasmid$species[tot_plasmid$species=="PKJK"]<-"pKJK5"
tot_plasmid$species[tot_plasmid$species=="RP4"]<-"RP4"


colours<-c("#D1C3DA","#8084D4","#C7E4D8","#D6D755","#6DE6D3","#DF5070","#D767C9","#91D794","#D78A47","#D89593","#72B5D3","#D99DDD","#D546E1","#7A46DE","#E2D4A0","#7CE45B","#0000FF","#000000")
orders<-c("Pseudomonadales","Flavobacteriales","Betaproteobacteriales","Aeromonadales","Rhizobiales","Xanthomonadales","Enterobacteriales","Alteromonadales","Sphingomonadales","Sphingobacteriales","Micrococcales","Cytophagales","Caulobacterales","Rhodobacterales","Clostridiales","Corynebacteriales","Others","unclassified")


tot_plasmid$new_lab<-names_freq[match(tot_plasmid$Var1, names_freq)]
tot_plasmid$new_lab<-ifelse( is.na(tot_plasmid$new_lab),"Others", tot_plasmid$new_lab   )
tot_plasmid<-tot_plasmid[tot_plasmid$sig>1,]
tot_plasmid$label[tot_plasmid$label=="Enterobacteriales"]<-"Enterobacterales"

ggplot(tot_plasmid, aes(fill=new_lab, x=factor(label) ,y=Freq))+
  geom_bar( stat = "identity")+
  facet_grid(~ species, , scales="free_x") +
  theme_bw() +
  xlab("kmers")+
  ylab("Frequency")+
  ylim(range(0,1))+
  theme(axis.text.x = element_text( size=13, angle=90,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )+
  scale_fill_manual(values=colours[match(sort(unique(tot_plasmid$new_lab)),orders)])

tot_plasmid
ggplot(tot_plasmid, aes( x=factor(label) ,y=sig))+
  geom_bar( stat = "identity")+
  facet_grid(~ species, , scales="free_x") +
  theme_bw() +
  xlab("kmers")+
  ylab("-log10(p-value)")+
  ylim(range(0,100))+
  theme(axis.text.x = element_text( size=13, angle=90,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )


input<-data.frame(list(label=c(rep("with", length(with)),rep("without", length(without))),value=c(df$PKJK[with],df$PKJK[without])))
boxplot(input$value ~ input$label, ylim=range(0,60))


#tree visualization
input<-data.frame(list(Name=c(paste0("seq",with),paste0("seq",without)), label=c(rep("with", length(with)),rep("without", length(without)))))
write_tsv(input, "/Users/moradigd/Documents/Plasmid/input_anno_PKJK.tsv")

#simulation visulization of data
library(gmodels)
library(tidyverse)
files<-dir("/Users/moradigd/Documents/Plasmid", pattern = "^results_simulation_")

files_total<-c()
for(i in 1:length(files)){
  file_tmp<-read_csv(paste0("/Users/moradigd/Documents/Plasmid/", files[i]),col_types = cols())
  file_tmp$alpha<-as.character(strsplit(files[i],"_")[[1]][3])
  file_tmp$mutation<-as.character(strsplit(files[i],"_")[[1]][7])
  file_tmp$size<-as.character(gsub(".fasta.csv","", strsplit(files[i],"_")[[1]][8]))
  files_total<-rbind(files_total, file_tmp)
  print(file_tmp)
}
colnames(files_total)<-c("index","cor","pvalue","sampleset","count","alpha","mutation","size")

files_total_sum<-files_total %>%
  #dplyr::filter(size=="2000") %>%
  dplyr::filter(sampleset=="train") %>%
  group_by(alpha,mutation, size,sampleset) %>%
  
  summarise(mean_cor=mean(cor), lowCI = ci(cor)[2],hiCI = ci(cor)[3], sd_correlation= ci(cor)[4])

files_total_sum$size <- factor(files_total_sum$size, levels=c("100", "500", "1000","2000"))


ggplot(files_total_sum, aes(fill=alpha, x=size ,y=mean_cor))+
  geom_bar( stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin= files_total_sum$lowCI, ymax=files_total_sum$hiCI), width=.2,
                position=position_dodge(.9))+
  facet_grid(~ mutation) +
  theme_bw() +
  ylab("")+
  ylim(range(0,1))+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )



files_total_sum<-files_total %>%
  dplyr::filter(alpha=="600") %>%
  dplyr::filter(sampleset=="train") %>%
  group_by(alpha,mutation, size,sampleset) %>%
  summarise(mean_cor=mean(cor), lowCI = ci(cor)[2],hiCI = ci(cor)[3], sd_correlation= ci(cor)[4])



ggplot(files_total_sum, aes(x=mutation, y=mean_cor,group=size,  color=size)) + 
  geom_line( size = 2, alpha = 0.75) +
  geom_point()+
  ylim(range(0,1))+
  theme_bw()+
  ylab("Correlation")+
  geom_errorbar(aes(ymin=files_total_sum$lowCI, ymax=files_total_sum$hiCI), width=.01,
                position=position_dodge(0.05), size = 0.35)+
  theme(axis.text.x = element_text( size=20, angle=90,hjust = 1),
        axis.text.y = element_text( size=20, hjust = 1),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        ))


#Tree for simulated sequences 
library(rncl)
library(phangorn)
library(ape)

dna<-read.dna("/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/seq_RP_aligned.fasta",format= "fasta")
distdna<-dist.dna(dna,model ="N",pairwise.deletion = TRUE,as.matrix = TRUE)
write.csv(distdna, file="/Users/moradigd/Documents/Applications/FastML.v3.11/programs/fastml/PKJK_dist-Core.csv",quote=FALSE, row.names = FALSE)
tre<-nj(distdna)
tre<-midpoint(tre)
write.tree(tre, file="/Users/moradigd/Documents/Plasmid/RP_sequence_alignment.fasta.tre")


#Comparison between BAPS and kmer informative groups 
baps_kmer_comp<-read_csv('/Users/moradigd/Documents/Plasmid/BAPS_kmer_comp.csv') %>%
  mutate(new_label=ifelse(Code=="0","0","1")) 

baps_kmer_comp$Plasmid[baps_kmer_comp$Plasmid=="PB"]<- "pB10"
baps_kmer_comp$Plasmid[baps_kmer_comp$Plasmid=="PKJK"]<- "pKJK5"
baps_kmer_comp$Plasmid[baps_kmer_comp$Plasmid=="RP"]<- "RP4"

factor(baps_kmer_comp$Code)

ggplot(baps_kmer_comp, aes(fill=factor(new_label), x=Model ,y=Permissivenss))+
  geom_boxplot(position="dodge", alpha=0.5, outlier.colour="transparent")+
  theme_bw() +
  facet_grid(~ Plasmid)+
  ylab("Permissiveness")+
  ylim(range(0,50))+
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )+
  scale_fill_manual(values=c("#0000FF", "#FF0000"))

#Tradis results 
library(tidyverse)
files<-list.files("/Users/moradigd/Documents/Tradis/new/surA", pattern = ".csv.essen.csv")

output<-matrix(0, nrow = length(files) , ncol = length(files)   )
for(i in 1:length(files)){
  for(j in 1:length(files)){
    file_1<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/surA/",files[i]),show_col_types = FALSE) 
    file_2<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/surA/",files[j]),show_col_types = FALSE) 
    output[i,j]<-length(intersect(file_1$locus_tag, file_2$locus_tag))/length(union(file_1$locus_tag, file_2$locus_tag))
  }
}


library(tidyverse)
files_ess<-list.files("/Users/moradigd/Documents/Tradis/new/bamE", pattern = ".csv.essen.csv")
files_ambig<-list.files("/Users/moradigd/Documents/Tradis/new/bamE", pattern = ".csv.ambig.csv")

output<-matrix(0, nrow = length(files_ess) , ncol = length(files_ess)   )
for(i in 1:length(files_ess)){
  for(j in 1:length(files_ess)){
    file_1_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/bamE/",files_ess[i]),show_col_types = FALSE) 
    file_2_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/bamE/",files_ess[j]),show_col_types = FALSE) 
    
    file_1_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/bamE/",files_ambig[i]),show_col_types = FALSE) 
    file_2_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/bamE/",files_ambig[j]),show_col_types = FALSE) 
    
    output[i,j]<-length(intersect( c(file_1_ess$locus_tag, file_1_ambig$locus_tag)   ,c(file_2_ess$locus_tag, file_2_ambig$locus_tag)   )  )/length(union(c(file_1_ess$locus_tag, file_1_ambig$locus_tag) , c(file_2_ess$locus_tag, file_2_ambig$locus_tag) ))
  }
}
colnames(output)<-paste0("run_", seq(dim(output)[1]))
row.names(output)<-paste0("run_", seq(dim(output)[1]))

library(corrplot)
corrplot(output, method = "number", type="upper",tl.col="black", tl.cex = 1, main="")





files_all<-list.files("/Users/moradigd/Documents/Tradis/new/bamE", pattern = ".csv.all.csv")
output<-matrix(0, nrow = length(files_all) , ncol = length(files_all)   )
for(i in 1:length(files_all)){
  for(j in 1:length(files_all)){
    file_1_all<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/bamE/",files_all[i]),show_col_types = FALSE) 
    file_2_all<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/bamE/",files_all[j]),show_col_types = FALSE) 
    output[i,j]<-cor(file_1_all$ins_index, file_2_all$ins_index)
  }
}
colnames(output)<-paste0("run_", seq(dim(output)[1]))
row.names(output)<-paste0("run_", seq(dim(output)[1]))

library(corrplot)
corrplot(output, method = "circle", type="upper",tl.col="black", tl.cex = 1, main="")


#wilde type
library(tidyverse)
gene<-"surA"
files_ess<-list.files(paste0("/Users/moradigd/Documents/Tradis/new/", gene), pattern = ".csv.essen.csv")
files_ess_wild<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/WT/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.essen.csv"))
files_ess_wild_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/WT/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.ambig.csv"))
files_ess_wild<-rbind(files_ess_wild, files_ess_wild_ambig)

unique_ess_forward<-c()
for(i in 1:length(files_ess)){
  file_input<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/",gene,"/",files_ess[i]),show_col_types = FALSE) 
  tmp_output<-file_input[which(!file_input$locus_tag %in% files_ess_wild$locus_tag),]
  tmp_output$file_tag<-files_ess[i]
  
  unique_ess_forward<-rbind(unique_ess_forward,tmp_output )
}

exclusive_genes_forward<-unique_ess_forward %>% 
  group_by(locus_tag, gene_name) %>% 
  summarise(n = n()) %>%
  filter(n==length(files_ess)) %>%
  select(-gene_name) %>%
  inner_join(unique_ess_forward, by=c("locus_tag"="locus_tag")) %>% 
  select(locus_tag, gene_name, start, end,  gene_length) %>%
  distinct_all(.)



#GO KEGG
library(edgeR)
keg<-kegga(exclusive_genes_forward$locus_tag, species.KEGG="eco")
go<-goana(exclusive_genes_forward$locus_tag, species.KEGG="eco")
topKEGG(keg)
topGO(go)

write_csv(topKEGG(keg), paste0("/Users/moradigd/Documents/Tradis/new/results/TopKEGG_",gene,".csv"))
write_csv(topKEGG(keg), paste0("/Users/moradigd/Documents/Tradis/new/results/TopGO_",gene,".csv"))
write_csv(exclusive_genes_forward, paste0("/Users/moradigd/Documents/Tradis/new/results/Exclusive_Essential_Gene_",gene,".csv"))

#comparative genes 
genes<-c("WT",	"bamB"	,"bamC"	,"bamE"	,"degP"	, "skp"  , "surA"	)

gene_extractor<-function(gene){
  files_ess_wild<-read_csv( paste0("/Users/moradigd/Documents/Tradis/new/",gene, "/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.essen.csv" ),show_col_types = FALSE)
  files_ess_wild_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/",gene, "/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.ambig.csv" ),show_col_types = FALSE)
  files_ess_wild<-rbind(files_ess_wild, files_ess_wild_ambig)
  return(files_ess_wild$locus_tag)
}

output<-matrix(0, nrow = length(genes) , ncol = length(genes)   )
for(i in 1:length(genes)){
  for(j in 1:length(genes)){
    output[i,j]<-length(intersect( gene_extractor(genes[i])   ,gene_extractor(genes[j])    )  )/length(union(gene_extractor(genes[i]) , gene_extractor(genes[j])  ))
  }
}
colnames(output)<-genes
row.names(output)<-genes


stats<-c()
for(i in 1:length(genes)){
  files_runs_tmp<-read_csv( paste0("/Users/moradigd/Documents/Tradis/new/",genes[i], "/fastqs.stats" ),show_col_types = FALSE)
  files_runs_tmp$gene<-genes[i]
  files_runs_tmp$sec<-paste0("run_", seq(dim(files_runs_tmp)[1]))
  stats<-rbind(stats, files_runs_tmp)
}
stats$info_reads<-stats$`% Matched`*stats$`% Mapped`*0.01


ggplot(stats, aes(fill=sec, x=gene,y=info_reads))+
  geom_bar( stat = "identity", position = "dodge")+
  theme_bw() +
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )+
  ylab("% Informative reads")+
  xlab("Genes")+ 
  ylim(range(0,100))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Sharing of essentiall genes plot 

library(tidyverse)
genes<-c("WT",	"bamB"	,"bamC"	,"bamE"	,"degP"	, "skp"  , "surA"	)
output<-matrix(0, nrow = length(genes) , ncol = length(genes)   )

for(i in 1:length(genes)){
  for(j in 1:length(genes)){
    file_1_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/",genes[i], "/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.essen.csv"  ),show_col_types = FALSE) 
    file_2_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/",genes[j], "/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.essen.csv"  ),show_col_types = FALSE) 
    
    file_1_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/",genes[i], "/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.ambig.csv"  ),show_col_types = FALSE) 
    file_2_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/new/",genes[j], "/total/merged.out_trimmed.fastq.NC_000913.3.tradis_gene_insert_sites.csv.ambig.csv"  ),show_col_types = FALSE) 
    
    output[i,j]<-length(intersect( c(file_1_ess$locus_tag, file_1_ambig$locus_tag)   ,c(file_2_ess$locus_tag, file_2_ambig$locus_tag)   )  )/length(union(c(file_1_ess$locus_tag, file_1_ambig$locus_tag) , c(file_2_ess$locus_tag, file_2_ambig$locus_tag) ))
    
  }
}
colnames(output)<-genes
row.names(output)<-genes

library(corrplot)
corrplot(output, method = "number", type="upper",tl.col="black", tl.cex = 1, main="")



#comparison with external 
library(tidyverse)
wt_1_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.all.csv")
wt_2_genes<-read_csv("/Users/moradigd/Documents/Tradis/essential_genes_ex.csv")


df<- data.frame(wt_1_genes$ins_index)
df$y<-wt_2_genes$`Insertion Index Score`[match(wt_1_genes$gene_name, wt_2_genes$Gene)]
colnames(df)<-c("Sample","External")
df<-df[complete.cases(df),]



ggplot(df, aes(x=Sample, y=External)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_text( size=13, angle=45,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )



plot(wt_1_genes$ins_index, wt_2_genes$`Insertion Index Score`[match(wt_1_genes$gene_name, wt_2_genes$Gene)])
cor(wt_1_genes$ins_index, wt_2_genes$`Insertion Index Score`[match(wt_1_genes$gene_name, wt_2_genes$Gene)],use = "pairwise.complete.obs")
cor(wt_1_genes$ins_index,wt_1_genes$ins_index)

wt_1_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.essen.csv")
x <- list(Sample=wt_1_genes$gene_name, External= wt_2_genes$Gene[wt_2_genes$Essential==TRUE])
library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)


wt_2_genes<-read_csv("/Users/moradigd/Documents/Tradis/essential_genes_ex.csv")
wt_ess_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.essen.csv")
wt_ambig_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.ambig.csv")

genes<-c(	"bamB"	,"bamC"	,"bamE"	,"degP"	, "skp"  , "surA"	)
df_out<-c()
for(i in 1:length(genes)){
  file_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.essen.csv"  )) 
  file_ess$gene<-genes[i]
  file_ess$status<-"essential"
  
  file_ess$essential_wt_external<-ifelse(file_ess$gene_name %in% wt_2_genes$Gene[wt_2_genes$Essential==TRUE], "shared", "not_shared" )
  file_ess$ambigious_wt_external<-ifelse(file_ess$gene_name %in% wt_2_genes$Gene[wt_2_genes$Unclear==TRUE], "shared", "not_shared" )
  
  file_ess$essential_wt_internal<-ifelse(file_ess$gene_name %in% wt_ess_genes$gene_name, "shared", "not_shared" )
  file_ess$ambigious_wt_internal<-ifelse(file_ess$gene_name %in% wt_ambig_genes$gene_name, "shared", "not_shared" )
  
  
  file_amb<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.ambig.csv"  )) 
  file_amb$gene<-genes[i]
  file_amb$status<-"ambig"
  file_amb$essential_wt_external<-ifelse(file_amb$gene_name %in% wt_2_genes$Gene[wt_2_genes$Essential==TRUE], "shared_with_external", "not_shared_with_external" )
  file_amb$ambigious_wt_external<-ifelse(file_amb$gene_name %in% wt_2_genes$Gene[wt_2_genes$Unclear==TRUE], "shared_with_external", "not_shared_with_external" )
  
  file_amb$essential_wt_internal<-ifelse(file_amb$gene_name %in% wt_ess_genes$gene_name, "shared", "not_shared" )
  file_amb$ambigious_wt_internal<-ifelse(file_amb$gene_name %in% wt_ambig_genes$gene_name, "shared", "not_shared" )
  
  df_out<-rbind(df_out, file_ess, file_amb)
}
View(df_out)
write_csv(df_out, "/Users/moradigd/Documents/Tradis/results_Danesh/results.csv")

#WT reverse new 
library(tidyverse)
wt_2_genes<-read_csv("/Users/moradigd/Documents/Tradis/essential_genes_ex.csv")
wt_ess_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.essen.csv")
wt_ambig_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.ambig.csv")

genes<-c(	"bamB"	,"bamC"	,"bamE"	,"degP"	, "skp"  , "surA"	)
result_tot<-c()
for(i in 1:6){
  result<-c()
  wt_tp<-wt_2_genes[wt_2_genes$`Non-essential`==FALSE & wt_2_genes$Unclear==FALSE,]
  result<-data.frame(cbind(result, wt_tp$Gene))
  result$status<-"external(Emily)_essential"
  
  label<-rep("not_shared", dim(wt_tp)[1])
  file_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.essen.csv"  )) 
  label[which(wt_tp$Gene %in% file_ess$gene_name)]<-"shared"
  result<-cbind(result,label)
  
  label<-rep("not_shared", dim(wt_tp)[1])
  file_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.ambig.csv"  )) 
  label[which(wt_tp$Gene %in% file_ambig$gene_name)]<-"shared"
  result<-cbind(result,label)
  
  colnames(result)<-c("genes","database_status", paste0("essential_", genes[i]), paste0("ambig_", genes[i]))
  
  result_amb<-c()
  wt_tp<-wt_2_genes[wt_2_genes$`Non-essential`==FALSE & wt_2_genes$Unclear==TRUE,]
  result_amb<-data.frame(cbind(result_amb, wt_tp$Gene))
  result_amb$status<-"external(Emily)_ambigious"
  
  label<-rep("not_shared", dim(wt_tp)[1])
  file_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.essen.csv"  )) 
  label[which(wt_tp$Gene %in% file_ess$gene_name)]<-"shared"
  result_amb<-cbind(result_amb,label)
  
  label<-rep("not_shared", dim(wt_tp)[1])
  file_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.ambig.csv"  )) 
  label[which(wt_tp$Gene %in% file_ambig$gene_name)]<-"shared"
  result_amb<-cbind(result_amb,label)
  
  colnames(result_amb)<-c("genes","database_status", paste0("essential_", genes[i]), paste0("ambig_", genes[i]))
  result_tot[[i]]<-rbind(result, result_amb)
}
result_tot<-cbind(result_tot[[1]], result_tot[[2]][,c(3,4)],result_tot[[3]][,c(3,4)],result_tot[[4]][,c(3,4)], result_tot[[5]][,c(3,4)], result_tot[[6]][,c(3,4)])
unique(result_tot$database_status)
tmp1<-result_tot

#wt with own 
wt_ess_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.essen.csv")
wt_ambig_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.ambig.csv")

genes<-c(	"bamB"	,"bamC"	,"bamE"	,"degP"	, "skp"  , "surA"	)
result_tot<-c()
for(i in 1:6){
  result<-c()
  wt_ess_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.essen.csv")
  
  wt_tp<-wt_ess_genes$gene_name
  result<-data.frame(cbind(result, wt_tp))
  result$status<-"Internal_essential"
  
  label<-rep("not_shared", length(wt_tp))
  file_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.essen.csv"  )) 
  label[which(wt_tp %in% file_ess$gene_name)]<-"shared"
  result<-cbind(result,label)
  
  label<-rep("not_shared", length(wt_tp))
  file_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.ambig.csv"  )) 
  label[which(wt_tp %in% file_ambig$gene_name)]<-"shared"
  result<-cbind(result,label)
  
  colnames(result)<-c("genes","database_status", paste0("essential_", genes[i]), paste0("ambig_", genes[i]))
  
  result_amb<-c()
  wt_ess_genes<-read_csv("/Users/moradigd/Documents/Tradis/results_Danesh/WT/tradis_gene_insert_sites.csv.ambig.csv")
  
  wt_tp<-wt_ess_genes$gene_name
  result_amb<-data.frame(cbind(result_amb, wt_tp))
  result_amb$status<-"Internal_ambigious"
  
  label<-rep("not_shared", length(wt_tp))
  file_ess<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.essen.csv"  )) 
  label[which(wt_tp %in% file_ess$gene_name)]<-"shared"
  result_amb<-cbind(result_amb,label)
  
  label<-rep("not_shared", length(wt_tp))
  file_ambig<-read_csv(paste0("/Users/moradigd/Documents/Tradis/results_Danesh/",genes[i], "/tradis_gene_insert_sites.csv.ambig.csv"  )) 
  label[which(wt_tp %in% file_ambig$gene_name)]<-"shared"
  result_amb<-cbind(result_amb,label)
  
  colnames(result_amb)<-c("genes","database_status", paste0("essential_", genes[i]), paste0("ambig_", genes[i]))
  result_tot[[i]]<-rbind(result, result_amb)
}
result_tot<-cbind(result_tot[[1]], result_tot[[2]][,c(3,4)],result_tot[[3]][,c(3,4)],result_tot[[4]][,c(3,4)], result_tot[[5]][,c(3,4)], result_tot[[6]][,c(3,4)])


write_csv(rbind(tmp1,result_tot), "/Users/moradigd/Documents/Tradis/wildtype_analysis.csv")



file_tags<-read_csv("/Users/moradigd/Documents/chemicalgenomic/morphology_app/gene_locus_file.csv")
gsub("\\.jpg","", file_tags$Tags)

paste0(input$gene,".jpg")


#Clsuters unsupervised learning
library(tidyverse)
library(edgeR)
library(circlize)

df_circularity<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_circularity.csv")
df_morphology<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_morphology.csv")
df_size<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colony_size.csv")
df_colour<-read_csv("/Users/moradigd/Documents/chemicalgenomic/df_clustering_similarity_colour.csv")

shared<-intersect(intersect(df_circularity$isolate, df_morphology$isolate), df_size$isolate)
df_morphology<-df_morphology[match(shared,df_morphology$isolate),]
df_size<-df_size[match(shared,df_size$isolate),]
df_colour<-df_colour[match(shared,df_size$isolate),]

genes<-c(
  "PA14_64050",
  "PA14_20700",
  "PA14_50060",
  "PA14_60870",
  "PA14_56790",
  'PA14_16500',
  "PA14_42220",
  "PA14_18550",
  "PA14_24510",
  "PA14_65540",
  'PA14_02110',
  "PA14_64050",
  'PA14_21190',
  "PA14_66320",
  "PA14_46030"
)
df_morphology$cluster[match(genes, df_morphology$isolate)]
df_size$cluster[match(genes, df_size$isolate)]
df_colour$cluster[match(genes, df_colour$isolate)]

#DolP paper
library("xlsx")
label<-c( "main","wt","yraP","aimA","amiB","amiC","envC","nlpD","ftsE","ftsX","dolPamiA","dolPamiB","dolPamiC","dolPenvC","dolPnlpD","dolPftsE","dolftsX","nlpDamiA","nlpDamiC")
data_wt <- read.xlsx("/Users/moradigd/Documents/dolP/TableS2.xlsx", 2) 

mean(data_wt$Total.length)

ttest<-c()
for(i in 3:length(label)){
  data_tmp <- read.xlsx("/Users/moradigd/Documents/dolP/TableS2.xlsx", i) 
  if(var.test(data_wt$Total.length, data_tmp$SHAPE.length)$p.value<0.05){
    ttest<-c(ttest, t.test(data_wt$Total.length, data_tmp$SHAPE.length, var.equal = FALSE)$p.value  )
  }
  else{
    ttest<-c(ttest,   t.test(data_wt$Total.length, data_tmp$SHAPE.length, var.equal = TRUE)$p.value )
  }
  print(var.test(data_wt$Total.length, data_tmp$SHAPE.length))
}
output<-data.frame(list(labels= label[3:length(label)], pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS2_pvalue_WTvsOthers.xlsx")

data_tmp_1 <- read.xlsx("/Users/moradigd/Documents/dolP/TableS2.xlsx", 11) 
data_tmp_2 <- read.xlsx("/Users/moradigd/Documents/dolP/TableS2.xlsx", 12) 
var.test(data_tmp_1$SHAPE.length, data_tmp_2$SHAPE.length)$p.value
t.test(data_tmp_1$SHAPE.length, data_tmp_2$SHAPE.length, var.equal = FALSE)$p.value
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS2_pvalue_dolPamiAvsdolPamiB.xlsx")


label<-c( "main","wt","yraP","bamB","bamC","bamE","amiA","yrap_amiA","bamBamiA","bamCamiA", "bamEamiA")
data_wt <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", 8) 
ttest<-c()
for(i in 9:length(label)){
  data_tmp <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", i) 
  if(var.test(data_wt$SHAPE.length, data_tmp$SHAPE.length)$p.value<0.05){
    ttest<-c(ttest, t.test(data_wt$SHAPE.length, data_tmp$SHAPE.length, var.equal = FALSE)$p.value  )
  }
  else{
    ttest<-c(ttest,   t.test(data_wt$SHAPE.length, data_tmp$SHAPE.length, var.equal = TRUE)$p.value )
  }
  print(var.test(data_wt$SHAPE.length, data_tmp$SHAPE.length))
}
output<-data.frame(list(labels= label[3:length(label)], pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS3_WTvsOthers.xlsx")



label<-c("bamBamiA","bamCamiA", "bamEamiA")
ttest<-c()
data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", 1) 
for(j in 8:10){
  res <- chisq.test(data_main[c(7,j),c(3,6)], correct=FALSE)$p.value
  ttest<-c(ttest,res)
}
output<-data.frame(list(labels= label, pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS3_septacell.xlsx")

colnames(data_main)[c(7),c(3,6)]
data_main[c(7,j),c(3,6)]



label<-c("bamBamiA","bamCamiA", "bamEamiA")
ttest<-c()
data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", 1) 
for(j in 8:10){
  res <- chisq.test(data_main[c(7,j),c(3,6)], correct=FALSE)$p.value
  ttest<-c(ttest,res)
}
output<-data.frame(list(labels= label, pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS3_septacell.xlsx")


label<-c("wt","yraP","aimA","amiB","amiC","envC","nlpD","ftsE","ftsX","dolPamiA","dolPamiB","dolPamiC","dolPenvC","dolPnlpD","dolPftsE","dolftsX","nlpDamiA","nlpDamiC")
data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS2.xlsx", 1) 

ttest<-c()
#data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", 1) 
for(j in 2:length(label)){
  res <- chisq.test(data_main[c(1,j),c(3,6)], correct=FALSE)$p.value
  ttest<-c(ttest,res)
}
output<-data.frame(list(labels= label[-1], pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS2_septacell.xlsx")




label<-c( "main","WT6.9","yraP6.9","yraPygeR6.9","ygeR6.9","WT5.2","yraP5.2","yraPygeR5.2","ygeR5.2")
ttest<-c()

set_1<-c(6,7,8,9)
set_2<-c(2,3,4,5)



for(i in 1:4){
  data_wt <- read.xlsx("/Users/moradigd/Documents/dolP/TableS4.xlsx", set_1[i]) 
  data_tmp <- read.xlsx("/Users/moradigd/Documents/dolP/TableS4.xlsx", set_2[i]) 
  if(var.test(data_wt$SHAPE.length, data_tmp$SHAPE.length)$p.value<0.05){
    ttest<-c(ttest, t.test(data_wt$SHAPE.length, data_tmp$SHAPE.length, var.equal = FALSE)$p.value  )
  }
  else{
    ttest<-c(ttest,   t.test(data_wt$SHAPE.length, data_tmp$SHAPE.length, var.equal = TRUE)$p.value )
  }
  print(var.test(data_wt$SHAPE.length, data_tmp$SHAPE.length))
}
output<-data.frame(list(labels= paste0(label[set_1],"_",label[set_2]), pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS4_WTvsOthers.xlsx")

ttest<-c()
#data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", 1) 
for(j in 1:4){
  res <- chisq.test(data_main[c(set_1[i],set_2[i]),c(3,6)], correct=FALSE)$p.value
  ttest<-c(ttest,res)
}
output<-data.frame(list(labels= paste0(label[set_1],"_",label[set_2]), pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS4_septacell.xlsx")



label<-c("bamBamiA","bamCamiA", "bamEamiA")
ttest<-c()
data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", 1) 
for(j in 8:10){
  res <- prop.test(x =c(data_main$Total.no..Of.septa[7],data_main$Total.no..Of.septa[j]) , n = c(data_main$No..Of.cells[7], data_main$No..Of.cells[j]) )
  ttest<-c(ttest,res$p.value)
}


t.test(c(1,1.00000001))



label<-c("wt","yraP","aimA","amiB","amiC","envC","nlpD","ftsE","ftsX","dolPamiA","dolPamiB","dolPamiC","dolPenvC","dolPnlpD","dolPftsE","dolftsX","nlpDamiA","nlpDamiC")
data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS2.xlsx", 1) 

ttest<-c()
#data_main <- read.xlsx("/Users/moradigd/Documents/dolP/TableS3.xlsx", 1) 
for(j in 2:length(label)){
  res <- chisq.test(data_main[c(1,j),c(3,6)], correct=FALSE)$p.value
  ttest<-c(ttest,res)
}
output<-data.frame(list(labels= label[-1], pvalue=ttest))
write_csv(output, "/Users/moradigd/Documents/dolP/pvalue_TableS2_septacell.xlsx")


library("xlsx")
data_1 <- read.xlsx("/Users/moradigd/Documents/dolP/dolPamiA_nlpDamiA.xlsx", 1) 
data_2 <- read.xlsx("/Users/moradigd/Documents/dolP/dolPamiA_nlpDamiA.xlsx", 3) 

chisq.test(data_main[c(1,2),c(2,3)], correct=FALSE)$p.value
t.test(data_1$SHAPE.length, data_2$SHAPE.length, var.equal = FALSE)$p.value
