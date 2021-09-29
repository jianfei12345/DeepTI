###loading package######
library(igraph)
library(dnet)
library(ggplot2)
library(data.table)
library(keras)

###loading PPI and preprocessing####
PPI<-fread('F:\\data\\string database\\9606.protein.links.v11.0.txt\\9606.protein.links.v11.0.txt')
protein_info<-fread('F:\\data\\string database\\9606.protein.info.v11.0.txt\\9606.protein.info.v11.0.txt',sep = '\t')
PPI$protein2<-protein_info$preferred_name[match(PPI$protein2,protein_info$protein_external_id)]
PPI$protein1<-protein_info$preferred_name[match(PPI$protein1,protein_info$protein_external_id)]
PPI<-PPI[,-3]
rm(protein_info)

###construction of PPI network####
PPI_igraph<-graph_from_data_frame(PPI,directed = F)
PPI_igraph<-simplify(PPI_igraph,remove.multiple = T,remove.loops = T,edge.attr.comb="median")
rm(PPI)

###perform dRWR algorithm###
PTmatrix <- dRWR(g=PPI_igraph, normalise="laplacian", restart=0.995,parallel=T,multicores = 6)
PPI_similarity <- as.matrix(PTmatrix)
rownames(PPI_similarity)<-rownames(PTmatrix)
colnames(PPI_similarity)<-colnames(PTmatrix)

###ppmi normalization###
PPI_similarity<-log2((PPI_similarity*sum(PPI_similarity))/(matrix(rowSums(PPI_similarity),ncol=1) %*% matrix(colSums(PPI_similarity),nrow=1)))
PPI_similarity[PPI_similarity<0]<-0

###loading label###
label <- fread('label.txt')

###preparation of dataset###
PPI_similarity <- PPI_similarity[match(label$V1,rownames(PPI_similarity)),]
DATA<-list()
DATA$feature<-PPI_similarity
DATA$labels<-label$V2

###construction of model###
model <- keras_model_sequential()
model %>% layer_dense(units = 100, input_shape = c(19339),activation = 'relu') %>% layer_dense(units = 16,  activation = 'relu') %>% layer_dense(units = 2,activation = 'sigmoid')
model %>% compile(
       optimizer = 'adam', 
       loss = 'sparse_categorical_crossentropy',
       metrics = c('accuracy')
   )
model %>% fit(DATA$train$feature[r,], DATA$train$label[r,], epochs = 20, verbose = 2)
model %>% evaluate(DATA$train$feature[-r,],DATA$train$label[-r,] , verbose = 0)

###construction of decision tree model###
library("mlr3")
task_imvigor210<-as_task_classif(temp,target = 'response')
learner<-lrn("classif.rpart", predict_type = "prob")
rr = resample(task_imvigor210, learner, rsmp("cv",folds=10), store_models = TRUE)
rr$aggregate(msr("classif.acc"))

###codes for figure 2###
library(wordcloud2)
wordcloud2(inhibit,size = 0.2)
wordcloud2(promote,size = 0.2)

###codes for figure 3###
#figure 3a: volcano plot#
library(ggplot2)
library(ggrepel)
ggplot(RESULT$fig1$DEG,aes(x=logFC,y=-log10(PValue)))
+geom_point(aes(col=regulation),alpha = 0.8, size = 1)
+scale_color_manual(values = c('Blue', 'gray30', 'red2'))
+labs(x="log2 fold change",  y="-log10 padj", title="")
+theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.position = c(0.8, 0.8)) + theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'))+geom_hline(yintercept=-log10(0.05),linetype=4, color = 'gray', size = 0.25) +  geom_vline(xintercept=c(-1,1),linetype=4, color = 'gray', size = 0.25)
+geom_text_repel(aes(x=logFC, y=-log10(PValue),label=rownames(RESULT$fig1$D)),size=2)
+theme(plot.title = element_text(hjust = 0.5)

#figure 3b: forestplot#
library(forestplot)
forestplot(
  as.matrix(data[,c(1:3,8:9,7)]),                    #文字
  data$V10,                                                     #HR
  data$V11,                                                     #low 95%
  data$V12,                                                      #up 95%
  graph.pos = 4,                                              #森林图位置
  is.summary = c(1,rep(0,24)),                       #第一行总结
  zero = 1,
  graphwidth = unit(90,"mm"),
  line.margin = unit(10,"mm"),
  lwd.xaxis = 2,
  lwd.zero = 2,
  lwd.ci = 2,
  ci.vertices = T,
  ci.vertices.height =0.2,
  title = "Multivariate analysis of factors that influence overall survival in 6,028 ESCC patients",
  col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
  cex=0.9, 
  lineheight = "auto",
  boxsize=0.5, colgap=unit(6,"mm"),
  txt_gp=fpTxtGp                                        #字体大小设置
  (
    label=gpar(cex=1.25),
    ticks=gpar(cex=1.1),
    xlab=gpar(cex = 1.2),
    title=gpar(cex = 1.2)
  ),
  hrzl_lines=list                                                  #行添加颜色
  (
    "2" = gpar(lwd=1, col="#99999922"),
    "7" = gpar(lwd=70, lineend="butt", columns=c(2:7), col="#99999922"), 
    "13" = gpar(lwd=60, lineend="butt", columns=c(2:7), col="#99999922"), 
    "21" = gpar(lwd=100, lineend="butt", columns=c(2:7), col="#99999922")
  )
)

#figure 3c-e: survival curve#
library(survminer)
ggsurvplot(survfit(Surv(A1_OS,event)~consensusclusterplus,
                   data=filter(RESULT$prognostic_value$preprocessed,consensusclusterplus!=4))
        ,data=filter(RESULT$prognostic_value$preprocessed,consensusclusterplus!=4),legend=c(0.7,0.8),legend.title="consensus group",pval = T,legend.labs=c("C1","C2","C3","C5"))

#figure 4: correlation heatmap matrix#
library(corrplot)
corrplot(res2$r,method = 'color',type = 'lower',number.cex = 0.7,cl.pos = 'n',tl.col = 1,tl.pos = 'full',addCoef.col = 1,p.mat = res2$P,sig.level = 0.3,insig = 'blank',addgrid.col = 1)

#figure 5: ROC curve#
library(pROC)
roc1<-roc(predictions[,3],predictions[,2],plot=TRUE,print.auc=TRUE,legacy.axes=T)
plot(roc1, col="blue")