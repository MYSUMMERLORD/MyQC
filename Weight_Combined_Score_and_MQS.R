library(ROCR)

user_parameter_vector = commandArgs(trailingOnly = TRUE)
quantile_all_file = user_parameter_vector[1]
outfile_folder = user_parameter_vector[2]
################Function#######################


function_AUC <- function(Score_vector,Lable_vector){
  digit_label = c()
  for(i in Lable_vector){
    digit_label = c(digit_label,switch (i,
      "GeneExpressionOutlier" = 0,
      "MainPopulationCell" = 1
    ))
  }
  if(!(max(digit_label)==min(digit_label))){
    pref = prediction(Score_vector,digit_label)
    perf = performance(pref,'auc')
    AUC = attributes(perf)$y.values[[1]]
  }else{
    AUC = 0
  }
  AUC
}

function_weight <- function(AUC){
  weight = (AUC - 0.5) / 0.5
  if(weight < 0){
    weight = 0
  }else{}
  weight
}

Stouffer.test <- function(p,w){ # p is a vector of p-values, here it is (1-quantile)
  if(missing(w) | max(w) == 0){
    w <- rep(0.0000000000000001,length(p))/length(p)
  }else{
    if(length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi_raw <- qnorm(1-p)
  Zi = c()
  for (i in Zi_raw) {
    if(i == Inf){
      Zi = c(Zi,8.5)
    }else if(i == (-Inf)){
      Zi = c(Zi,8.5)
    }else{
      Zi = c(Zi,i)
    }
  }
  Z <- sum(w * Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  if(Z == Inf){
    Z = 8.5
  }else if(Z==(-Inf)){
    Z= -8.5
  }else{}
  return (Z)
}
###############################################
data = read.table(quantile_all_file,header = T,sep = '\t',check.names=F)

expression_type = as.vector(data$Type)
quantile_mapped_reads = as.vector(data$Quantile.MappedReads)
quantile_mapping_rate = as.vector(data$Quantile.MappingRate)
quantile_library_complexity = as.vector(data$Quantile.LibraryComplexity)
quantile_gene_detected = as.vector(data$Quantile.GeneDetected)
min_quantile_score = pmin(quantile_mapped_reads,quantile_mapping_rate,quantile_library_complexity,quantile_gene_detected)

AUC_mapped_reads = function_AUC(quantile_mapped_reads,expression_type)
AUC_mapping_rate = function_AUC(quantile_mapping_rate,expression_type)
AUC_library_complexity = function_AUC(quantile_library_complexity,expression_type)
AUC_gene_detected = function_AUC(quantile_gene_detected,expression_type)

weight_mapped_reads = function_weight(AUC_mapped_reads)
weight_mapping_rate = function_weight(AUC_mapping_rate)
weight_library_complexity = function_weight(AUC_library_complexity)
weight_gene_detected = function_weight(AUC_gene_detected)

weight_vector = c(weight_mapped_reads,weight_mapping_rate,weight_library_complexity,weight_gene_detected)
header = c('Weight.MappedReads','Weight.MappingRate','Weight.LibraryComplexity','Weight.GeneDetected')
outmatrix_weight = rbind(header,weight_vector)
outfile_name = paste(outfile_folder,'Weighting_Factors.txt',sep = '')
write.table(x = outmatrix_weight,file = outfile_name,sep = '\t',quote = F,row.names = F,col.names = F)
names(weight_vector) = c('Weight.MappedReads','Weight.MappingRate','Weight.LibraryComplexity','Weight.GeneDetected')
fig_outfile_name = paste(outfile_folder,'Weighting_Factors.pdf',sep='')
pdf(fig_outfile_name,useDingbats = FALSE)
user_mar_offset = c(0.3,2,-0.2,0) #[Bottom,Left,Top,Right]
par(mar=par()$mar+user_mar_offset)
barplot(weight_vector,col = 'blue',xlab = '',ylab = "Importance",ylim = c(0,1),cex.lab=1.4)
dev.off()

weighted_combined_quality_score_vector = c()
for(i in 1:length(quantile_mapped_reads)){
  this_quantile_vector = c(quantile_mapped_reads[i],quantile_mapping_rate[i],quantile_library_complexity[i],quantile_gene_detected[i])
  this_score = Stouffer.test((1-this_quantile_vector),weight_vector)
  weighted_combined_quality_score_vector = c(weighted_combined_quality_score_vector,this_score)
}

outmatrix_score = cbind(data,min_quantile_score,weighted_combined_quality_score_vector)
header = c(colnames(data),'MinQuantileScore','WeightedCombinedQualityScore')
colnames(outmatrix_score) = header
outfile_name = paste(outfile_folder,'Samples_MQS_and_WeightedCombinedQualityScore.txt',sep = '')

write.table(x = outmatrix_score,file = outfile_name,sep = '\t',quote = F,row.names = F,col.names = T)
