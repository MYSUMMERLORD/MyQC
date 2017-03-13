user_parameter_vector = commandArgs(trailingOnly = TRUE)
file_samples_MQS_and_WeightedCombinedQualityScore = user_parameter_vector[1]
outfile_folder = user_parameter_vector[2]

data = read.table(file_samples_MQS_and_WeightedCombinedQualityScore,sep = '\t',header = T,check.names = F)

type = as.vector(data$Type)
MQS = as.vector(data$MinQuantileScore)
WCQS = as.vector(data$WeightedCombinedQualityScore)
total_MainPopulationCell = length(type[type=='MainPopulationCell'])
total_GeneExpressionOutlier = length(type[type=='GeneExpressionOutlier'])

pool_MQS = c(0,sort(MQS))
WCQS_offset = sort(WCQS)[2]-sort(WCQS)[1]
add_WCQS = sort(WCQS)[1]-WCQS_offset
pool_WCQS = c(add_WCQS,sort(WCQS)) # i don't understand why add a add_WCQS?

MQS_cutoff_or = c()
WCQS_cutoff_or = c()

N_GeneExpressionOutlier_AND_EqLower = c()
N_MainPopulationCell_AND_EqLower = c()

for(i in pool_MQS){
  for(j in pool_WCQS){
    this_N_GeneExpressionOutlier_AND_EqLower = length(type[type=='GeneExpressionOutlier' & (MQS<=i | WCQS<=j)])
    this_N_MainPopulationCell_AND_EqLower = length(type[type=='MainPopulationCell' & (MQS<=i | WCQS<=j)])
    N_GeneExpressionOutlier_AND_EqLower = c(N_GeneExpressionOutlier_AND_EqLower,this_N_GeneExpressionOutlier_AND_EqLower)
    N_MainPopulationCell_AND_EqLower = c(N_MainPopulationCell_AND_EqLower,this_N_MainPopulationCell_AND_EqLower)
    MQS_cutoff_or = c(MQS_cutoff_or,i)
    WCQS_cutoff_or = c(WCQS_cutoff_or,j)
  }
}

Fraction_Artifacts = N_GeneExpressionOutlier_AND_EqLower/total_GeneExpressionOutlier
False_Positive_Rate = N_MainPopulationCell_AND_EqLower/total_MainPopulationCell

out_matrix = cbind(MQS_cutoff_or,WCQS_cutoff_or,N_GeneExpressionOutlier_AND_EqLower,Fraction_Artifacts,N_MainPopulationCell_AND_EqLower,False_Positive_Rate)
header = c('MQS_OR','wCQS_OR','Artifacts','Fraction_Artifacts','False_Positive','False_Positive_Rate')
out_matrix = rbind(header,out_matrix)

outfile_name = paste(outfile_folder,'False_Positive_Rate_Table_With_all_combinations.txt',sep='')
write.table(x=out_matrix,file = outfile_name,sep = '\t',quote = F,row.names = F,col.names = F)




