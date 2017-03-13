user_parameter_vector = commandArgs(trailingOnly = TRUE)
distinct_PValue_file = user_parameter_vector[1]
distinct_PValue_mi_cutoff = as.numeric(user_parameter_vector[2])
distinct_PValue_pearson_cutoff = as.numeric(user_parameter_vector[3])
distinct_PValue_spearman_cutoff = as.numeric(user_parameter_vector[4])
logical_tag = user_parameter_vector[5]
sample_seqQC_file = user_parameter_vector[6]
outfile_folder = user_parameter_vector[7]

######### Function ##########

function_quantile <- function(Target_element, Control_vector)
{
  N_Positives = length(Control_vector[Target_element>=Control_vector])
  N_Negatives = length(Control_vector[Target_element<Control_vector])
  Quantile = N_Positives/(N_Positives+N_Negatives)
  Quantile
}

#############################
data = read.table(distinct_PValue_file,sep = '\t',header = T,check.names = F)
sample_name_vactor = as.vector(data[,1])
distinct_mi_p_value_vector = as.vector(data[,2])
distinct_pearson_p_value_vector = as.vector(data[,3])
distinct_spearman_p_value_vector = as.vector(data[,4])
sample_type_vector = c()
if(logical_tag == 'AND'){
  for(i in 1:length(sample_name_vactor)){
    if(distinct_mi_p_value_vector[i]<distinct_PValue_mi_cutoff & distinct_pearson_p_value_vector[i]<distinct_PValue_pearson_cutoff & distinct_spearman_p_value_vector[i]<distinct_PValue_spearman_cutoff){
      this_type = 'GeneExpressionOutlier'
    }else{
      this_type = 'MainPopulationCell'
    }
    sample_type_vector=c(sample_type_vector,this_type)
  }
}else if(logical_tag == 'OR'){
  for(i in 1:length(sample_name_vactor)){
    if(distinct_mi_p_value_vector[i]<distinct_PValue_mi_cutoff | distinct_pearson_p_value_vector[i]<distinct_PValue_pearson_cutoff | distinct_spearman_p_value_vector[i]<distinct_PValue_spearman_cutoff){
      this_type = 'GeneExpressionOutlier'
    }else{
      this_type = 'MainPopulationCell'
    }
    sample_type_vector=c(sample_type_vector,this_type)
  }
}else{
  print('Error......')
}

dataframe_sample_to_type = data.frame(cbind(sample_name_vactor,sample_type_vector))
colnames(dataframe_sample_to_type) = c('Sample.Name','Type')

data = read.table(sample_seqQC_file,sep = '\t',header = T,check.names = F)
merged_sample_type_to_qc = merge(dataframe_sample_to_type,data,by.x="Sample.Name",by.y="Sample.Name")

sample_name = as.vector(merged_sample_type_to_qc[,1])
sample_type = as.vector(merged_sample_type_to_qc[,2])
total_reads = as.vector(merged_sample_type_to_qc[,3])
mapped_reads = as.vector(merged_sample_type_to_qc[,4])
mapping_rate = as.vector(merged_sample_type_to_qc[,5])
library_complexity = as.vector(merged_sample_type_to_qc[,6])
gene_detected = as.vector(merged_sample_type_to_qc[,7])

total_reads_quantile_in_all_vector=c()
mapped_reads_quantile_in_all_vector = c()
mapping_rate_quantile_in_all_vector = c()
library_complexity_quantile_in_all_vector = c()
gene_detected_quantile_in_all_vector = c()

for (i in 1:length(sample_name)) {
  this_total_reads_quantile = function_quantile(total_reads[i],total_reads)
  this_mapped_reads_quantile = function_quantile(mapped_reads[i],mapped_reads)
  this_mapping_rate_quantile = function_quantile(mapping_rate[i],mapping_rate)
  this_library_complexity_quantile = function_quantile(library_complexity[i],library_complexity)
  this_gene_detected_quantile = function_quantile(gene_detected[i],gene_detected)
  total_reads_quantile_in_all_vector=c(total_reads_quantile_in_all_vector,this_total_reads_quantile)
  mapped_reads_quantile_in_all_vector = c(mapped_reads_quantile_in_all_vector,this_mapped_reads_quantile)
  mapping_rate_quantile_in_all_vector = c(mapping_rate_quantile_in_all_vector,this_mapping_rate_quantile)
  library_complexity_quantile_in_all_vector = c(library_complexity_quantile_in_all_vector,this_library_complexity_quantile)
  gene_detected_quantile_in_all_vector = c(gene_detected_quantile_in_all_vector,this_gene_detected_quantile)
  
}

out_matrix = cbind(as.matrix(merged_sample_type_to_qc),total_reads_quantile_in_all_vector,mapped_reads_quantile_in_all_vector,mapping_rate_quantile_in_all_vector,library_complexity_quantile_in_all_vector,gene_detected_quantile_in_all_vector)

header = c('Quantile.TotalReads','Quantile.MappedReads','Quantile.MappingRate','Quantile.LibraryComplexity','Quantile.GeneDetected')
header = c(colnames(merged_sample_type_to_qc),header)

out_matrix = rbind(header,out_matrix)

outfile_name = paste(outfile_folder,'Samples_Seq_Quality_in_all_Quantile.txt',sep = '')
write.table(x = out_matrix,file = outfile_name,sep = '\t',quote = F,row.names = F,col.names = F)
