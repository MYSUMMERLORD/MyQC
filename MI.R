user_parameter_vector = commandArgs(trailingOnly=TRUE)
expression_matrix_file = user_parameter_vector[1]
outfile_folder = user_parameter_vector[2]

########## Function ##################

library(infotheo)

function_Cor_One_Vs_Others_Pearson <- function(One_vector,Others_Matrix)
{
  cor_vector=as.vector(cor(One_vector,Others_Matrix,method='pearson'))
  cor_vector
}

function_Cor_One_Vs_Others_Spearman <- function(One_vector,Others_Matrix)
{
  cor_vector=as.vector(cor(One_vector,Others_Matrix,method='spearman'))
  cor_vector
}


function_Cor_Pairwise_to_vector_Pearson <- function(Input_Matrix)
{
  cor_matrix=cor(Input_Matrix,method='pearson')
  pairwise_cor_vector=cor_matrix[upper.tri(cor_matrix,diag=F)]
  pairwise_cor_vector
}

function_Cor_Pairwise_to_vector_Spearman <- function(Input_Matrix)
{
  cor_matrix=cor(Input_Matrix,method='spearman')
  pairwise_cor_vector=cor_matrix[upper.tri(cor_matrix,diag=F)]
  pairwise_cor_vector
}

function_Cor_One_Vs_Others_MI <- function(i,MI_matrix)
{
  cor_vector <- c()
  cor_vector <- MI_matrix[i,-c(i)]
  #for(i in 1:dim(Others_Matrix)[2]){
  #  mi <- mutinformation(One_vector,Others_Matrix[,i])
  #  cor_vector <- c(cor_vector,mi)
  #}
  cor_vector
}

function_Cor_Pairwise_to_vector_MI <- function(MI_matrix)
{
  pairwise_cor_vector <- c()
  for(i in 1:dim(MI_matrix)[2]){
    pairwise_cor_vector <- c(pairwise_cor_vector,MI_matrix[i,-c(i)])
  }
  #for(i in 1:dim(Input_Matrix)[2]){
  #  for(j in 1:dim(Input_Matrix)[2]){
  #    if(i!=j){
  #      mi <- mutinformation(Input_Matrix[,i],Input_Matrix[,j])
  #      pairwise_cor_vector <- c(pairwise_cor_vector,mi)
  #    }
  #  }
  #}
  pairwise_cor_vector
}

#expression_matrix_file="/home/liulj/single_cell/data/tpm.csv"

data = read.table(expression_matrix_file,sep='\t',header=TRUE,check.names=F)
inmatrix = as.matrix(data[2:dim(data)[2]])
sample_name_vector= colnames(inmatrix)

TPM_matrix <- inmatrix

P_value_Spearman_vector=c()


for (i in 1:dim(TPM_matrix)[2]) {
  One_vector=TPM_matrix[,i]
  Others_Matrix=TPM_matrix[,-i]
  a=function_Cor_One_Vs_Others_Spearman(One_vector,Others_Matrix)
  b=function_Cor_Pairwise_to_vector_Spearman(Others_Matrix)
  P_value_Spearman=wilcox.test(a,b,alternative='less')$p.value
  P_value_Spearman_vector=c(P_value_Spearman_vector,P_value_Spearman)
}
P_value_Pearson_vector=c()
for (i in 1:dim(TPM_matrix)[2]) {
  One_vector=TPM_matrix[,i]
  Others_Matrix=TPM_matrix[,-i]
  a=function_Cor_One_Vs_Others_Pearson(One_vector,Others_Matrix)
  b=function_Cor_Pairwise_to_vector_Pearson(Others_Matrix)
  P_value_Pearson=wilcox.test(a,b,alternative='less')$p.value
  P_value_Pearson_vector=c(P_value_Pearson_vector,P_value_Pearson)
}


discdata <- discretize(t(inmatrix), "equalwidth") 

TPM_matrix <- t(discdata)
P_value_MI_vector=c()

sampleNum <- dim(TPM_matrix)[2]
MI_matrix <- array(1:sampleNum*sampleNum,c(sampleNum,sampleNum))

for(i in 1:(sampleNum-1)){
  for (j in (i+1):sampleNum) {
    mi = entropy(TPM_matrix[,i])-condentropy(TPM_matrix[,i],TPM_matrix[,j])
    MI_matrix[i,j]=mi
    MI_matrix[j,i]=mi
  }
}

for (i in 1:sampleNum) {
  #One_vector=TPM_matrix[,i]
  #Others_Matrix=TPM_matrix[,-i]
  a=function_Cor_One_Vs_Others_MI(i,MI_matrix)
  b=function_Cor_Pairwise_to_vector_MI(MI_matrix[-c(i),-c(i)])
  P_value_MI=wilcox.test(a,b,alternative='less')$p.value
  P_value_MI_vector=c(P_value_MI_vector,P_value_MI)
}

Out_Matrix = cbind(sample_name_vector,P_value_MI_vector,P_value_Pearson_vector,P_value_Spearman_vector)

header = c('Sample','Distinct.PValue.MI','Distinct.PValue.Pearson','Distinct.PValue.Spearman')
out_matrix = rbind(header,Out_Matrix)
outfile_name = paste(outfile_folder,'Distinct_PValue.txt',sep='')
#outfile_name = "/home/liulj/single_cell/data/MI2.txt"
write.table(x=out_matrix,file=outfile_name,sep='\t',quote=F,row.names=F,col.names=F)
