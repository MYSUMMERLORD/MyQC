library('stats')

User_parameter_vector=commandArgs(trailingOnly=TRUE)
TPM_Matrix_file=User_parameter_vector[1]
SinQC_information_file=User_parameter_vector[2]
Outfile_folder=User_parameter_vector[3]

###
data=read.csv(TPM_Matrix_file,sep='\t', header=T, check.names=F)

data_matrix=as.matrix(data[1:dim(data)[1],3:dim(data)[2]])
rownames(data_matrix)=data$Gene.ID

PCA=prcomp(data_matrix)
component=PCA$rotation
variance=PCA$sdev


PCA_out1_file=paste(Outfile_folder,'PCA_variance.txt',sep='')
write.table(x=cbind(colnames(component),variance), file=PCA_out1_file, quote=F, sep='\t', row.names=F, col.names=T)

header_component=c('Sample',colnames(component))
component_out=cbind(rownames(component),component)

component_out=rbind(header_component,component_out)

PCA_out2_file=paste(Outfile_folder,'PCA_Samples.txt',sep='')
write.table(x=component_out, file=PCA_out2_file, quote=F, sep='\t', row.names=F, col.names=F)


Data_frame_component=data.frame(component)
Selected_component=cbind(Data_frame_component$PC1,Data_frame_component$PC2,Data_frame_component$PC3)
colnames(Selected_component)=c('PC1','PC2','PC3')
rownames(Selected_component)=rownames(Data_frame_component)
Selected_component=data.frame(Selected_component)

##
data=read.table(SinQC_information_file,sep='\t',header=T, check.names=F)
QC=data.frame(data$QC)
rownames(QC)=data$Sample.Name

Merged_QC_and_PC=merge(QC, Selected_component, by="row.names")
Finial_QC=as.vector(Merged_QC_and_PC[,2])

color_vector=c()
for (i in 1:length(Finial_QC)) {
  if(Finial_QC[i]=='PASS') {
    color_vector=c(color_vector,'blue')
  } else {color_vector=c(color_vector,'red') }
}

Fig_outfile_name=paste(Outfile_folder,'PCA_Plot.pdf',sep='')
pdf(Fig_outfile_name,useDingbats = FALSE)
pairs(~PC1+PC2+PC3,Merged_QC_and_PC, main = "Principal Component Analysis (PCA)", pch=20, col=color_vector,cex=1)
par(xpd=TRUE)

legend(x='topright',c('PASS','Artifact'), col=c('blue','red'), pch=20, horiz=T, bty='n')

dev.off()


