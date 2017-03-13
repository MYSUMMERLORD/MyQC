user_parameter_vector = commandArgs(trailingOnly = TRUE)
file_format_FPR = user_parameter_vector[1]
outfile_folder = user_parameter_vector[2]

data = read.table(file_format_FPR,sep = '\t',header = T,check.names = F)

FPR = as.vector(data$False_Positive_Rate)
Artifact_Rate = as.vector(data$Fraction_Artifacts)
Fig_outfile_name = paste(outfile_folder,'False_Positive_Rate_Vs_Artifact_Rate_Plot.pdf',sep = '')

pdf(Fig_outfile_name,useDingbats = FALSE)
user_mar_offset = c(0.3,2,-0.2,0) #[Bottom,Left,Top,Right]
par(mar=par()$mar+user_mar_offset)
barplot(FPR,Artifact_Rate,xlab = 'False Positive Rate',ylab = "Artifact Rate in Gene Expression Outliers",xlim=c(0,1), ylim = c(0,1),type='1',lwd=1.5,cex.lab=1.4)
dev.off()