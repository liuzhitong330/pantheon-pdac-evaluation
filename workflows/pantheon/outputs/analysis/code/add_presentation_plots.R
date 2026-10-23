options(warn=1)
suppressPackageStartupMessages(library(edgeR))
o<-Sys.getenv("ANALYSIS_OUTPUT","outputs/analysis")
a<-readRDS(file.path(o,"primary7_min10_no_wu16","analysis_objects.rds"))
p<-read.delim(file.path(o,"primary_all_genes_annotated.tsv"),check.names=FALSE)
examples<-read.delim(file.path(o,"primary_example_genes.tsv"),check.names=FALSE)
m<-a$metadata;logc<-cpm(a$y,log=TRUE,prior.count=.5)
pdf(file.path(o,"primary_expression_and_robustness.pdf"),width=11,height=8.5)
par(mfrow=c(2,3),mar=c(6,4,3,1))
for(g in examples$ensembl_id) {
 symbol<-examples$gene_symbol[examples$ensembl_id==g]
 plot(c(.8,2.2),range(logc[g,])+c(-.5,.5),type="n",xaxt="n",xlab="",ylab="TMM log2 CPM (prior count 0.5)",main=symbol)
 axis(1,at=1:2,labels=c("Normal","Malignant"))
 for(i in seq_along(levels(m$patient))) {
  pat<-levels(m$patient)[i];pos<-which(m$patient==pat)
  lines(as.integer(m$status[pos]),logc[g,pos],col=i,type="b",pch=16)
 }
}
par(mfrow=c(1,1),mar=c(5,5,4,2))
plot.new();legend("center",legend=levels(m$patient),col=seq_along(levels(m$patient)),lty=1,pch=16,title="Patient colors in paired expression plots")
mat<-logc[examples$ensembl_id,,drop=FALSE];z<-t(scale(t(mat)));z[!is.finite(z)]<-0
image(seq_len(ncol(z)),seq_len(nrow(z)),t(z[nrow(z):1,,drop=FALSE]),col=colorRampPalette(c("navy","white","firebrick"))(101),axes=FALSE,xlab="",ylab="",main="Primary examples: gene-wise standardized logCPM")
axis(1,at=seq_len(ncol(z)),labels=m$pseudobulk_id,las=2,cex.axis=.6);axis(2,at=seq_len(nrow(z)),labels=rev(examples$gene_symbol),las=2,cex.axis=.7)
par(mfrow=c(2,3),mar=c(4,4,3,1))
spec<-read.delim(file.path(o,"cohort_specifications.tsv"),check.names=FALSE)
for(id in spec$analysis_id[spec$role %in% c("qc_sensitivity","prespecified_sensitivity")]) {
 s<-read.delim(file.path(o,id,"all_genes_results.tsv"));idx<-match(p$ensembl_id,s$ensembl_id);keep<-is.finite(p$logFC)&is.finite(s$logFC[idx])
 plot(p$logFC[keep],s$logFC[idx][keep],pch=16,cex=.25,col=ifelse(p$FDR_lt_0.05[keep],"firebrick","grey70"),xlab="Primary7 log2FC",ylab="Sensitivity log2FC",main=id,cex.main=.75);abline(0,1,col="blue",lty=2)
}
par(mfrow=c(1,2))
for(which in c("primary7","original8")) {
 z<-read.delim(file.path(o,paste0("influence_",which,"_gene_summary.tsv")),check.names=FALSE)
 keep<-z$parent_FDR<.05
 plot(z$parent_logFC[keep],z$max_abs_loo_effect_change[keep],pch=16,cex=.4,xlab="Parent log2FC",ylab="Largest absolute LOO effect change",main=paste(which,"influence: parent discoveries"))
}
dev.off()
# Exact MDS and dispersion/offset summaries are also available in per-analysis tables and RDS.
cat("Paired expression, gene-wise heatmap, sensitivity effect comparison and distinct LOO plots saved.\n")
