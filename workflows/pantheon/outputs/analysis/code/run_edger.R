options(warn=1)
suppressPackageStartupMessages(library(edgeR))
o <- Sys.getenv("ANALYSIS_OUTPUT","outputs/analysis"); prep <- Sys.getenv("PREP_OUTPUT","outputs/preparation")
dir.create(o,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x, file, row.names=FALSE) write.table(x,file,sep="\t",quote=FALSE,row.names=row.names,col.names=if(row.names) NA else TRUE,na="NA")
stopifnot(as.character(getRversion())=="4.4.2",as.character(packageVersion("edgeR"))=="4.4.0",as.character(packageVersion("limma"))=="3.62.0",as.character(packageVersion("statmod"))=="1.5.0",as.character(packageVersion("BiocVersion"))=="3.20.0")
writeLines(capture.output(sessionInfo()),file.path(o,"sessionInfo.txt"))
ip<-installed.packages();write_tsv(ip[,c("Package","Version","Built"),drop=FALSE],file.path(o,"package_manifest.tsv"))
writeLines(capture.output(lapply(list(filterByExpr.default=edgeR::filterByExpr.default,calcNormFactors.default=edgeR::calcNormFactors.default,estimateDisp.default=edgeR::estimateDisp.default,glmQLFit.default=edgeR::glmQLFit.default,glmFit.default=edgeR::glmFit.default,glmQLFTest=edgeR::glmQLFTest),formals)),file.path(o,"function_defaults.txt"))
features <- read.delim(file.path(prep,"feature_mapping.tsv"),check.names=FALSE)
read_counts <- function(f) {x<-as.matrix(read.delim(gzfile(file.path(prep,f)),row.names=1,check.names=FALSE));storage.mode(x)<-"double";stopifnot(identical(rownames(x),features$ensembl_id),all(is.finite(x)),all(x>=0),all(x==floor(x)));x}
counts <- read_counts("raw_pseudobulk_counts_all.tsv.gz")
subcounts <- read_counts("raw_subtype_pseudobulk_counts_all.tsv.gz")
meta <- read.delim(file.path(prep,"pseudobulk_metadata_all.tsv"),check.names=FALSE)
submeta <- read.delim(file.path(prep,"subtype_pseudobulk_metadata_all.tsv"),check.names=FALSE)
stopifnot(identical(colnames(counts),meta$pseudobulk_id),identical(colnames(subcounts),submeta$pseudobulk_id),all(colSums(counts)==meta$raw_library_umi),all(colSums(subcounts)==submeta$raw_library_umi))
spec <- read.delim(file.path(o,"cohort_specifications.tsv"),check.names=FALSE)
cache <- new.env();summaries<-list();warnings_all<-list()
for (ii in seq_len(nrow(spec))) {
 s<-spec[ii,]; id<-s$analysis_id; path<-file.path(o,id);dir.create(path,showWarnings=FALSE)
 cat("\nSTART",id,format(Sys.time(),tz="UTC"),"\n");flush.console()
 warn<-character()
 result<-withCallingHandlers({
 samples<-strsplit(s$samples,";",fixed=TRUE)[[1]]
 if(s$source=="classical") {m<-submeta[submeta$sample_id %in% samples & submeta$group %in% c("Normal Epithelial","Malignant - Classical"),];m$status<-ifelse(m$group=="Normal Epithelial","Normal","Malignant");xall<-subcounts} else {m<-meta[meta$sample_id %in% samples,];m$status<-m$group;xall<-counts}
 m<-m[order(match(m$sample_id,samples),match(m$status,c("Normal","Malignant"))),]
 m$patient<-factor(m$patient_id,levels=unique(m$patient_id));m$status<-factor(m$status,levels=c("Normal","Malignant"))
 x<-xall[,m$pseudobulk_id,drop=FALSE]
 stopifnot(identical(colnames(x),m$pseudobulk_id),nrow(m)==2*length(samples),all(table(m$patient,m$status)==1),all(colSums(x)==m$raw_library_umi))
 design<-model.matrix(~patient+status,data=m);rownames(design)<-m$pseudobulk_id
 stopifnot(qr(design)$rank==ncol(design),identical(levels(m$status),c("Normal","Malignant")),"statusMalignant" %in% colnames(design))
 contrast<-setNames(rep(0,ncol(design)),colnames(design));contrast["statusMalignant"]<-1
 for(pat in levels(m$patient)) {a<-which(m$patient==pat & m$status=="Malignant");b<-which(m$patient==pat & m$status=="Normal");stopifnot(all(design[a,]-design[b,]==contrast))}
 write_tsv(m,file.path(path,"metadata.tsv"));write_tsv(design,file.path(path,"design.tsv"),TRUE);write_tsv(data.frame(coefficient=names(contrast),weight=contrast),file.path(path,"contrast.tsv"))
 y0<-DGEList(counts=x,group=m$status,genes=features)
 if(s$filter_mode=="parent_mask") {parent<-get(s$parent,envir=cache);keep<-parent$keep;keep<-keep & rowSums(x)>0;filter_note<-"Parent gene mask; all-zero reduced-cohort genes flagged not estimable"} else {keep<-filterByExpr(y0,group=m$status,min.count=10,min.total.count=15,large.n=10,min.prop=.7);filter_note<-"Independent prespecified filterByExpr by status"}
 stopifnot(any(keep),!anyNA(keep))
 mask<-data.frame(ensembl_id=rownames(x),keep=keep,raw_total_umi=rowSums(x),status=ifelse(keep,"tested",ifelse(s$filter_mode=="parent_mask" & rowSums(x)==0,"not_estimable_all_zero_in_reduced_cohort","excluded_abundance_filter_or_parent_mask")))
 if(s$role=="influence_diagnostic")mask$status[keep]<-"effect_diagnostic_only_no_test"
 write_tsv(mask,file.path(path,"filtering_mask.tsv"))
 minsize<-min(table(m$status));if(minsize>10)minsize<-10+(minsize-10)*.7
 write_tsv(data.frame(filter_mode=s$filter_mode,note=filter_note,min_count=10,min_total_count=15,large_n=10,min_prop=.7,median_raw_library_size=median(colSums(x)),derived_CPM_cutoff=10/median(colSums(x))*1e6,independent_min_libraries=minsize,n_genes_retained=sum(keep),parent=s$parent),file.path(path,"filter_parameters.tsv"))
 y<-y0[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method="TMM",logratioTrim=.3,sumTrim=.05,doWeighting=TRUE,Acutoff=-1e10)
 norm<-data.frame(pseudobulk_id=m$pseudobulk_id,raw_library_size=colSums(x),retained_library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors,effective_library_size=y$samples$lib.size*y$samples$norm.factors)
 norm$log_offset<-log(norm$effective_library_size);write_tsv(norm,file.path(path,"normalization_offsets.tsv"))
 stopifnot(all(is.finite(norm$log_offset)),all(norm$effective_library_size>0))
 y<-estimateDisp(y,design,robust=TRUE)
 fit<-glmQLFit(y,design,robust=TRUE,abundance.trend=TRUE,legacy=FALSE)
 stopifnot(identical(rownames(fit$coefficients),rownames(y)),all(is.finite(fit$coefficients)),length(fit$s2.post)==nrow(y),all(is.finite(fit$s2.post)),all(fit$s2.post>0),!any(fit$failed))
 writeLines(capture.output(str(fit,max.level=1)),file.path(path,"fit_structure.txt"))
 genepar<-data.frame(ensembl_id=rownames(y))
 for(nm in c("dispersion","var.post","var.prior","df.prior","s2.prior","s2.post")) if(!is.null(fit[[nm]]))genepar[[nm]]<-rep(fit[[nm]],length.out=nrow(y))
 for(nm in c("df.residual","df.residual.zeros","df.residual.adj","deviance")) if(!is.null(fit[[nm]]))genepar[[nm]]<-rep(fit[[nm]],length.out=nrow(y))
 write_tsv(genepar,file.path(path,"gene_model_parameters.tsv"))
 logc<-cpm(y,log=TRUE,prior.count=.5)
 write_tsv(logc,gzfile(file.path(path,"normalized_logCPM.tsv.gz")),TRUE)
 pairdiff<-sapply(levels(m$patient),function(p)logc[,m$patient==p & m$status=="Malignant"]-logc[,m$patient==p & m$status=="Normal"])
 rownames(pairdiff)<-rownames(y);write_tsv(pairdiff,file.path(path,"patient_logCPM_differences.tsv"),TRUE)
 out<-features;out$analysis_id<-id;out$fdr_family<-s$fdr_family;out$n_pairs<-length(samples);out$test_status<-mask$status
 out$logFC<-NA_real_;out$logCPM<-NA_real_;out$F<-NA_real_;out$PValue<-NA_real_;out$FDR<-NA_real_
 if(s$role!="influence_diagnostic") {
  test<-glmQLFTest(fit,coef=which(colnames(design)=="statusMalignant"),poisson.bound=TRUE)
  tab<-test$table;stopifnot(all(is.finite(tab$PValue)),all(tab$PValue>=0 & tab$PValue<=1))
  tab$FDR<-p.adjust(tab$PValue,"BH")
  out[keep,c("logFC","logCPM","F","PValue","FDR")]<-tab[,c("logFC","logCPM","F","PValue","FDR")]
  if(!is.null(test$df.total)) write_tsv(data.frame(ensembl_id=rownames(y),QL_test_df_total=rep(test$df.total,length.out=nrow(y))),file.path(path,"test_df.tsv"))
 } else {test<-NULL;out$logFC[keep]<-fit$coefficients[,"statusMalignant"]/log(2);out$logCPM[keep]<-fit$AveLogCPM}
 out$FDR_lt_0.05<-!is.na(out$FDR)&out$FDR<.05
 out$abs_log2FC_ge_1<-!is.na(out$logFC)&abs(out$logFC)>=1
 out$paired_direction_agreement<-NA_real_;out$paired_direction_agreement[keep]<-rowMeans(sign(pairdiff)==sign(out$logFC[keep]))
 write_tsv(out,file.path(path,"all_genes_results.tsv"))
 write_tsv(out[keep,],file.path(path,"tested_genes_results.tsv"))
 if(s$role!="influence_diagnostic") {
  sig<-out[out$FDR_lt_0.05,];sig<-sig[order(sig$FDR,-abs(sig$logFC)),];write_tsv(sig,file.path(path,"significant_FDR05.tsv"))
 }
 saveRDS(list(specification=s,counts=x,metadata=m,design=design,contrast=contrast,filter=keep,y=y,fit=fit,test=test,results=out),file.path(path,"analysis_objects.rds"),compress="gzip")
 pdf(file.path(path,"diagnostics.pdf"),width=10,height=8)
 par(mfrow=c(2,2))
 barplot(norm$raw_library_size/1e6,names.arg=m$pseudobulk_id,las=2,cex.names=.5,main=paste(id,"raw UMIs (millions)"))
 barplot(norm$TMM_factor,names.arg=m$pseudobulk_id,las=2,cex.names=.5,main="TMM composition factors")
 plotBCV(y,main="NB dispersion (BCV)");plotQLDisp(fit,main="Robust v4 QL dispersion")
 par(mfrow=c(1,1))
 md<-plotMDS(y,labels=paste(m$sample_id,as.character(m$status)),col=ifelse(m$status=="Malignant","firebrick","steelblue"),main=paste(id,"MDS"),top=500,gene.selection="pairwise")
 write_tsv(data.frame(pseudobulk_id=m$pseudobulk_id,MDS1=md$x,MDS2=md$y),file.path(path,"MDS_coordinates.tsv"))
 if(!is.null(test)) {plotMD(test,status=as.integer(sign(out$logFC[keep])*out$FDR_lt_0.05[keep]),main=paste(id,"Malignant minus Normal"));abline(h=c(-1,1),lty=2,col="grey")}
 dev.off()
 assign(id,list(keep=keep,results=out),envir=cache)
 info<-data.frame(analysis_id=id,role=s$role,fdr_family=s$fdr_family,contains_wu16="wu16" %in% samples,n_pairs=length(samples),n_cells_Malignant=sum(m$n_cells[m$status=="Malignant"]),n_cells_Normal=sum(m$n_cells[m$status=="Normal"]),n_features=nrow(x),n_retained=sum(keep),rank=qr(design)$rank,nominal_residual_df=nrow(m)-ncol(design),n_FDR05=sum(out$FDR_lt_0.05),n_up=sum(out$FDR_lt_0.05 & out$logFC>0),n_down=sum(out$FDR_lt_0.05 & out$logFC<0),n_FDR05_absFC1=sum(out$FDR_lt_0.05 & out$abs_log2FC_ge_1),common_NB_dispersion=y$common.dispersion,min_effective_library=min(norm$effective_library_size),max_effective_library=max(norm$effective_library_size))
 info
 },warning=function(w){msg<-conditionMessage(w);warn<<-c(warn,msg);cat("WARNING",id,msg,"\n");invokeRestart("muffleWarning")})
 writeLines(if(length(warn))warn else "No warnings.",file.path(path,"warnings.txt"))
 if(s$role=="influence_diagnostic") result[,c("n_FDR05","n_up","n_down","n_FDR05_absFC1")]<-NA
 warnings_all[[id]]<-data.frame(analysis_id=id,warning=if(length(warn))warn else "NONE")
 summaries[[id]]<-result;write_tsv(do.call(rbind,summaries),file.path(o,"analysis_summary.tsv"));write_tsv(do.call(rbind,warnings_all),file.path(o,"warnings_all.tsv"))
 cat("COMPLETE",id,"genes",result$n_retained,"FDR05",result$n_FDR05,"\n");flush.console()
}
cat("All analyses completed. Primary is primary7_min10_no_wu16; sensitivities are separate families, LOO effects only.\n")
