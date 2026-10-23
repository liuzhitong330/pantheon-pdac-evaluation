options(warn=1)
suppressPackageStartupMessages(library(edgeR))
o<-Sys.getenv("ANALYSIS_OUTPUT","outputs/analysis")
spec<-read.delim(file.path(o,"cohort_specifications.tsv"),check.names=FALSE)
write_tsv<-function(x,p)write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
rows<-list()
for(i in seq_len(nrow(spec))) {
 id<-spec$analysis_id[i];p<-file.path(o,id);a<-readRDS(file.path(p,"analysis_objects.rds"));f<-a$fit
 stopifnot(length(f$s2.post)==nrow(a$y),all(is.finite(f$s2.post)),all(f$s2.post>0),all(is.finite(f$coefficients)),!any(f$failed),all(is.finite(f$fitted.values)),all(f$fitted.values>=0),all(is.finite(f$df.residual.adj)),all(f$df.residual.adj>=0),all(as.matrix(f$offset)==matrix(log(a$y$samples$lib.size*a$y$samples$norm.factors),nrow(a$y),ncol(a$y),byrow=TRUE)))
 # Save raw count input per analysis in a portable table in addition to RDS.
 write.table(a$counts,gzfile(file.path(p,"raw_counts.tsv.gz")),sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)
 pars<-read.delim(file.path(p,"gene_model_parameters.tsv"),check.names=FALSE)
 pars$NB_common_estimate<-a$y$common.dispersion;pars$NB_trended_estimate<-a$y$trended.dispersion;pars$NB_tagwise_estimate<-a$y$tagwise.dispersion
 pars$df_residual_adj_zero<-f$df.residual.adj<1e-8
 write_tsv(pars,file.path(p,"gene_model_parameters.tsv"))
 # These fields are named s2.post/prior in edgeR v4, not var.post/prior.
 # A zero adjusted residual df does not imply no QL information: moderated test uses prior df.
 res<-a$results;res$QL_df_residual_adj<-NA_real_;res$QL_df_prior<-NA_real_;res$QL_s2_post<-NA_real_;res$QL_df_total<-NA_real_
 res$QL_df_residual_adj[a$filter]<-f$df.residual.adj;res$QL_df_prior[a$filter]<-rep(f$df.prior,length.out=nrow(a$y));res$QL_s2_post[a$filter]<-f$s2.post
 if(!is.null(a$test))res$QL_df_total[a$filter]<-a$test$df.total
 res$flag_near_zero_adjusted_residual_df<-!is.na(res$QL_df_residual_adj)&res$QL_df_residual_adj<.01
 res$flag_zero_adjusted_residual_df<-!is.na(res$QL_df_residual_adj)&res$QL_df_residual_adj<1e-8
 write_tsv(res,file.path(p,"all_genes_results.tsv"));write_tsv(res[a$filter,],file.path(p,"tested_genes_results.tsv"))
 if(!is.null(a$test))write_tsv(res[order(res$FDR,na.last=TRUE),][which(res$FDR[order(res$FDR,na.last=TRUE)]<.05),],file.path(p,"significant_FDR05.tsv"))
 a$results<-res;saveRDS(a,file.path(p,"analysis_objects.rds"),compress="gzip")
 writeLines(capture.output(lapply(as.list(f),function(z)list(class=class(z),dim=dim(z),length=length(z)))),file.path(p,"fit_structure.txt"))
 rows[[id]]<-data.frame(analysis_id=id,n_failed=sum(f$failed),n_genes=nrow(a$y),n_nonfinite_coefficients=sum(!is.finite(f$coefficients)),min_s2post=min(f$s2.post),max_s2post=max(f$s2.post),min_df_residual_adj=min(f$df.residual.adj),n_zero_df_residual_adj=sum(f$df.residual.adj<1e-8),min_df_prior=min(f$df.prior),n_sig_zero_df=if(is.null(a$test))NA else sum(res$FDR_lt_0.05 & res$flag_zero_adjusted_residual_df),fit_NB_dispersion_min=min(f$dispersion),fit_NB_dispersion_max=max(f$dispersion))
}
write_tsv(do.call(rbind,rows),file.path(o,"fit_audit.tsv"))
sumtab<-read.delim(file.path(o,"analysis_summary.tsv"),check.names=FALSE)
sumtab[sumtab$role=="influence_diagnostic",c("n_FDR05","n_up","n_down","n_FDR05_absFC1")]<-NA
write_tsv(sumtab,file.path(o,"analysis_summary.tsv"))
cat("All fits audited: actual v4 s2.post/prior, convergence, coefficients, offsets, residual df and input count export verified.\n")
print(do.call(rbind,rows))
