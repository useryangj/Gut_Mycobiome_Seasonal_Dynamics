
GENUS_FILE   <- file.path("C:/Users/xuzy0/Desktop/xzy/split/chouping/genus_table1.txt")
DIET_FILE    <- file.path( "C:/Users/xuzy0/Desktop/xzy/每月食物结构.xlsx")
CLIMATE_FILE <- file.path("C:/Users/xuzy0/Desktop/xzy/每月温度湿度.xlsx")
OUT_DIR <- file.path("Monthly_RDA_final1")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
N_PERM <- 9999
SEED <- 20261106
MONTHS <- c("Jul","Aug","Sep","Oct","Nov","Dec",
            "Jan","Feb","Mar","Apr","May","Jun")
DIET_VARS <- c("Young leaves","Mature leaves","Fruits","Flowers",
               "Invertebrate","Seeds","Roots","Macro-fungi",
               "Barks","Fallen leaves","Others")
ENV_VARS <- c(DIET_VARS,"Temperature","Humidity")
FALLBACK_COMPONENTS <- c("Mature leaves","Roots","Barks","Fallen leaves")
CORE_PREVALENCE <- 0.50
CORE_MEAN_ABUND <- 0.001
HEATMAP_MAX_GENERA <- 20
pkgs <- c("vegan","readxl","dplyr","tibble","ggplot2","pheatmap","ggrepel")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Please install packages: ",
                         paste(missing,collapse=", "),
                         "\ninstall.packages(c(",
                         paste(sprintf('"%s"',missing),collapse=", "),"))")
suppressPackageStartupMessages({
  library(vegan); library(readxl); library(dplyr)
  library(tibble); library(ggplot2); library(pheatmap); library(ggrepel)
})
norm_month <- function(x) {
  z <- trimws(as.character(x))
  z <- sub("^2021-01-01.*","Jan",z)
  z <- sub("^202[0-9]-01-01.*","Jan",z)
  replacements <- c("Sept"="Sep","September"="Sep",
    "January"="Jan","February"="Feb","March"="Mar","April"="Apr",
    "June"="Jun","July"="Jul","August"="Aug",
    "October"="Oct","November"="Nov","December"="Dec")
  hit <- z %in% names(replacements)
  z[hit] <- unname(replacements[z[hit]])
  z
}
as_numeric_table <- function(x) {
  x[] <- lapply(x, function(a) suppressWarnings(as.numeric(as.character(a))))
  x
}
read_monthly <- function(path, label, expected) {
  d <- as.data.frame(readxl::read_excel(path),check.names=FALSE)
  names(d)[1] <- "Month"
  d$Month <- norm_month(d$Month)
  d <- d[!is.na(d$Month) & d$Month != "",,drop=FALSE]
  if (anyDuplicated(d$Month)) stop(label,": duplicated month names")
  if (!setequal(d$Month, MONTHS)) {
    stop(label,": expected Jul-Jun. Found: ",paste(d$Month,collapse=", "))
  }
  if (!all(expected %in% names(d)))
    stop(label,": missing column(s): ",
         paste(setdiff(expected,names(d)),collapse=", "))
  d <- d[match(MONTHS,d$Month),c("Month",expected),drop=FALSE]
  d[expected] <- as_numeric_table(d[expected])
  if (anyNA(d[expected])) stop(label,": missing/non-numeric values in predictors")
  d
}
test_rda <- function(model, seed, permutations=N_PERM, by=NULL) {
  set.seed(seed)
  if (is.null(by))
    vegan::anova.cca(model, permutations=permutations)
  else
    vegan::anova.cca(model, permutations=permutations,by=by)
}
r2_values <- function(model) {
  s <- vegan::RsquareAdj(model)
  c(R2=unname(s$r.squared), AdjR2=unname(s$adj.r.squared))
}
extract_test <- function(a) {
  data.frame(F=as.numeric(a$F[1]),
             P=as.numeric(a$`Pr(>F)`[1]))
}
adjusted_r2_for <- function(y,x) {
  stopifnot(nrow(y)==nrow(x))
  vegan::RsquareAdj(vegan::rda(y, as.matrix(x)))$adj.r.squared
}
safe_vif <- function(model) {
  v <- vegan::vif.cca(model)
  data.frame(Predictor=names(v),VIF=unname(as.numeric(v)))
}
save_pdf <- function(file, plot, w, h) {
  ggplot2::ggsave(filename=file.path(OUT_DIR,file),
                  plot=plot,width=w,height=h,units="in",device="pdf")
}
if (!all(file.exists(c(GENUS_FILE,DIET_FILE,CLIMATE_FILE))))
  stop("Input file missing. Set BASE_DIR correctly.")
g0 <- read.table(GENUS_FILE,header=TRUE,row.names=1,sep="\t",
                 check.names=FALSE,quote="",comment.char="",fill=FALSE)
counts <- t(as.matrix(g0))
storage.mode(counts) <- "numeric"
if (anyNA(counts) || any(counts<0)) stop("Genus count matrix has NA/negative values.")
if (anyDuplicated(rownames(counts))) stop("Duplicated sample ID")
if (anyDuplicated(colnames(counts))) stop("Duplicated taxon label")
sample_month <- norm_month(sub(".*\\.","",rownames(counts)))
if (!setequal(unique(sample_month),MONTHS))
  stop("Some sample months are missing/invalid: ",
       paste(unique(sample_month),collapse=", "))
depth <- rowSums(counts)
if (any(depth<=0)) stop("All-zero sample(s) in genus count matrix")
diet <- read_monthly(DIET_FILE,"Diet",DIET_VARS)
climate <- read_monthly(CLIMATE_FILE,"Climate",c("Temperature","Humidity"))
env <- cbind(diet,climate[,c("Temperature","Humidity"),drop=FALSE])
stopifnot(identical(env$Month,MONTHS))
if (any(vapply(env[ENV_VARS],sd,numeric(1))==0))
  stop("Constant monthly predictor(s) detected.")
diet_sums <- rowSums(env[DIET_VARS])
if (any(abs(diet_sums-100)>1)) warning("Some monthly diet totals differ from 100% by >1 percentage point.")
sample_n <- as.integer(table(factor(sample_month,levels=MONTHS)))
season <- c("Summer","Summer","Summer","Autumn","Autumn","Winter",
            "Winter","Winter","Winter","Spring","Spring","Spring")
qc <- data.frame(Month=MONTHS,Season=season,
                 N_samples=sample_n,
                 Diet_sum_pct=round(diet_sums,3),
                 Temperature=env$Temperature,Humidity=env$Humidity)
write.csv(qc,file.path(OUT_DIR,"01_monthly_QC.csv"),row.names=FALSE)
sample_relative <- sweep(counts,1,depth,"/")
monthly_relative <- vapply(MONTHS,function(m)
  colMeans(sample_relative[sample_month==m,,drop=FALSE]),numeric(ncol(counts)))
monthly_relative <- t(monthly_relative)
rownames(monthly_relative) <- MONTHS
colnames(monthly_relative) <- colnames(counts)
stopifnot(all(abs(rowSums(monthly_relative)-1)<1e-8))
Y <- vegan::decostand(monthly_relative,method="hellinger")
rownames(Y) <- MONTHS
E <- as.data.frame(env[ENV_VARS],check.names=FALSE)
rownames(E) <- MONTHS
single <- lapply(seq_along(ENV_VARS),function(i) {
  name <- ENV_VARS[i]
  x <- scale(as.numeric(E[[name]]))[,1]
  mod <- vegan::rda(Y ~ x)
  a <- test_rda(mod,SEED+i)
  cbind(data.frame(Variable=name,N_months=12),
        as.data.frame(as.list(r2_values(mod))),
        extract_test(a))
})
single <- dplyr::bind_rows(single)
single$BH_q_exploratory <- p.adjust(single$P,method="BH")
write.csv(single,file.path(OUT_DIR,"02_single_factor_RDA.csv"),row.names=FALSE)
single$Variable <- factor(single$Variable,levels=rev(ENV_VARS))
p_single <- ggplot(single,aes(x=Variable,y=100*R2)) +
  geom_col(fill="#347F91",width=.7) +
  geom_text(aes(label=sprintf("%.1f%%",100*R2)),hjust=-.15,size=3) +
  coord_flip(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0,.22))) +
  labs(x=NULL,y="Unadjusted explained variation (%)",
       title="Monthly fungal community: single-factor RDA",
       subtitle="13 separate exploratory models; n = 12 months") +
  theme_classic(base_size=11) + theme(plot.title=element_text(face="bold"))
save_pdf("Figure_1_single_factor_RDA.pdf",p_single,9,6)
fallback <- rowSums(E[FALLBACK_COMPONENTS])
clim <- scale(E[,c("Temperature","Humidity"),drop=FALSE])
climate_pca <- prcomp(clim,center=FALSE,scale.=FALSE)
climate_pc1 <- as.numeric(climate_pca$x[,1])
if (cor(climate_pc1,E$Temperature)<0) climate_pc1 <- -climate_pc1
climate_loading <- climate_pca$rotation[,1]
if (cor(as.numeric(climate_pca$x[,1]),E$Temperature)<0)
  climate_loading <- -climate_loading
climate_variance <- summary(climate_pca)$importance["Proportion of Variance",1]
X0 <- data.frame(Macro_fungi=E[["Macro-fungi"]],
                 Fallback_foods=fallback,Climate_PC1=climate_pc1,
                 row.names=MONTHS)
if (any(vapply(X0,sd,numeric(1))==0))
  stop("Reduced model has a constant predictor")
X <- as.data.frame(scale(X0))
rownames(X) <- MONTHS
model <- vegan::rda(Y ~ Macro_fungi + Fallback_foods + Climate_PC1,data=X)
overall <- extract_test(test_rda(model,SEED+100))
marg <- test_rda(model,SEED+101,by="margin")
marg_tbl <- data.frame(Predictor=rownames(marg),
                       F=as.numeric(marg$F),
                       P=as.numeric(marg$`Pr(>F)`))
marg_tbl <- marg_tbl[marg_tbl$Predictor!="Residual",,drop=FALSE]
marg_tbl$BH_q_exploratory <- p.adjust(marg_tbl$P,method="BH")
vif <- safe_vif(model)
summary_model <- data.frame(
  N_months=12,N_predictors=3,
  R2=unname(r2_values(model)["R2"]),
  AdjR2=unname(r2_values(model)["AdjR2"]),
  F_overall=overall$F,P_overall_exploratory=overall$P,
  Climate_PC1_variance=climate_variance,
  Temp_loading=unname(climate_loading["Temperature"]),
  Humidity_loading=unname(climate_loading["Humidity"]),
  Diet_climate_correlation=cor(X$Fallback_foods,X$Climate_PC1),
  N_free_permutations=N_PERM
)
write.csv(summary_model,file.path(OUT_DIR,"03_reduced_RDA_model.csv"),row.names=FALSE)
write.csv(marg_tbl,file.path(OUT_DIR,"04_reduced_RDA_marginal.csv"),row.names=FALSE)
write.csv(vif,file.path(OUT_DIR,"05_reduced_RDA_VIF.csv"),row.names=FALSE)
sites <- as.data.frame(scores(model,display="sites",choices=1:2,scaling=2))
names(sites)[1:2] <- c("RDA1","RDA2"); sites$Month <- rownames(sites)
arrows <- as.data.frame(scores(model,display="bp",choices=1:2,scaling=2))
names(arrows)[1:2] <- c("RDA1","RDA2"); arrows$Variable <- rownames(arrows)
scale_factor <- min(.8*max(abs(sites$RDA1))/max(abs(arrows$RDA1)),
                    .8*max(abs(sites$RDA2))/max(abs(arrows$RDA2)))
if (!is.finite(scale_factor)) scale_factor <- 1
arrows$RDA1 <- arrows$RDA1*scale_factor
arrows$RDA2 <- arrows$RDA2*scale_factor
eig <- vegan::eigenvals(model,model="constrained")
axis_pct <- 100*eig/sum(eig)
p_joint <- ggplot(sites,aes(RDA1,RDA2)) +
  geom_hline(yintercept=0,linewidth=.25,color="grey60") +
  geom_vline(xintercept=0,linewidth=.25,color="grey60") +
  geom_point(size=2.7,color="#347F91") +
  ggrepel::geom_text_repel(aes(label=Month),size=3.6,seed=SEED) +
  geom_segment(data=arrows,aes(x=0,y=0,xend=RDA1,yend=RDA2),
               inherit.aes=FALSE,linewidth=.7,color="#B26C41",
               arrow=grid::arrow(length=grid::unit(.16,"cm"))) +
  ggrepel::geom_text_repel(data=arrows,
                           aes(x=RDA1,y=RDA2,label=Variable),
                           inherit.aes=FALSE,size=3.7,seed=SEED+1) +
  labs(title="Reduced monthly RDA",
       subtitle="Macro-fungi + fallback foods + climate PC1 (n=12)",
       x=sprintf("RDA1 (%.1f%% of constrained inertia)",axis_pct[1]),
       y=sprintf("RDA2 (%.1f%% of constrained inertia)",axis_pct[2])) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"))
save_pdf("Figure_2_reduced_RDA.pdf",p_joint,8,6)
Xd <- X[,c("Macro_fungi","Fallback_foods"),drop=FALSE]
Xc <- X[,"Climate_PC1",drop=FALSE]
vp <- vegan::varpart(Y,Xd,Xc)
r_d <- adjusted_r2_for(Y,Xd)
r_c <- adjusted_r2_for(Y,Xc)
r_all <- unname(r2_values(model)["AdjR2"])
frac <- c(Diet_unique=r_all-r_c,Climate_unique=r_all-r_d,
          Shared=r_d+r_c-r_all,Unexplained=1-r_all)
partition <- data.frame(Fraction=names(frac),Adjusted_R2=as.numeric(frac))
write.csv(partition,file.path(OUT_DIR,"06_variance_partition.csv"),row.names=FALSE)
partition$Fraction <- factor(partition$Fraction,
  levels=c("Diet_unique","Climate_unique","Shared","Unexplained"))
p_vp <- ggplot(partition,aes(Fraction,100*Adjusted_R2,fill=Fraction)) +
  geom_col(width=.62,show.legend=FALSE) +
  geom_hline(yintercept=0,linewidth=.25) +
  geom_text(aes(label=sprintf("%.1f%%",100*Adjusted_R2)),
            vjust=ifelse(partition$Adjusted_R2>=0,-.4,1.4),size=3.5) +
  scale_fill_manual(values=c(Diet_unique="#347F91",
                             Climate_unique="#D6A25A",
                             Shared="#8D76AA",
                             Unexplained="#B8C5C6")) +
  labs(x=NULL,y="Adjusted R² fraction (%)",
       title="Exploratory variance partitioning",
       subtitle="Diet (2 predictors) vs climate PC1 (1 predictor)") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=15,hjust=1),
        plot.title=element_text(face="bold"))
save_pdf("Figure_3_variance_partition.pdf",p_vp,8,5)
loo <- lapply(seq_len(12),function(i) {
  yi <- Y[-i,,drop=FALSE]
  xi <- X[-i,,drop=FALSE]
  mi <- vegan::rda(yi ~ Macro_fungi + Fallback_foods + Climate_PC1,data=xi)
  data.frame(Left_out_month=MONTHS[i],
             R2=unname(r2_values(mi)["R2"]),
             AdjR2=unname(r2_values(mi)["AdjR2"]),
             Max_VIF=max(vegan::vif.cca(mi)))
})
loo <- dplyr::bind_rows(loo)
write.csv(loo,file.path(OUT_DIR,"07_leave_one_month_out.csv"),row.names=FALSE)
shift_tbl <- do.call(rbind,lapply(0:11,function(k) {
  idx <- ((seq_len(12)-1+k) %% 12)+1
  xs <- X[idx,,drop=FALSE]
  ms <- vegan::rda(Y,as.matrix(xs))
  data.frame(Shift_months=k,AdjR2=unname(r2_values(ms)["AdjR2"]),
             R2=unname(r2_values(ms)["R2"]))
}))
write.csv(shift_tbl,file.path(OUT_DIR,"08_circular_shift_sensitivity.csv"),
          row.names=FALSE)
rda_res <- residuals(model)
res_pc1 <- stats::prcomp(rda_res,center=TRUE,scale.=FALSE)$x[,1]
lag1_resid <- cor(res_pc1[-1],res_pc1[-length(res_pc1)])
lag1_clim <- cor(X$Climate_PC1[-1],X$Climate_PC1[-12])
lag1_fallback <- cor(X$Fallback_foods[-1],X$Fallback_foods[-12])
prevalence <- colMeans(counts>0)
mean_rel <- colMeans(sample_relative)
core_all <- names(prevalence)[prevalence>=CORE_PREVALENCE &
                              mean_rel>=CORE_MEAN_ABUND]
if (length(core_all)>=2) {
  rank_order <- order(-prevalence[core_all],-mean_rel[core_all])
  core_plot <- core_all[rank_order][seq_len(min(length(core_all),HEATMAP_MAX_GENERA))]
  rho <- pval <- matrix(NA_real_,nrow=length(core_plot),ncol=length(ENV_VARS),
                       dimnames=list(core_plot,ENV_VARS))
  for (i in seq_along(core_plot)) for(j in seq_along(ENV_VARS)) {
    y <- monthly_relative[,core_plot[i]]
    x <- E[[ENV_VARS[j]]]
    if (sd(y)>0 && sd(x)>0) {
      s <- suppressWarnings(cor.test(y,x,method="spearman",exact=FALSE))
      rho[i,j] <- unname(s$estimate); pval[i,j] <- s$p.value
    }
  }
  q <- matrix(p.adjust(as.vector(pval),method="BH"),
              nrow=nrow(pval),ncol=ncol(pval),dimnames=dimnames(pval))
  heat_results <- expand.grid(Genus=rownames(rho),
                              Variable=colnames(rho),stringsAsFactors=FALSE)
  heat_results$Rho <- as.vector(rho)
  heat_results$P_exploratory <- as.vector(pval)
  heat_results$BH_q_exploratory <- as.vector(q)
  write.csv(heat_results,file.path(OUT_DIR,"09_core_genera_correlations.csv"),
            row.names=FALSE)
  stars <- matrix(ifelse(is.finite(q) & q<0.05,"*",""),
                  nrow=nrow(q),dimnames=dimnames(q))
  pheatmap::pheatmap(rho,cluster_cols=FALSE,cluster_rows=TRUE,
    display_numbers=stars,number_color="black",
    color=colorRampPalette(c("#4575b4","white","#d73027"))(101),
    breaks=seq(-1,1,length.out=102),border_color=NA,
    main="Prevalent fungal taxa: monthly correlations (exploratory)",
    filename=file.path(OUT_DIR,"Figure_4_genus_correlation_heatmap.pdf"),
    width=11,height=max(5,3+.24*nrow(rho)),angle_col=45)
} else {
  warning("Fewer than two core taxa meet the preset thresholds; heatmap skipped.")
}
audit <- c(
  "Monthly community / diet / climate analysis (EXPLORATORY)",
  "Samples are fecal collections, not identified macaque individuals.",
  "Single observed Jul-Jun annual cycle; month n=12, temporal dependence possible.",
  paste("Total fecal samples:",nrow(counts)),
  paste("Taxonomic entries in count table:",ncol(counts)),
  paste("Monthly sample sizes:",paste(paste0(MONTHS,":",sample_n),collapse="; ")),
  paste("Seasons defined as meteorological: ",
        paste(paste(names(tapply(sample_n,season,sum)),
                    as.integer(tapply(sample_n,season,sum)),sep="="),
              collapse="; ")),
  paste("Reduced model R2:",round(summary_model$R2,4)),
  paste("Reduced model adjusted R2:",round(summary_model$AdjR2,4)),
  paste("Exploratory unrestricted-permutation P:",round(summary_model$P_overall_exploratory,6)),
  paste("Climate PC1 variance explained:",round(100*climate_variance,1),"%"),
  paste("Fallback vs climate PC1 correlation:",round(summary_model$Diet_climate_correlation,4)),
  paste("VIF:",paste(paste(vif$Predictor,round(vif$VIF,3),sep="="),collapse="; ")),
  paste("Leave-one-month-out adjusted R2 range:",
        paste(round(range(loo$AdjR2),4),collapse=" to ")),
  paste("Lag-1 RDA residual PC1 correlation:",round(lag1_resid,4)),
  paste("Lag-1 climate PC1 correlation:",round(lag1_clim,4)),
  paste("Lag-1 fallback food correlation:",round(lag1_fallback,4)),
  paste("Circular shift observed AdjR2:",round(shift_tbl$AdjR2[1],4)),
  paste("Circular shift other max AdjR2:",round(max(shift_tbl$AdjR2[-1]),4)),
  paste("Core taxa meeting prevalence and abundance thresholds:",length(core_all)),
  "INTERPRETATION: All P and BH q refer to unrestricted exploratory permutations.",
  "BH adjustment does not repair temporal dependence or unmeasured repeat sampling.",
  "Single-factor R2 are overlapping associations, not independent contributions.",
  "Variance fractions are adjusted-R2 associations, not causal effects.",
  "Circular-shift comparisons are diagnostic; not confirmatory P-values.",
  "No inference about within-individual longitudinal change is supported.",
  "",
  "Session information:",capture.output(sessionInfo())
)
writeLines(audit,file.path(OUT_DIR,"10_analysis_readme.txt"))
cat(paste(audit[1:20],collapse="\n"),"\n\nOutput folder:",OUT_DIR,"\n")