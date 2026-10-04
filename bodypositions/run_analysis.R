run_analysis <- function(dvtype=c("Physiology","Self-Reports","Behavior"),
                         prereg=c(1,2), poseposture=c(1,2), culture=c(1,2),
                         coverstory=c(1,2), checked=c(1,2),
                         year=c(1982,2023), outliers=c(-6,6),
                         data_path="bp_full.csv", outdir=".") {

  ppbase <- read.csv(data_path, stringsAsFactors=FALSE)
  ppbase$year <- as.numeric(gsub("[^0-9]","",ppbase$reference))
  ppbase$include <- ppbase$include == "TRUE" | ppbase$include == TRUE

  pp <- ppbase[ppbase$include == TRUE, ]
  pp <- pp[!is.na(pp$include), ]
  pp <- pp[as.numeric(pp$se) != 0, ]
  pp$esg <- as.numeric(pp$g)
  pp <- pp[!is.na(as.numeric(pp$culture)), ]
  pp <- pp[!is.na(as.numeric(pp$pose_posture)), ]
  pp <- pp[!is.na(as.numeric(pp$coverstory)), ]
  pp$checked <- 1
  pp$reference_studyno <- paste(pp$reference, ", Study ", pp$study_no, sep="")
  pp$vi <- metafor::escalc(measure="SMD", yi=pp$esg, sei=as.numeric(pp$se), ni=as.numeric(pp$nhighlow), data=pp)$vi
  pp$zg <- as.numeric(scale(pp$esg))

  # filters
  if (!("Physiology" %in% dvtype)) pp <- pp[pp$dvtype != "Physiology", ]
  if (!("Self-Reports" %in% dvtype)) pp <- pp[pp$dvtype != "Self-Reports", ]
  if (!("Behavior" %in% dvtype)) pp <- pp[pp$dvtype != "Behavior", ]
  if (!(1 %in% culture)) pp <- pp[pp$culture != 0, ]
  if (!(2 %in% culture)) pp <- pp[pp$culture != 1, ]
  if (!(1 %in% poseposture)) pp <- pp[pp$pose_posture != 0, ]
  if (!(2 %in% poseposture)) pp <- pp[pp$pose_posture != 1, ]
  if (!(1 %in% coverstory)) pp <- pp[pp$coverstory != 0, ]
  if (!(2 %in% coverstory)) pp <- pp[pp$coverstory != 1, ]
  pp <- pp[pp$year >= year[1], ]
  pp <- pp[pp$year <= year[2], ]
  pp <- pp[pp$zg >= outliers[1], ]
  pp <- pp[pp$zg <= outliers[2], ]
  if (!(1 %in% prereg)) pp <- pp[as.numeric(pp$preregistered) != 1, ]
  if (!(2 %in% prereg)) pp <- pp[as.numeric(pp$preregistered) != 0, ]
  if (!(1 %in% checked)) pp <- pp[pp$checked != 0, ]
  if (!(2 %in% checked)) pp <- pp[pp$checked != 1, ]
  pp <- pp[!is.na(pp$esg), ]
  pp$esg <- as.numeric(pp$esg); pp$se <- as.numeric(pp$se); pp$vi <- as.numeric(pp$vi)
  pp$id2 <- as.numeric(pp$id2); pp$id1 <- as.numeric(pp$id1)
  pp$nhighlow <- as.numeric(pp$nhighlow); pp$nhigh <- as.numeric(pp$nhigh); pp$nlow <- as.numeric(pp$nlow); pp$ntotal <- as.numeric(pp$ntotal)
  pp$year <- as.numeric(pp$year)

  # sample overview
  nstudies <- length(unique(pp$reference))
  nindependent <- length(unique(pp$reference_studyno))
  neffects <- length(pp$reference)
  nsample <- sum(aggregate(nhighlow ~ id2, data=transform(pp, nhighlow=as.numeric(nhighlow)), FUN=function(x) min(x))$nhighlow)
  sample <- data.frame(Articles=nstudies, Studies=nindependent, Effects=neffects, Samplesize=round(nsample,0))

  # rma.mv total
  metapp_total <- tryCatch(metafor::rma.mv(yi=esg, V=vi, random=~1|id2/id1, tdist=TRUE, data=pp, method="ML"), error=function(e) NULL)
  het <- data.frame(Sigma2_Level1=metapp_total$sigma2[1], Sigma2_Level2=metapp_total$sigma2[2],
                    Tau=metapp_total$tau2, Q=round(metapp_total$QE,2),
                    Q_p=ifelse(round(metapp_total$QEp,3)==0, "< .001", round(metapp_total$QEp,3)))

  # aggregated
  pptemp <- pp[!is.na(pp$se) & !is.na(pp$esg), ]
  pp_agg <- aggregate(esg ~ reference_studyno, data=pptemp, FUN="mean")
  pp_agg$se <- aggregate(se ~ reference_studyno, data=pptemp, FUN="mean")$se
  pp_agg$ntotal <- aggregate(ntotal ~ reference_studyno, data=pptemp, FUN="min")$ntotal
  pp_agg$reference <- aggregate(reference ~ reference_studyno, data=pptemp, FUN="min")$reference
  pp_agg$nhighlow <- aggregate(nhighlow ~ reference_studyno, data=pptemp, FUN="min")$nhighlow
  pp_agg$vi <- metafor::escalc(measure="SMD", yi=pp_agg$esg, sei=pp_agg$se, ni=pp_agg$nhighlow, data=pp_agg)$vi

  # metagen (for trim-and-fill / funnel / eggers)
  metapp_metagen_mean <- tryCatch(metagen(TE=esg, seTE=se, data=pp_agg, studlab=pp_agg$reference_studyno, fixed=TRUE, random=TRUE, method.tau="ML", hakn=TRUE, prediction=TRUE, sm="SMD"), error=function(e) NULL)

  # 9 models
  taf_agg <- tryCatch(metafor::trimfill(metapp_metagen_mean), error=function(e) NULL)
  punif_agg <- tryCatch(puniform::puni_star(yi=pp_agg$esg, vi=pp_agg$se^2, alpha=.05, side="right", method="ML", boot=FALSE), error=function(e) NULL)
  pet_agg <- lm(esg~sqrt(vi), data=pp_agg, weights=1/vi)
  peese_agg <- lm(esg~vi, data=pp_agg, weights=1/vi)
  rve <- tryCatch(robumeta::robu(esg~1, data=pp, studynum=reference_studyno, var.eff.size=se*sqrt(nhighlow), small=FALSE), error=function(e) NULL)
  hvsm <- tryCatch(weightr::weightfunct(pp_agg$esg, pp_agg$vi, steps=c(0.025,1), fe=FALSE), error=function(e) NULL)

  pppos <- pp[!is.na(pp$p_position), ]
  pp_first <- pppos[(pppos$p_position==1 | pppos$p_position==3), ]
  pp_last  <- pppos[(pppos$p_position==2 | pppos$p_position==3), ]
  pc_first <- tryCatch(pcurve(data.frame(TE=pp_first$esg, seTE=pp_first$se, studlab=pp_first$reference, n=pp_first$nhighlow, id=pp_first$id2, nhighlow=pp_first$nhighlow), effect.estimation=TRUE, N=pp_first$nhighlow, dmin=0, dmax=1), error=function(e) NULL)
  pc_last <- tryCatch(pcurve(data.frame(TE=pp_last$esg, seTE=pp_last$se, studlab=pp_last$reference, n=pp_last$nhighlow, id=pp_last$id2, nhighlow=pp_last$nhighlow), effect.estimation=TRUE, N=pp_last$nhighlow, dmin=0, dmax=1), error=function(e) NULL)

  estimates <- data.frame(Model=as.character(), g=as.numeric(), LCL=as.numeric(), UCL=as.numeric(), k=as.numeric(), stringsAsFactors=FALSE)
  estimates[1,] <- c("Random-Effects Multilevel Model", metapp_total$b, metapp_total$ci.lb, metapp_total$ci.ub, metapp_total$k)
  if (!is.null(rve)) estimates[2,] <- c("Robust Variance Estimation", rve$reg_table$b.r, rve$reg_table$CI.L, rve$reg_table$CI.U, length(rve$k)) else estimates[2,] <- c("Robust Variance Estimation", NA, NA, NA, NA)
  if (!is.null(taf_agg)) estimates[3,] <- c("Trim-and-fill", taf_agg$TE.random, taf_agg$lower.random, taf_agg$upper.random, taf_agg$k) else estimates[3,] <- c("Trim-and-fill", NA, NA, NA, NA)
  if (!is.null(punif_agg)) estimates[4,] <- c("P-uniform star", punif_agg$est, punif_agg$ci.lb, punif_agg$ci.ub, punif_agg$k) else estimates[4,] <- c("P-uniform star", NA, NA, NA, NA)
  if (!is.null(hvsm)) estimates[5,] <- c("Hedges-Vevea Selection Model", hvsm$output_adj$par[2], hvsm$ci.lb_adj[2], hvsm$ci.ub_adj[2], hvsm$k) else estimates[5,] <- c("Hedges-Vevea Selection Model", NA, NA, NA, NA)
  if (!is.null(pc_first)) estimates[6,] <- c("P-Curve (first value)", tryCatch(pc_first$dEstimate, error=function(e) NA), NA, NA, tryCatch(pc_first$kAnalyzed, error=function(e) 0)) else estimates[6,] <- c("P-Curve (first value)", NA, NA, NA, 0)
  if (!is.null(pc_last)) estimates[7,] <- c("P-Curve (last value)", tryCatch(pc_last$dEstimate, error=function(e) NA), NA, NA, tryCatch(pc_last$kAnalyzed, error=function(e) 0)) else estimates[7,] <- c("P-Curve (last value)", NA, NA, NA, 0)
  estimates[8,] <- c("Precision Effect Test", as.numeric(pet_agg$coefficients[1]), confint(pet_agg)[1,1], confint(pet_agg)[1,2], pet_agg$df+2)
  estimates[9,] <- c("Precision Effect Estimate using Standard Error", as.numeric(peese_agg$coefficients[1]), confint(peese_agg)[1,1], confint(peese_agg)[1,2], peese_agg$df+2)
  estimates$g <- as.numeric(estimates$g); estimates$LCL <- as.numeric(estimates$LCL); estimates$UCL <- as.numeric(estimates$UCL); estimates$k <- as.numeric(estimates$k)
  estimates[,2:4] <- round(estimates[,2:4], 3)

  # Egger's test on metagen
  eggers <- tryCatch(meta::metabias(metapp_metagen_mean, k.min=3, method="linreg"), error=function(e) NULL)
  eggers_table <- if (!is.null(eggers)) { d <- data.frame(Intercept=eggers$estimate[1], Tau2=eggers$tau, t=eggers$statistic, p=ifelse(round(eggers$p.value,3)==0, "< .001", round(eggers$p.value,3))); row.names(d) <- NULL; d } else data.frame(Intercept=NA, Tau2=NA, t=NA, p=NA)

  res <- list(sample=sample, estimates=estimates, heterogeneity=het, eggers=eggers_table)
  jsonlite::write_json(res, file.path(outdir, "results.json"), auto_unbox=TRUE, na="null")

  # ===== plots =====
  # forest estimates
  estp <- estimates; estp$LCL[is.nan(estp$LCL)] <- NA; estp$UCL[is.nan(estp$UCL)] <- NA
  estp$Model <- factor(estp$Model, levels=rev(estp$Model))
  png(file.path(outdir,"plot_forest_estimates.png"), width=880, height=660, res=110)
  print(ggplot(estp, aes(x=Model, y=g)) + geom_point(stat="identity") + geom_abline(slope=0, intercept=0, linetype=2) + ylab("Hedges's g") + geom_errorbar(aes(ymin=LCL, ymax=UCL), stat="identity", width=0.2) + theme_bw() + coord_flip() + theme(text=element_text(size=14)))
  dev.off()

  # funnel
  rma_agg <- tryCatch(metafor::rma(yi=esg, vi=vi, data=pp_agg, method="REML"), error=function(e) NULL)
  png(file.path(outdir,"plot_funnel.png"), width=900, height=700, res=110)
  if (!is.null(rma_agg)) metafor::funnel(rma_agg, xlab="Hedges's g", studlab=FALSE, contour=.95, col.contour="light grey") else plot.new()
  dev.off()

  # pcurve
  pppcurve <- data.frame(TE=pp$esg, seTE=pp$se, studlab=pp$reference, n=pp$nhighlow, id=pp$id2, nhighlow=pp$nhighlow)
  pppcurve_first <- pppcurve[!duplicated(pppcurve$id), ]
  pc_plot <- tryCatch(pcurve(pppcurve_first, effect.estimation=FALSE, N=pppcurve_first$nhighlow, dmin=0, dmax=1), error=function(e) NULL)
  png(file.path(outdir,"plot_pcurve.png"), width=770, height=605, res=110)
  if (!is.null(pc_plot)) {
    pd <- pc_plot$PlotData; names(pd) <- c("bin","observed","power33","flat"); pd$bin <- factor(pd$bin, levels=pd$bin)
    print(ggplot(pd, aes(x=bin, y=observed)) + geom_col(fill="#4472c4", width=0.7) + geom_hline(yintercept=20, linetype="dashed", color="red", size=1) + geom_line(aes(x=as.numeric(bin), y=power33, group=1), color="#70ad47", size=1) + geom_point(aes(x=as.numeric(bin), y=power33), color="#70ad47") + scale_y_continuous(limits=c(0,60)) + labs(x="p-value", y="Percentage of p-values") + theme_bw())
  } else plot.new()
  dev.off()

  # zcurve
  pppcurve$z <- abs(pppcurve$TE / pppcurve$seTE)
  pppcurve_first_z <- pppcurve[!duplicated(pppcurve$id), ]
  zc <- tryCatch(zcurve(pppcurve_first_z$z, bootstrap=0), error=function(e) NULL)
  png(file.path(outdir,"plot_zcurve.png"), width=900, height=700, res=110)
  if (!is.null(zc)) plot.zcurve(zc, annotation=TRUE, main="") else plot.new()
  dev.off()

  # violin
  ppv <- pp[!is.na(pp$zg), ]
  gm <- mean(ppv$esg); gsd <- sd(ppv$esg)
  png(file.path(outdir,"plot_violin.png"), width=660, height=550, res=110)
  print(ggplot(data=ppv, aes(x=1, y=esg)) + xlab("") + ylab("Effect size") + geom_violin(fill=rgb(100/255,180/255,1,.5)) + theme_bw() + scale_y_continuous(name="Hedges's g", sec.axis=sec_axis(~./gsd-gm, name="Standardized g")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + geom_boxplot(width=.25) + theme(text=element_text(size=10)))
  dev.off()

  # forest studies (robu) - big
  ppf <- pp; ppf$reference_studyno <- paste(ppf$reference, ", Study ", ppf$study_no, sep="")
  rve_f <- tryCatch(robumeta::robu(esg~1, data=ppf, studynum=reference_studyno, var.eff.size=se, small=FALSE), error=function(e) NULL)
  png(file.path(outdir,"plot_forest_studies.png"), width=1100, height=min(max(500, 25*nrow(ppf)), 9000), res=110)
  if (!is.null(rve_f)) robumeta::forest.robu(rve_f, es.lab="dv", study.lab="reference_studyno", "Effect size"=esg) else plot.new()
  dev.off()

  list(ok=TRUE, k=nrow(pp))
}