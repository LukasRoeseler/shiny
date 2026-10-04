parse_ts <- function(s){
  s <- tolower(trimws(as.character(s)))
  fam <- substr(s, 1, 1)
  m <- regmatches(s, regexpr("=([-0-9.eE+]+)\\s*$", s))
  val <- if (length(m) && nzchar(m)) as.numeric(sub("=","",m)) else NA
  list(fam=fam, val=val)
}
effect_r <- function(fam, tv, n){
  if (is.na(tv) || is.na(n) || n<=0) return(NA)
  if (fam=='t') return(tv/sqrt(tv^2 + n - 2))
  if (fam=='f') return(sqrt(tv)/sqrt(tv + n - 2))
  if (fam=='c'){ d <- 2*sqrt(tv/(n-tv)); return(d/sqrt(d^2+4)) }
  if (fam=='r') return(tv)
  if (fam=='z') return(tv/sqrt(tv^2 + n))
  NA
}

tb_run <- function(data_path="toolbox_data_JCP.csv", outdir="."){
  raw <- read.csv(data_path, stringsAsFactors=FALSE, sep=";", check.names=FALSE)
  if (!("test.statistic" %in% colnames(raw))) raw <- read.csv(data_path, stringsAsFactors=FALSE, check.names=FALSE)
  if (nrow(raw)==0) return(list(ok=FALSE, msg="empty data"))
  ts <- raw$test.statistic; n <- as.numeric(raw$sample.size)
  # compute effect size, se, z, p
  TE <- rep(NA, length(ts)); seTE <- rep(NA, length(ts)); z <- rep(NA, length(ts))
  for (i in seq_along(ts)){
    if (is.na(ts[i]) || trimws(as.character(ts[i]))=="" ) next
    p <- parse_ts(ts[i])
    if (is.na(p$val) || is.na(n[i])) next
    r <- effect_r(p$fam, p$val, n[i])
    if (is.na(r)) next
    TE[i] <- r
    seTE[i] <- sqrt((1 - abs(r)^2)/(n[i]-2))
    if (p$fam=='z') z[i] <- p$val else z[i] <- abs(r/seTE[i])
  }
  ok <- !is.na(TE) & !is.na(seTE) & !is.na(z)
  dat <- data.frame(TE=TE[ok], seTE=seTE[ok], N=n[ok], z=z[ok],
                    studlab=if ("study.label" %in% colnames(raw)) raw$study.label[ok] else paste0("Study ",seq_len(sum(ok))),
                    stringsAsFactors=FALSE)
  k <- nrow(dat)

  # p-curve
  pc <- tryCatch(pcurve(data.frame(TE=dat$TE, seTE=dat$seTE, studlab=dat$studlab, n=dat$N, id=seq_len(k), nhighlow=dat$N),
                        effect.estimation=TRUE, N=dat$N, dmin=0, dmax=1), error=function(e) NULL)
  pc_d <- if (!is.null(pc)) tryCatch(pc$dEstimate, error=function(e) NA) else NA
  pc_k <- if (!is.null(pc)) tryCatch(pc$kAnalyzed, error=function(e) k) else k

  # z-curve
  zc <- tryCatch(zcurve(dat$z, bootstrap=0), error=function(e) NULL)
  zc_coef <- if (!is.null(zc)) as.list(zc$coefficients) else list(ERR=NA, EDR=NA, Z0=NA)
  zc_n <- if (!is.null(zc)) zc$N_obs else k
  zc_nsig <- if (!is.null(zc)) zc$N_sig else sum(dat$z > 1.96, na.rm=TRUE)

  res <- list(k=k, pcurve_d=pc_d, pcurve_k=pc_k, zc_coef=zc_coef, zc_n=zc_n, zc_nsig=zc_nsig)
  jsonlite::write_json(res, file.path(outdir, "tb_results.json"), auto_unbox=TRUE, na="null")

  # p-curve plot
  png(file.path(outdir,"tb_pcurve.png"), width=770, height=605, res=110)
  if (!is.null(pc)){
    pd <- pc$PlotData; names(pd) <- c("bin","observed","power33","flat"); pd$bin <- factor(pd$bin, levels=pd$bin)
    print(ggplot(pd, aes(x=bin, y=observed)) + geom_col(fill="#4472c4", width=0.7) +
          geom_hline(yintercept=20, linetype="dashed", color="red", size=1) +
          geom_line(aes(x=as.numeric(bin), y=power33, group=1), color="#70ad47", size=1) +
          geom_point(aes(x=as.numeric(bin), y=power33), color="#70ad47") +
          scale_y_continuous(limits=c(0, max(c(60, max(pd$observed, na.rm=TRUE)+5)))) +
          labs(x="p-value", y="Percentage of p-values") + theme_bw())
  } else plot.new()
  dev.off()

  # z-curve plot
  png(file.path(outdir,"tb_zcurve.png"), width=900, height=700, res=110)
  if (!is.null(zc)) plot.zcurve(zc, annotation=TRUE, main="") else plot.new()
  dev.off()

  list(ok=TRUE, k=k)
}