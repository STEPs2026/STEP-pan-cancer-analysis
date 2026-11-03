# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
rm(list=ls())
options(stringsAsFactors=FALSE, timeout=900)
set.seed(2026)

TARGET_REGEX <- "^Kinget_2024_KIRC_"

PRIMARY_SCORE <- "STEP_UP_resid_CAF3_z"
RAW_UP_SCORE <- "STEP_UP_z"
CAF3_SCORE <- "CAF3_z"
DOWN_SCORE <- "STEP_DOWN_z"

PRIMARY_LABEL <- "STEP_UP CAF3-adjusted"
RAW_UP_LABEL <- "STEP_UP raw"
CAF3_LABEL <- "CAF3"
DOWN_LABEL <- "STEP_DOWN raw"

COLOR_HIGH <- "#D73027"
COLOR_LOW <- "#4575B4"
COLOR_NONRESP <- "#D73027"
COLOR_RESP <- "#4575B4"

ROC_COLORS <- c(
  "STEP_UP raw"="#D73027",
  "CAF3"="#4DAF4A",
  "STEP_UP CAF3-adjusted"="#377EB8",
  "STEP_DOWN raw"="#984EA3"
)

RESP4_COLORS <- c(
  "CR"="#2166AC",
  "PR"="#67A9CF",
  "SD"="#FDAE61",
  "PD"="#D73027"
)

MIN_RESPONSE_PER_GROUP <- 3
MIN_SURVIVAL_N <- 10
MIN_SURVIVAL_EVENTS <- 5

KINGET_CORRECT_EVENT_CODING <- TRUE

SUMMARY_WIDTH <- 22
SUMMARY_HEIGHT <- 11.5
DPI <- 300

required_pkgs <- c(
  "ggplot2","dplyr","tidyr","survival","survminer",
  "pROC","patchwork","gridExtra","scales"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    install.packages(p, repos="https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(pROC)
  library(patchwork)
  library(gridExtra)
  library(scales)
})

find_analysis_dir <- function(start_dir) {
  if (grepl("^PRECOG_STEP_FORMAL_ANALYSIS_v1$", basename(start_dir), ignore.case=TRUE)) {
    return(start_dir)
  }
  dirs <- list.dirs(start_dir, recursive=FALSE, full.names=TRUE)
  hit <- dirs[grepl("^PRECOG_STEP_FORMAL_ANALYSIS_v1$", basename(dirs), ignore.case=TRUE)]
  if (length(hit)>=1) return(hit[1])
  NULL
}

analysis_dir <- find_analysis_dir(getwd())

if (is.null(analysis_dir)) {
  if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
    selected <- rstudioapi::selectDirectory(
      caption="请选择 PRECOG_STEP_FORMAL_ANALYSIS_v1 文件夹"
    )
    if (!is.null(selected) && selected!="") analysis_dir <- selected
  }
}

if (is.null(analysis_dir) || !dir.exists(analysis_dir)) {
  stop("Could not locate PRECOG_STEP_FORMAL_ANALYSIS_v1.")
}

analysis_dir <- normalizePath(analysis_dir, mustWork=TRUE)

score_dir <- file.path(
  analysis_dir,
  "01_MAIN_HIGH_COVERAGE",
  "01_SCORES"
)

if (!dir.exists(score_dir)) stop("MAIN patient-score folder not found:\n", score_dir)

out_root <- file.path(
  analysis_dir,
  "07_PRECOG_VERSION_B_SUMMARY_TEST_KINGET_KIRC_FROM_V1"
)

for (d in c("00_QC","01_PANELS","02_SUMMARY_FIGURES","03_SOURCE_DATA")) {
  dir.create(file.path(out_root,d), recursive=TRUE, showWarnings=FALSE)
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+","_",x)
  substr(x,1,220)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_z <- function(x) {
  x <- safe_numeric(x)
  out <- rep(NA_real_,length(x))
  ok <- is.finite(x)
  if (sum(ok)<2) return(out)
  s <- sd(x[ok],na.rm=TRUE)
  if (!is.finite(s) || s<=1e-12) return(out)
  out[ok] <- (x[ok]-mean(x[ok],na.rm=TRUE))/s
  out
}

p_fmt <- function(p) {
  if (length(p) == 0) return(character(0))

  p <- suppressWarnings(as.numeric(p))

  vapply(
    p,
    function(x) {
      if (is.na(x) || !is.finite(x)) return("p = NA")
      if (x < 0.001) return("p < 0.001")
      paste0("p = ", sprintf("%.3f", x))
    },
    character(1)
  )
}

make_group <- function(x) {
  x <- safe_numeric(x)
  out <- rep(NA_character_,length(x))
  ok <- is.finite(x)
  if (sum(ok)<4) return(factor(out,levels=c("Low","High")))
  cut <- median(x[ok],na.rm=TRUE)
  out[ok & x<=cut] <- "Low"
  out[ok & x>cut] <- "High"
  factor(out,levels=c("Low","High"))
}

theme_pub <- function(base_size=10.5) {
  theme_classic(base_size=base_size)+
    theme(
      plot.title=element_text(size=9.5,face="bold",hjust=.5),
      plot.subtitle=element_text(size=8.2,hjust=.5),
      axis.title=element_text(size=9.5),
      axis.text=element_text(size=8.7),
      legend.title=element_text(size=8.5),
      legend.text=element_text(size=8.2)
    )
}

save_plot <- function(p,path_no_ext,width=6,height=5) {
  ggsave(paste0(path_no_ext,".pdf"),p,width=width,height=height,limitsize=FALSE)
  ggsave(paste0(path_no_ext,".png"),p,width=width,height=height,dpi=DPI,limitsize=FALSE)
  print(p)
}

placeholder_panel <- function(title,msg) {
  ggplot()+
    annotate("text",x=.5,y=.55,label=msg,size=4,lineheight=1.1)+
    xlim(0,1)+ylim(0,1)+labs(title=title)+theme_void()+
    theme(
      plot.title=element_text(size=9.5,face="bold",hjust=.5),
      panel.border=element_rect(fill=NA,color="grey75",linewidth=.4)
    )
}

pretty_dataset <- function(dataset) {
  tx <- "ICI"
  if (grepl("CTLA4_PD1",dataset,ignore.case=TRUE)) tx <- "CTLA4 + PD1"
  else if (grepl("_PD1_",dataset,ignore.case=TRUE)) tx <- "PD1"
  else if (grepl("_PDL1_",dataset,ignore.case=TRUE)) tx <- "PD-L1"
  paste0("Kinget 2024 | KIRC | ",tx," | PRECOG")
}

score_files <- list.files(score_dir,pattern="_patient_scores\\.csv$",full.names=TRUE)

score_index <- list()

for (f in score_files) {
  tmp <- tryCatch(
    read.csv(f,stringsAsFactors=FALSE,check.names=FALSE,nrows=5),
    error=function(e) NULL
  )
  if (is.null(tmp) || !"Dataset"%in%colnames(tmp) || nrow(tmp)==0) next
  ds <- unique(tmp$Dataset)
  ds <- ds[!is.na(ds) & ds!=""]
  if (length(ds)==1) score_index[[ds]] <- f
}

target_datasets <- names(score_index)[
  grepl(TARGET_REGEX,names(score_index),ignore.case=TRUE)
]

if (length(target_datasets)==0) {
  stop("No Kinget KIRC PRECOG patient-score tables found.")
}

parse_response4 <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_,length(y))
  out[y %in% c("CR","COMPLETE RESPONSE")] <- "CR"
  out[y %in% c("PR","PARTIAL RESPONSE")] <- "PR"
  out[y %in% c("SD","STABLE DISEASE")] <- "SD"
  out[y %in% c("PD","PROGRESSIVE DISEASE")] <- "PD"
  factor(out,levels=c("CR","PR","SD","PD"),ordered=TRUE)
}

detect_response4 <- function(d) {
  candidates <- grep(
    "best.*response|response.*category|recist|bestresponse|best_response|(^|_)bor($|_)",
    colnames(d),ignore.case=TRUE,value=TRUE
  )
  empty <- factor(rep(NA_character_,nrow(d)),levels=c("CR","PR","SD","PD"),ordered=TRUE)
  if (length(candidates)==0) return(list(column=NA_character_,value=empty))
  n_valid <- vapply(candidates,function(nm)sum(!is.na(parse_response4(d[[nm]]))),numeric(1))
  if (max(n_valid)==0) return(list(column=NA_character_,value=empty))
  best <- candidates[which.max(n_valid)]
  list(column=best,value=parse_response4(d[[best]]))
}

reconstruct_summary_events <- function(d,dataset) {
  d$Summary_PFS_event <- if ("PFS_event"%in%names(d)) safe_numeric(d$PFS_event) else NA_real_
  d$Summary_OS_event <- if ("OS_event"%in%names(d)) safe_numeric(d$OS_event) else NA_real_

  pfs_rule <- "PRECOG v1.1 PFS_event"
  os_rule <- "PRECOG v1.1 OS_event"

  if (KINGET_CORRECT_EVENT_CODING && grepl("Kinget_2024_KIRC",dataset,ignore.case=TRUE)) {
    if ("PFS_Censor"%in%names(d)) {
      raw <- safe_numeric(d$PFS_Censor)
      valid <- !is.na(raw) & raw%in%c(0,1)
      if (sum(valid)>0) {
        d$Summary_PFS_event[valid] <- raw[valid]
        pfs_rule <- "Kinget override: raw PFS_Censor used directly as event (1=event)"
      }
    }
    if ("OS_Censor"%in%names(d)) {
      raw <- safe_numeric(d$OS_Censor)
      valid <- !is.na(raw) & raw%in%c(0,1)
      if (sum(valid)>0) {
        d$Summary_OS_event[valid] <- raw[valid]
        os_rule <- "Kinget override: raw OS_Censor used directly as event (1=event)"
      }
    }
  }

  qc <- data.frame(
    Dataset=dataset,
    N=nrow(d),
    Old_v11_PFS_events=if("PFS_event"%in%names(d))sum(safe_numeric(d$PFS_event)==1,na.rm=TRUE) else NA,
    Summary_PFS_events=sum(d$Summary_PFS_event==1,na.rm=TRUE),
    PFS_rule=pfs_rule,
    Old_v11_OS_events=if("OS_event"%in%names(d))sum(safe_numeric(d$OS_event)==1,na.rm=TRUE) else NA,
    Summary_OS_events=sum(d$Summary_OS_event==1,na.rm=TRUE),
    OS_rule=os_rule
  )

  list(data=d,qc=qc)
}

fit_response_model <- function(d,score_var,score_label) {
  use <- d[complete.cases(d$NonResponder_binary,d[[score_var]]),,drop=FALSE]
  if (
    nrow(use)<2*MIN_RESPONSE_PER_GROUP ||
    length(unique(use$NonResponder_binary))<2
  ) return(NULL)

  md <- data.frame(
    y=safe_numeric(use$NonResponder_binary),
    score=safe_z(use[[score_var]])
  )
  md <- md[complete.cases(md),,drop=FALSE]

  fit <- tryCatch(
    glm(y~score,family=binomial(),data=md),
    error=function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  co <- summary(fit)$coefficients["score",]
  b <- co["Estimate"]
  se <- co["Std. Error"]
  p <- co["Pr(>|z|)"]

  data.frame(
    Endpoint="Response",
    Score=score_label,
    Score_variable=score_var,
    N=nrow(md),
    Events=sum(md$y==1),
    Effect=exp(b),
    CI_low=exp(b-1.96*se),
    CI_high=exp(b+1.96*se),
    P=p,
    Effect_name="OR",
    Direction="OR>1 = more non-response"
  )
}

fit_survival_model <- function(d,endpoint,score_var,score_label) {
  time_var <- if(endpoint=="OS") "OS_time" else "PFS_time"
  event_var <- if(endpoint=="OS") "Summary_OS_event" else "Summary_PFS_event"

  if (!all(c(time_var,event_var,score_var)%in%names(d))) return(NULL)

  use <- d[
    complete.cases(d[[time_var]],d[[event_var]],d[[score_var]]) &
      safe_numeric(d[[time_var]])>0,
    ,drop=FALSE
  ]
  events <- sum(safe_numeric(use[[event_var]])==1,na.rm=TRUE)

  if (nrow(use)<MIN_SURVIVAL_N || events<MIN_SURVIVAL_EVENTS) return(NULL)

  md <- data.frame(
    time=safe_numeric(use[[time_var]]),
    event=safe_numeric(use[[event_var]]),
    score=safe_z(use[[score_var]])
  )

  fit <- tryCatch(
    coxph(Surv(time,event)~score,data=md),
    error=function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  co <- summary(fit)$coef["score",]
  b <- co["coef"]
  se <- co["se(coef)"]
  p <- co["Pr(>|z|)"]

  data.frame(
    Endpoint=endpoint,
    Score=score_label,
    Score_variable=score_var,
    N=nrow(md),
    Events=events,
    Effect=exp(b),
    CI_low=exp(b-1.96*se),
    CI_high=exp(b+1.96*se),
    P=p,
    Effect_name="HR",
    Direction="HR>1 = worse survival"
  )
}

make_response_box <- function(d,score_var,score_label,panel_title) {
  use <- d[complete.cases(d$Responder_binary,d[[score_var]]),,drop=FALSE]
  if (nrow(use)<2*MIN_RESPONSE_PER_GROUP) {
    return(placeholder_panel(panel_title,"Response data unavailable"))
  }

  use$Response <- factor(
    ifelse(safe_numeric(use$Responder_binary)==1,"Responder","Non-responder"),
    levels=c("Responder","Non-responder")
  )

  if (length(unique(use$Response))<2) {
    return(placeholder_panel(panel_title,"Only one response class available"))
  }

  wilcox_p <- tryCatch(
    wilcox.test(use[[score_var]]~use$Response)$p.value,
    error=function(e) NA_real_
  )

  mr <- fit_response_model(d,score_var,score_label)

  model_text <- ""
  if (!is.null(mr)) {
    model_text <- paste0(
      "\nContinuous OR=",sprintf("%.2f",mr$Effect),
      " [",sprintf("%.2f",mr$CI_low),", ",sprintf("%.2f",mr$CI_high),"], ",
      p_fmt(mr$P)
    )
  }

  nn <- table(use$Response)
  labs_x <- setNames(
    paste0(names(nn),"\n(n=",as.integer(nn),")"),
    names(nn)
  )

  ggplot(use,aes(Response,.data[[score_var]],fill=Response))+
    geom_violin(trim=FALSE,alpha=.28,color="grey35")+
    geom_boxplot(width=.20,outlier.shape=NA,alpha=.82)+
    geom_jitter(width=.09,size=1.4,alpha=.68)+
    scale_fill_manual(values=c("Responder"=COLOR_RESP,"Non-responder"=COLOR_NONRESP))+
    scale_x_discrete(labels=labs_x)+
    labs(
      title=panel_title,
      subtitle=paste0("Wilcoxon ",p_fmt(wilcox_p),model_text),
      x=NULL,
      y=paste0(score_label,"\n(within-cohort z-score)")
    )+
    theme_pub()+
    theme(legend.position="none")
}

nonresponse_or_high_vs_low <- function(group,nonresponse) {
  use <- data.frame(group=group,y=safe_numeric(nonresponse))
  use <- use[complete.cases(use),,drop=FALSE]
  use$group <- factor(use$group,levels=c("Low","High"))
  tb <- table(use$group,use$y)

  if (!all(dim(tb)==c(2,2))) {
    return(list(p=NA,OR=NA,CI_low=NA,CI_high=NA))
  }

  fp <- tryCatch(fisher.test(tb)$p.value,error=function(e)NA_real_)

  a <- tb["High","1"]; b <- tb["High","0"]
  c0 <- tb["Low","1"]; d0 <- tb["Low","0"]

  OR <- ((a+.5)/(b+.5))/((c0+.5)/(d0+.5))
  se <- sqrt(1/(a+.5)+1/(b+.5)+1/(c0+.5)+1/(d0+.5))

  list(
    p=fp,
    OR=OR,
    CI_low=exp(log(OR)-1.96*se),
    CI_high=exp(log(OR)+1.96*se)
  )
}

make_response_rate <- function(d,score_var,panel_title) {
  use <- d[complete.cases(d$Responder_binary,d[[score_var]]),,drop=FALSE]
  if (nrow(use)<2*MIN_RESPONSE_PER_GROUP) {
    return(placeholder_panel(panel_title,"Response data unavailable"))
  }

  use$ScoreGroup <- make_group(use[[score_var]])
  use <- use[complete.cases(use$ScoreGroup,use$Responder_binary),,drop=FALSE]

  rate_df <- use %>%
    group_by(ScoreGroup) %>%
    summarise(
      N=n(),
      Responders=sum(safe_numeric(Responder_binary)==1),
      Response_rate=Responders/N,
      .groups="drop"
    )

  st <- nonresponse_or_high_vs_low(use$ScoreGroup,use$NonResponder_binary)

  subtitle <- paste0(
    "Fisher ",p_fmt(st$p),
    "; Non-response OR High vs Low=",
    ifelse(is.finite(st$OR),sprintf("%.2f",st$OR),"NA")
  )

  if (is.finite(st$CI_low)) {
    subtitle <- paste0(
      subtitle," [",sprintf("%.2f",st$CI_low),", ",sprintf("%.2f",st$CI_high),"]"
    )
  }

  ggplot(rate_df,aes(ScoreGroup,Response_rate*100,fill=ScoreGroup))+
    geom_col(width=.62)+
    geom_text(
      aes(label=paste0(Responders,"/",N,"\n",sprintf("%.1f%%",Response_rate*100))),
      vjust=-.18,size=3.5
    )+
    scale_fill_manual(values=c("Low"=COLOR_LOW,"High"=COLOR_HIGH))+
    coord_cartesian(ylim=c(0,110),clip="off")+
    labs(
      title=panel_title,
      subtitle=subtitle,
      x=paste0(PRIMARY_LABEL," group"),
      y="ICI response rate (%)"
    )+
    theme_pub()+
    theme(legend.position="none")
}

make_response_composition <- function(d,score_var,panel_title) {
  use <- d[complete.cases(d$Responder_binary,d[[score_var]]),,drop=FALSE]
  if (nrow(use)<2*MIN_RESPONSE_PER_GROUP) {
    return(placeholder_panel(panel_title,"Response data unavailable"))
  }

  use$ScoreGroup <- make_group(use[[score_var]])
  use$Response <- factor(
    ifelse(safe_numeric(use$Responder_binary)==1,"Responder","Non-responder"),
    levels=c("Responder","Non-responder")
  )

  use <- use[complete.cases(use$ScoreGroup,use$Response),,drop=FALSE]
  tb <- table(use$ScoreGroup,use$Response)

  if (nrow(tb)<2 || ncol(tb)<2) {
    return(placeholder_panel(panel_title,"Both response classes are required"))
  }

  fp <- tryCatch(fisher.test(tb)$p.value,error=function(e)NA_real_)

  pd <- as.data.frame(tb)
  colnames(pd) <- c("ScoreGroup","Response","N")
  pd <- pd %>%
    group_by(ScoreGroup) %>%
    mutate(
      Total=sum(N),
      Percent=N/Total*100,
      Label=ifelse(N>0,paste0(N,"\n",sprintf("%.1f%%",Percent)),"")
    ) %>%
    ungroup()

  ggplot(pd,aes(ScoreGroup,Percent,fill=Response))+
    geom_col(width=.68)+
    geom_text(aes(label=Label),position=position_stack(vjust=.5),size=3.3)+
    scale_fill_manual(values=c("Responder"=COLOR_RESP,"Non-responder"=COLOR_NONRESP))+
    labs(
      title=panel_title,
      subtitle=paste0("Fisher ",p_fmt(fp)),
      x=paste0(PRIMARY_LABEL," group"),
      y="Patients (%)",
      fill=NULL
    )+
    theme_pub()
}

make_multi_roc <- function(d,panel_title) {
  score_map <- c(
    "STEP_UP raw"=RAW_UP_SCORE,
    "CAF3"=CAF3_SCORE,
    "STEP_UP CAF3-adjusted"=PRIMARY_SCORE,
    "STEP_DOWN raw"=DOWN_SCORE
  )

  rocs <- list()
  auc_txt <- c()

  for (nm in names(score_map)) {
    sc <- score_map[[nm]]
    use <- d[complete.cases(d$NonResponder_binary,d[[sc]]),,drop=FALSE]
    if (
      nrow(use)<2*MIN_RESPONSE_PER_GROUP ||
      length(unique(use$NonResponder_binary))<2
    ) next

    rr <- tryCatch(
      pROC::roc(
        response=safe_numeric(use$NonResponder_binary),
        predictor=safe_numeric(use[[sc]]),
        levels=c(0,1),
        direction="<",
        quiet=TRUE
      ),
      error=function(e) NULL
    )

    if (is.null(rr)) next
    rocs[[nm]] <- rr
    auc_txt <- c(
      auc_txt,
      paste0(nm,"=",sprintf("%.3f",as.numeric(pROC::auc(rr))))
    )
  }

  if (length(rocs)==0) return(placeholder_panel(panel_title,"ROC unavailable"))

  pROC::ggroc(rocs,legacy.axes=TRUE,linewidth=.85)+
    geom_abline(slope=1,intercept=0,linetype=2,linewidth=.5)+
    scale_color_manual(values=ROC_COLORS[names(rocs)])+
    coord_equal()+
    labs(
      title=panel_title,
      subtitle=paste(auc_txt,collapse="; "),
      x="1 - Specificity",
      y="Sensitivity",
      color=NULL
    )+
    theme_pub()
}

make_response4_panel <- function(d,response4,panel_title) {
  use <- data.frame(Score=safe_numeric(d[[PRIMARY_SCORE]]),Response4=response4)
  use <- use[complete.cases(use),,drop=FALSE]

  if (nrow(use)<10 || length(unique(use$Response4))<3) return(NULL)

  kw <- tryCatch(kruskal.test(Score~Response4,data=use)$p.value,error=function(e)NA_real_)
  nn <- table(factor(use$Response4,levels=c("CR","PR","SD","PD")))
  labs_x <- setNames(paste0(names(nn),"\n(n=",as.integer(nn),")"),names(nn))

  ggplot(use,aes(Response4,Score,fill=Response4))+
    geom_violin(trim=FALSE,alpha=.32)+
    geom_boxplot(width=.20,outlier.shape=NA)+
    geom_jitter(width=.08,size=1.3,alpha=.65)+
    scale_fill_manual(values=RESP4_COLORS,drop=FALSE)+
    scale_x_discrete(labels=labs_x,drop=FALSE)+
    labs(
      title=panel_title,
      subtitle=paste0("Kruskal-Wallis ",p_fmt(kw)),
      x="Best response",
      y=PRIMARY_LABEL
    )+
    theme_pub()+
    theme(legend.position="none")
}

make_disease_control_panel <- function(d,response4,panel_title) {
  dc <- ifelse(
    response4%in%c("CR","PR","SD"),1,
    ifelse(response4=="PD",0,NA)
  )

  temp <- d
  temp$DC <- dc
  temp$Responder_binary <- temp$DC
  temp$NonResponder_binary <- ifelse(is.na(temp$DC),NA,1-temp$DC)

  make_response_box(
    temp,
    PRIMARY_SCORE,
    PRIMARY_LABEL,
    panel_title
  )
}

endpoint_usable <- function(d,endpoint) {
  time_var <- if(endpoint=="OS") "OS_time" else "PFS_time"
  event_var <- if(endpoint=="OS") "Summary_OS_event" else "Summary_PFS_event"

  if (!all(c(time_var,event_var)%in%names(d))) return(FALSE)

  use <- complete.cases(d[[time_var]],d[[event_var]]) &
    safe_numeric(d[[time_var]])>0

  n <- sum(use)
  ev <- sum(safe_numeric(d[[event_var]][use])==1,na.rm=TRUE)

  n>=MIN_SURVIVAL_N && ev>=MIN_SURVIVAL_EVENTS
}

make_km_panel <- function(d,endpoint,score_var,score_label,panel_title) {
  time_var <- if(endpoint=="OS") "OS_time" else "PFS_time"
  event_var <- if(endpoint=="OS") "Summary_OS_event" else "Summary_PFS_event"

  if (!all(c(time_var,event_var,score_var)%in%names(d))) return(NULL)

  use <- d[
    complete.cases(d[[time_var]],d[[event_var]],d[[score_var]]) &
      safe_numeric(d[[time_var]])>0,
    ,drop=FALSE
  ]

  events <- sum(safe_numeric(use[[event_var]])==1,na.rm=TRUE)

  if (nrow(use)<MIN_SURVIVAL_N || events<MIN_SURVIVAL_EVENTS) return(NULL)

  use$ScoreGroup <- make_group(use[[score_var]])
  use <- use[complete.cases(use$ScoreGroup),,drop=FALSE]

  md <- data.frame(
    time=safe_numeric(use[[time_var]]),
    event=safe_numeric(use[[event_var]]),
    group=factor(use$ScoreGroup,levels=c("Low","High"))
  )

  fit <- survfit(Surv(time,event)~group,data=md)

  lp <- tryCatch({
    sdif <- survdiff(Surv(time,event)~group,data=md)
    1-pchisq(sdif$chisq,df=length(sdif$n)-1)
  },error=function(e)NA_real_)

  cr <- fit_survival_model(d,endpoint,score_var,score_label)

  ctxt <- ""
  if (!is.null(cr)) {
    ctxt <- paste0(
      "Continuous HR=",sprintf("%.2f",cr$Effect),
      " [",sprintf("%.2f",cr$CI_low),", ",sprintf("%.2f",cr$CI_high),"], ",
      p_fmt(cr$P)
    )
  }

  g <- survminer::ggsurvplot(
    fit,
    data=md,
    risk.table=TRUE,
    pval=paste0("Log-rank ",p_fmt(lp)),
    conf.int=FALSE,
    censor=TRUE,
    palette=c(COLOR_LOW,COLOR_HIGH),
    legend.title="STEP group",
    legend.labs=c("Low","High"),
    xlab=paste0(endpoint," time (days)"),
    ylab=paste0(endpoint," probability"),
    title=panel_title,
    risk.table.height=.25,
    ggtheme=theme_classic(base_size=9.5),
    tables.theme=survminer::theme_cleantable(base_size=8.2)
  )

  g$plot <- g$plot+
    annotate(
      "text",x=Inf,y=.16,hjust=1.05,vjust=0,
      label=paste0(ctxt,"\nN=",nrow(md),"; events=",events),
      size=2.8
    )+
    theme(
      plot.title=element_text(size=9.5,face="bold",hjust=.5),
      legend.title=element_text(size=7.5),
      legend.text=element_text(size=7.5)
    )

  grob <- gridExtra::arrangeGrob(
    g$plot,g$table,ncol=1,heights=c(3,1)
  )

  patchwork::wrap_elements(full=grob)
}

make_effect_forest <- function(d,preferred_endpoint,panel_title) {
  score_map <- c(
    "STEP_UP raw"=RAW_UP_SCORE,
    "CAF3"=CAF3_SCORE,
    "STEP_UP CAF3-adjusted"=PRIMARY_SCORE,
    "STEP_DOWN raw"=DOWN_SCORE
  )

  rows <- list()
  k <- 1

  for (nm in names(score_map)) {
    sc <- score_map[[nm]]

    rr <- fit_response_model(d,sc,nm)
    if (!is.null(rr)) {
      rows[[k]] <- rr
      k <- k+1
    }

    sr <- fit_survival_model(d,preferred_endpoint,sc,nm)
    if (!is.null(sr)) {
      rows[[k]] <- sr
      k <- k+1
    }
  }

  res <- if(length(rows)) bind_rows(rows) else data.frame()

  if (!nrow(res)) return(placeholder_panel(panel_title,"Effect estimates unavailable"))

  res$Label <- paste0(res$Endpoint," | ",res$Score)
  res$Label <- factor(res$Label,levels=rev(unique(res$Label)))
  res$Txt <- paste0(
    res$Effect_name,"=",sprintf("%.2f",res$Effect),"\n",p_fmt(res$P)
  )

  ggplot(res,aes(Label,Effect,ymin=CI_low,ymax=CI_high,shape=Endpoint))+
    geom_hline(yintercept=1,linetype=2,linewidth=.5)+
    geom_errorbar(width=.12)+
    geom_point(size=2.3)+
    geom_text(aes(label=Txt),hjust=-.08,size=2.6)+
    scale_y_log10()+
    coord_flip(clip="off")+
    labs(
      title=panel_title,
      subtitle=paste0("Response OR>1 worse; ",preferred_endpoint," HR>1 worse"),
      x=NULL,
      y="OR / HR (log scale)",
      shape=NULL
    )+
    theme_pub()+
    theme(
      axis.text.y=element_text(size=6.8),
      plot.margin=margin(5.5,45,5.5,5.5)
    )
}

existing_cols <- function(d, cols) {
  unique(cols[cols %in% names(d)])
}

export_ABE_source_data <- function(d, dataset, preferred_endpoint, detailed_n, source_dir) {
  dir.create(source_dir, recursive=TRUE, showWarnings=FALSE)

  id_cols <- existing_cols(d, c("Dataset","Patient_ID","Sample_ID","geo_accession"))

  if (!is.na(preferred_endpoint)) {
    time_var <- if(preferred_endpoint=="OS") "OS_time" else "PFS_time"
    event_var <- if(preferred_endpoint=="OS") "Summary_OS_event" else "Summary_PFS_event"
    old_event_var <- if(preferred_endpoint=="OS") "OS_event" else "PFS_event"
    raw_time_var <- if(preferred_endpoint=="OS") "OS_Days" else "PFS_Days"
    raw_censor_var <- if(preferred_endpoint=="OS") "OS_Censor" else "PFS_Censor"

    useA <- d[
      complete.cases(d[[time_var]], d[[event_var]], d[[PRIMARY_SCORE]]) &
        safe_numeric(d[[time_var]]) > 0,
      ,drop=FALSE
    ]

    if (nrow(useA)>0) {
      useA$ScoreGroup_exact <- make_group(useA[[PRIMARY_SCORE]])
      cutoffA <- median(safe_numeric(useA[[PRIMARY_SCORE]]), na.rm=TRUE)

      a_cols <- existing_cols(
        useA,
        c(
          id_cols,
          raw_time_var, raw_censor_var,
          time_var, old_event_var, event_var,
          PRIMARY_SCORE,
          "ScoreGroup_exact"
        )
      )

      write.csv(
        useA[,a_cols,drop=FALSE],
        file.path(source_dir, paste0("A_",preferred_endpoint,"_KM_patient_level_source_data.csv")),
        row.names=FALSE
      )

      mdA <- data.frame(
        time=safe_numeric(useA[[time_var]]),
        event=safe_numeric(useA[[event_var]]),
        group=factor(useA$ScoreGroup_exact,levels=c("Low","High"))
      )

      lp <- tryCatch({
        sdif <- survdiff(Surv(time,event)~group,data=mdA)
        1-pchisq(sdif$chisq,df=length(sdif$n)-1)
      }, error=function(e) NA_real_)

      coxA <- fit_survival_model(d, preferred_endpoint, PRIMARY_SCORE, PRIMARY_LABEL)

      event_rule <- if(
        KINGET_CORRECT_EVENT_CODING &&
          grepl("Kinget_2024_KIRC",dataset,ignore.case=TRUE) &&
          raw_censor_var %in% names(d) &&
          sum(!is.na(safe_numeric(d[[raw_censor_var]])))>0
      ) {
        paste0("Kinget override: raw ",raw_censor_var," used directly as event (1=event)")
      } else {
        paste0("Existing ",event_var," used")
      }

      statA <- data.frame(
        Dataset=dataset,
        Panel="A",
        Endpoint=preferred_endpoint,
        Score_variable=PRIMARY_SCORE,
        Median_cutoff=cutoffA,
        N=nrow(mdA),
        Events=sum(mdA$event==1,na.rm=TRUE),
        Low_N=sum(mdA$group=="Low",na.rm=TRUE),
        High_N=sum(mdA$group=="High",na.rm=TRUE),
        Log_rank_P=lp,
        Continuous_HR=if(!is.null(coxA)) coxA$Effect else NA_real_,
        Continuous_HR_CI_low=if(!is.null(coxA)) coxA$CI_low else NA_real_,
        Continuous_HR_CI_high=if(!is.null(coxA)) coxA$CI_high else NA_real_,
        Continuous_Cox_P=if(!is.null(coxA)) coxA$P else NA_real_,
        Event_rule=event_rule,
        stringsAsFactors=FALSE
      )

      write.csv(
        statA,
        file.path(source_dir, paste0("A_",preferred_endpoint,"_KM_statistics.csv")),
        row.names=FALSE
      )
    }
  }

  useB <- d[complete.cases(d$Responder_binary,d[[PRIMARY_SCORE]]),,drop=FALSE]

  if (nrow(useB)>0) {
    useB$Response_exact <- factor(
      ifelse(safe_numeric(useB$Responder_binary)==1,"Responder","Non-responder"),
      levels=c("Responder","Non-responder")
    )

    b_cols <- existing_cols(
      useB,
      c(
        id_cols,
        "Responder", "Responder_binary", "NonResponder_binary",
        PRIMARY_SCORE,
        "Response_exact"
      )
    )

    write.csv(
      useB[,b_cols,drop=FALSE],
      file.path(source_dir,"B_Response_box_patient_level_source_data.csv"),
      row.names=FALSE
    )

    wilcoxB <- tryCatch(
      wilcox.test(useB[[PRIMARY_SCORE]]~useB$Response_exact)$p.value,
      error=function(e) NA_real_
    )
    rrB <- fit_response_model(d,PRIMARY_SCORE,PRIMARY_LABEL)

    statB <- data.frame(
      Dataset=dataset,
      Panel="B",
      Score_variable=PRIMARY_SCORE,
      N=nrow(useB),
      Responder_N=sum(useB$Response_exact=="Responder",na.rm=TRUE),
      NonResponder_N=sum(useB$Response_exact=="Non-responder",na.rm=TRUE),
      Wilcoxon_P=wilcoxB,
      Continuous_nonresponse_OR=if(!is.null(rrB)) rrB$Effect else NA_real_,
      Continuous_OR_CI_low=if(!is.null(rrB)) rrB$CI_low else NA_real_,
      Continuous_OR_CI_high=if(!is.null(rrB)) rrB$CI_high else NA_real_,
      Continuous_logistic_P=if(!is.null(rrB)) rrB$P else NA_real_,
      stringsAsFactors=FALSE
    )

    write.csv(
      statB,
      file.path(source_dir,"B_Response_box_statistics.csv"),
      row.names=FALSE
    )
  }

  valid_r4 <- !is.na(d$Response4_summary)
  has_detailed_E <- detailed_n>=10 && length(unique(d$Response4_summary[valid_r4]))>=3

  if (has_detailed_E) {
    useE <- d[complete.cases(d[[PRIMARY_SCORE]],d$Response4_summary),,drop=FALSE]
    useE$Response4_exact <- d$Response4_summary[complete.cases(d[[PRIMARY_SCORE]],d$Response4_summary)]

    e_cols <- existing_cols(
      useE,
      c(id_cols,PRIMARY_SCORE,"Response4_summary","Response4_exact")
    )

    write.csv(
      useE[,e_cols,drop=FALSE],
      file.path(source_dir,"E_CR_PR_SD_PD_patient_level_source_data.csv"),
      row.names=FALSE
    )

    kwE <- tryCatch(
      kruskal.test(safe_numeric(useE[[PRIMARY_SCORE]])~useE$Response4_exact)$p.value,
      error=function(e) NA_real_
    )

    summaryE <- as.data.frame(table(useE$Response4_exact))
    colnames(summaryE) <- c("Response4","N")
    summaryE$Percent <- summaryE$N/sum(summaryE$N)*100
    summaryE$Kruskal_Wallis_P <- kwE

    write.csv(
      summaryE,
      file.path(source_dir,"E_CR_PR_SD_PD_summary_and_statistics.csv"),
      row.names=FALSE
    )

  } else {
    useE <- d[complete.cases(d$Responder_binary,d[[PRIMARY_SCORE]]),,drop=FALSE]
    useE$ScoreGroup_exact <- make_group(useE[[PRIMARY_SCORE]])
    useE$Response_exact <- factor(
      ifelse(safe_numeric(useE$Responder_binary)==1,"Responder","Non-responder"),
      levels=c("Responder","Non-responder")
    )
    useE <- useE[complete.cases(useE$ScoreGroup_exact,useE$Response_exact),,drop=FALSE]

    e_cols <- existing_cols(
      useE,
      c(
        id_cols,
        "Responder", "Responder_binary", "NonResponder_binary",
        PRIMARY_SCORE,
        "ScoreGroup_exact", "Response_exact"
      )
    )

    write.csv(
      useE[,e_cols,drop=FALSE],
      file.path(source_dir,"E_Response_composition_patient_level_source_data.csv"),
      row.names=FALSE
    )

    tbE <- table(useE$ScoreGroup_exact,useE$Response_exact)
    fisherE <- tryCatch(fisher.test(tbE)$p.value,error=function(e) NA_real_)

    summaryE <- as.data.frame(tbE)
    colnames(summaryE) <- c("STEP_group","Response","N")
    summaryE <- summaryE %>%
      group_by(STEP_group) %>%
      mutate(
        Total=sum(N),
        Percent=N/Total*100
      ) %>%
      ungroup()
    summaryE$Fisher_P <- fisherE

    write.csv(
      summaryE,
      file.path(source_dir,"E_Response_composition_summary_and_statistics.csv"),
      row.names=FALSE
    )
  }
}

event_qc_list <- list()
availability_list <- list()
model_list <- list()

for (dataset in target_datasets) {

  message("\nBuilding Version-B figure: ",dataset)

  d <- read.csv(
    score_index[[dataset]],
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  required_cols <- c(
    PRIMARY_SCORE,RAW_UP_SCORE,CAF3_SCORE,DOWN_SCORE,
    "Responder_binary","NonResponder_binary"
  )

  missing_cols <- setdiff(required_cols,names(d))

  if (length(missing_cols)) {
    warning(dataset,": missing columns: ",paste(missing_cols,collapse=", "))
    next
  }

  ec <- reconstruct_summary_events(d,dataset)
  d <- ec$data
  event_qc_list[[dataset]] <- ec$qc

  r4 <- detect_response4(d)
  d$Response4_summary <- r4$value
  detailed_n <- sum(!is.na(d$Response4_summary))

  pfs_ok <- endpoint_usable(d,"PFS")
  os_ok <- endpoint_usable(d,"OS")

  preferred_endpoint <- if(pfs_ok) "PFS" else if(os_ok) "OS" else NA_character_
  secondary_endpoint <- if(pfs_ok && os_ok) "OS" else NA_character_

  panel_A <- if(!is.na(preferred_endpoint)) {
    make_km_panel(
      d,preferred_endpoint,PRIMARY_SCORE,PRIMARY_LABEL,
      paste0(PRIMARY_LABEL," – ",preferred_endpoint)
    )
  } else NULL

  if (is.null(panel_A)) {
    panel_A <- placeholder_panel("Survival","No usable PRECOG OS/PFS endpoint")
  }

  panel_B <- make_response_box(
    d,PRIMARY_SCORE,PRIMARY_LABEL,
    paste0(PRIMARY_LABEL," vs ICI response")
  )

  panel_C <- make_response_rate(
    d,PRIMARY_SCORE,
    "ICI response rate by STEP group"
  )

  panel_D <- make_multi_roc(
    d,
    "Prediction of ICI non-response"
  )

  if (
    detailed_n>=10 &&
    length(unique(d$Response4_summary[!is.na(d$Response4_summary)]))>=3
  ) {
    panel_E <- make_response4_panel(
      d,d$Response4_summary,
      "CR / PR / SD / PD"
    )
  } else {
    panel_E <- make_response_composition(
      d,PRIMARY_SCORE,
      "Responder / Non-responder composition"
    )
  }

  if (detailed_n>=10) {
    panel_F <- make_disease_control_panel(
      d,d$Response4_summary,
      "Disease control: CR/PR/SD vs PD"
    )
  } else {
    panel_F <- make_response_box(
      d,RAW_UP_SCORE,RAW_UP_LABEL,
      "Raw STEP_UP vs ICI response"
    )
  }

  panel_G <- make_response_box(
    d,DOWN_SCORE,DOWN_LABEL,
    "STEP_DOWN raw vs ICI response"
  )

  if (!is.na(secondary_endpoint)) {
    panel_H <- make_km_panel(
      d,secondary_endpoint,PRIMARY_SCORE,PRIMARY_LABEL,
      paste0(PRIMARY_LABEL," – ",secondary_endpoint)
    )
    if (is.null(panel_H)) {
      panel_H <- make_effect_forest(d,preferred_endpoint,"Response / survival effect summary")
    }
  } else {
    panel_H <- make_effect_forest(
      d,
      ifelse(is.na(preferred_endpoint),"PFS",preferred_endpoint),
      "Response / survival effect summary"
    )
  }

  ds_safe <- safe_filename(dataset)
  source_dir <- file.path(out_root,"03_SOURCE_DATA",ds_safe)
  export_ABE_source_data(
    d=d,
    dataset=dataset,
    preferred_endpoint=preferred_endpoint,
    detailed_n=detailed_n,
    source_dir=source_dir
  )

  panel_dir <- file.path(out_root,"01_PANELS",ds_safe)
  dir.create(panel_dir,recursive=TRUE,showWarnings=FALSE)

  save_plot(
    panel_A,
    file.path(panel_dir,"A_Survival_KM"),
    7.2,
    6.5
  )

  std_panels <- list(
    B_Response_box=panel_B,
    C_Response_rate=panel_C,
    D_ROC=panel_D,
    E_Response_detail_or_composition=panel_E,
    F_DiseaseControl_or_RawUP=panel_F,
    G_DOWN_response=panel_G,
    H_Secondary_or_forest=panel_H
  )

  for (nm in names(std_panels)) {
    obj <- std_panels[[nm]]
    if (inherits(obj,"ggplot")) {
      save_plot(obj,file.path(panel_dir,nm),6.2,5.2)
    }
  }

  summary_plot <- (
    panel_A | panel_B | panel_C | panel_D
  ) / (
    panel_E | panel_F | panel_G | panel_H
  ) +
    plot_layout(
      widths=c(1.25,1,1,1),
      heights=c(1,1)
    ) +
    plot_annotation(
      title=pretty_dataset(dataset),
      subtitle=paste0(
        "PRECOG v1 source | Version-B summary | primary score: ",PRIMARY_LABEL,
        " | High = red, Low = blue",
        ifelse(
          KINGET_CORRECT_EVENT_CODING &&
            grepl("Kinget_2024_KIRC",dataset,ignore.case=TRUE),
          " | Kinget survival event coding corrected",
          ""
        )
      ),
      tag_levels="A",
      theme=theme(
        plot.title=element_text(size=17,face="bold",hjust=.5),
        plot.subtitle=element_text(size=10,hjust=.5),
        plot.tag=element_text(size=13,face="bold")
      )
    )

  save_plot(
    summary_plot,
    file.path(
      out_root,"02_SUMMARY_FIGURES",
      paste0(ds_safe,"__PRECOG_VersionB_2x4_summary")
    ),
    SUMMARY_WIDTH,
    SUMMARY_HEIGHT
  )

  response_n <- sum(!is.na(d$Responder_binary))
  responders <- sum(safe_numeric(d$Responder_binary)==1,na.rm=TRUE)
  nonresponders <- sum(safe_numeric(d$Responder_binary)==0,na.rm=TRUE)

  availability_list[[dataset]] <- data.frame(
    Dataset=dataset,
    N_patient_score=nrow(d),
    Response_N=response_n,
    Responder_N=responders,
    NonResponder_N=nonresponders,
    Detailed_response_column=r4$column,
    Detailed_response_N=detailed_n,
    PFS_usable=pfs_ok,
    PFS_N=sum(complete.cases(d$PFS_time,d$Summary_PFS_event)&safe_numeric(d$PFS_time)>0),
    PFS_events=sum(d$Summary_PFS_event==1,na.rm=TRUE),
    OS_usable=os_ok,
    OS_N=sum(complete.cases(d$OS_time,d$Summary_OS_event)&safe_numeric(d$OS_time)>0),
    OS_events=sum(d$Summary_OS_event==1,na.rm=TRUE),
    Preferred_survival_panel=preferred_endpoint,
    Panel_E_content=ifelse(detailed_n>=10,"CR/PR/SD/PD","Binary response composition"),
    Panel_F_content=ifelse(detailed_n>=10,"CR/PR/SD vs PD","Raw STEP_UP response box"),
    Panel_G_content="Raw STEP_DOWN response box",
    Panel_H_content=ifelse(
      !is.na(secondary_endpoint),
      paste0(secondary_endpoint," KM"),
      "Response + survival OR/HR forest"
    )
  )

  score_map <- c(
    "STEP_UP raw"=RAW_UP_SCORE,
    "CAF3"=CAF3_SCORE,
    "STEP_UP CAF3-adjusted"=PRIMARY_SCORE,
    "STEP_DOWN raw"=DOWN_SCORE
  )

  for (snm in names(score_map)) {
    sc <- score_map[[snm]]

    rr <- fit_response_model(d,sc,snm)
    if (!is.null(rr)) {
      rr$Dataset <- dataset
      model_list[[length(model_list)+1]] <- rr
    }

    for (ep in c("PFS","OS")) {
      sr <- fit_survival_model(d,ep,sc,snm)
      if (!is.null(sr)) {
        sr$Dataset <- dataset
        model_list[[length(model_list)+1]] <- sr
      }
    }
  }
}

event_qc_df <- if(length(event_qc_list)) bind_rows(event_qc_list) else data.frame()
availability_df <- if(length(availability_list)) bind_rows(availability_list) else data.frame()
model_df <- if(length(model_list)) bind_rows(model_list) else data.frame()

write.csv(
  event_qc_df,
  file.path(out_root,"00_QC","00_Kinget_PRECOG_event_coding_QC.csv"),
  row.names=FALSE
)

write.csv(
  availability_df,
  file.path(out_root,"00_QC","01_Kinget_PRECOG_VersionB_panel_availability.csv"),
  row.names=FALSE
)

write.csv(
  model_df,
  file.path(out_root,"00_QC","02_Kinget_PRECOG_continuous_OR_HR_results.csv"),
  row.names=FALSE
)

if (nrow(model_df)>0) {

  primary_model <- model_df %>%
    filter(Score=="STEP_UP CAF3-adjusted")

  if (nrow(primary_model)>0) {

    primary_model$Dataset_short <- ifelse(
      grepl("CTLA4_PD1",primary_model$Dataset,ignore.case=TRUE),
      "Kinget CTLA4+PD1",
      "Kinget PD1"
    )

    primary_model$Label <- paste0(
      primary_model$Dataset_short," | ",primary_model$Endpoint
    )

    primary_model$Label <- factor(
      primary_model$Label,
      levels=rev(unique(primary_model$Label))
    )

    primary_model$Txt <- paste0(
      primary_model$Effect_name,"=",sprintf("%.2f",primary_model$Effect),
      " [",sprintf("%.2f",primary_model$CI_low),", ",sprintf("%.2f",primary_model$CI_high),"]\n",
      p_fmt(primary_model$P)
    )

    pc <- ggplot(
      primary_model,
      aes(Label,Effect,ymin=CI_low,ymax=CI_high,shape=Endpoint)
    )+
      geom_hline(yintercept=1,linetype=2)+
      geom_errorbar(width=.12)+
      geom_point(size=2.7)+
      geom_text(aes(label=Txt),hjust=-.08,size=3)+
      scale_y_log10()+
      coord_flip(clip="off")+
      labs(
        title=paste0("Kinget KIRC PRECOG: ",PRIMARY_LABEL),
        subtitle="Response OR>1 = more non-response; survival HR>1 = worse outcome",
        x=NULL,
        y="OR / HR (log scale)",
        shape=NULL
      )+
      theme_pub()+
      theme(plot.margin=margin(5.5,100,5.5,5.5))

    save_plot(
      pc,
      file.path(
        out_root,"02_SUMMARY_FIGURES",
        "Kinget_KIRC_PRECOG_primary_adjusted_effect_comparison"
      ),
      9,
      max(5,2.8+.45*nrow(primary_model))
    )
  }
}

manifest <- data.frame(
  Panel=LETTERS[1:8],
  Intended_content=c(
    "Preferred PFS/OS KM: adjusted STEP High vs Low + risk table + log-rank P + continuous HR/P",
    "response plots_box: adjusted STEP in Responder vs Non-responder",
    "response plots_rate: response rate in adjusted STEP Low vs High + n/% + Fisher P + non-response OR",
    "ROC: raw UP / CAF3 / adjusted UP / raw DOWN",
    "CR/PR/SD/PD if genuinely available; otherwise responder/non-responder stacked composition",
    "CR/PR/SD vs PD if detailed response exists; otherwise raw STEP_UP response box",
    "Reciprocal-support panel: raw STEP_DOWN response box",
    "Second survival KM if both PFS/OS usable; otherwise response + survival OR/HR forest"
  )
)

write.csv(
  manifest,
  file.path(out_root,"00_QC","99_VersionB_panel_manifest.csv"),
  row.names=FALSE
)

writeLines(
  c(
    "PRECOG VERSION-B SUMMARY TEST: KINGET KIRC | SOURCE = PRECOG_STEP_FORMAL_ANALYSIS_v1 | v1.1 fix1",
    "",
    "v2.0: retains the vector-safe p_fmt fix and adds individual Panel A + exact A/B/E source-data exports.",
    "",
    "Primary score: STEP_UP_resid_CAF3_z",
    "Response OR > 1 = more non-response.",
    "Survival HR > 1 = worse outcome.",
    "",
    "High = red; Low = blue.",
    "Non-responder = red; Responder = blue.",
    "",
    "Kinget survival QC:",
    "This script corrects Kinget event coding from raw PRECOG PFS_Censor/OS_Censor when present.",
    "The correction is limited to datasets whose names contain Kinget_2024_KIRC.",
    "",
    "Detailed CR/PR/SD/PD is never inferred from binary Responder.",
    "If unavailable, scientifically honest fallback panels are used.",
    "",
    "Check first:",
    "00_QC/00_Kinget_PRECOG_event_coding_QC.csv",
    "00_QC/01_Kinget_PRECOG_VersionB_panel_availability.csv",
    "00_QC/02_Kinget_PRECOG_continuous_OR_HR_results.csv",
    "02_SUMMARY_FIGURES/",
    "03_SOURCE_DATA/<dataset>/  (exact manuscript source data for Panels A/B/E)"
  ),
  file.path(out_root,"README_PRECOG_VersionB_Kinget_test.txt")
)

cat("\n==================================================================================================\n")
cat("PRECOG VERSION-B KINGET KIRC SUMMARY COMPLETE | SOURCE = PRECOG v1\n")
cat("==================================================================================================\n")
cat("Analysis directory:\n",analysis_dir,"\n\n",sep="")
cat("Output directory:\n",out_root,"\n\n",sep="")
cat("Datasets processed:\n")
for (ds in target_datasets) cat("  - ",ds,"\n",sep="")
cat("\nCHECK FIRST:\n")
cat("  00_QC/00_Kinget_PRECOG_event_coding_QC.csv\n")
cat("  00_QC/01_Kinget_PRECOG_VersionB_panel_availability.csv\n")
cat("  00_QC/02_Kinget_PRECOG_continuous_OR_HR_results.csv\n")
cat("  01_PANELS/<dataset>/A_Survival_KM.pdf + .png\n")
cat("  02_SUMMARY_FIGURES/*__PRECOG_VersionB_2x4_summary.pdf\n")
cat("  03_SOURCE_DATA/<dataset>/A_* + B_* + E_* source-data CSVs\n")
cat("==================================================================================================\n")
