library(broom.mixed)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)
library(ggpubr)
library(parallel)
library(readxl)
library(scales)
library(sjPlot)
library(survival)
library(survminer)
library(writexl)

options(mc.cores = detectCores() - 1) 

# Working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-HYPO)")

# Output directory
outdir <- paste0("results/analyses_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# ------------------- Data importation and preprocessing -------------------- #

# Import data
fileName <- "data-raw/Final_data_for_statistical_analyses_08.02.2022.xlsx"
dta <- read_xlsx(fileName, sheet = "data")

# Remove columns which contain the units
dta <- dta[!sapply(names(dta), grepl, pattern = " unit$")]

# Unused variables (`HPP is true` == HPP)
with(dta, table(`HPP is true`, HPP, useNA = "ifany"))
dta$HPP[is.na(dta$HPP)] <- 0
dta <- dta[!(names(dta) %in% c("HPP is true", "HPP symptoms"))]

# Rename variables
names(dta) <- gsub("( |-)", "_", names(dta))

# Convert character variables which contain numeric values to numeric
for (j in 1:ncol(dta)) {
  x <- dta[[j]]
  if (all(is.na(x) | grepl("^(0|1)$", x))) {
    dta[[j]] <- as.logical(as.numeric(x))
  } else if (all(is.na(x) | grepl("^(-)?[0-9]+$", x))) {
    dta[[j]] <- as.integer(x)
  } else if (all(is.na(x) | grepl("^(-)?[0-9]+(\\.[0-9]+)?$", x))) {
    dta[[j]] <- as.numeric(x)
  }
}
rm(j, x)

# Recode HPP_timepoint and Nadir_timepoint
for (v in c("HPP_timepoint", "Nadir_timepoint")) {
  x <- dta[[v]]
  x <- sub("M", "/12", sub("Y", "", sub("^.+_", "", x)))
  dta[[v]] <- unname(sapply(x, function(z) eval(parse(text = z))))
}
rm(v, x)

# Recode Gender (factor to logical)
names(dta)[names(dta) == "Gender"] <- "Male"
dta$Male <- dta$Male == "M"

# Follow-up time
dta$Date_BS <- as.Date(dta$Date_BS, format = "%d/%m/%Y")
dta$FU_time <- as.numeric(as.Date("2021-08-07") - dta$Date_BS) / 365.2425
dta$FU_time <- round(dta$FU_time * 4) / 4

# New variable: HPP_post_nadir
dta$HPP_post_nadir <- dta$HPP_timepoint > dta$Nadir_timepoint
with(dta, table(HPP, HPP_post_nadir, useNA = "ifany"))

# -------------------------- Descriptive analyses --------------------------- #

# Descriptive analyses
X <- list(
  bin = names(dta)[sapply(dta, class) == "logical"],
  num = names(dta)[sapply(dta, class) %in% c("integer", "numeric")]
)
descr_tbl <- list(
  binary = do.call(rbind, mclapply(X[["bin"]], function(x) {
    d <- rbind(data.frame(x = dta[[x]], grp = ifelse(dta$HPP, "HPP", "noHPP")),
               data.frame(x = dta[[x]], grp = "all"))
    d <- na.omit(d)
    d$grp <- factor(d$grp, c("all", "noHPP", "HPP"))
    N <- function(x) sum(!is.na(x))
    Merge <- function(u, v) merge(u, v, by = "grp")
    fcts <- list(N = N, n = sum, p = mean)
    r <- Reduce(Merge, lapply(names(fcts), function(z) {
      r <- aggregate(x ~ grp, d, fcts[[z]], drop = FALSE)
      names(r)[2] <- z
      return(r)
    }))
    r <- unlist(lapply(r$grp, function(z) {
      w <- r[r$grp == z, names(r) != "grp"]
      names(w) <- paste(z, names(w), sep = ".")
      return(w)
    }))
    v <- names(r)[grep("\\.(N|n)$", names(r))]
    r[v][is.na(r[v])] <- 0
    if (x %in% c("T2D_po", "HPP", "HPP_post_nadir")) {
      pv1 <- NA
    } else {
      pv1 <- chisq.test(table(d$x, d$grp))$p.value
    }
    if (x %in% c("HPP", "HPP_post_nadir")) {
      pv2 <- NA
    } else {
      pv2 <- fisher.test(table(d$x, d$grp))$p.value
    }
    cbind(variable = x, as.data.frame(t(r)), chisq.test.pvalue = pv1,
          fisher.test.pvalue = pv2)
  })),
  numeric = do.call(rbind, mclapply(X[["num"]], function(x) {
    d <- na.omit(dta[c("HPP", x)])
    u <- d[[x]]
    u0 <- d[!d$HPP, x, drop = TRUE]
    u1 <- d[d$HPP, x, drop = TRUE]
    if (length(u0) >=1 & length(u1) >= 1) {
      t.test.pval <- t.test(u0, u1)$p.value
      wilcox.test.pval <- wilcox.test(u0, u1, exact = FALSE)$p.value
    } else {
      t.test.pval <- NA
      wilcox.test.pval <- NA
    }
    data.frame(
      variable = x,
      all.n = length(u),
      all.mean = mean(u),
      all.sd = sd(u),
      all.min = if (length(u) > 0) min(u) else NA,
      all.q25 = quantile(u, .25)[[1]],
      all.median = median(u),
      all.q75 = quantile(u, .75)[[1]],
      all.max = if (length(u) > 0) max(u) else NA,
      noHPP.n = length(u0),
      noHPP.mean = mean(u0),
      noHPP.sd = sd(u0),
      noHPP.min = if (length(u0) > 0) min(u0) else NA,
      noHPP.q25 = quantile(u0, .25)[[1]],
      noHPP.median = median(u0),
      noHPP.q75 = quantile(u0, .75)[[1]],
      noHPP.max = if (length(u0) > 0) max(u0) else NA,
      HPP.n = length(u1),
      HPP.mean = mean(u1),
      HPP.sd = sd(u1),
      HPP.min = if (length(u1) > 0) min(u1) else NA,
      HPP.q25 = quantile(u1, .25)[[1]],
      HPP.median = median(u1),
      HPP.q75 = quantile(u1, .75)[[1]],
      HPP.max = if (length(u1) > 0) max(u1) else NA,
      t.test.pval = t.test.pval,
      wilcox.test.pval = wilcox.test.pval
    )
  }))
)
write_xlsx(descr_tbl, file.path(outdir, "descriptive_analyses.xlsx"))

# Descriptive analyses - Figures
descr_figs <- unlist(recursive = FALSE, mclapply(X[["num"]], function(x) {
  d <- na.omit(dta[c("HPP", x)])
  d$HPP <- factor(ifelse(d$HPP, "HPP", "No HPP"), c("No HPP", "HPP"))
  d$HPP <- droplevels(d$HPP)
  fig1 <- ggplot(d, aes_string(x = "HPP", y = x)) +
    geom_boxplot() +
    theme_bw() +
    labs(x = "", title = x)
  fig2 <- ggplot(d, aes_string(x = x, fill = "HPP")) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = c("#5e81ac", "#8fbcbb")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(title = x)
  list(fig1, fig2)
}))
pdf(file.path(outdir, "descriptive_analyses.pdf"))
for(fig in descr_figs) {
  print(fig)
}
dev.off()

# -------------------------- Logistic regressions --------------------------- #

# Explanatory variables
X <- names(dta)
X <- X[!(X %in% c("Subject_ID", "HPP", "Date_BS", "FU_time"))]
X <- X[!sapply(X, function(x) all(is.na(dta[!dta$HPP, x])))]
X <- X[X != "T2D_po"]  # T2D_po -> only one case

# Univariable regression - Odds ratio table
uv_reg_tbl <- do.call(rbind, lapply(X, function(x) {
  fml <- as.formula(paste("HPP ~", x))
  fit <- glm(fml, family = binomial, data = dta)
  or <- exp(cbind(odds.ratio = coef(fit), suppressMessages(confint(fit))))
  or <- cbind(or, p.value = coef(summary(fit))[, 4])[-1, , drop = FALSE]
  rownames(or) <- NULL
  cbind(data.frame(variable = x, n = nrow(fit$model)), or)
}))
write_xlsx(uv_reg_tbl, file.path(outdir, "univariable_regressions.xlsx"))

# Univariable regression - Graphical summary
d <- t(sapply(setNames(X, X), function(u) {
  y <- dta$HPP
  x <- dta[[u]]
  if (!is.logical(x)) x <- x / sd(x, na.rm = TRUE)
  fit <- glm(y ~ x, family = binomial)
  exp(cbind(or = coef(fit), suppressMessages(confint(fit))))[-1, ]
}))
d <- cbind(data.frame(x = rownames(d)), d)
d$x <- factor(d$x, rev(X))
uv_reg_fig <- ggplot(data = d, aes(x = or, y = x)) +
    geom_point() +
    geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`), height = 0.6) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(x = "Odds ratio", y = "", title = "Odds ratios graphical summary",
         subtitle = "Standardized explanatory variables")
tiff(file.path(outdir, "univariable_regressions.tiff"),  width = 4800,
     height = 4800, res = 600, compression = "zip")
print(uv_reg_fig)
dev.off()
rm(d)

# Univariable regression - Figures of predicted values
uv_reg_figs <- mclapply(X, function(x) {
  fml <- as.formula(paste("HPP ~", x))
  fit <- glm(fml, family = binomial, data = dta)
  d <- cbind(fit$model, prob = predict(fit, type = "response"))
  if (is.logical(d[[x]])) {
    fig <- ggplot(unique(d[c(x, "prob")]), aes_string(x = x, y = "prob")) +
       geom_point() +
       scale_y_continuous(limits = 0:1) +
       theme_bw()
  } else {
    d$HPP2 <- factor(ifelse(d$HPP, "HPP", "No HPP"), c("No HPP", "HPP"))
    fig <- ggplot(d, aes_string(x = x)) +
      geom_point(aes(y = as.numeric(HPP), color = HPP2), shape = 3) +
      geom_point(aes(y = prob, color = HPP2)) +
      geom_line(aes(y = prob), color = "grey70") +
      scale_color_manual(values = c("#5e81ac", "#bf616a")) +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank())
  }
  fig + labs(y = "Predicted probability of HPP", title = x)
})
pdf(file.path(outdir, "univariable_regressions.pdf"))
for(fig in uv_reg_figs) {
  print(fig)
}
dev.off()
rm(fig)

# Multivariable analyses
X <- c("weightPO", "BMI_po", "glyc0", "glyc120", "Insul0", "Insul120",
       "HbA1c", "HOMA_IR", "HOMA_B", "Matsuda", "Age_BS", "Male", "HTApo",
       "DyslipidemiaPO", "HyperuricemiePO", "HMG", "Steatose")
X <- list(X, c(X, "FMTotpc", "FMTissAndpc", "FMI", "VAT", "RAG"))
names(X) <- paste0("with", c("out", ""), "_DXA")
mv_reg_tbl <- mclapply(X, function(x) {
  fml <- as.formula(paste("HPP ~", paste(x, collapse = " + ")))
  fit1 <- glm(fml, family = binomial, data = na.omit(dta[c("HPP", x)]))
  fit2 <- stats::step(fit1, trace = FALSE)
  fit3 <- glm(formula(fit2), family = binomial, data = dta)
  or <- exp(cbind(odds.ratio = coef(fit3), suppressMessages(confint(fit3))))
  or <- cbind(or, p.value = coef(summary(fit3))[, 4])
  z <- data.frame(n = c(nrow(fit3$model), rep(NA, nrow(or) - 1)),
                  coef_name = rownames(or))
  rownames(or) <- NULL
  cbind(z, or)
})
names(mv_reg_tbl) <- c("model1", "model2")
write_xlsx(mv_reg_tbl, file.path(outdir, "multivariable_regressions.xlsx"))
rm(X)

# ------------------------------ HPP incidence ------------------------------ #

# Survival analysis
dta$HPP_survtime <- with(dta, ifelse(HPP, HPP_timepoint, FU_time))
HPP_surv_fit <- survfit(Surv(HPP_survtime, HPP) ~ 1, data = dta)

# Cumulative incidence table
HPP_surv_tbl <- c("time", "n.risk", "n.event", "surv", "lower", "upper") %>%
 sapply(function(z) HPP_surv_fit[[z]]) %>%
 as_tibble() %>%
 mutate(time_inc = time - lag(time, default = 0),
        incidence.rate = cumsum(n.event) / cumsum(n.risk * time_inc),
        cumulative.incidence = 1 - surv,
        lower = 1 - upper, lower = 1 - upper) %>%
 select(time, n.risk, n.event, cumulative.incidence, lower, upper,
        incidence.rate)
write_xlsx(HPP_surv_tbl, file.path(outdir, "HPP_incidence.xlsx"))

# Cumulative incidence figure
HPP_surv_fig <- ggsurvplot(fit = HPP_surv_fit, fun = "event",
                           xlab = "Years after bariatric surgery",
                           ylab = "Cumulative incidence", legend = "none")
tiff(filename = file.path(outdir, "HPP_incidence.tiff"),
     height = 3600, width = 5400, res = 1024, compression = "zip")
print(HPP_surv_fig)
dev.off()

# HPP events
HPP_events_fig <- dta %>%
  mutate(HPP = factor(HPP, c(FALSE, TRUE), c("NoHPP", "HPP")),
         Subject_ID = factor(Subject_ID, Subject_ID[order(HPP_survtime)])) %>%
  ggplot(aes(Subject_ID, HPP_survtime)) + 
  geom_bar(stat = "identity", width = 0.2) + 
  geom_point(aes(color = HPP, shape = HPP)) +
  coord_flip() +
  scale_color_manual(values = c("#5e81ac", "#bf616a")) +
  theme_minimal() + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 6)) +
  labs(y = "Years after bariatric surgery")
tiff(filename = file.path(outdir, "HPP_events.tiff"), height = 5400,
     width = 8100, res = 256, compression = "zip")
print(HPP_events_fig)
dev.off()

# ---------------------------- EBMIL post Nadir ----------------------------- #

# Longitudinal data
v <- grep("EBMIL_", names(dta), value = TRUE)
lg <- reshape(as.data.frame(dta[c("Subject_ID", v)]), varying = v,
              v.names = "EBMIL", timevar = "Time", idvar = "Subject_ID",
              times = sub("EBMIL_", "", v), direction = "long")
lg <- lg[!is.na(lg$EBMIL), ]
v <- c("Subject_ID", "HPP", "HPP_timepoint", "Nadir_timepoint",
       "BMI_po", "BMI_nadir", "Age_BS", "Male")
lg <- merge(lg, dta[v], by = "Subject_ID")
lg$Time0 <- lg$Time
lg$Time <- as.numeric(ifelse(lg$Time == "HPP", lg$HPP_timepoint, ifelse(
  lg$Time == "nadir", lg$Nadir_timepoint, sub("Y", "", lg$Time))))
lg$Time_post_nadir <- lg$Time - lg$Nadir_timepoint
lg$HPP_time_post_nadir <- lg$HPP_timepoint - lg$Nadir_timepoint
lg$HPP[lg$HPP & lg$Time_post_nadir < lg$HPP_time_post_nadir] <- FALSE
lg <- lg[lg$Time_post_nadir >= 0, ]
lg$HPP <- factor(lg$HPP, c(FALSE, TRUE), c("NoHPP", "HPP"))
lg$Male <- factor(lg$Male, c(FALSE, TRUE), c("Female", "Male"))
rm(v)

# Regressions
K <- c(without_covariables = 1, with_covariables = 2)
ebmil_reg <- mclapply(K, function(k) {
  fml <- EBMIL ~ HPP * Time_post_nadir
  fml <- list(fml, update(fml, . ~ . + BMI_po + BMI_nadir + Age_BS + Male))
  fml <- lapply(fml, function(z) update(z, . ~ . + (1 | Subject_ID)))[[k]]
  fit <- lmer(fml, data = lg)
  summary(fit)
  tbl <- tidy(fit, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(
      term = sub("HPPHPP", "HPP", term),
      term = sub("MaleMale", "Male", term),
      term = sub(":", " x ", term),
      term = sub("sd__", "SD ", term),
    )
  sttl <- paste(c("Without", "With")[k], "covariables")
  figs <- list()
  dgp <- plot_model(fit, type = "diag")
  dgp[[2]] <- dgp[[2]]$Subject_ID
  figs$diag <- ggarrange(plotlist = dgp, nrow = 2, ncol = 2) %>%
    annotate_figure(top = text_grob(paste0("Diagnostic plots\n", sttl),
                                    face = "bold", size = 16)) %>%
    suppressMessages()
  figs$pred1 <- augment(fit) %>%
    group_by(Subject_ID) %>%
    filter(n() > 1) %>%
    ggplot(aes(Time_post_nadir, EBMIL)) +
    geom_point(aes(colour = HPP)) +
    scale_color_manual(values = c("#5e81ac", "#bf616a")) +
    geom_line(aes(y = .fitted)) +
    facet_wrap(~ Subject_ID) +
    labs(title = "Individual predictions", subtitle = sttl) +
    theme(legend.position = "bottom", legend.title = element_blank())
  figs$pred2  <- ggpredict(fit, c("Time_post_nadir", "HPP")) %>%
    as_tibble() %>%
    select(x, predicted, conf.low, conf.high, group) %>%
    ggplot(aes(x = x, y = predicted)) +
    geom_line(aes(colour = group)) +
    scale_color_manual(values = c("#5e81ac", "#bf616a")) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
                alpha = 0.3, show.legend = FALSE) +
    labs(x = "Time_post_nadir", y = "EBMIL",
         title = "Fixed effects predictions", subtitle = sttl) +
    theme(legend.position = "bottom", legend.title = element_blank())
  list(fit = fit, tbl = tbl, figs = figs)
})
rm(K)

# Any singular fit ?
b <- sapply(ebmil_reg, function(r) isSingular(r$fit))
if (any(unlist(b))) warning("Singular fit")
rm(b)

# Export coefficient tables
write_xlsx(lapply(ebmil_reg, function(r) r$tbl),
           file.path(outdir, "EBMIL_regression_coefficients.xlsx"))

# Export results - Figures
o <- file.path(outdir, "EBMIL_regression_figures")
if (!dir.exists(o)) dir.create(o)
for (k in 1:2) {
  s <- c("nocov", "cov")[k]
  tiff(filename = file.path(o, paste0("diagnostic_plots_", s, ".tiff")),
       height = 3600, width = 5400, res = 384, compression = "zip")
  print(ebmil_reg[[k]]$figs$diag)
  dev.off()
  tiff(filename = file.path(o, paste0("individual_predictions_", s, ".tiff")),
       height = 7200, width = 10800, res = 512, compression = "zip")
  print(ebmil_reg[[k]]$figs$pred1)
  dev.off()
  tiff(
    filename = file.path(o, paste0("fixed_effects_predictions_", s, ".tiff")),
    height = 3600, width = 5400, res = 1024, compression = "zip")
  print(ebmil_reg[[k]]$figs$pred2)
  dev.off()
}
rm(k, o, s)

# ------------------------------ Session info ------------------------------- #

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
