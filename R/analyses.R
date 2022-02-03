library(ggplot2)
library(parallel)
library(readxl)
library(scales)
library(writexl)

options(mc.cores = detectCores() - 1) 

# Working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-HYPO)")

# Output directory
outdir <- paste0("results/analyses_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# ------------------- Data importation and preprocessing -------------------- #

# Import data
fileName <- "data-raw/Final_data_for_statistical_analyses_24.01.2022.xlsx"
dta <- read_xlsx(fileName, sheet = "data")

# Remove columns which contain the units
dta <- dta[!sapply(names(dta), grepl, pattern = " unit$")]

# Unused variables (`HPP is true` == HPP)
dta <- dta[!(names(dta) %in% c("HPP is true", "Date BS", "HPP symptoms"))]

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
      r <- aggregate(x ~ grp, d, fcts[[z]])
      names(r)[2] <- z
      return(r)
    }))
    r <- unlist(lapply(r$grp, function(z) {
      w <- r[r$grp == z, names(r) != "grp"]
      names(w) <- paste(z, names(w), sep = ".")
      return(w)
    }))
    if (x %in% c("T2D_po", "HPP")) {
      pv1 <- NA
    } else {
      pv1 <- chisq.test(table(d$x, d$grp))$p.value
    }
    if (x == "HPP") {
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
X <- names(dta)[!(names(dta) %in% c("Subject_ID", "HPP"))]
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

