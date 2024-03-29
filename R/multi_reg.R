library(readxl)
library(writexl)

# Working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-HYPO)")

# Output directory
outdir <- paste0("results/analyses_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Import data
fileName <- "data-raw/Final_data_for_statistical_analyses_08.02.2022.xlsx"
dta <- read_xlsx(fileName, sheet = "data")

# Rename variables
names(dta) <- gsub("( |-)", "_", names(dta))

# Recode HPP
dta$HPP[is.na(dta$HPP)] <- 0

# Select variable and keep complete cases
dta <- na.omit(dta[c("HPP", "glyc0", "glyc120", "Age_BS")])

# character -> numeric
dta <- as.data.frame(lapply(dta, as.numeric))

# Multivariable analysis
fml <- HPP ~ glyc0 + glyc120 + Age_BS
fit <- glm(HPP ~ glyc0 + glyc120 + Age_BS, family = binomial, data = dta)
or <- exp(cbind(odds.ratio = coef(fit), suppressMessages(confint(fit))))
or <- cbind(or, p.value = coef(summary(fit))[, 4])
McFadden_PseudoR2 <- as.numeric(1 - logLik(fit) / logLik(update(fit, . ~ 1)))
McFadden_PseudoR2

# Univariate analyses on the complete cases
or0 <- do.call(rbind, lapply(c("glyc0", "glyc120", "Age_BS"), function(x) {
  fit0 <- glm(formula(paste("HPP ~", x)), family = binomial, data = fit$model)
  or <- exp(cbind(odds.ratio = coef(fit0), suppressMessages(confint(fit0))))
  or <- cbind(or, p.value = coef(summary(fit0))[, 4])
  or[2, , drop = FALSE]
}))

# Unify and export tables
t0 <- or0
colnames(t0) <- paste(colnames(t0), "(univ)")
t0 <- cbind(data.frame(variable = rownames(t0)), t0)
t1 <- data.frame(variable = rownames(or),
                 n = c(nrow(fit$model), rep(NA, nrow(or) - 1)))
t1 <- cbind(t1, or)
or_tbl <- merge(t0, t1, by = "variable", all = TRUE)
or_tbl <- or_tbl[, c(1, 6, 2:5, 7:10)]
write_xlsx(or_tbl, file.path(outdir, "multivariable_regression.xlsx"))
