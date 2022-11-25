library(readxl)

# Working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-HYPO)")

if (FALSE) {
# Output directory
outdir <- paste0("results/analyses_", format(Sys.Date(), "%Y%m%d"))
#outdir <- "results/analyses_dev"
if (!dir.exists(outdir)) dir.create(outdir)
}

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

# Multivariable analyses
fml <- HPP ~ glyc0 + glyc120 + Age_BS
fit <- glm(HPP ~ glyc0 + glyc120 + Age_BS, family = binomial, data = dta)
or <- exp(cbind(odds.ratio = coef(fit), suppressMessages(confint(fit))))
or <- cbind(or, p.value = coef(summary(fit))[, 4])
McFadden_PseudoR2 <- as.numeric(1 - logLik(fit) / logLik(update(fit, . ~ 1)))
McFadden_PseudoR2
