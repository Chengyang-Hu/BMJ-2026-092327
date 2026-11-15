# BMJ-2026-092327: Esomeprazole and postoperative nausea and vomiting
# Run this script section by section in RStudio.
# 1. Working directory and packages --------------------------------------------
getwd()
setwd("/Users/cy.hu/Downloads/")  # Replace "." with the folder containing BMJ_data.xlsx.
# Example for Windows: setwd("D:/BMJ_PONV")
# Example for macOS:   setwd("/Users/yourname/Documents/BMJ_PONV")

# Run the following line once if these packages are not installed:
# install.packages(c("readxl", "tableone", "sandwich", "lmtest"))
library(readxl)
library(tableone)
library(sandwich)
library(lmtest)

dir.create("BMJ_results", showWarnings = FALSE)

# 2. Read and check the data ----------------------------------------------------
data <- as.data.frame(read_excel("BMJ_data.xlsx", sheet = "data"))
dim(data)
names(data)
str(data)

# group: 0 = control; 1 = esomeprazole. Binary outcomes: 0 = no; 1 = yes.
# itt_flag and pp_flag: 1 = included; 0 = excluded.
table(data$group, useNA = "ifany")
table(data$itt_flag, useNA = "ifany")
table(data$pp_flag, useNA = "ifany")
colSums(is.na(data))
stopifnot(all(data$group %in% c(0, 1)))
stopifnot(all(data$itt_flag %in% c(0, 1)))
stopifnot(all(data$pp_flag %in% c(0, 1)))
stopifnot(!anyDuplicated(data$no))

# Retain the original allocation. Do not remove participants because of outcomes.
itt <- subset(data, itt_flag == 1)
pp  <- subset(data, itt_flag == 1 & pp_flag == 1)
table(itt$group)
table(pp$group)

# A labelled group is used for descriptive tables only.
# Numeric group (0/1) is retained for all regression models.
itt$treatment <- factor(itt$group, levels = c(1, 0),
                        labels = c("Esomeprazole", "Control"))

# Check the current Apfel coding: most participants have score 3.
table(itt$apfel_score, itt$treatment)

# These two variables are binary 24-hour events, NOT times in hours.
table(itt$passing_flatus_within_24h_after_surgery, itt$treatment)
table(itt$defecation_within_24h_after_surgery, itt$treatment)

# Confirm that the two supplied overall PONV fields agree.
table(itt$ponv_occurrence, itt$ponv_0_24h)
# The recorded ponv_0_24h field is the primary endpoint.
# No missing values are imputed. Models use complete cases for their own variables.

# 3. Baseline characteristics ---------------------------------------------------
baseline_variables <- c(
  "age", "height", "weight", "bmi", "hypertension", "diabetes",
  "cardiac_disease", "stroke", "asthma", "migraine_history",
  "previous_chemotherapy", "alcohol_use", "smoking_history",
  "motion_sickness_history", "previous_ponv_history",
  "asa_physical_status", "apfel_score")
baseline_categorical <- c(
  "hypertension", "diabetes", "cardiac_disease", "stroke", "asthma",
  "migraine_history", "previous_chemotherapy", "alcohol_use",
  "smoking_history", "motion_sickness_history", "previous_ponv_history",
  "asa_physical_status", "apfel_score")

baseline_table <- CreateTableOne(vars = baseline_variables, strata = "treatment",
  data = itt, factorVars = baseline_categorical, test = FALSE)
print(baseline_table, showAllLevels = TRUE, contDigits = 1, catDigits = 1)
baseline_export <- print(baseline_table, showAllLevels = TRUE,
                         printToggle = FALSE, noSpaces = TRUE, contDigits = 1, catDigits = 1)
write.csv(baseline_export, "BMJ_results/Table1_baseline.csv")
# Baseline data are descriptive; P values are not needed to assess randomisation.

# 4. Intraoperative characteristics ---------------------------------------------
intra_variables <- c("duration_of_anaesthesia", "duration_of_surgery",
  "propofol_dose", "remifentanil_dose", "sufentanil_dose",
  "intraoperative_fluid_administration")
intra_table <- CreateTableOne(vars = intra_variables, strata = "treatment", data = itt)
intra_export <- print(intra_table, nonnormal = intra_variables,
                      printToggle = FALSE, noSpaces = TRUE, contDigits = 1, catDigits = 1)
intra_export
write.csv(intra_export, "BMJ_results/Table1_intraoperative.csv")
# These variables are summarised as median (IQR). With two arms, the rank-based
# test used by tableone is Kruskal-Wallis, asymptotically equivalent to the
# uncorrected two-sided Mann-Whitney test. 


# 5. Table 2: Primary and key secondary PONV outcomes ----------------------------
# Run Sections 1-4 first. Each table is a printable data frame of formatted text.
# Percentages and median (IQR): 1 decimal. RR and CI: 2 decimals.
# P values: 3 decimals, or <0.001. 
# All calculations use unrounded values. 
# Binary comparisons: Pearson chi-square without continuity correction;
# Fisher's exact test if any expected cell is below 5.
# RR confidence intervals: log-Wald. Zero-event RR intervals are shown as NE;
# no continuity correction is silently added.
# Adjusted RR: modified Poisson regression with HC0 robust standard errors.
# Continuous comparisons: two-sided Wilcoxon rank-sum, continuity correction.
# Missing outcomes are excluded only from the corresponding analysis.
options(width = 240, max.print = 10000)
n_eso <- sum(itt$group == 1)
n_control <- sum(itt$group == 0)
arm_names <- c(sprintf("Esomeprazole group (n=%d)", n_eso),
               sprintf("Control group (n=%d)", n_control))


# PONV within 24 h
d <- itt[complete.cases(itt[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_ponv_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_ponv_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_ponv_0_24h <- glm(ponv_0_24h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_ponv_0_24h, vcov. = vcovHC(fit_ponv_0_24h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_ponv_0_24h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_ponv_0_24h <- c("PONV within 24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Nausea within 24 h
d <- itt[complete.cases(itt[, c("group", "nausea_0_24h")]), ]
e1 <- sum(d$nausea_0_24h[d$group == 1] == 1)
e0 <- sum(d$nausea_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_nausea_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_nausea_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_nausea_0_24h <- glm(nausea_0_24h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_nausea_0_24h, vcov. = vcovHC(fit_nausea_0_24h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_nausea_0_24h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_nausea_0_24h <- c("Nausea within 24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Vomiting within 24 h
d <- itt[complete.cases(itt[, c("group", "vomiting_0_24h")]), ]
e1 <- sum(d$vomiting_0_24h[d$group == 1] == 1)
e0 <- sum(d$vomiting_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_vomiting_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_vomiting_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_vomiting_0_24h <- glm(vomiting_0_24h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_vomiting_0_24h, vcov. = vcovHC(fit_vomiting_0_24h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_vomiting_0_24h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_vomiting_0_24h <- c("Vomiting within 24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Clinically important PONV within 24 h
d <- itt[complete.cases(itt[, c("group", "clinically_important_ponv_0_24h")]), ]
e1 <- sum(d$clinically_important_ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$clinically_important_ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_clinically_important_ponv_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_clinically_important_ponv_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_clinically_important_ponv_0_24h <- glm(clinically_important_ponv_0_24h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_clinically_important_ponv_0_24h, vcov. = vcovHC(fit_clinically_important_ponv_0_24h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_clinically_important_ponv_0_24h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_clinically_important_ponv_0_24h <- c("Clinically important PONV within 24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Rescue antiemetic use
d <- itt[complete.cases(itt[, c("group", "rescue_antiemetic_use")]), ]
e1 <- sum(d$rescue_antiemetic_use[d$group == 1] == 1)
e0 <- sum(d$rescue_antiemetic_use[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_rescue_antiemetic_use <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_rescue_antiemetic_use <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_rescue_antiemetic_use <- glm(rescue_antiemetic_use ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_rescue_antiemetic_use, vcov. = vcovHC(fit_rescue_antiemetic_use, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_rescue_antiemetic_use <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_rescue_antiemetic_use <- c("Rescue antiemetic use",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


Table2 <- as.data.frame(rbind(
  c("Primary outcome", "", "", "", "", "", ""),
  row_ponv_0_24h,
  c("Key secondary outcomes", "", "", "", "", "", ""),
  row_nausea_0_24h,
  row_vomiting_0_24h,
  row_clinically_important_ponv_0_24h,
  row_rescue_antiemetic_use), stringsAsFactors = FALSE)

names(Table2) <- c("Outcome", arm_names, "Risk ratio (95% CI)", "P value", "Adjusted risk ratio (95% CI)", "P value")
cat("\n", "Table 2 | Primary and key secondary PONV-related outcomes in the intention-to-treat population. Values are numbers (percentages) unless stated otherwise.", "\n", sep = "")
print(Table2, row.names = FALSE, right = FALSE)
write.csv(Table2, "BMJ_results/Table2.csv", row.names = FALSE, fileEncoding = "UTF-8")


cat("Adjusted for age, BMI, diabetes, hypertension, previous PONV, motion sickness, and migraine.\n")
# Absolute risk reduction and NNT for the primary and clinically important outcomes.


r <- raw_ponv_0_24h
arr <- unname(r["e0"]/r["n0"] - r["e1"]/r["n1"])
arr_se <- sqrt((r["e0"]/r["n0"])*(1-r["e0"]/r["n0"])/r["n0"] +
               (r["e1"]/r["n1"])*(1-r["e1"]/r["n1"])/r["n1"])
cat(sprintf("ponv_0_24h: ARR %.1f percentage points (95%% CI %.1f to %.1f); NNT %.1f.\n",
  100*arr, 100*(arr-qnorm(0.975)*arr_se), 100*(arr+qnorm(0.975)*arr_se), 1/arr))


r <- raw_clinically_important_ponv_0_24h
arr <- unname(r["e0"]/r["n0"] - r["e1"]/r["n1"])
arr_se <- sqrt((r["e0"]/r["n0"])*(1-r["e0"]/r["n0"])/r["n0"] +
               (r["e1"]/r["n1"])*(1-r["e1"]/r["n1"])/r["n1"])
cat(sprintf("clinically_important_ponv_0_24h: ARR %.1f percentage points (95%% CI %.1f to %.1f); NNT %.1f.\n",
  100*arr, 100*(arr-qnorm(0.975)*arr_se), 100*(arr+qnorm(0.975)*arr_se), 1/arr))


# 6. Table 3: Recovery, pain, QoR-15, and intraoperative safety ------------------

# Median (IQR) QoR-15 preoperative
x1 <- itt$qor15_pre_total[itt$group == 1]
x0 <- itt$qor15_pre_total[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_pre_total <- c("Median (IQR) QoR-15 preoperative",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Median (IQR) QoR-15 postoperative
x1 <- itt$qor15_post_total[itt$group == 1]
x0 <- itt$qor15_post_total[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_post_total <- c("Median (IQR) QoR-15 postoperative",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Median (IQR) VAS at 1 h
x1 <- itt$pain_vas_score_at_1h[itt$group == 1]
x0 <- itt$pain_vas_score_at_1h[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_pain_vas_score_at_1h <- c("Median (IQR) VAS at 1 h",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Median (IQR) VAS at 6 h
x1 <- itt$pain_vas_score_at_6h[itt$group == 1]
x0 <- itt$pain_vas_score_at_6h[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_pain_vas_score_at_6h <- c("Median (IQR) VAS at 6 h",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Median (IQR) VAS at 24 h
x1 <- itt$pain_vas_score_at_24h[itt$group == 1]
x0 <- itt$pain_vas_score_at_24h[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_pain_vas_score_at_24h <- c("Median (IQR) VAS at 24 h",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Emergence agitation
d <- itt[complete.cases(itt[, c("group", "emergence_agitation")]), ]
e1 <- sum(d$emergence_agitation[d$group == 1] == 1)
e0 <- sum(d$emergence_agitation[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_emergence_agitation <- p_raw


row_emergence_agitation <- c("Emergence agitation",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Median (IQR) Ramsay sedation score
x1 <- itt$ramsay_sedation_score[itt$group == 1]
x0 <- itt$ramsay_sedation_score[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_ramsay_sedation_score <- c("Median (IQR) Ramsay sedation score",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Median (IQR) PACU stay, min
x1 <- itt$pacu_stay[itt$group == 1]
x0 <- itt$pacu_stay[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_pacu_stay <- c("Median (IQR) PACU stay, min",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Median (IQR) length of hospital stay, days
x1 <- itt$length_of_hospital_stay[itt$group == 1]
x0 <- itt$length_of_hospital_stay[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_length_of_hospital_stay <- c("Median (IQR) length of hospital stay, days",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Passing flatus within 24 h
d <- itt[complete.cases(itt[, c("group", "passing_flatus_within_24h_after_surgery")]), ]
e1 <- sum(d$passing_flatus_within_24h_after_surgery[d$group == 1] == 1)
e0 <- sum(d$passing_flatus_within_24h_after_surgery[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_passing_flatus_within_24h_after_surgery <- p_raw


row_passing_flatus_within_24h_after_surgery <- c("Passing flatus within 24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Defecation within 24 h
d <- itt[complete.cases(itt[, c("group", "defecation_within_24h_after_surgery")]), ]
e1 <- sum(d$defecation_within_24h_after_surgery[d$group == 1] == 1)
e0 <- sum(d$defecation_within_24h_after_surgery[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_defecation_within_24h_after_surgery <- p_raw


row_defecation_within_24h_after_surgery <- c("Defecation within 24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Rescue analgesic use
d <- itt[complete.cases(itt[, c("group", "rescue_analgesic_use")]), ]
e1 <- sum(d$rescue_analgesic_use[d$group == 1] == 1)
e0 <- sum(d$rescue_analgesic_use[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_rescue_analgesic_use <- p_raw


row_rescue_analgesic_use <- c("Rescue analgesic use",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Intraoperative hypotension
d <- itt[complete.cases(itt[, c("group", "intraoperative_hypotension")]), ]
e1 <- sum(d$intraoperative_hypotension[d$group == 1] == 1)
e0 <- sum(d$intraoperative_hypotension[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_intraoperative_hypotension <- p_raw


row_intraoperative_hypotension <- c("Intraoperative hypotension",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Intraoperative hypertension
d <- itt[complete.cases(itt[, c("group", "intraoperative_hypertension")]), ]
e1 <- sum(d$intraoperative_hypertension[d$group == 1] == 1)
e0 <- sum(d$intraoperative_hypertension[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_intraoperative_hypertension <- p_raw


row_intraoperative_hypertension <- c("Intraoperative hypertension",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Intraoperative bradycardia
d <- itt[complete.cases(itt[, c("group", "intraoperative_bradycardia")]), ]
e1 <- sum(d$intraoperative_bradycardia[d$group == 1] == 1)
e0 <- sum(d$intraoperative_bradycardia[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_intraoperative_bradycardia <- p_raw


row_intraoperative_bradycardia <- c("Intraoperative bradycardia",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Intraoperative tachycardia
d <- itt[complete.cases(itt[, c("group", "intraoperative_tachycardia")]), ]
e1 <- sum(d$intraoperative_tachycardia[d$group == 1] == 1)
e0 <- sum(d$intraoperative_tachycardia[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_intraoperative_tachycardia <- p_raw


row_intraoperative_tachycardia <- c("Intraoperative tachycardia",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Intraoperative allergic reaction
d <- itt[complete.cases(itt[, c("group", "intraoperative_allergic_reaction")]), ]
e1 <- sum(d$intraoperative_allergic_reaction[d$group == 1] == 1)
e0 <- sum(d$intraoperative_allergic_reaction[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_intraoperative_allergic_reaction <- p_raw


row_intraoperative_allergic_reaction <- c("Intraoperative allergic reaction",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


Table3 <- as.data.frame(rbind(
  c("QoR-15", "", "", ""),
  row_qor15_pre_total,
  row_qor15_post_total,
  c("Pain and recovery outcomes", "", "", ""),
  row_pain_vas_score_at_1h,
  row_pain_vas_score_at_6h,
  row_pain_vas_score_at_24h,
  row_emergence_agitation,
  row_ramsay_sedation_score,
  row_pacu_stay,
  row_length_of_hospital_stay,
  row_passing_flatus_within_24h_after_surgery,
  row_defecation_within_24h_after_surgery,
  row_rescue_analgesic_use,
  c("Intraoperative safety outcomes", "", "", ""),
  row_intraoperative_hypotension,
  row_intraoperative_hypertension,
  row_intraoperative_bradycardia,
  row_intraoperative_tachycardia,
  row_intraoperative_allergic_reaction), stringsAsFactors = FALSE)

names(Table3) <- c("Outcome", arm_names, "P value")
cat("\n", "Table 3 | Recovery, pain, quality of recovery, and intraoperative safety outcomes. Values are number (percentage) or median (IQR) unless stated otherwise.", "\n", sep = "")
print(Table3, row.names = FALSE, right = FALSE)
write.csv(Table3, "BMJ_results/Table3.csv", row.names = FALSE, fileEncoding = "UTF-8")


# 7. Table 4: Subgroup analyses of PONV within 24 h -----------------------------
# Interactions use unadjusted robust Poisson models on the RR scale.
# Although the interaction P values differed from those obtained using the 
# original logistic regression models, the interpretation remained unchanged: 
# neither approach provided evidence of treatment-by-subgroup interaction.
# Surgery duration is measured after randomisation; this subgroup is exploratory.

# Overall
d <- itt[complete.cases(itt[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_overall <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_overall <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_overall <- c("Overall",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


itt$bmi_subgroup <- factor(ifelse(itt$bmi < 25, "Low", "High"), levels = c("Low", "High"))
fit_interaction_bmi <- glm(ponv_0_24h ~ group * bmi_subgroup,
  family = poisson(link = "log"), data = itt)
interaction_test <- coeftest(fit_interaction_bmi,
  vcov. = vcovHC(fit_interaction_bmi, type = "HC0"))
p_interaction_bmi <- interaction_test["group:bmi_subgroupHigh", 4]


sub_data <- subset(itt, bmi_subgroup == "Low")

# <25
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_bmi_Low <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_bmi_Low <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_bmi_Low <- c("<25",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


sub_data <- subset(itt, bmi_subgroup == "High")

# ≥25
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_bmi_High <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_bmi_High <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_bmi_High <- c("≥25",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


itt$age_subgroup <- factor(ifelse(itt$age < 50, "Low", "High"), levels = c("Low", "High"))
fit_interaction_age <- glm(ponv_0_24h ~ group * age_subgroup,
  family = poisson(link = "log"), data = itt)
interaction_test <- coeftest(fit_interaction_age,
  vcov. = vcovHC(fit_interaction_age, type = "HC0"))
p_interaction_age <- interaction_test["group:age_subgroupHigh", 4]


sub_data <- subset(itt, age_subgroup == "Low")

# <50
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_age_Low <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_age_Low <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_age_Low <- c("<50",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


sub_data <- subset(itt, age_subgroup == "High")

# ≥50
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_age_High <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_age_High <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_age_High <- c("≥50",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


itt$asa_subgroup <- factor(ifelse(itt$asa_physical_status == 1, "Low", "High"), levels = c("Low", "High"))
fit_interaction_asa <- glm(ponv_0_24h ~ group * asa_subgroup,
  family = poisson(link = "log"), data = itt)
interaction_test <- coeftest(fit_interaction_asa,
  vcov. = vcovHC(fit_interaction_asa, type = "HC0"))
p_interaction_asa <- interaction_test["group:asa_subgroupHigh", 4]


sub_data <- subset(itt, asa_subgroup == "Low")

# I
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_asa_Low <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_asa_Low <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_asa_Low <- c("I",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


sub_data <- subset(itt, asa_subgroup == "High")

# II–III
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_asa_High <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_asa_High <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_asa_High <- c("II–III",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


itt$surgery_subgroup <- factor(ifelse(itt$duration_of_surgery < 72, "Low", "High"), levels = c("Low", "High"))
fit_interaction_surgery <- glm(ponv_0_24h ~ group * surgery_subgroup,
  family = poisson(link = "log"), data = itt)
interaction_test <- coeftest(fit_interaction_surgery,
  vcov. = vcovHC(fit_interaction_surgery, type = "HC0"))
p_interaction_surgery <- interaction_test["group:surgery_subgroupHigh", 4]


sub_data <- subset(itt, surgery_subgroup == "Low")

# <72
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_surgery_Low <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_surgery_Low <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_surgery_Low <- c("<72",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


sub_data <- subset(itt, surgery_subgroup == "High")

# ≥72
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_surgery_High <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_surgery_High <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_surgery_High <- c("≥72",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


Table4 <- as.data.frame(rbind(
  c(row_overall, ""),
  c("BMI, kg/m²", "", "", "", "", ifelse(is.na(p_interaction_bmi), "NA", ifelse(p_interaction_bmi < 0.001, "<0.001", sprintf("%.3f", p_interaction_bmi)))),
  c(row_bmi_Low, ""),
  c(row_bmi_High, ""),
  c("Age, years", "", "", "", "", ifelse(is.na(p_interaction_age), "NA", ifelse(p_interaction_age < 0.001, "<0.001", sprintf("%.3f", p_interaction_age)))),
  c(row_age_Low, ""),
  c(row_age_High, ""),
  c("ASA physical status", "", "", "", "", ifelse(is.na(p_interaction_asa), "NA", ifelse(p_interaction_asa < 0.001, "<0.001", sprintf("%.3f", p_interaction_asa)))),
  c(row_asa_Low, ""),
  c(row_asa_High, ""),
  c("Duration of surgery, min", "", "", "", "", ifelse(is.na(p_interaction_surgery), "NA", ifelse(p_interaction_surgery < 0.001, "<0.001", sprintf("%.3f", p_interaction_surgery)))),
  c(row_surgery_Low, ""),
  c(row_surgery_High, "")), stringsAsFactors = FALSE)

names(Table4) <- c("Subgroup", "Esomeprazole group", "Control group", "Risk ratio (95% CI)", "P value", "P for interaction")
cat("\n", "Table 4 | Subgroup analyses of the primary outcome of PONV within 24 h. Values are events/participants (percentages) unless stated otherwise.", "\n", sep = "")
print(Table4, row.names = FALSE, right = FALSE)
write.csv(Table4, "BMJ_results/Table4.csv", row.names = FALSE, fileEncoding = "UTF-8")


# 8. Supplementary table S1: Complete PONV outcomes ----------------------------

# PONV 0–1 h
d <- itt[complete.cases(itt[, c("group", "ponv_0_1h")]), ]
e1 <- sum(d$ponv_0_1h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_1h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_ponv_0_1h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_ponv_0_1h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_ponv_0_1h <- glm(ponv_0_1h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_ponv_0_1h, vcov. = vcovHC(fit_ponv_0_1h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_ponv_0_1h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_ponv_0_1h <- c("PONV 0–1 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# PONV 1–6 h
d <- itt[complete.cases(itt[, c("group", "ponv_1_6h")]), ]
e1 <- sum(d$ponv_1_6h[d$group == 1] == 1)
e0 <- sum(d$ponv_1_6h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_ponv_1_6h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_ponv_1_6h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_ponv_1_6h <- glm(ponv_1_6h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_ponv_1_6h, vcov. = vcovHC(fit_ponv_1_6h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_ponv_1_6h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_ponv_1_6h <- c("PONV 1–6 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# PONV 6–24 h
d <- itt[complete.cases(itt[, c("group", "ponv_6_24h")]), ]
e1 <- sum(d$ponv_6_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_6_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_ponv_6_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_ponv_6_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_ponv_6_24h <- glm(ponv_6_24h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_ponv_6_24h, vcov. = vcovHC(fit_ponv_6_24h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_ponv_6_24h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_ponv_6_24h <- c("PONV 6–24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Nausea 0–1 h
d <- itt[complete.cases(itt[, c("group", "nausea_0_1h")]), ]
e1 <- sum(d$nausea_0_1h[d$group == 1] == 1)
e0 <- sum(d$nausea_0_1h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_nausea_0_1h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_nausea_0_1h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_nausea_0_1h <- glm(nausea_0_1h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_nausea_0_1h, vcov. = vcovHC(fit_nausea_0_1h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_nausea_0_1h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_nausea_0_1h <- c("Nausea 0–1 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Nausea 1–6 h
d <- itt[complete.cases(itt[, c("group", "nausea_1_6h")]), ]
e1 <- sum(d$nausea_1_6h[d$group == 1] == 1)
e0 <- sum(d$nausea_1_6h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_nausea_1_6h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_nausea_1_6h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_nausea_1_6h <- glm(nausea_1_6h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_nausea_1_6h, vcov. = vcovHC(fit_nausea_1_6h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_nausea_1_6h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_nausea_1_6h <- c("Nausea 1–6 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Nausea 6–24 h
d <- itt[complete.cases(itt[, c("group", "nausea_6_24h")]), ]
e1 <- sum(d$nausea_6_24h[d$group == 1] == 1)
e0 <- sum(d$nausea_6_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_nausea_6_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_nausea_6_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_nausea_6_24h <- glm(nausea_6_24h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_nausea_6_24h, vcov. = vcovHC(fit_nausea_6_24h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_nausea_6_24h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_nausea_6_24h <- c("Nausea 6–24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Vomiting 0–1 h
d <- itt[complete.cases(itt[, c("group", "vomiting_0_1h")]), ]
e1 <- sum(d$vomiting_0_1h[d$group == 1] == 1)
e0 <- sum(d$vomiting_0_1h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_vomiting_0_1h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_vomiting_0_1h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_vomiting_0_1h <- glm(vomiting_0_1h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_vomiting_0_1h, vcov. = vcovHC(fit_vomiting_0_1h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_vomiting_0_1h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_vomiting_0_1h <- c("Vomiting 0–1 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Vomiting 1–6 h
d <- itt[complete.cases(itt[, c("group", "vomiting_1_6h")]), ]
e1 <- sum(d$vomiting_1_6h[d$group == 1] == 1)
e0 <- sum(d$vomiting_1_6h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_vomiting_1_6h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_vomiting_1_6h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_vomiting_1_6h <- glm(vomiting_1_6h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_vomiting_1_6h, vcov. = vcovHC(fit_vomiting_1_6h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_vomiting_1_6h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_vomiting_1_6h <- c("Vomiting 1–6 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


# Vomiting 6–24 h
d <- itt[complete.cases(itt[, c("group", "vomiting_6_24h")]), ]
e1 <- sum(d$vomiting_6_24h[d$group == 1] == 1)
e0 <- sum(d$vomiting_6_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_vomiting_6_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_vomiting_6_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

fit_vomiting_6_24h <- glm(vomiting_6_24h ~ group + age + bmi + diabetes + hypertension +
  previous_ponv_history + motion_sickness_history + migraine_history,
  family = poisson(link = "log"), data = itt)
robust <- coeftest(fit_vomiting_6_24h, vcov. = vcovHC(fit_vomiting_6_24h, type = "HC0"))
beta <- robust["group", 1]
se <- robust["group", 2]
p_adj <- robust["group", 4]
adj_vomiting_6_24h <- c(rr=exp(beta), lower=exp(beta-qnorm(0.975)*se),
  upper=exp(beta+qnorm(0.975)*se), p=p_adj)
adj_text <- sprintf("%.2f (%.2f to %.2f)", exp(beta),
                    exp(beta-qnorm(0.975)*se), exp(beta+qnorm(0.975)*se))

row_vomiting_6_24h <- c("Vomiting 6–24 h",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))),
  adj_text,
  ifelse(is.na(p_adj), "NA", ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))))


TableS1 <- as.data.frame(rbind(
  row_ponv_0_1h,
  row_ponv_1_6h,
  row_ponv_6_24h,
  row_ponv_0_24h,
  row_nausea_0_1h,
  row_nausea_1_6h,
  row_nausea_6_24h,
  row_nausea_0_24h,
  row_vomiting_0_1h,
  row_vomiting_1_6h,
  row_vomiting_6_24h,
  row_vomiting_0_24h,
  row_clinically_important_ponv_0_24h,
  row_rescue_antiemetic_use), stringsAsFactors = FALSE)

names(TableS1) <- c("Outcome", arm_names, "Risk ratio (95% CI)", "P value", "Adjusted risk ratio (95% CI)", "P value")
cat("\n", "Supplementary table S1 | Complete PONV-related outcomes in the intention-to-treat population.", "\n", sep = "")
print(TableS1, row.names = FALSE, right = FALSE)
write.csv(TableS1, "BMJ_results/TableS1.csv", row.names = FALSE, fileEncoding = "UTF-8")


# 9. Supplementary table S2: QoR-15 domain scores ------------------------------
# Use supplied domain scores; do not reverse-score them again.

# Pain (preoperative)
x1 <- itt$qor15_pre_pain[itt$group == 1]
x0 <- itt$qor15_pre_pain[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_pre_pain <- c("Pain (preoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Pain (postoperative)
x1 <- itt$qor15_post_pain[itt$group == 1]
x0 <- itt$qor15_post_pain[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_post_pain <- c("Pain (postoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Physical comfort (preoperative)
x1 <- itt$qor15_pre_physical_comfort[itt$group == 1]
x0 <- itt$qor15_pre_physical_comfort[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_pre_physical_comfort <- c("Physical comfort (preoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Physical comfort (postoperative)
x1 <- itt$qor15_post_physical_comfort[itt$group == 1]
x0 <- itt$qor15_post_physical_comfort[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_post_physical_comfort <- c("Physical comfort (postoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Physical independence (preoperative)
x1 <- itt$qor15_pre_physical_independence[itt$group == 1]
x0 <- itt$qor15_pre_physical_independence[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_pre_physical_independence <- c("Physical independence (preoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Physical independence (postoperative)
x1 <- itt$qor15_post_physical_independence[itt$group == 1]
x0 <- itt$qor15_post_physical_independence[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_post_physical_independence <- c("Physical independence (postoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Psychological support (preoperative)
x1 <- itt$qor15_pre_psychological_support[itt$group == 1]
x0 <- itt$qor15_pre_psychological_support[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_pre_psychological_support <- c("Psychological support (preoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Psychological support (postoperative)
x1 <- itt$qor15_post_psychological_support[itt$group == 1]
x0 <- itt$qor15_post_psychological_support[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_post_psychological_support <- c("Psychological support (postoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Emotional state (preoperative)
x1 <- itt$qor15_pre_emotional_state[itt$group == 1]
x0 <- itt$qor15_pre_emotional_state[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_pre_emotional_state <- c("Emotional state (preoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Emotional state (postoperative)
x1 <- itt$qor15_post_emotional_state[itt$group == 1]
x0 <- itt$qor15_post_emotional_state[itt$group == 0]
q1 <- quantile(x1, c(0.5, 0.25, 0.75), na.rm = TRUE)
q0 <- quantile(x0, c(0.5, 0.25, 0.75), na.rm = TRUE)
p_raw <- wilcox.test(x1, x0, exact = FALSE, correct = TRUE)$p.value
row_qor15_post_emotional_state <- c("Emotional state (postoperative)",
  sprintf("%.1f (%.1f–%.1f)", q1[1], q1[2], q1[3]),
  sprintf("%.1f (%.1f–%.1f)", q0[1], q0[2], q0[3]), ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


TableS2 <- as.data.frame(rbind(
  row_qor15_pre_pain,
  row_qor15_post_pain,
  row_qor15_pre_physical_comfort,
  row_qor15_post_physical_comfort,
  row_qor15_pre_physical_independence,
  row_qor15_post_physical_independence,
  row_qor15_pre_psychological_support,
  row_qor15_post_psychological_support,
  row_qor15_pre_emotional_state,
  row_qor15_post_emotional_state), stringsAsFactors = FALSE)

names(TableS2) <- c("Variable", arm_names, "P value")
cat("\n", "Supplementary table S2 | QoR-15 domain scores. Values are median (IQR).", "\n", sep = "")
print(TableS2, row.names = FALSE, right = FALSE)
write.csv(TableS2, "BMJ_results/TableS2.csv", row.names = FALSE, fileEncoding = "UTF-8")


# 10. Supplementary table S3: Postoperative complications and adverse events ----

# Insomnia
d <- itt[complete.cases(itt[, c("group", "insomnia")]), ]
e1 <- sum(d$insomnia[d$group == 1] == 1)
e0 <- sum(d$insomnia[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_insomnia <- p_raw


row_insomnia <- c("Insomnia",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Electrolyte disturbance
d <- itt[complete.cases(itt[, c("group", "electrolyte_disturbance")]), ]
e1 <- sum(d$electrolyte_disturbance[d$group == 1] == 1)
e0 <- sum(d$electrolyte_disturbance[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_electrolyte_disturbance <- p_raw


row_electrolyte_disturbance <- c("Electrolyte disturbance",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Peripheral oedema
d <- itt[complete.cases(itt[, c("group", "peripheral_oedema")]), ]
e1 <- sum(d$peripheral_oedema[d$group == 1] == 1)
e0 <- sum(d$peripheral_oedema[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_peripheral_oedema <- p_raw


row_peripheral_oedema <- c("Peripheral oedema",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Abdominal pain
d <- itt[complete.cases(itt[, c("group", "abdominal_pain")]), ]
e1 <- sum(d$abdominal_pain[d$group == 1] == 1)
e0 <- sum(d$abdominal_pain[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_abdominal_pain <- p_raw


row_abdominal_pain <- c("Abdominal pain",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Abdominal distension
d <- itt[complete.cases(itt[, c("group", "abdominal_distension")]), ]
e1 <- sum(d$abdominal_distension[d$group == 1] == 1)
e0 <- sum(d$abdominal_distension[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_abdominal_distension <- p_raw


row_abdominal_distension <- c("Abdominal distension",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Constipation
d <- itt[complete.cases(itt[, c("group", "constipation")]), ]
e1 <- sum(d$constipation[d$group == 1] == 1)
e0 <- sum(d$constipation[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_constipation <- p_raw


row_constipation <- c("Constipation",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Increased bowel movements
d <- itt[complete.cases(itt[, c("group", "increased_bowel_movements")]), ]
e1 <- sum(d$increased_bowel_movements[d$group == 1] == 1)
e0 <- sum(d$increased_bowel_movements[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_increased_bowel_movements <- p_raw


row_increased_bowel_movements <- c("Increased bowel movements",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Dry mouth
d <- itt[complete.cases(itt[, c("group", "dry_mouth")]), ]
e1 <- sum(d$dry_mouth[d$group == 1] == 1)
e0 <- sum(d$dry_mouth[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_dry_mouth <- p_raw


row_dry_mouth <- c("Dry mouth",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Cough
d <- itt[complete.cases(itt[, c("group", "cough")]), ]
e1 <- sum(d$cough[d$group == 1] == 1)
e0 <- sum(d$cough[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_cough <- p_raw


row_cough <- c("Cough",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Sputum production
d <- itt[complete.cases(itt[, c("group", "sputum_production")]), ]
e1 <- sum(d$sputum_production[d$group == 1] == 1)
e0 <- sum(d$sputum_production[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_sputum_production <- p_raw


row_sputum_production <- c("Sputum production",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Sore throat
d <- itt[complete.cases(itt[, c("group", "sore_throat")]), ]
e1 <- sum(d$sore_throat[d$group == 1] == 1)
e0 <- sum(d$sore_throat[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_sore_throat <- p_raw


row_sore_throat <- c("Sore throat",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Hoarseness
d <- itt[complete.cases(itt[, c("group", "hoarseness")]), ]
e1 <- sum(d$hoarseness[d$group == 1] == 1)
e0 <- sum(d$hoarseness[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_hoarseness <- p_raw


row_hoarseness <- c("Hoarseness",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Postoperative hypotension
d <- itt[complete.cases(itt[, c("group", "postoperative_hypotension")]), ]
e1 <- sum(d$postoperative_hypotension[d$group == 1] == 1)
e0 <- sum(d$postoperative_hypotension[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_postoperative_hypotension <- p_raw


row_postoperative_hypotension <- c("Postoperative hypotension",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Postoperative hypertension
d <- itt[complete.cases(itt[, c("group", "postoperative_hypertension")]), ]
e1 <- sum(d$postoperative_hypertension[d$group == 1] == 1)
e0 <- sum(d$postoperative_hypertension[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_postoperative_hypertension <- p_raw


row_postoperative_hypertension <- c("Postoperative hypertension",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Postoperative tachycardia
d <- itt[complete.cases(itt[, c("group", "postoperative_tachycardia")]), ]
e1 <- sum(d$postoperative_tachycardia[d$group == 1] == 1)
e0 <- sum(d$postoperative_tachycardia[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_postoperative_tachycardia <- p_raw


row_postoperative_tachycardia <- c("Postoperative tachycardia",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Chest tightness
d <- itt[complete.cases(itt[, c("group", "chest_tightness")]), ]
e1 <- sum(d$chest_tightness[d$group == 1] == 1)
e0 <- sum(d$chest_tightness[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_chest_tightness <- p_raw


row_chest_tightness <- c("Chest tightness",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Chest pain
d <- itt[complete.cases(itt[, c("group", "chest_pain")]), ]
e1 <- sum(d$chest_pain[d$group == 1] == 1)
e0 <- sum(d$chest_pain[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_chest_pain <- p_raw


row_chest_pain <- c("Chest pain",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Palpitations
d <- itt[complete.cases(itt[, c("group", "palpitations")]), ]
e1 <- sum(d$palpitations[d$group == 1] == 1)
e0 <- sum(d$palpitations[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_palpitations <- p_raw


row_palpitations <- c("Palpitations",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Dizziness
d <- itt[complete.cases(itt[, c("group", "dizziness")]), ]
e1 <- sum(d$dizziness[d$group == 1] == 1)
e0 <- sum(d$dizziness[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_dizziness <- p_raw


row_dizziness <- c("Dizziness",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Headache
d <- itt[complete.cases(itt[, c("group", "headache")]), ]
e1 <- sum(d$headache[d$group == 1] == 1)
e0 <- sum(d$headache[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_headache <- p_raw


row_headache <- c("Headache",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Urinary frequency urgency
d <- itt[complete.cases(itt[, c("group", "urinary_frequency_urgency")]), ]
e1 <- sum(d$urinary_frequency_urgency[d$group == 1] == 1)
e0 <- sum(d$urinary_frequency_urgency[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_urinary_frequency_urgency <- p_raw


row_urinary_frequency_urgency <- c("Urinary frequency urgency",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Dysuria
d <- itt[complete.cases(itt[, c("group", "dysuria")]), ]
e1 <- sum(d$dysuria[d$group == 1] == 1)
e0 <- sum(d$dysuria[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_dysuria <- p_raw


row_dysuria <- c("Dysuria",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Urinary tract infection
d <- itt[complete.cases(itt[, c("group", "urinary_tract_infection")]), ]
e1 <- sum(d$urinary_tract_infection[d$group == 1] == 1)
e0 <- sum(d$urinary_tract_infection[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_urinary_tract_infection <- p_raw


row_urinary_tract_infection <- c("Urinary tract infection",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Fever
d <- itt[complete.cases(itt[, c("group", "fever")]), ]
e1 <- sum(d$fever[d$group == 1] == 1)
e0 <- sum(d$fever[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_fever <- p_raw


row_fever <- c("Fever",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Chills
d <- itt[complete.cases(itt[, c("group", "chills")]), ]
e1 <- sum(d$chills[d$group == 1] == 1)
e0 <- sum(d$chills[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_chills <- p_raw


row_chills <- c("Chills",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Postoperative bleeding
d <- itt[complete.cases(itt[, c("group", "postoperative_bleeding")]), ]
e1 <- sum(d$postoperative_bleeding[d$group == 1] == 1)
e0 <- sum(d$postoperative_bleeding[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_postoperative_bleeding <- p_raw


row_postoperative_bleeding <- c("Postoperative bleeding",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


TableS3 <- as.data.frame(rbind(
  row_insomnia,
  row_electrolyte_disturbance,
  row_peripheral_oedema,
  row_abdominal_pain,
  row_abdominal_distension,
  row_constipation,
  row_increased_bowel_movements,
  row_dry_mouth,
  row_cough,
  row_sputum_production,
  row_sore_throat,
  row_hoarseness,
  row_postoperative_hypotension,
  row_postoperative_hypertension,
  row_postoperative_tachycardia,
  row_chest_tightness,
  row_chest_pain,
  row_palpitations,
  row_dizziness,
  row_headache,
  row_urinary_frequency_urgency,
  row_dysuria,
  row_urinary_tract_infection,
  row_fever,
  row_chills,
  row_postoperative_bleeding), stringsAsFactors = FALSE)

names(TableS3) <- c("Variable", arm_names, "P value")
cat("\n", "Supplementary table S3 | Postoperative complications and adverse events. Values are numbers (percentages).", "\n", sep = "")
print(TableS3, row.names = FALSE, right = FALSE)
write.csv(TableS3, "BMJ_results/TableS3.csv", row.names = FALSE, fileEncoding = "UTF-8")


# 11. Supplementary table S4: Per-protocol analysis -----------------------------

# PONV within 24 h
d <- pp[complete.cases(pp[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_pp_ponv_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_pp_ponv_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_pp_ponv_0_24h <- c("PONV within 24 h",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Nausea within 24 h
d <- pp[complete.cases(pp[, c("group", "nausea_0_24h")]), ]
e1 <- sum(d$nausea_0_24h[d$group == 1] == 1)
e0 <- sum(d$nausea_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_pp_nausea_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_pp_nausea_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_pp_nausea_0_24h <- c("Nausea within 24 h",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Vomiting within 24 h
d <- pp[complete.cases(pp[, c("group", "vomiting_0_24h")]), ]
e1 <- sum(d$vomiting_0_24h[d$group == 1] == 1)
e0 <- sum(d$vomiting_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_pp_vomiting_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_pp_vomiting_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_pp_vomiting_0_24h <- c("Vomiting within 24 h",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Clinically important PONV within 24 h
d <- pp[complete.cases(pp[, c("group", "clinically_important_ponv_0_24h")]), ]
e1 <- sum(d$clinically_important_ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$clinically_important_ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_pp_clinically_important_ponv_0_24h <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_pp_clinically_important_ponv_0_24h <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_pp_clinically_important_ponv_0_24h <- c("Clinically important PONV within 24 h",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


# Rescue antiemetic use
d <- pp[complete.cases(pp[, c("group", "rescue_antiemetic_use")]), ]
e1 <- sum(d$rescue_antiemetic_use[d$group == 1] == 1)
e0 <- sum(d$rescue_antiemetic_use[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_pp_rescue_antiemetic_use <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_pp_rescue_antiemetic_use <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_pp_rescue_antiemetic_use <- c("Rescue antiemetic use",
  sprintf("%d/%d (%.1f)", e1, n1, 100*e1/n1),
  sprintf("%d/%d (%.1f)", e0, n0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


TableS4 <- as.data.frame(rbind(
  row_pp_ponv_0_24h,
  row_pp_nausea_0_24h,
  row_pp_vomiting_0_24h,
  row_pp_clinically_important_ponv_0_24h,
  row_pp_rescue_antiemetic_use), stringsAsFactors = FALSE)

names(TableS4) <- c("Outcome", sprintf("Esomeprazole group (n=%d)", sum(pp$group == 1)), sprintf("Control group (n=%d)", sum(pp$group == 0)), "Risk ratio (95% CI)", "P value")
cat("\n", "Supplementary table S4 | Per-protocol analysis of primary and key secondary outcomes. Values are events/participants (percentages).", "\n", sep = "")
print(TableS4, row.names = FALSE, right = FALSE)
write.csv(TableS4, "BMJ_results/TableS4.csv", row.names = FALSE, fileEncoding = "UTF-8")


# 12. Supplementary table S5: Multiplicity adjustment ---------------------------
# Use the full-precision P values for all four PONV time windows.
p_windows <- c(p_ponv_0_1h, p_ponv_1_6h, p_ponv_6_24h, p_ponv_0_24h)
p_bonferroni <- p.adjust(p_windows, method = "bonferroni")
p_fdr <- p.adjust(p_windows, method = "BH")


TableS5 <- data.frame(
  "PONV time window" = c("0–1 h", "1–6 h", "6–24 h", "0–24 h (primary)"),
  "Original P value" = ifelse(is.na(p_windows), "NA", ifelse(p_windows < 0.001, "<0.001", sprintf("%.3f", p_windows))),
  "Bonferroni adjusted P value" = ifelse(is.na(p_bonferroni), "NA", ifelse(p_bonferroni < 0.001, "<0.001", sprintf("%.3f", p_bonferroni))),
  "FDR adjusted P value" = ifelse(is.na(p_fdr), "NA", ifelse(p_fdr < 0.001, "<0.001", sprintf("%.3f", p_fdr))), check.names = FALSE)
cat("\nSupplementary table S5 | Multiplicity adjustment for PONV across postoperative time windows.\n")
print(TableS5, row.names = FALSE, right = FALSE)
write.csv(TableS5, "BMJ_results/TableS5.csv", row.names = FALSE, fileEncoding = "UTF-8")

# 13. Supplementary table S6: Additional sensitivity analyses -------------------

sub_data <- subset(itt, apfel_score == 3)

# Apfel score 3 subgroup
d <- sub_data[complete.cases(sub_data[, c("group", "ponv_0_24h")]), ]
e1 <- sum(d$ponv_0_24h[d$group == 1] == 1)
e0 <- sum(d$ponv_0_24h[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_apfel3 <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_apfel3 <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_apfel3 <- c("Apfel score 3 subgroup",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


itt$ponv_or_rescue <- as.integer(itt$ponv_0_24h == 1 | itt$rescue_antiemetic_use == 1)

# Composite endpoint: PONV or rescue antiemetic
d <- itt[complete.cases(itt[, c("group", "ponv_or_rescue")]), ]
e1 <- sum(d$ponv_or_rescue[d$group == 1] == 1)
e0 <- sum(d$ponv_or_rescue[d$group == 0] == 1)
n1 <- sum(d$group == 1)
n0 <- sum(d$group == 0)
cells <- matrix(c(e1, n1-e1, e0, n0-e0), nrow = 2, byrow = TRUE)
p_raw <- 1
if (all(colSums(cells) > 0)) {
  chi <- suppressWarnings(chisq.test(cells, correct = FALSE))
  p_raw <- if (any(chi$expected < 5)) fisher.test(cells)$p.value else chi$p.value
}
p_composite <- p_raw


rr <- lower <- upper <- NA_real_
if (e1 > 0 && e0 > 0) {
  rr <- (e1/n1)/(e0/n0)
  se_log <- sqrt(1/e1 - 1/n1 + 1/e0 - 1/n0)
  lower <- exp(log(rr) - qnorm(0.975)*se_log)
  upper <- exp(log(rr) + qnorm(0.975)*se_log)
}
rr_text <- if (is.na(rr)) "NE" else sprintf("%.2f (%.2f to %.2f)", rr, lower, upper)


raw_composite <- c(e1=e1, n1=n1, e0=e0, n0=n0, rr=rr, lower=lower, upper=upper, p=p_raw)

row_composite <- c("Composite endpoint: PONV or rescue antiemetic",
  sprintf("%d (%.1f)", e1, 100*e1/n1),
  sprintf("%d (%.1f)", e0, 100*e0/n0),
  rr_text,
  ifelse(is.na(p_raw), "NA", ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw))))


r <- adj_ponv_0_24h
s6_adjusted <- sprintf("Adjusted RR %.3f (95%% CI %.3f to %.3f); P=%s",
  r["rr"], r["lower"], r["upper"], ifelse(is.na(r["p"]), "NA", ifelse(r["p"] < 0.001, "<0.001", sprintf("%.3f", r["p"]))))

r <- raw_apfel3
s6_apfel <- sprintf("RR %.3f (95%% CI %.3f to %.3f); P=%s",
  r["rr"], r["lower"], r["upper"], ifelse(is.na(r["p"]), "NA", ifelse(r["p"] < 0.001, "<0.001", sprintf("%.3f", r["p"]))))

r <- raw_composite
s6_composite <- sprintf("RR %.3f (95%% CI %.3f to %.3f); P=%s",
  r["rr"], r["lower"], r["upper"], ifelse(is.na(r["p"]), "NA", ifelse(r["p"] < 0.001, "<0.001", sprintf("%.3f", r["p"]))))

TableS6 <- as.data.frame(rbind(
  c("Adjusted Poisson regression with robust variance", s6_adjusted),
  c("Apfel score 3 subgroup", s6_apfel),
  c("Composite endpoint: PONV or rescue antiemetic", s6_composite)), stringsAsFactors = FALSE)

names(TableS6) <- c("Analysis", "Result")
cat("\n", "Supplementary table S6 | Additional sensitivity and risk-stratified analyses.", "\n", sep = "")
print(TableS6, row.names = FALSE, right = FALSE)
write.csv(TableS6, "BMJ_results/TableS6.csv", row.names = FALSE, fileEncoding = "UTF-8")


# 14. Print any table again and record the software environment -----------------
# Examples (run after the corresponding table has been created):
# print(Table2, row.names = FALSE, right = FALSE)
# print(Table3, row.names = FALSE, right = FALSE)
# print(Table4, row.names = FALSE, right = FALSE)
# print(TableS1, row.names = FALSE, right = FALSE)
# print(TableS2, row.names = FALSE, right = FALSE)
# print(TableS3, row.names = FALSE, right = FALSE)
# print(TableS4, row.names = FALSE, right = FALSE)
# print(TableS5, row.names = FALSE, right = FALSE)
# print(TableS6, row.names = FALSE, right = FALSE)
# Formatted table cells are text; use raw_* and adj_* objects for calculations.
capture.output(sessionInfo(), file = "BMJ_results/sessionInfo.txt")

