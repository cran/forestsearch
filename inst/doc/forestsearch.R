## ----setup--------------------------------------------------------------------
#| message: false
#| warning: false
knitr::opts_chunk$set(comment = "#>", fig.align = "center")
library(forestsearch)
library(survival)


## ----data---------------------------------------------------------------------
df <- survival::gbsg
df <- within(df, {
  id          <- seq_len(nrow(df))              # the fit's subject identifier
  time_months <- rfstime / 30.4375              # recurrence-free survival, days to months
  grade3      <- ifelse(grade == "3", 1, 0)     # grade 3 against grade 1 or 2
})
confounders <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
stopifnot(nrow(df) == 686, sum(df$status) == 299)

whole <- coxph(Surv(time_months, status) ~ hormon, data = df)
whole_hr <- c(est = unname(exp(coef(whole))), exp(confint(whole))[1, ])
c(patients = nrow(df), hormonal_therapy = sum(df$hormon == 1),
  none = sum(df$hormon == 0), events = sum(df$status))
sprintf("Whole-trial hazard ratio %.3f (%.3f, %.3f), hormonal therapy vs none",
        whole_hr[["est"]], whole_hr[[2]], whole_hr[[3]])


## ----fit-fs-------------------------------------------------------------------
#| message: false
#| results: hide
t0 <- proc.time()[["elapsed"]]
fit_fs <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  # the candidate family
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  # the identifier, its thresholds and the selection rule
  subgroup_method = "consistency", consistency_method = "resample",
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  stop_threshold = NULL,                      # evaluate every candidate; no early stop
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  quiet = TRUE, parallel_args = list(plan = "sequential"))
fit_sec_fs <- proc.time()[["elapsed"]] - t0


## ----fs-selected--------------------------------------------------------------
in_H <- fit_fs$df.est$treat.recommend[match(df$id, fit_fs$df.est$id)] == 0
stopifnot(!anyNA(in_H))
cat(sprintf("Identified subgroup: %s\n", paste(fit_fs$sg.harm, collapse = " & ")))
cat(sprintf("Patients in it: %d of %d. Fit: %.1f s\n", sum(in_H), nrow(df), fit_sec_fs))


## ----fs-funnel----------------------------------------------------------------
counts <- fit_fs$find.grps$filter_counts
funnel <- data.frame(
  stage = c("enumerated (at most two cuts, size-eligible)",
            "fittable (ten events in each arm, Cox fit obtained)",
            "clearing the screening threshold, HR > 1.0",
            "passing the consistency screen, p_cons >= 0.90"),
  candidates = c(counts$n_evaluated, counts$n_passed_cox, counts$n_passed_hr,
                 fit_fs$grp.consistency$n_passed))
knitr::kable(funnel, row.names = FALSE, format.args = list(big.mark = ","))


## ----fs-candidates------------------------------------------------------------
cand <- as.data.frame(fit_fs$grp.consistency$out_sg$result)
cand$candidate <- ifelse(cand$M.2 == "", cand$M.1, paste(cand$M.1, cand$M.2, sep = " & "))
cand <- cand[order(-cand$N), ]
knitr::kable(cand[, c("candidate", "N", "hr", "Pcons")], digits = 3,
             row.names = FALSE, col.names = c("candidate", "N", "HR", "p_cons"),
             caption = "The candidates passing the consistency screen, largest first.")
band <- (1 - 0.20) * max(cand$hr)
within_band <- cand$hr >= band
sel <- cand[within_band, ][which.max(cand$N[within_band]), ]
cat(sprintf("Largest HR over the family: %.2f (N %d). Band: HR >= %.2f.\n",
            max(cand$hr), cand$N[which.max(cand$hr)], band))
cat(sprintf("max N(eps) takes the largest candidate within the band: N %d, HR %.2f.\n",
            sel$N, sel$hr))
cat(sprintf("max N alone would take N %d (HR %.2f); max beta the N %d.\n",
            max(cand$N), cand$hr[which.max(cand$N)], cand$N[which.max(cand$hr)]))
stopifnot(sel$N == sum(in_H))   # the rule, re-applied, picks the fit's own selection


## ----unadjusted---------------------------------------------------------------
hr_ci <- function(d) {
  f <- coxph(Surv(time_months, status) ~ hormon, data = d)
  ci <- exp(confint(f))
  sprintf("%.3f (%.3f, %.3f)", exp(coef(f)), ci[1, 1], ci[1, 2])
}
knitr::kable(
  data.frame(set = c("identified region", "complement"),
             N = c(sum(in_H), sum(!in_H)),
             events = c(sum(df$status[in_H]), sum(df$status[!in_H])),
             unadjusted = c(hr_ci(df[in_H, ]), hr_ci(df[!in_H, ]))),
  col.names = c("", "N", "events", "unadjusted HR (95% CI)"), row.names = FALSE,
  caption = "The identified region and its complement: unadjusted estimates, before any correction for the selection.")


## ----fit-fs-inference---------------------------------------------------------
#| eval: false
# fit_fs_inf <- forestsearch(
#   df.analysis = df,
#   outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
#   confounders.name = confounders,
#   is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
#   conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
#   n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
#   subgroup_method = "consistency", consistency_method = "resample",
#   hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
#   stop_threshold = NULL,
#   sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
#   # selection-adjusted inference
#   mr_inference = TRUE,
#   mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
#                            field_scale_complement = "selected"),
#   quiet = TRUE, parallel_args = list(plan = "sequential"))
# mr_estimates_table(fit_fs_inf, est.scale = "hr", digits = 3)


## ----session------------------------------------------------------------------
sessionInfo()

