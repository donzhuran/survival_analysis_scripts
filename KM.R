install.packages(c("survival", "survminer"))
library("survival")
library("survminer")
View(lung)
head(lung)
?lung

#The function survfit() can be used to compute kaplan-Meier survival estimate.
fit <- survfit(Surv(time, status) ~ sex, data = lung) # calculating survival between time and event (status); comparing data grouped by sex
print(fit)
summary(fit)

res.sum <- surv_summary(fit)
head(res.sum)

# Access to the sort summary table
summary(fit)$table
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
head(d)


# Statistics
#Log-Rank test comparing survival curves: survdiff()
surv_diff <- survdiff(Surv(time, status) ~ sex, data = lung)
surv_diff # there is statistically significant difference

# Visualization
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval=TRUE, conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           xlim = c(0, 600))

ggsurvplot(fit,
           pval=TRUE, conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           fun = "event")

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


surv_object <- Surv(time, status, data = lung)
fit <- survfit(surv_object~sex)
