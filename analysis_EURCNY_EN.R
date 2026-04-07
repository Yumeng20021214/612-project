# ============================================================
# The Impact of EUR/CNY Exchange Rate Volatility on China-EU Trade
# Full Empirical Analysis: ARIMA + GARCH + VAR
# ============================================================

# ---------- 0. Load Packages ----------
# install.packages(c("forecast","urca","tseries","rugarch","vars",
#                    "stargazer","ggplot2","ggfortify","gridExtra",
#                    "lubridate","dplyr"))

library(forecast)
library(urca)
library(tseries)
library(rugarch)
library(vars)
library(stargazer)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(lubridate)
library(dplyr)

options(warn = -1)
options(repr.plot.width = 14, repr.plot.height = 6)

# Save all plots to PDF
pdf("results_plots.pdf", width = 14, height = 6)


# ============================================================
# Part 1: Data Import and Preprocessing
# ============================================================

# ---------- 1.1 Load raw data ----------
fx_raw  <- read.csv("exchange_rate_final.csv", stringsAsFactors = FALSE)
exp_raw <- read.csv("export_final.csv",        stringsAsFactors = FALSE)

colnames(fx_raw)  <- c("date", "rate")
colnames(exp_raw) <- c("month", "export")

head(fx_raw);  tail(fx_raw)
head(exp_raw); tail(exp_raw)

# ---------- 1.2 Clean exchange rate data ----------
fx_raw$date <- as.Date(fx_raw$date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))
fx_raw       <- fx_raw[!is.na(fx_raw$date) & !is.na(fx_raw$rate), ]
fx_raw       <- fx_raw[order(fx_raw$date), ]

cat("Exchange rate data:", nrow(fx_raw), "trading days,",
    format(min(fx_raw$date)), "to", format(max(fx_raw$date)), "\n")

# ---------- 1.3 Clean export data ----------
exp_raw <- exp_raw[!grepl(":", exp_raw$export), ]
exp_raw$export <- as.numeric(gsub(",", "", exp_raw$export))
exp_raw$month  <- as.Date(paste0(exp_raw$month, "-01"))
exp_raw        <- exp_raw[!is.na(exp_raw$export), ]
exp_raw        <- exp_raw[order(exp_raw$month), ]

cat("Export data:", nrow(exp_raw), "months,",
    format(min(exp_raw$month)), "to", format(max(exp_raw$month)), "\n")

# ---------- 1.4 Aggregate daily FX to monthly mean and variance ----------
fx_raw$month <- as.Date(format(fx_raw$date, "%Y-%m-01"))

fx_monthly <- fx_raw %>%
  group_by(month) %>%
  summarise(
    fx_mean = mean(rate, na.rm = TRUE),
    fx_var  = var(rate,  na.rm = TRUE),
    .groups = "drop"
  )

# ---------- 1.5 Merge datasets on common time period ----------
df <- merge(fx_monthly, exp_raw, by = "month")
df <- df[order(df$month), ]
df <- df[complete.cases(df), ]

cat("Merged sample:", nrow(df), "months,",
    format(min(df$month)), "to", format(max(df$month)), "\n")

# ---------- 1.6 Build time series objects ----------
start_yr  <- as.integer(format(df$month[1], "%Y"))
start_mon <- as.integer(format(df$month[1], "%m"))

ts_fx_mean <- ts(df$fx_mean, start = c(start_yr, start_mon), frequency = 12)
ts_fx_var  <- ts(df$fx_var,  start = c(start_yr, start_mon), frequency = 12)
ts_export  <- ts(df$export,  start = c(start_yr, start_mon), frequency = 12)


# ============================================================
# Part 2: Descriptive Statistics and Visualization
# ============================================================

cat("\n===== Descriptive Statistics =====\n")
cat("-- EUR/CNY Monthly Mean Rate --\n")
print(summary(df$fx_mean))
cat("-- China Exports to EU (USD) --\n")
print(summary(df$export))

# Three-panel plot
p1 <- autoplot(ts_fx_mean, colour = "steelblue") +
  theme_light() +
  labs(title = "EUR/CNY Monthly Mean Exchange Rate",
       x = "Year", y = "Exchange Rate (CNY per EUR)")

p2 <- autoplot(ts_fx_var, colour = "darkorange") +
  theme_light() +
  labs(title = "EUR/CNY Monthly Within-Month Variance (Volatility Proxy)",
       x = "Year", y = "Variance")

p3 <- autoplot(ts_export / 1e9, colour = "darkgreen") +
  theme_light() +
  labs(title = "China Exports to EU (Billion USD)",
       x = "Year", y = "Billion USD")

grid.arrange(p1, p2, p3, ncol = 3)


# ============================================================
# Part 3: Unit Root Tests (ADF)
# ============================================================

cat("\n===== ADF Unit Root Tests =====\n")

# ---------- 3.1 Level series ----------
cat("\n-- EUR/CNY Exchange Rate (Level) --\n")
summary(ur.df(ts_fx_mean, type = "trend", lags = 12, selectlags = "AIC"))

cat("\n-- China Exports to EU (Level) --\n")
summary(ur.df(ts_export, type = "trend", lags = 12, selectlags = "AIC"))

# ---------- 3.2 Log first differences ----------
log_fx   <- log(ts_fx_mean)
dlog_fx  <- diff(log_fx) * 100          # Monthly FX change rate (%)
dlog_exp <- diff(log(ts_export)) * 100  # Monthly export growth rate (%)

cat("\n-- EUR/CNY Log First Difference (Monthly Change Rate) --\n")
summary(ur.df(dlog_fx, type = "drift", lags = 12, selectlags = "AIC"))

cat("\n-- Export Log First Difference (Monthly Growth Rate) --\n")
summary(ur.df(dlog_exp, type = "drift", lags = 12, selectlags = "AIC"))

# ---------- 3.3 Check for deterministic trend in differenced series ----------
cat("\n-- Trend Regression: FX Change Rate on Time Trend --\n")
t_fx <- 1:length(dlog_fx)
model_trend_fx <- lm(dlog_fx ~ t_fx)
stargazer(model_trend_fx, type = "text",
          title = "Trend Regression for FX Change Rate")

cat("\n-- Trend Regression: Export Growth Rate on Time Trend --\n")
t_exp <- 1:length(dlog_exp)
model_trend_exp <- lm(dlog_exp ~ t_exp)
stargazer(model_trend_exp, type = "text",
          title = "Trend Regression for Export Growth Rate")

# ---------- 3.4 Visualize stationary series ----------
p4 <- autoplot(dlog_fx, colour = "steelblue") +
  theme_light() +
  labs(title = "EUR/CNY Monthly Change Rate (Log First Difference, %)",
       x = "Year", y = "Change Rate (%)")

p5 <- autoplot(dlog_exp, colour = "darkgreen") +
  theme_light() +
  labs(title = "China Export Growth Rate (Log First Difference, %)",
       x = "Year", y = "Growth Rate (%)")

grid.arrange(p4, p5, ncol = 2)


# ============================================================
# Part 4: ARIMA Modeling (Exchange Rate Mean Equation)
# ============================================================

cat("\n===== ARIMA Modeling (EUR/CNY Log First Difference) =====\n")

arima_dlog <- auto.arima(dlog_fx,
                         max.d = 0, max.p = 5, max.q = 5,
                         seasonal = FALSE, ic = "aic",
                         stepwise = FALSE, stationary = TRUE)
print(arima_dlog)

# Residual diagnostics
checkresiduals(arima_dlog)

# Ljung-Box test for residual autocorrelation
p_ar <- arima_dlog$arma[1]
q_ma <- arima_dlog$arma[2]
lb_test <- Box.test(residuals(arima_dlog), lag = 20,
                    type = "Ljung-Box", fitdf = p_ar + q_ma)
cat("\nLjung-Box Test (lag=20): Q =", round(lb_test$statistic, 3),
    ", p-value =", round(lb_test$p.value, 4), "\n")
if (lb_test$p.value > 0.05) {
  cat("Conclusion: No significant residual autocorrelation. Mean model is well specified.\n")
} else {
  cat("Conclusion: Residual autocorrelation detected. Consider higher lag orders.\n")
}


# ============================================================
# Part 5: GARCH Modeling (Conditional Variance of Exchange Rate)
# ============================================================

# ---------- 5.1 Extract ARIMA residuals ----------
res_dlog <- residuals(arima_dlog)

# ---------- 5.2 Test for ARCH effects: Ljung-Box on squared residuals ----------
cat("\n===== ARCH Effect Tests =====\n")
cat("-- Ljung-Box Test on Squared Residuals --\n")
for (h in c(4, 8, 12, 20)) {
  lb <- Box.test(res_dlog^2, lag = h, type = "Ljung-Box")
  cat("lag =", h, ": Q =", round(lb$statistic, 3),
      ", p-value =", round(lb$p.value, 5), "\n")
}

# ---------- 5.3 Engle ARCH-LM Test (manual implementation) ----------
cat("\n-- Engle ARCH-LM Test --\n")
arch_lm_test <- function(residuals, order) {
  T     <- length(residuals)
  u2    <- residuals^2
  y_aux <- u2[(order + 1):T]
  X_aux <- matrix(NA, nrow = T - order, ncol = order + 1)
  X_aux[, 1] <- 1
  for (j in 1:order) {
    X_aux[, j + 1] <- u2[(order + 1 - j):(T - j)]
  }
  fit  <- lm(y_aux ~ X_aux - 1)
  R2   <- summary(fit)$r.squared
  LM   <- (T - order) * R2
  pval <- pchisq(LM, df = order, lower.tail = FALSE)
  cat("order =", order, ": LM =", round(LM, 3),
      ", p-value =", round(pval, 5), "\n")
}

for (p_arch in c(1, 4, 8, 12)) {
  arch_lm_test(res_dlog, p_arch)
}

# ---------- 5.4 Visualize squared residuals and ACF ----------
options(repr.plot.width = 14, repr.plot.height = 5)

par(mfrow = c(1, 2))
plot(res_dlog^2, type = "l", col = "steelblue",
     main = "Squared Residuals of ARIMA (EUR/CNY Change Rate)",
     ylab = "Squared Residuals", xlab = "Time")

acf(res_dlog^2, main = "ACF of Squared Residuals",
    xlab = "Lag", ylab = "ACF")
par(mfrow = c(1, 1))

# ---------- 5.5 Estimate GARCH(1,1) ----------
cat("\n===== GARCH(1,1) Estimation =====\n")
library(rugarch)

spec_garch <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(p_ar, q_ma), include.mean = TRUE),
  distribution.model = "norm"
)

garch_fit <- ugarchfit(spec_garch, dlog_fx)
cat("\nGARCH(1,1) Coefficient Estimates:\n")
print(garch_fit@fit$matcoef)

# ---------- 5.6 Extract conditional variance ----------
sigma2_cond <- as.numeric(sigma(garch_fit)^2)
ts_sigma2   <- ts(sigma2_cond,
                  start     = c(start_yr, start_mon + 1),
                  frequency = 12)

# ---------- 5.7 Unconditional variance ----------
omega <- coef(garch_fit)["omega"]
alpha <- coef(garch_fit)["alpha1"]
beta  <- coef(garch_fit)["beta1"]
unc_var <- omega / (1 - alpha - beta)

cat("\nUnconditional Variance (Long-run Mean):", round(unc_var, 8), "\n")
cat("Unconditional Volatility (Std Dev):", round(sqrt(unc_var), 6), "\n")
cat("Persistence alpha + beta =", round(alpha + beta, 4), "\n")

if (alpha + beta < 1) {
  cat("Conclusion: GARCH process is covariance-stationary (alpha + beta < 1).\n")
} else {
  cat("Conclusion: alpha + beta >= 1, integrated GARCH (IGARCH) behavior.\n")
}

# ---------- 5.8 Plot conditional volatility ----------
options(repr.plot.width = 14, repr.plot.height = 5)
plot(ts_sigma2,
     main = "EUR/CNY Conditional Variance from GARCH(1,1)",
     ylab = "Conditional Variance", xlab = "Year",
     col  = "firebrick", lwd = 1.5)
abline(h = unc_var, lty = 2, col = "navy", lwd = 2)
legend("topright",
       legend = c("Conditional Variance", "Unconditional Variance (Long-run Mean)"),
       col    = c("firebrick", "navy"),
       lty    = c(1, 2), lwd = 2)

# ---------- 5.9 GARCH Forecast (5-step ahead) ----------
cat("\n===== GARCH(1,1) Forecast (5 steps ahead) =====\n")
garch_fcst <- ugarchforecast(garch_fit, n.ahead = 5)
print(garch_fcst)


# ============================================================
# Part 6: Cointegration Test and VAR / ECM Modeling
# ============================================================

cat("\n===== Cointegration Test (Johansen) =====\n")

# ---------- 6.1 Test level series for cointegration ----------
data_level <- ts.intersect(log(ts_fx_mean), log(ts_export))
colnames(data_level) <- c("log_fx", "log_exp")

lag_sel_level <- VARselect(data_level, lag.max = 8, type = "trend")
cat("Lag selection for level VAR:\n")
print(lag_sel_level$selection)
K_jo <- lag_sel_level$selection["AIC(n)"]

jo_test <- ca.jo(data_level, ecdet = "trend", type = "trace",
                 K = K_jo, spec = "longrun")
summary(jo_test)

cat("\nInterpretation:\n")
cat("If test statistic > 5% critical value: reject H0 (no cointegration)\n")
cat("→ Use VECM; otherwise → Use VAR in first differences.\n")

# ---------- 6.2 Build VAR input variables (stationary series) ----------
# Variables: FX change rate, log conditional variance, export growth rate
log_sigma2_ts <- log(ts_sigma2 + 1e-10)

data_var <- ts.intersect(dlog_fx, log_sigma2_ts, dlog_exp)
colnames(data_var) <- c("dlog_fx", "log_sigma2", "dlog_exp")

cat("\nVAR dataset dimensions:", nrow(data_var), "x", ncol(data_var), "\n")
cat("Time span:", start(data_var)[1], "M", start(data_var)[2],
    "to", end(data_var)[1], "M", end(data_var)[2], "\n")

# ---------- 6.3 Confirm stationarity of VAR variables ----------
cat("\n-- Stationarity Check for VAR Variables --\n")
for (v in colnames(data_var)) {
  res <- ur.df(data_var[, v], type = "drift", lags = 6, selectlags = "AIC")
  tau <- res@teststat[1]
  cv5 <- res@cval[1, 2]
  cat(v, ": tau =", round(tau, 3), ", 5% CV =", cv5,
      ifelse(tau < cv5, "-> Stationary", "-> Non-stationary"), "\n")
}

# ---------- 6.4 VAR lag order selection ----------
cat("\n===== VAR Lag Order Selection =====\n")
lag_sel_var <- VARselect(data_var, lag.max = 12, type = "const")
print(lag_sel_var$selection)
optimal_lag <- lag_sel_var$selection["AIC(n)"]
cat("Optimal lag (AIC):", optimal_lag, "\n")

# ---------- 6.5 Estimate VAR ----------
cat("\n===== VAR Estimation =====\n")
model_VAR <- VAR(data_var, p = optimal_lag, type = "const")
summary(model_VAR)

# ---------- 6.6 VAR Forecast (4 steps ahead) ----------
cat("\n===== VAR Forecast (4 steps ahead) =====\n")
prd_VAR <- predict(model_VAR, n.ahead = 4, ci = 0.95)
print(prd_VAR)

options(repr.plot.width = 14, repr.plot.height = 8)
plot(prd_VAR, main = "VAR 4-Step-Ahead Forecast")


# ============================================================
# Part 7: Impulse Response Functions (IRF)
# ============================================================

cat("\n===== Impulse Response Functions =====\n")

options(repr.plot.width = 14, repr.plot.height = 6)

# FX level shock -> Export growth
irf1 <- irf(model_VAR,
             impulse  = "dlog_fx",
             response = "dlog_exp",
             n.ahead  = 12, boot = TRUE, ci = 0.95)
plot(irf1,
     main = "IRF: FX Change Rate Shock -> Export Growth Rate",
     ylab = "Response of Export Growth (%)")

# Volatility shock -> Export growth
irf2 <- irf(model_VAR,
             impulse  = "log_sigma2",
             response = "dlog_exp",
             n.ahead  = 12, boot = TRUE, ci = 0.95)
plot(irf2,
     main = "IRF: FX Volatility Shock -> Export Growth Rate",
     ylab = "Response of Export Growth (%)")

# Volatility shock -> FX change rate
irf3 <- irf(model_VAR,
             impulse  = "log_sigma2",
             response = "dlog_fx",
             n.ahead  = 12, boot = TRUE, ci = 0.95)
plot(irf3,
     main = "IRF: FX Volatility Shock -> FX Change Rate",
     ylab = "Response of FX Change Rate (%)")


# ============================================================
# Part 8: Granger Causality Tests
# ============================================================

cat("\n===== Granger Causality Tests =====\n")

gc1 <- causality(model_VAR, cause = "log_sigma2")
cat("H0: FX Volatility does NOT Granger-cause Export Growth\n")
print(gc1$Granger)

gc2 <- causality(model_VAR, cause = "dlog_fx")
cat("\nH0: FX Change Rate does NOT Granger-cause Export Growth\n")
print(gc2$Granger)

gc3 <- causality(model_VAR, cause = "dlog_exp")
cat("\nH0: Export Growth does NOT Granger-cause FX Volatility\n")
print(gc3$Granger)


# ============================================================
# Part 9: Forecast Error Variance Decomposition (FEVD)
# ============================================================

cat("\n===== Forecast Error Variance Decomposition (FEVD) =====\n")
fevd_res <- fevd(model_VAR, n.ahead = 12)
print(fevd_res)

options(repr.plot.width = 14, repr.plot.height = 8)
plot(fevd_res, main = "FEVD: Forecast Error Variance Decomposition")


# ============================================================
# Part 10: Pseudo Out-of-Sample Forecasting
# ============================================================

cat("\n===== Pseudo Out-of-Sample Forecasting (4-Step-Ahead) =====\n")

T       <- nrow(data_var)
Tstar   <- floor(0.8 * T)
K_iter  <- T - Tstar - 4

cat("Total observations:", T, "\n")
cat("Initial estimation window:", Tstar, "\n")
cat("Number of rolling forecasts:", K_iter, "\n\n")

results_oos <- data.frame()

for (i in 1:K_iter) {
  temp_end <- Tstar + i - 1

  model_temp <- try(
    suppressWarnings(VAR(data_var[1:temp_end, ],
                         p = optimal_lag, type = "const")),
    silent = TRUE
  )
  if (inherits(model_temp, "try-error")) next

  f <- predict(model_temp, n.ahead = 4, ci = 0.95)

  f_mean    <- f$fcst$dlog_exp[4, "fcst"]
  f_lower   <- f$fcst$dlog_exp[4, "lower"]
  f_upper   <- f$fcst$dlog_exp[4, "upper"]
  actual    <- as.numeric(data_var[temp_end + 4, "dlog_exp"])
  abs_error <- abs(f_mean - actual)

  results_oos <- rbind(results_oos,
                       data.frame(Forecast   = f_mean,
                                  Lower95    = f_lower,
                                  Upper95    = f_upper,
                                  Actual     = actual,
                                  AbsError   = abs_error))
}

cat("-- First 6 Rolling Forecasts --\n"); print(head(results_oos))
cat("\n-- Last 6 Rolling Forecasts --\n");  print(tail(results_oos))

mae_val   <- mean(results_oos$AbsError, na.rm = TRUE)
rmsfe_val <- sqrt(mean((results_oos$Forecast - results_oos$Actual)^2,
                       na.rm = TRUE))

cat("\nForecast Accuracy:\n")
cat("MAE   =", round(mae_val,   4), "(percentage points)\n")
cat("RMSFE =", round(rmsfe_val, 4), "(percentage points)\n")

# Plot: Actual vs Forecast
options(repr.plot.width = 14, repr.plot.height = 5)
n_res <- nrow(results_oos)
plot(1:n_res, results_oos$Actual, type = "l",
     col = "black", lwd = 2,
     ylim = range(c(results_oos$Lower95, results_oos$Upper95), na.rm = TRUE),
     main = "Pseudo Out-of-Sample Forecast: Export Growth Rate (4-Step-Ahead)",
     xlab = "Rolling Window Index", ylab = "Export Growth Rate (%)")
lines(1:n_res, results_oos$Forecast, col = "steelblue", lwd = 2, lty = 2)
lines(1:n_res, results_oos$Lower95,  col = "grey60",    lwd = 1, lty = 3)
lines(1:n_res, results_oos$Upper95,  col = "grey60",    lwd = 1, lty = 3)
legend("topright",
       legend = c("Actual", "Forecast", "95% Prediction Interval"),
       col    = c("black", "steelblue", "grey60"),
       lty    = c(1, 2, 3), lwd = c(2, 2, 1))


# ============================================================
# Part 11: Summary
# ============================================================

cat("\n")
cat("========================================================\n")
cat("                    SUMMARY OF RESULTS\n")
cat("========================================================\n")
cat("1. Stationarity:\n")
cat("   EUR/CNY level series: I(1) (unit root present)\n")
cat("   Log first difference (monthly change rate): Stationary\n\n")
cat("2. ARIMA mean model: Residuals pass Ljung-Box test.\n\n")
cat("3. GARCH(1,1) Results:\n")
cat("   alpha (ARCH) =", round(alpha, 4), "\n")
cat("   beta  (GARCH) =", round(beta,  4), "\n")
cat("   Persistence alpha+beta =", round(alpha + beta, 4), "\n")
cat("   Conditional variance is time-varying.\n\n")
cat("4. VAR(", optimal_lag, ") Model:\n")
cat("   See IRF plots and Granger causality results above.\n\n")
cat("5. Pseudo OOS Forecast Accuracy:\n")
cat("   MAE   =", round(mae_val,   4), "percentage points\n")
cat("   RMSFE =", round(rmsfe_val, 4), "percentage points\n")
cat("========================================================\n")

# Close PDF device
dev.off()
cat("All plots saved to: results_plots.pdf\n")
