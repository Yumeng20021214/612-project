# 0. Load packages
library(readr)
library(dplyr)
library(ggplot2)
library(tseries)
install.packages("FinTS")
library(FinTS)
install.packages("forecast")
library(forecast)
install.packages("lubridate")
library(lubridate)
install.packages("rugarch")
library(rugarch)
library(tseries)
install.packages("vars")
library(vars)

# 1. Import data：exchange rate(daily)
exchange_rate <- read_csv("Downloads/exchange_rate_final.csv")
colnames(exchange_rate) <- c("date", "rate")
str(exchange_rate)
exchange_rate$date <- as.Date(exchange_rate$date, format = "%m/%d/%Y")

# 2. Plot exchange rate level: nonstationary

ggplot(exchange_rate, aes(x = date, y = rate)) +
  geom_line() +
  labs(title = "EUR/CNY Exchange Rate Level",
       x = "Date",
       y = "Exchange Rate") +
  theme_minimal()

# 3. Construct log return
# 收益率r_t = log(e_t) - log(e_{t-1})

exchange_rate <- exchange_rate %>%
  mutate(
    log_rate = log(rate),
    return = log_rate - lag(log_rate)
  )

# 删除第一行 NA
exchange_rate_ret <- exchange_rate %>%
  filter(!is.na(return))

# 4. Plot log return：return大致在0附近波动，没有明显趋势但不同时段波动幅度不同

ggplot(exchange_rate_ret, aes(x = date, y = return)) +
  geom_line() +
  labs(title = "Log Returns of EUR/CNY",
       x = "Date",
       y = "Log Return") +
  theme_minimal()

# 5. Plot squared return
# 把return平方，去掉正负方向只保留波动强弱，用来看 volatility clustering

exchange_rate_ret <- exchange_rate_ret %>%
  mutate(return_sq = return^2)

ggplot(exchange_rate_ret, aes(x = date, y = return_sq)) +
  geom_line() +
  labs(title = "Squared Log Returns of EUR/CNY",
       x = "Date",
       y = "Squared Log Return") +
  theme_minimal()

# 6. ACF/PACF plots

# ACF of return：return本身无明显自相关，均值过程接近白噪声
acf(exchange_rate_ret$return,
    main = "ACF of Log Returns")

# ACF of squared return：很多滞后项显著，且缓慢衰减
acf(exchange_rate_ret$return_sq,
    main = "ACF of Squared Log Returns")

#说明汇率收益率在均值上不可预测，但其波动是可预测的

# PACF
pacf(exchange_rate_ret$return)

# 7. ADF tests：检验平稳性

# 对汇率水平做ADF：non-stationary
adf_level <- adf.test(exchange_rate$rate)
print(adf_level)

# 对收益率做ADF: stationary
adf_return <- adf.test(exchange_rate_ret$return)
print(adf_return)

# 8. ARIMA: 自动选择出了ARIMA(0,0,0) with zero mean
#说明 return 序列在均值方程上基本就是白噪声，没有显著的 AR 或 MA 结构。
auto.arima(exchange_rate_ret$return)
checkresiduals(auto.arima(exchange_rate_ret$return))

# 9. ARCH effect test

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0)),  # 因为return基本是白噪声
  distribution.model = "norm"
)
fit <- ugarchfit(spec = spec, data = exchange_rate_ret$return)
fit
#alpha1 = 0.032936 表示短期冲击对当前波动的影响
#beta1 = 0.963593 表示过去波动对当前波动的持续影响
#alpha+beta接近1 说明汇率波动具有非常强的持续性，一旦波动上升，会在较长时间内维持高位

# 10. 提取Volatility
volatility <- sigma(fit)^2 #GARCH估计出来的条件方差
vol_df <- data.frame(
  date = exchange_rate_ret$date,
  vol = volatility
)

ggplot(vol_df, aes(x = date, y = vol)) +
  geom_line() +
  labs(title = "Conditional Variance (Exchange Rate Volatility)",
       x = "Date",
       y = "Volatility") +
  theme_minimal()

# 11. import trade data
export <- read_csv(file.choose())
str(export)
export <- export %>%
  mutate(Month = as.Date(paste0(Month, "-01")),
         Export = as.numeric(Export))

# 将rate整合成月度数据
vol_df <- data.frame(
  date = exchange_rate_ret$date,
  vol = volatility
)
# 对汇率取月平均
ex_monthly <- exchange_rate %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(rate = mean(rate, na.rm = TRUE))

# 对波动取月平均
vol_monthly <- vol_df %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(vol = mean(vol, na.rm = TRUE))

export <- export %>%
  rename(month = Month)

str(data_var)

#merge data
data_var <- ex_monthly %>%
  inner_join(vol_monthly, by = "month") %>%
  inner_join(export, by = "month")
data_var <- data_var %>%
  mutate(
    log_export = log(Export)
  )

# 12. 检查平稳性

data_var <- na.omit(data_var)
adf.test(data_var$rate)
adf.test(data_var$vol)
adf.test(data_var$log_export)

# 12. 对rate做差分
data_var <- data_var %>%
  mutate(
    d_rate = c(NA, diff(rate))
  ) %>%
  na.omit()
adf.test(data_var$d_rate)

# 13. 构建VAR模型

ts_data <- ts(data_var[, c("d_rate", "vol", "log_export")])

#选择lag：2
VARselect(ts_data)
model_var <- VAR(ts_data, p = 2, type = "const")

# 14. Granger因果检验：汇率波动会影响出口

causality(model_var, cause = "vol")

# 15. IRF：怎么影响，影响多久
irf_result <- irf(model_var,
                  impulse = "vol",
                  response = "log_export",
                  n.ahead = 12)

plot(irf_result)