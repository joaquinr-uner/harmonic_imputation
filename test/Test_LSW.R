

library(forecast)
library(forecastLSW)
# Get the data points in form of a R vector.

rain <- c(987,1025,978,774,1563,569,1456,789,1479,566,1563,1698)


# Convert it to a time series object.

rain_ts <- ts(rain,start = c(2020,1),frequency = 12)


# Print the timeseries data.

print(rain_ts)


plot(rain_ts, main = "Time series data")


summary(rain_ts)

plot(rain_ts)

abline(reg=lm(rain_ts~time(rain_ts)))

model <- arima(rain_ts, order = c(1,0,0))


model


model <- arima(rain_ts, order = c(2,0,0))

model

pacf(ts(diff(log10(rain_ts))),main='PACF for rain_ts')


library(lpacf)

forecast <- forecastlpacf(rain_ts,h=5,lag.max=6,filter.number=1,family="DaubExPhase",
                          forecast.type='extend')

pred <- forecast$mean

plot(pred)
rlapcf <- lpacf(rain_ts,lag.max=2)

