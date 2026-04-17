# -------------------------------------------------------------------------
# SCRIPT: Rainfall Data Cleaning & Mann-Kendall Trend Analysis
# STUDY AREA: Dhaka & Sylhet, Bangladesh (1981-2025)
# -------------------------------------------------------------------------
# SCRIPT: Rainfall Data Cleaning & Trend Analysis (Corrected)
# -------------------------------------------------------------------------

# install.packages(c("dplyr", "Kendall", "trend"))

library(dplyr)  
library(Kendall)  
library(trend)    

data <- read.csv(file.choose(), header = TRUE)

# Data cleaning

clean_data <- data %>%
  select(year, ADM2_NAME, mean) %>%
  
  mutate(year = as.integer(year)) %>%
  filter(!is.na(mean))


dhaka_data <- clean_data %>%
  filter(ADM2_NAME == "Dhaka") %>%
  arrange(year)

sylhet_data <- clean_data %>%
  filter(ADM2_NAME == "Sylhet") %>%
  arrange(year)

get_summary_df <- function(data, name) {
  data.frame(
    Region = name,
    Mean = round(mean(data), 2),
    Median = round(median(data), 2),
    Std_Dev = round(sd(data), 2),
    Skewness = round(skewness(data), 2),
    Kurtosis = round(kurtosis(data), 2),
    CV = round(sd(data)/mean(data), 2)
  )
}

# Combine Dhaka and Sylhet statistics into one table
summary_table <- rbind(
  get_summary_df(dhaka_data$mean, "Dhaka"),
  get_summary_df(sylhet_data$mean, "Sylhet")
)

# 1. View the table in the console
print(summary_table)

# 2. Save as CSV file
write.csv(summary_table, "Rainfall_Summary_Statistics.csv", row.names = FALSE)

cat("\nSummary Statistics saved as 'Rainfall_Summary_Statistics.csv'!\n")





#Mann Kendall Test
cat("\n--- Mann-Kendall Test Result: Dhaka ---\n")
mk_dhaka <- MannKendall(dhaka_data$mean)
print(summary(mk_dhaka))

cat("\n--- Mann-Kendall Test Result: Sylhet ---\n")
mk_sylhet <- MannKendall(sylhet_data$mean)
print(summary(mk_sylhet))

#  (Sen's Slope)
cat("\n--- Sen's Slope Estimation ---\n")
sens_dhaka <- sens.slope(dhaka_data$mean)
sens_sylhet <- sens.slope(sylhet_data$mean)

print("Sen's Slope (Dhaka):")
print(sens_dhaka)

print("Sen's Slope (Sylhet):")
print(sens_sylhet)


# -------------------------------------------------------------------------
# 11. Creating and Saving Trend Analysis Summary Table
# -------------------------------------------------------------------------


trend_results <- data.frame(
  Region = c("Dhaka", "Sylhet"),
  
  # Mann-Kendall P-value
  P_Value = c(round(as.numeric(mk_dhaka$sl), 4), 
              round(as.numeric(mk_sylhet$sl), 4)),
  
  # Kendall's Tau 
  Tau = c(round(as.numeric(mk_dhaka$tau), 3), 
          round(as.numeric(mk_sylhet$tau), 3)),
  
  # Sen's Slope (Rainfall change in mm/year)
  Sens_Slope_mm_per_year = c(round(as.numeric(sens_dhaka$estimates), 3), 
                             round(as.numeric(sens_sylhet$estimates), 3)),
  
  # Significant test
  Significance = c(ifelse(mk_dhaka$sl < 0.05, "Significant", "Non-Significant"),
                   ifelse(mk_sylhet$sl < 0.05, "Significant", "Non-Significant"))
)


cat("\n--- Rainfall Trend Analysis Summary Table ---\n")
print(trend_results)


write.csv(trend_results, "Rainfall_Trend_Analysis_Results.csv", row.names = FALSE)

cat("\nSuccess! The table has been saved as 'Rainfall_Trend_Analysis_Results.csv' in your folder.\n")
getwd

#Graphical Representation
par(mfrow=c(1, 2), mar=c(4.5, 4.5, 3, 1.5))
year_ticks <- seq(1980, 2025, length.out = 5)

# Dhaka Plot
plot(dhaka_data$year, dhaka_data$mean, type="o", pch=19, col="blue",
     main="Rainfall Trend: Dhaka", xlab="Year", ylab="Rainfall (mm)",
     xaxt="n")
axis(1, at = year_ticks, labels = round(year_ticks))
abline(lm(dhaka_data$mean ~ dhaka_data$year), col="red", lwd=2, lty=2)

# Sylhet Plot
plot(sylhet_data$year, sylhet_data$mean, type="o", pch=19, col="darkgreen",
     main="Rainfall Trend: Sylhet", xlab="Year", ylab="Rainfall (mm)",
     xaxt="n")
axis(1, at = year_ticks, labels = round(year_ticks))
abline(lm(sylhet_data$mean ~ sylhet_data$year), col="red", lwd=2, lty=2)
getwd

# -------------------------------------------------------------------------
# ARIMA Model Fitting & Diagnostics (The Check)
# -------------------------------------------------------------------------
library(forecast) # For ARIMA modeling and forecasting
# Creating Time Series objects (Starting 1981, Annual)
ts_dhaka <- ts(dhaka_data$mean, start = 1981, frequency = 1)
ts_sylhet <- ts(sylhet_data$mean, start = 1981, frequency = 1)

# Fitting the best ARIMA model
model_arima_dhaka <- auto.arima(ts_dhaka, stepwise = FALSE, approximation = FALSE)
model_arima_sylhet <- auto.arima(ts_sylhet, stepwise = FALSE, approximation = FALSE)

# Residual Check 
cat("\n--- Checking Residuals for Dhaka Model ---\n")
checkresiduals(model_arima_dhaka)
png("Dhaka_Model_Diagnostics.png", width = 1000, height = 700, res = 150)


checkresiduals(model_arima_dhaka)


dev.off()
cat("\n--- Checking Residuals for Sylhet Model ---\n")
checkresiduals(model_arima_sylhet)
# Print Accuracy Metrics (RMSE, MAPE, etc.)
cat("\n--- Accuracy Metrics ---\n")
print(accuracy(model_arima_dhaka))
print(accuracy(model_arima_sylhet))

#  15-Year Rainfall Forecasting (2026 - 2040)
# -------------------------------------------------------------------------
forecast_dhaka_arima <- forecast(model_arima_dhaka, h = 15)
forecast_sylhet_arima <- forecast(model_arima_sylhet, h = 15)

# STEP 7: Visualization
# -------------------------------------------------------------------------
par(mfrow=c(1, 2), mar=c(4.5, 4.5, 3, 1.5))

# Plot Dhaka Forecast
plot(forecast_dhaka_arima, main="15-Year Rain Forecast: Dhaka", 
     xlab="Year", ylab="Rainfall (mm)", col="blue", lwd=2)

# Plot Sylhet Forecast
plot(forecast_sylhet_arima, main="15-Year Rain Forecast: Sylhet", 
     xlab="Year", ylab="Rainfall (mm)", col="darkgreen", lwd=2)

# Reset plot window
par(mfrow=c(1,1))

# STEP 8: Save Forecasted Values to CSV
# -------------------------------------------------------------------------
arima_output <- data.frame(
  Year = 2026:2040,
  Dhaka_Pred_mm = round(as.numeric(forecast_dhaka_arima$mean), 2),
  Sylhet_Pred_mm = round(as.numeric(forecast_sylhet_arima$mean), 2)
)
write.csv(arima_output, "ARIMA_15Year_Rainfall_Forecast.csv", row.names = FALSE)

cat("\nSuccess! Analysis complete and results saved to 'ARIMA_15Year_Rainfall_Forecast.csv'\n")


Rainfall Anomaly Detection (Standardized Anomaly Index - SAI)
# -------------------------------------------------------------------------

# The Standardized Anomaly Index (SAI) identifies years with 
# exceptionally high or low rainfall compared to the long-term average.

# Calculate anomalies for Dhaka and Sylhet
# Formula: (Annual Value - Mean) / Standard Deviation
dhaka_data <- dhaka_data %>%
  mutate(anomaly = (mean - mean(mean)) / sd(mean))

sylhet_data <- sylhet_data %>%
  mutate(anomaly = (mean - mean(mean)) / sd(mean))

# Visualization: Creating Anomaly Bar Charts
# Setting up the plot area for two stacked charts
par(mfrow=c(2,1), mar=c(4,4,3,1))

# --- Plot 1: Dhaka Rainfall Anomaly ---
barplot(dhaka_data$anomaly, 
        names.arg = dhaka_data$year, 
        col = ifelse(dhaka_data$anomaly > 0, "#3498db", "#e74c3c"), # Blue for wet, Red for dry
        main = "Rainfall Anomaly Index: Dhaka (1981-2025)",
        ylab = "Standardized Index", 
        border = NA,
        las = 2,      # Rotate axis labels for better readability
        cex.names = 0.7)
abline(h = 0, lwd = 1.5, lty = 2) # Adding a dashed baseline at zero

# --- Plot 2: Sylhet Rainfall Anomaly ---
barplot(sylhet_data$anomaly, 
        names.arg = sylhet_data$year, 
        col = ifelse(sylhet_data$anomaly > 0, "#27ae60", "#f39c12"), # Green for wet, Orange for dry
        main = "Rainfall Anomaly Index: Sylhet (1981-2025)",
        ylab = "Standardized Index", 
        border = NA,
        las = 2, 
        cex.names = 0.7)
abline(h = 0, lwd = 1.5, lty = 2)

# Resetting the graphical parameters
par(mfrow=c(1,1))


