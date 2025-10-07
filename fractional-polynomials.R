# Load necessary libraries
library(readr)
library(mfp)
library(ggplot2)
library(stringr)
source("blood_test_trends/frac_polynomials/modelling/config.R")

base_data_path <- "/data/WIPH-CanDetect/clean_data/Ben/blood_test_trends/frac_polynomials/"
output_plot_dir <- "/data/WIPH-CanDetect/results/Ben/blood_test_trends/frac_polynomials"

run_single_analysis <- function(test, cancer, control, interval) {
    cat("Running analysis for:", test, ",", cancer, ",", control, ",", interval, "\n")
    file_path <- paste0(base_data_path, test, "/", cancer, "/", test, "_", cancer, "_vs_", control, "_", interval, "_updated.csv")

    # Import data
    tryCatch({
        current_df <- read_csv(file_path, show_col_types = FALSE)
        cat("Loaded:", file_path, "\n")
        }, error = function(e) {
            warning(paste("Could not load", file_path, ":", e$message))
            })

    # Establish cancer status as a binary factor
    current_df$cancer_status <- factor(current_df$cancer_status, levels = c(0, 1))

    # Remove missing data
    current_df = na.omit(current_df)

    # Build fp model with up to two terms
    tryCatch({
        model <- mfp(cancer_status ~ fp(value, df = 4), data = current_df, family = binomial(link = "logit"))
        cat("Fractional polynomial model built successfully.\n")

        # Print model summary
        print(paste("Summary of model for", test, "_", cancer, "_vs_", control, interval, ":"))
        print(summary(model))

        }, error = function(e) {
            warning(paste("  Error building model for", test, "_", cancer, "_vs_", control, interval, ":", e$message))
            })

    # Check if model was built successfully
    if (is.null(model)) {
        cat("Skipping analysis for", test, ",", cancer, ",", control, ",", interval, ", due to model building issues.\n")
    }

    # Calculate predicted probabilities and odds on 100 bins
    n_bins <- 100
    current_df$bin <- cut(current_df$value,
    breaks = unique(quantile(current_df$value, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)),
    include.lowest = TRUE)

    # Take mean value for each bin
    bin_means <- aggregate(value ~ bin, data = current_df, FUN = mean, na.rm = TRUE)
    smooth_df <- data.frame(value = bin_means$value)
    cat("Calculated mean values for", n_bins, "bins.\n")

    pred <- predict(model, newdata = smooth_df, type = "link", se.fit = TRUE)
    smooth_df$logit_fit <- pred$fit
    smooth_df$logit_se <- pred$se.fit

    # 95% CI on logit scale
    smooth_df$logit_lower <- smooth_df$logit_fit - 1.96 * smooth_df$logit_se
    smooth_df$logit_upper <- smooth_df$logit_fit + 1.96 * smooth_df$logit_se

    # Transform to probability scale
    smooth_df$predicted_prob <- plogis(smooth_df$logit_fit)
    smooth_df$predicted_prob_lower <- plogis(smooth_df$logit_lower)
    smooth_df$predicted_prob_upper <- plogis(smooth_df$logit_upper)
    cat("Predicted probabilities and odds calculated.\n")

    # Convert to odds
    smooth_df$predicted_odds <- smooth_df$predicted_prob / (1 - smooth_df$predicted_prob)
    smooth_df$predicted_odds_lower <- smooth_df$predicted_prob_lower / (1 - smooth_df$predicted_prob_lower)
    smooth_df$predicted_odds_upper <- smooth_df$predicted_prob_upper / (1 - smooth_df$predicted_prob_upper)

    # Mean values
    mean_value <- mean(current_df$value, na.rm = TRUE)
    cat("Mean value for", test, "_", cancer, "_vs_", control, interval, ":", mean_value, "\n")

    # Create dataframe for prediction at the mean
    mean_df <- data.frame(value = mean_value)

    # Predict probability at the mean
    mean_prob <- predict(model, newdata = mean_df, type = "response")

    # Calc odds at the mean
    mean_odds <- mean_prob / (1 - mean_prob)
    if (mean_prob == 1) mean_odds <- Inf
    if (mean_prob == 0) mean_odds <- 0

    # Calculate Odds Ratio (OR) for each observation relative to the mean odds
    # Avoid division by zero if mean_odds is 0
    if (mean_odds == 0) {
        smooth_df$odds_ratio_vs_mean <- Inf # If mean odds is 0, any non-zero odds will have infinite OR
        warning(paste0("Mean odds for ", test, "_", cancer, "_vs_", control, interval, " is 0. Odds ratios set to Inf where predicted_odds > 0."))
    } else {
        smooth_df$odds_ratio_vs_mean <- smooth_df$predicted_odds / mean_odds
        smooth_df$odds_ratio_lower <- smooth_df$predicted_odds_lower / mean_odds
        smooth_df$odds_ratio_upper <- smooth_df$predicted_odds_upper / mean_odds
    }
    cat("Odds and Odds Ratios calculated.\n")

    # Store final dataframes
    output_path <- paste0(base_data_path, test, "/", cancer, "/", test, "_", cancer, "_vs_", control, "_", interval, "_results.csv")
    write_csv(smooth_df, output_path)
    cat("Results saved to:", output_path, "\n")

    min_value <- NULL
    max_value <- NULL

    tryCatch({
        min_value <- blood_test_dict[[test]]$min
        max_value <- blood_test_dict[[test]]$max
    }, error = function(e) {
        warning(paste("Error in retrieving min/max values for", test, ":", e$message))
    })

    # Plotting
    odds_ratio_plot <- ggplot(smooth_df, aes(x = value, y = odds_ratio_vs_mean)) +
    geom_ribbon(aes(ymin = odds_ratio_lower, ymax = odds_ratio_upper), fill = "blue", alpha = 0.2) +
    geom_line(linewidth = 1, color = "blue") +
    #geom_point(alpha = 0.5, color = "blue") +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(xlim = c(min_value, max_value)) +
    labs(
    title = paste0("Odds Ratio - ", blood_test_figure_titles[test], " (", str_to_title(cancer), " vs. ",str_to_title(control), "): ", interval),
    x = paste0(blood_test_figure_titles[test], " ", units_named_vector[test]),
    y = "Odds Ratio (vs. mean)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

    # Define plot filename and path
    plot_filename <- paste0(test, "_", cancer, "_vs_", control, "_", interval, "_odds_ratio_plot.png")
    #plot_path <- file.path(output_plot_dir, test, cancer, plot_filename) temp ignore
    plot_path <- file.path(output_plot_dir, plot_filename)

    # Save the plot
    ggsave(plot_path, plot = odds_ratio_plot, width = 10, height = 7, dpi=400, device = "png")
    cat("  Odds Ratio plot for ", interval, " saved to:", plot_path, "\n")
}

# Loop through combinations of blood tests, cancers, and controls
for (test in blood_tests) {
  for (cancer in cancers) {
    for (control in controls) {
        run_single_analysis(test, cancer, control, "int1") # Cancer vs control 1 month to 2 years
        run_single_analysis(test, cancer, control, "int2") # Cancer vs control 2 to 5 years
    }
  }
}
cat("All analysis and plots complete")