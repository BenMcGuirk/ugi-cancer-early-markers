import pandas as pd
import piecewise_regression as pr
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from blood_test_trends.joinpoint_regression.preprocessing.config import BLOOD_TESTS, BLOOD_TEST_RANGES, group_titles
from blood_test_trends.joinpoint_regression.modelling.config import *
from blood_test_trends.joinpoint_regression.preprocessing.preprocessing import groups
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

np.random.seed(42)

def fit_and_plot_data(test, data_type, plot_title_suffix, y_label, data_key, output_filename_suffix):
    """
    Fits models and generates a single plot with all groups overlaid (one line per group).
    """
    data_and_models = {}

    # Load data for the current test
    for group in groups:
        values_path = f'/data/WIPH-CanDetect/clean_data/Ben/blood_test_trends/joinpoint_regression/{test}/{group}/all_values.csv'
        values_df = pd.read_csv(values_path)

        ci_path = f'/data/WIPH-CanDetect/clean_data/Ben/blood_test_trends/joinpoint_regression/{test}/{group}/all_values_ci.csv'
        ci_values = pd.read_csv(ci_path)

        prop_path = f'/data/WIPH-CanDetect/clean_data/Ben/blood_test_trends/joinpoint_regression/{test}/{group}/proportions.csv'
        patient_proportions = pd.read_csv(prop_path)

        abnormal_prop_path = f'/data/WIPH-CanDetect/clean_data/Ben/blood_test_trends/joinpoint_regression/{test}/{group}/abnormal_proportions.csv'
        abnormal_proportions = pd.read_csv(abnormal_prop_path)
        logger.info(f"Loaded data for {test} and {group} for {data_type}")

        data_and_models[group] = {
            'values': values_df,
            'ci': ci_values,
            'proportions': patient_proportions,
            'abnormal_proportions': abnormal_proportions,
        }

    logger.info(f"Loaded all data for {test} ({data_type}) and saved to dictionary.")

    # Fit models for each group
    for group in groups:
        x = data_and_models[group][data_key]['months_pre_diagnosis'].values
        y = data_and_models[group][data_key]['proportion' if data_key != 'values' else 'value'].values

        ms = pr.ModelSelection(x, y, max_breakpoints=3, min_distance_between_breakpoints=0.1, verbose=False)
        filtered_summaries = [s for s in ms.model_summaries if s['n_breakpoints'] != 0]
        valid_bic_summaries = [entry for entry in filtered_summaries if entry['bic'] is not None]

        if not valid_bic_summaries: # If no valid models, fit linear regression
            fitted_model = LinearRegression()
            fitted_model.fit(x.reshape(-1, 1), y)
            y_pred_data = fitted_model.predict(x.reshape(-1, 1))
            x_plot = x
            y_pred_plot = y_pred_data
        else:
            lowest_bic_entry = min(valid_bic_summaries, key=lambda x: x['bic'])
            optimal_n_breakpoints = lowest_bic_entry['n_breakpoints']

            for i in range(optimal_n_breakpoints-1): # To match random initialisation state of ModelSelection (which loops from 0 to max breakpoints)
                pr.Fit(x, y, n_breakpoints=i+1)
            fitted_model = pr.Fit(x, y, n_breakpoints=optimal_n_breakpoints, min_distance_between_breakpoints=0.1)
            y_pred_data = fitted_model.predict(x)
            logger.info(f"Fitted {data_type} model for {test} and {group} with {optimal_n_breakpoints} breakpoints")
            logger.info(f"Model summary for {test}, {group}, {data_type}:")
            logger.info(fitted_model.summary())
            davies_results_dict = fitted_model.get_results()

            if davies_results_dict.get('davies') > 0.05: # If Davies test p-value > 0.05, fit linear regression
                fitted_model = LinearRegression()
                fitted_model.fit(x.reshape(-1, 1), y)
                y_pred_data = fitted_model.predict(x.reshape(-1, 1))
                x_plot = x
                y_pred_plot = y_pred_data
            else:
                x_plot = np.linspace(x.min(), x.max(), 100)
                y_pred_plot = fitted_model.predict(x_plot)

        # Save fitted model and predictions
        data_and_models[group]['model_' + data_type] = fitted_model
        data_and_models[group]['x_plot_' + data_type] = x_plot
        data_and_models[group]['y_pred_plot_' + data_type] = y_pred_plot

    # One plot with all groups
    colours = plt.cm.tab10.colors
    group_colors = {group: colours[i % len(colours)] for i, group in enumerate(groups)}
    # Override the last group's color to black
    last_group = groups[-1]
    group_colors[last_group] = "black"
    fig, ax = plt.subplots(figsize=(12, 8))
    xticks = np.arange(-58, -1, 6)

    for i, group in enumerate(groups):
        group_title = group_titles[group]['title']
        group_data = data_and_models[group]

        x = group_data[data_key]['months_pre_diagnosis']
        y = group_data[data_key]['proportion' if data_key != 'values' else 'value']
        x_plot_line = group_data['x_plot_' + data_type]
        y_pred_plot_line = group_data['y_pred_plot_' + data_type]

        colour = group_colors[group]
        if group == 'general_controls' or group == 'benign_GI_controls':
            linestyle = ':' # Dotted line for control groups
        else:
            linestyle = '-' # Solid line for cancer groups
        ax.plot(x_plot_line, y_pred_plot_line, linewidth=2, label=group_title, color=colour, linestyle=linestyle)

    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks])
    ax.set_xlabel("Months pre diagnosis")
    ax.set_ylabel(y_label)
    ax.set_title(f"{BLOOD_TEST_RANGES[test]['title']} - {plot_title_suffix}", fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(title="Groups")

    plt.tight_layout()
    output_dir = f"/data/WIPH-CanDetect/results/Ben/blood_test_trends/joinpoint_regression_grouped_plots/"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f"{output_dir}{test}_{output_filename_suffix}_all_groups.png", dpi=400)
    logger.info(f"Saved combined {data_type} plot for {test}")
    plt.close(fig)

def fit_models():
    """
    Fit models to the blood test data for each test and group.
    """
    tests = BLOOD_TESTS
    for test in tests:
        units = BLOOD_TEST_RANGES[test]['unit']

        # Call the function for proportions
        #fit_and_plot_data(
        #    test=test,
        #    data_type='proportions',
        #    plot_title_suffix='Proportion of Patients with a Test',
        #    y_label='Proportion',
        #    data_key='proportions',
        #    output_filename_suffix='proportions'
        #)

        # Call function for abnormal proportions
        fit_and_plot_data(
            test=test,
            data_type='abnormal_proportions',
            plot_title_suffix='Proportion of Abnormal Tests',
            y_label='Abnormal Proportion',
            data_key='abnormal_proportions',
            output_filename_suffix='abnormal_proportions'
        )

        # Call function for values
        fit_and_plot_data(
            test=test,
            data_type='values',
            plot_title_suffix='Mean Test Values',
            y_label=units,
            data_key='values',
            output_filename_suffix='values'
        )

if __name__ == "__main__":
    fit_models()
