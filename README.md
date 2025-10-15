# Blood tests and body mass index as early markers of upper-gastrointestinal cancer
[![python](https://img.shields.io/badge/python-3.12.6-blue)](https://www.python.org)
[![R](https://img.shields.io/badge/R-4.5.1-blue)](https://www.r-project.org/)
[![piecewise_regression](https://img.shields.io/badge/piecewise_regression-1.5.0-blue)](https://piecewise-regression.readthedocs.io/en/latest/)
[![mfp](https://img.shields.io/badge/mfp-1.5.5.1-blue)](https://cran.r-project.org/web/packages/mfp/index.html)
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/BenMcGuirk/ugi-cancer-early-markers/blob/main/LICENSE)

## Introduction
This repository contains analysis code used for the paper "Tests, comorbidities, medications and symptoms as early markers for upper gastrointestinal cancer: a longitudinal matched case control study". Temporal trends of blood tests and BMI in the five years pre UGI cancer diagnosis were analysed using joinpoint regression. Temporal associations of blood tests and BMI were analysed using fractional polynomials with up to two terms.

Joinpoint regression was performed using the piecewise regression python package: https://piecewise-regression.readthedocs.io/en/latest/  

Fractional polynomials were built using the "mfp" R package: https://cran.r-project.org/web/packages/mfp/index.html

## Usage
⚠️ Note: The analysis was conducted on CPRD Aurum, a restricted dataset. The data cannot be shared.
The code is provided for transparency and documentation purposes, however, researchers with access to CPRD Aurum may use this code to reproduce the analyses, with preprocessing steps outlined below.

The CSV file paths in the scripts point to local directories. 
Researchers using this code should replace these paths with their own file locations.

## Preprocessing
CPRD raw data consists of flat text files.
Patient text files are used to identify populations, and observation files are used to identify medical events.

Preprocessing steps differ slightly between joinpoint regression and fractional polynomials.  
Once populations have been identified, index dates assigned and test data extracted:
### Joinpoint regression
1. Global preprocessing (performed on entire dataset)
- Remove duplicate data
- Calculate how many months pre index the test occured (using 'obsdate' column)
- Filter data to only include data from 1 month pre index to 61 months pre index (5 year period)*
2. Per test preprocessing**
- Loop over each test and filter using each test's codelist
- Remove implausible values (choose ranges with clinical guidance)
3. Per test and group preprocessing (for each test and group)
- Ensure patients only have one test result per 3 month interval
- Calculate mean test value in each 3 month interval 
- Calculate the proportion of patients with an abnormal result. Reference ranges provided by CPRD in 'numrangehigh' and 'numrangelow' columns, as reference ranges can vary by lab. For tests missing ref range, choose normal ranges with clinical guidance.  
Final output of preprocessing for mean test values should be a dataframe with two columns, months pre index (x-axis) and mean value (y-axis).  
Final output of preprocessing for abnormal proportions is the same but with proportion instead of mean value for y-axis.

### Fractional polynomials
1. Global preprocessing
- As above, except no need to group into 3 month intervals
2. Per test preprocessing**
- Loop over each test and filter using each test's codelist
- Remove implausible values (choose ranges with clinical guidance)
- Add 'cancer_status' as binary variable
- For each group split data into two intervals (1 month-2 years pre diagnosis, 2-5 years pre diagnosis)
- Ensure patients only have one test result per interval
- Concatenate and shuffle all combinations of case vs control comparisons e.g. pancreatic vs general controls (interval 1), oesophageal vs benign controls (interval 2) etc  
Final output of preprocessing for fractional polynomials is a dataframe for each comparison of cases and controls, and each interval, per test. Each dataframe has just two columns, value (x axis) and cancer status. Odds ratios are calculated with respect to odds at the mean test value.

\* For joinpoint regression we group our data into 3 month intervals starting with 1-4 months. This means that the last interval is 58-61 months due to shifting the 5-year period 1 month backwards. If you have enough data, 1 month intervals may be appropriate.  
** Preprocessing steps for BMI and NLR differ to the other tests:
Extra BMI data is calculated by identifying height and weight measurements taken on the same day for patients without an explicit BMI result.
NLR data calculated using neutrophil count and lymphocyte count measurements taken on the same day for each patient.

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
