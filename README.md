# Blood tests and body mass index as early markers of upper-gastrointestinal cancer

## Introduction
This repository contains analysis code used for the paper "Tests, comorbidities, medications and symptoms as early markers for upper gastrointestinal cancer: a longitudinal matched case control study". Temporal trends of blood tests and BMI in the five years pre UGI cancer diagnosis were analysed using joinpoint regression. Temporal associations of blood tests and BMI were analysed using fractional polynomials with up to two terms.

Joinpoint regression was performed using the piecewise regression python package: https://piecewise-regression.readthedocs.io/en/latest/  

Fractional polynomials were built using the "mfp" R package: https://cran.r-project.org/web/packages/mfp/index.html

## Usage
⚠️ Note: The analysis was conducted on CPRD Aurum, a restricted dataset. The data cannot be shared.
The code is provided for transparency and documentation purposes only. 
Researchers with access to CPRD Aurum may use this code to reproduce the analyses.

The CSV file paths in the scripts point to local directories. 
Researchers using this code should replace these paths with their own file locations. 
The underlying data (CPRD Aurum) cannot be shared.

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
