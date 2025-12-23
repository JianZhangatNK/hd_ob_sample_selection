ILLUSTRATION: These are codes and data for our working paper entitled"Has the Dual Economic Transformation Truly Eliminate Hukou Discrimination?—A Re-examination Based on High-Dimensional Sample Selection Model with Oaxaca-Blinder Decomposition"
Originally, the title in Chinese is "二元经济转型是否真正消除了户籍歧视？——基于高维样本选择模型Oaxaca-Blinder分解的再考察".

(I) Main codes
simulation_for_hd_ob_ss_par.R: R code for simulation. It provides comparison between the model with sample selection and without sample selection.

Empirical_income_hd_ob_ss.R: R code for empirical study (Outcome Y is income). Hukou discrimination in labor income.
datacl_income.R: Data clearing code for "Empirical_income_hd_ob_ss.R".

Empirical_GJ_hd_ob_ss.R: R code for empirical study (Outcome Y is "Good Job" indicator). Hukou discrimination in labor income.
datacl_GJ.R: Data clearing code for "Empirical_income_hd_ob_ss.R".


(II) Source functions

source_hdobss.R: R codes for related functions (created by ourselves)
AutoML.rar： Functions for Automatic Debias Machine Learning (created by Chernozhukov et al., 2022;  https://doi.org/10.3982/ECTA18515)

(III) Data:

HDOBcfps2018_dataV2.csv :  CFPS2018
HDOBcfps2022_dataV2.xls :   CFPS2022
