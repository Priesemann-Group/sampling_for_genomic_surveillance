#libraries to use
import pandas as pd
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind

def define_significance_value(p_value):

    summary = ""
    if p_value > 0.05:
        summary = "NS"
    elif p_value > 0.01 and p_value <= 0.05:
        summary = "*"
    elif p_value > 0.001 and p_value <= 0.01:
        summary = "**"
    elif p_value > 0.001 and p_value <= 0.001:
        summary = "***"
    else:
        summary = "****"

    return summary

#read datasets ignored columns
df_adaptativo = pd.read_csv("comparative_data\\Res_Adaptive.csv")
df_constant = pd.read_csv("comparative_data\\Res_Constant.csv")

#adding columns name
name_columns = ["column_{}".format(i+1) for i in range(len(df_adaptativo.columns))]
df_adaptativo.columns = name_columns
df_constant.columns = name_columns

matrix_data = []

#start comparison for each columns
for column in name_columns:

    #get values of the columns using the proposed name
    values_adaptativo = df_adaptativo[column]
    values_constante = df_constant[column]

    #apply different statistical test, all returns the statistic value and the p-value
    u_test_result = mannwhitneyu(values_adaptativo, values_constante, alternative='two-sided')
    kolmogorov = ks_2samp(values_adaptativo, values_constante)
    different_variance = ttest_ind(values_adaptativo, values_constante, equal_var=False)

    #get the number of * by definitions
    sig_u_test = define_significance_value(u_test_result[1])
    sig_k_test = define_significance_value(kolmogorov[1])
    sig_d_test = define_significance_value(different_variance[1])

    #create row with different element
    row = [column, u_test_result[0], u_test_result[1], sig_u_test, kolmogorov[0], kolmogorov[1], sig_k_test, different_variance[0], different_variance[1], sig_d_test]

    #adding row to matrix
    matrix_data.append(row)

#creating df to export results in a csv format
df_export = pd.DataFrame(matrix_data, columns=['column_compare', 'statistic-u-test', 'p-value-u-test', 'significance-level-u-test', 'statistic-ks-test', 'p-value-ks-test', 'significance-level-ks-test', 'statistic-t-test', 'p-value-t-test', 'significance-level-t-test'])

df_export.to_csv("comparative_data\\results_test.csv", index=False)