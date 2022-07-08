%%%Script to extract statistical results generated with python routines%%%

clear all
close all
clc

%% Defining empty structures for the extraction loop


Kindice = [25 50 100];
str_base = "C:\Users\scontreras\Dropbox\Trabajo\CEBIB\Proyectos pendientes\Paper muestreo\mFiles\Fig_1_Statistical_analysis\results_new_process\exp_";
str_end = "\results_comparison.csv";
Cell_Tab_u = cell(16,6);
Cell_Tab_ks = cell(16,6);
Cell_Tab_t = cell(16,6);
for expno = 1:4
    for i = 1:3
        P.K = Kindice(i);
        id_fila = 4*(expno-1)+1;
        id_col_num = 2*(i-1)+1;
        id_col_ast = 2*i;
        [pvalueutest, significancelevelutest, pvaluekstest, significancelevelkstest, pvaluettest, significancelevelttest] = importfile1(strcat(str_base,num2str(expno),"_k_",num2str(P.K),str_end), [2, Inf]);
        for j = 1:4
            id_fila = 4*(expno-1)+j;
            Cell_Tab_u{id_fila,id_col_num} = pvalueutest(j);
            Cell_Tab_u{id_fila,id_col_ast} = significancelevelutest(j);
            Cell_Tab_ks{id_fila,id_col_num} = pvaluekstest(j);
            Cell_Tab_ks{id_fila,id_col_ast} = significancelevelkstest(j);
            Cell_Tab_t{id_fila,id_col_num} = pvaluettest(j);
            Cell_Tab_t{id_fila,id_col_ast} = significancelevelttest(j);
        end
    end
end
T = cell2table(Cell_Tab_u,...
    "VariableNames",["p-value (K25)" "sig. level (K25)" "p-value (K50)" "sig. level(K50)" "p-value (K100)" "sig. level(K100)"]);
writetable(T,'Summary_Statistical_Results_utest.csv')
T = cell2table(Cell_Tab_ks,...
    "VariableNames",["p-value (K25)" "sig. level (K25)" "p-value (K50)" "sig. level(K50)" "p-value (K100)" "sig. level(K100)"]);
writetable(T,'Summary_Statistical_Results_kstest.csv')
T = cell2table(Cell_Tab_t,...
    "VariableNames",["p-value (K25)" "sig. level (K25)" "p-value (K50)" "sig. level(K50)" "p-value (K100)" "sig. level(K100)"]);
writetable(T,'Summary_Statistical_Results_ttest.csv')
