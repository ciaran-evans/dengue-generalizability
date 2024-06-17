# Assessing the generalizability of dengue classifiers

This repository contains the code and data to reproduce the results in "Assessing generalizability of a dengue classifier across multiple
datasets" (Lu, Li, and Evans). 

## File list

### data: folder containing the five different datasets analyzed in the paper

* `data/set-01/dengue-data-01.xls`: the raw data from [Tuan *et al.* (2015)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0003638)
* `data/set-02/dengue-data-02.xlsx`: the raw data from [Gasem *et al.* (2020)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0007927)
* `data/set-03/dengue-data-03.xlsx`: the raw data from [Saito *et al.* (2022)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0010414)
* `data/set-04/dengue-data-04.xlsx`: the raw data from [Park *et al.* (2018)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0006799)
* `data/set-05/dengue-data-05.xls`: the raw data from [Ngim *et al.* (2021)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0009445)

### code: folder containing all code to reproduce the analyses, figures, and tables in the paper

* `analysis.R`: the main analysis for the manuscript, covering data import and cleaning with the main explanatory variables (Age, WBC, and PLT) for the five different datasets; logistic regression modeling and in-sample performance; and assessment of generalizability between datasets. Covers Figures 1--5 and Table 2 in the main manuscript, Supplementary Tables 1--3, and Supplementary Figures 1--2
* `alternative_subset_1.R`: assessment of in-sample performance and generalizability for logistic regression with the variables in Alternative Subset 1 (see Table 1 in the manuscript for a description of the variables in each subset). Produces Supplementary Table 4 
* `alternative_subset_2.R`: assessment of in-sample performance and generalizability for logistic regression with the variables in Alternative Subset 2 (see Table 1 in the manuscript for a description of the variables in each subset). Produces Supplementary Table 5
* `alternative_subset_3.R`: assessment of in-sample performance and generalizability for logistic regression with the variables in Alternative Subset 3 (see Table 1 in the manuscript for a description of the variables in each subset). Produces Supplementary Table 6
* `model_comparison_original.R`: comparison of in-sample performance and generalizability for logistic regression, SVMs, and CART classification trees with the original explanatory variables (Age, WBC, and PLT). Produces Supplementary Figure 3
* `model_comparison_alternative_1.R`: comparison of in-sample performance and generalizability for logistic regression, SVMs, and CART classification trees with the variables in Alternative Subset 1 (see Table 1 in the manuscript for a description of the variables in each subset). Produces Supplementary Figure 4
* `model_comparison_alternative_2.R`: comparison of in-sample performance and generalizability for logistic regression, SVMs, and CART classification trees with the variables in Alternative Subset 2 (see Table 1 in the manuscript for a description of the variables in each subset). Produces Supplementary Figure 5
* `model_comparison_alternative_3.R`: comparison of in-sample performance and generalizability for logistic regression, SVMs, and CART classification trees with the variables in Alternative Subset 3 (see Table 1 in the manuscript for a description of the variables in each subset). Produces Supplementary Figure 6


### images: folder containing pdf images for the figures in the main manuscript and supplementary material

(Image names match numbers of the figures in the manuscript and supplement)

* `figure_1.pdf`
* `figure_2.pdf`
* `figure_3.pdf`
* `figure_4.pdf`
* `figure_5.pdf`
* `supplementary_figure_1.pdf`
* `supplementary_figure_2.pdf`
* `supplementary_figure_3.pdf`
* `supplementary_figure_4.pdf`
* `supplementary_figure_5.pdf`
* `supplementary_figure_6.pdf`


