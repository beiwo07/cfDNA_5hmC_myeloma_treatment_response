# Identifying 5-hydroxymethylcytosine (5hmC) signatures in circulating cell-free DNA and treatment response in multiple myeloma with a machine learning approach (2023-2025)

# Background and importance 

Despite recent therapeutic advancements in multiple myeloma (MM), many patients do not achieve an adequate response to initial drug regimens, even with novel agents. Currently, there is no reliable method to predict treatment response in front-line setting of MM. Genome-wide 5-Hydroxymethylcytosine (5hmC) alterations in circulating cell-free DNA (cfDNA) have been associated with survival outcomes of MM and are predictive of treatment response in several human cancers. The role of cfDNA-based genome-wide 5hmC in predicting treatment response in MM has not been assessed.

# Method 

The current study included 258 patients newly diagnosed with MM between 2010-2017. We profiled genome-wide 5hmC in cfDNA extracted from peripheral blood plasma collected at the time of diagnosis. The best response achieved within the line of initial therapy was assessed and abstracted from patient’s EHRs. We also collected baseline demographic and key clinical variables from EHRs. Achieving a complete response (CR) or stringent CR (sCR) was used as the threshold to dichotomize treatment response into two groups: >=CR and <CR. We divided patients into a training (n=159) and a testing set (n=87). In the training set, differentially modified 5hmC regions were identified comparing >=CR and <CR (p-value<0.05, log2 fold change>0.1), followed by logistic regression models with the elastic net regularization for feature selection through 10-fold cross validation. The optimal “alpha” and “lambda” were identified to maximize model performance measured by the AUC-ROC. This selection step was repeated for 100 times in the training set, and 5hmC-modified genes appearing across all iterations were identified as signature genes. The differential 5hmC enrichment levels of the signature genes were used to calculate a wp-score for each patient. Lastly, we assessed the association between the wp-score and response among patients in the independent testing set. 

<img src="plots/5hmc_response_design.drawio.png" width="600">

## **Citation**  

> Wang, B., Derman, B. A., Zhang, Z., Kowitwanich, K., Cursio, J. F., Appelbaum, D., ... & Chiu, B. (2024). Identifying 5-Hydroxymethylcytosine Signatures in Circulating Cell-Free DNA and Treatment Response in Multiple Myeloma with a Multi-Step Machine Learning Approach. Blood, 144, 3336.






