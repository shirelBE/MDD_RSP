clear 
import delimited using "..\data_files\MDD_RSP_Final121Sample_Demog.csv", varnames(1)
rename patno participant
tabstat age, stat (n mean sd min max)
tab gender
tab religion
tab income
tab education //many values
tabstat education, stat(mean sd min max)
*generating ordinal edu variable:
gen edu_cat=1
drop edu_cat

recode education (0/11 = 1) (12 = 2) (13/20 = 3) (21/30 = 4), gen(edu_cat)
tabulate education edu_cat
tab edu_cat
