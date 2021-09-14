/*
Classic EM EXP was conducted at the Psychotherapy Lab - University of Haifa. 
the patients were diagnosed as clinicallly depressed (different levels of MDD), using Hemilton questionnaire, MINI, and BDI (score: 19+)
The psychotherapy treatment includes 20 sessions. 16 weekly session and then 4 followups once a month.
T1= part1= before the treatment, intake+1st session; T2= part2= 1st followup (not included in the current analysis).  
Dep level was measured at the beggining of each session.
*/

/*
///////////////////////////////////
////Preparing Data_set + FILTERS//////
//////////////////////////////////
clear
use "..\data_files\taskANDdep_UnFiltered.dta"

*Drop non-depressed subjs according to HRSD
tab subnum if totalhrsd17<14 
drop if totalhrsd17<14  

*Drop non-depressed subjs according to BDI-II
list subnum if bdi<19
drop if bdi<19
tab  cont

merge m:m subnum using "..\data_files\RSP_MDD_MINI_data.dta" 
drop if _merge==2
drop if mdd_mini==0 
tab  cont
drop _merge

*** Create handy variables ***
by subnum, sort:egen mean_correct=mean(correct) // create mean correct per subj
gen rt =output_rt*1000
drop output_rt

by subnum, sort:egen mean_rt=mean (rt) // create mean rt per sub

sort subnum cont trialnum
*/

clear
global apply_filters 1

*mean_correct by condition (cont)
tabstat mean_corr, stat(mean sd n) by(cont) 
hist mean_correct, freq 
hist mean_correct, freq by (cont) 
graph box mean_correct, over(cont)

by cont, sort: stem mean_correct // 
hist rt if mean_c<.5, freq by (cont) // 
tab subnum if mean_c<0.5
drop if mean_correct < 0.5 // 
tabstat mean_correct, stat(mean sd n) by(cont) //
graph box mean_correct, over(cont)

*Note! We are using different filter in this depressed population (.50, regular filter is >.85)

if $apply_filters{
**********************************************************
********************  EM  FILTERS  ***********************
**********************************************************

*1. filter by %correct (per subj we dropped earier- <.50%, now per cont)
*we drop subjects with low percent of correct answers per cond (outliers), so that we can see the varience of subjs in each cond 
*, then, we will calculate the sd's of RT.
by cont, sort:egen mean_correct_cont=mean (mean_correct)
by cont, sort:egen SDcorrect_cont= sd(mean_correct)
gen z_correct=(mean_correct-mean_correct_cont)/SDcorrect
gen include=1 if z_correct<-2 | z_correct>2
tab subnum if include==1
drop if include==1 // 
drop include  


*2. filter by rt
*drop observations when participant didnt press at all (and was coded as incorrect)
drop if rt < 200  // drop very fast trials    
drop if rt > 800 // drop very slow trials. Note! because pop is MDD, regular filter of > 700 was changed; response time window was set as 850 msec as in all our previous EM experiments

graph box mean_rt, over(cont)


by cont, sort: egen mean_rt_cont = mean (mean_rt) // create mean of RT by cont
by cont, sort: egen sd_rt_cont = sd (rt) // create SD of RT by cont
gen z_rt = (rt - mean_rt_cont)/sd_rt_cont // create Z scores of RT by cont 
gen include=1 if z_rt<-2 | z_rt>2
drop if include==1 //
drop include 

drop z_correct 
drop z_rt

*3. drop incorrect trials 
drop if correct==0 //

sort cont subnum trialnum
tabstat subnum, stat(mean sd n) by (cont) // their obs.
tabstat subnum, stat(n) by (cont)
tab subnum if cont==0 // 
tab subnum if cont==1 // 
tabulate subnum cont //
}
save "..\data_files\MDD_RSP_final_filters_applied_$apply_filters.dta", replace

***************************************************
*******DESCRIPTIVES*****


*LOOKING AT THE DiSTRIBUTION OF BDI BY CONDITION
destring bdi,replace
hist bdi
hist bdi, by(cont)
tabstat bdi, statistics (min max mean range sd) by(cont)
graph box bdi, over(cont)


*LOOKING AT THE DiSTRIBUTION OF HRSD BY CONDITION
destring totalhrsd17,replace
hist totalhrsd17
hist totalhrsd17, percent by(cont)
tabstat totalhrsd17, statistics (mean min max range sd) by(cont)
graph box totalhrsd17, over(cont)
hist cbdi, percent by(cont)
hist c_hrsd, percent by(cont)

*LOOKING AT ACCURACY BY CONDITION 
tabstat mean_correct, stat(mean sd n) by(cont)  
scatter mean_correct mean_rt 
scatter mean_correct totalhrsd17 
graph box mean_correct, over(cont)
hist mean_correct, by (cont)
 

*LOOKING AT THE DISTRIBUTION OF RT BY CONDITION
hist rt, by(cont) //
by cont, sort:tabstat rt, stat(mean sd n range)    //

hist mean_rt, freq 
hist mean_rt, freq by (cont) // 
hist mean_rt, percent by (cont) //




*********LOOKING AT DYNAMICS*******
local w 20
sort subnum trialnum
gen smoothed_rt = 0
forvalues i = 1/`w'{
	by subnum, sort: replace smoothed_rt = smoothed_rt + rt[_n-`i'] ///
		if trialnum > trialnum[_n-`i'] 
		
}
replace smoothed_rt = smoothed_rt / `w'
replace smoothed_rt = . if smoothed_rt == 0

by trialnum cont , sort: egen smoothed_rt_mean = mean(smoothed_rt)


scatter smoothed_rt trialnum if trialnum > `w' ///
, by(cont) msize(tiny) ||   connected smoothed_rt_mean trialnum if trialnum >`w' ///
 , by(cont) 



*rt X trialnum
scatter rt trialnum if rt<=800 & cont==1, msize(tiny) mcolor(orange) || lowess rt trialnum if rt<=800 & cont==1 || scatter rt trialnum if rt<=800 & cont==0, msize(tiny) mcolor(green) || lowess rt trialnum if rt<=800 & cont==0

*Looking at first 50 trials:
scatter rt trialnum if rt<=800 & cont==1 & trialnum<=50 , msize(tiny) || lowess rt trialnum if rt<=800 & cont==1  & trialnum<=50 || scatter rt trialnum if rt<=800 & cont==0  & trialnum<=50, msize(tiny) || lowess rt trialnum if rt<=800 & cont==0  & trialnum<=50
*Looking at first 20 trials
scatter rt trialnum if rt<=800 & cont==1 & trialnum<=20 , msize(tiny) || lowess rt trialnum if rt<=800 & cont==1  & trialnum<=20 || scatter rt trialnum if rt<=800 & cont==0  & trialnum<=20, msize(tiny) || lowess rt trialnum if rt<=800 & cont==0  & trialnum<=20
** looks like performance trend in both conds is similar in first 10 trials, i'll take a look on the trials where learning suppose to happen, maybe one source of the differences. // nope, nothing suspicious...
*Looking at trials 10=>T<=20 
scatter rt trialnum if rt<=800 & cont==1 & trialnum<=20 & trialnum>=10, msize(tiny) || lowess rt trialnum if rt<=800 & cont==1  & trialnum<=20  & trialnum>=10 || scatter rt trialnum if rt<=800 & cont==0  & trialnum<=20 & trialnum>=10, msize(tiny) || lowess rt trialnum if rt<=800 & cont==0  & trialnum<=20 & trialnum>=10
*Looking at last 50 trials, last 10 trials
scatter rt trialnum if rt<=800 & cont==1 & trialnum>=130 , msize(tiny) || lowess rt trialnum if rt<=800 & cont==1  & trialnum>=130 || scatter rt trialnum if rt<=800 & cont==0  & trialnum>=130, msize(tiny) || lowess rt trialnum if rt<=800 & cont==0  & trialnum>=130
*Last 10 trials
scatter rt trialnum if rt<=800 & cont==1 & trialnum>=170 , msize(tiny) || lowess rt trialnum if rt<=800 & cont==1  & trialnum>=170 || scatter rt trialnum if rt<=800 & cont==0  & trialnum>=170, msize(tiny) || lowess rt trialnum if rt<=800 & cont==0  & trialnum>=170

collapse totalhrsd17 hrsd8 bdi cont rt  mean_correct mean_correct_cont SDcorrect_cont mean_rt_cont mean_rt sd_rt_cont med mdd_mini, by (subnum)


reg mean_rt c.mean_correct##i.cont
regress, beta
margins cont, at(mean_correct=( .9 (3) .99 )) 
marginsplot, addplot(scatter mean_rt mean_correct if cont==0, msymb(Oh)  ///
		 xlabel( .9 (3) .99 ) xmtick( .9 (3) .99 )  || ///
	                 scatter mean_rt mean_correct if cont==1, msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Effectiveness Feedback"))) ///
					  title("{bf:RT as Predicted by Effectiveness Feedback & Accuracy level}") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Accuracy Level }") 
graph export "model_figure.png", replace width(3000)
cohend mean_rt cont


*1st MODEL: PREDICTORS: CENTERED HRSD , CONTROL, INT.
egen mean_hrsd= mean(totalhrsd17) // 20.5
gen c_hrsd=totalhrsd17-mean_hrsd
tabulate c_hrsd // range: -7 - 11

reg mean_rt c.c_hrsd##i.cont // 
regress, beta
*estat esize
estat vif

margins cont, at(c_hrsd=( -7 (3) 11 )) 
				 
contrast cont#c.c_hrsd, effects


*2nd MODEL : NOW WITH LOG_RT
gen logRT= log(mean_rt) 

reg logRT c.c_hrsd##i.cont, vce(cluster sub) // 
regress, beta
*estat esize
estat vif

margins cont, at(c_hrsd=( -7 (3) 11 )) 
marginsplot, addplot(scatter logRT c_hrsd if cont==0, msymb(Oh)  ///
		 xlabel( -7 (3) 11 ) xmtick( -7 (3) 11 )  || ///
	                 scatter logRT c_hrsd if cont==1, msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Control Feedback"))) ///
					  title("{bf:RT as Predicted by Control Feedback & Depression level}") ytitle("{bf:logRT (ms)}") ///
					 xtitle("{bf:Depression Severity (HRSD Score)}") 
graph export "model_figure_logRT.png", replace width(3000)



*BASING ONLY ON CURRENT DATA POINTS
*for 1st MODEL
reg mean_rt c.c_hrsd##i.cont if c_hrsd>-9 & c_hrsd<9 , vce(cluster sub) // 
regress, beta
*estat esize
estat vif
margins cont, at(c_hrsd=( -8 (3) 9 )) 
marginsplot, addplot(scatter mean_rt c_hrsd if cont==0, msymb(Oh)  ///
		 xlabel(-9(3)9) xmtick(-9(3)9)  || ///
	                 scatter mean_rt c_hrsd if cont==1, jitter(1) msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Effectiveness Feedback"))) ///
					  title("{bf:RT as Predicted by Effectiveness Feedback & Depression level}") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Depression Severity (HRSD Score)}")

					 
reg mean_rt c.c_hrsd##i.cont if c_hrsd>-9 & c_hrsd<9 , vce(cluster sub)
margins cont,at(c_hrsd=(-9(3)9)) 
marginsplot
marginsplot, xdimension(c_hrsd)recast(scatter) plotopts(mlabel(_margin)) 
					 
					 
					 
*for 2nd MODEL
reg logRT c.c_hrsd##i.cont if c_hrsd>-9 & c_hrsd<9 // 
regress, beta
*estat esize
estat vif
margins cont, at(c_hrsd=( -9 (3) 9 )) 
marginsplot, addplot(scatter logRT c_hrsd if cont==0, msymb(Oh)  ///
		 xlabel(-13(3)13) xmtick(-9(3)9)  || ///
	                 scatter logRT c_hrsd if cont==1, msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Control Feedback"))) ///
					  title("{bf:RT as Predicted by Control Feedback & Depression level}") ytitle("{bf:logRT (ms)}") ///
					 xtitle("{bf:Depression Severity (HRSD Score)}") 



margins, dydx(c_hrsd)  at(cont=(1 0) )
marginsplot, yline(0)
marginsplot, addplot(scatter mean_rt c_hrsd if cont==0, msymb(Oh))



mixed mean_rt c.c_hrsd##i.cont if c_hrsd>-9 & c_hrsd<7 || subnum: // 

*3rd REG MODEL, PREDICTORS: BDI (centered) ,CONDITION
*centering BDI 
egen mean_bdi=mean(bdi)
gen cbdi=bdi-mean_bdi
tabstat cbdi, statistics (min max mean range sd)

reg mean_rt c.cbdi##i.cont 
regress, beta
*estat esize  
margins cont, at(cbdi=( -18 (4) 31 )) 
marginsplot, addplot(scatter mean_rt cbdi if cont==0, msymb(Oh)  ///
		 xlabel(-18(4) 31) xmtick(-18 (4) 31)  || ///
	                 scatter mean_rt cbdi if cont==1, msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Effectiveness Feedback"))) ///
					  title("{bf:RT as Predicted by Effectiveness Feedback & Depression level}") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Depression Severity (BDI-II Score)}") 
graph export "model_figure.png", replace width(3000)
					 
/*
*Range restriction // same pattern					 
margins cont, at(cbdi=( -18 (4) 14 )) 
marginsplot, addplot(scatter mean_rt cbdi if cont==0, msymb(Oh)  ///
		 xlabel(-18 (4) 14) xmtick(-18 (4) 14)  || ///
	                 scatter mean_rt cbdi if cont==1, msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Effectiveness Feedback"))) ///
					  title("{bf:RT as Predicted by Effectiveness Feedback & Depression level}") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Depression Severity (BDI-II Score)}") 
graph export "model_figure.png", replace width(3000)
					 					 
*/*[ DEP Threshold: BDI >19 HRSD>14]

scatter mean_rt bdi if cont==1 & bdi>=19, msize(small) || scatter mean_rt bdi if cont==0 & bdi>=19 , msize(small) || lowess mean_rt bdi if cont==1 & bdi>=19 || lowess mean_rt bdi if cont==0 & bdi>=19

scatter mean_rt bdi if cont==1 & bdi>=19, msize(small) || scatter mean_rt bdi if cont==0 & bdi>=19, msize(small)|| lfit mean_rt bdi if cont==1 & bdi>=19 || lfit mean_rt bdi if cont==0 & bdi>=19



reg mean_rt c.cbdi##i.cont
margins cont,at(cbdi=(14(5)40)) 
marginsplot
marginsplot, xdimension(cbdi)recast(scatter) plotopts(mlabel(_margin)) 

pwcorr cbdi c_hrsd, sig
pwcorr cbdi c_hrsd if cont==1, sig
pwcorr cbdi c_hrsd if cont==0, sig

twoway (scatter cbdi c_hrsd if cont==1) || (scatter cbdi c_hrsd if cont==0) 

scatter cbdi c_hrsd if cont==1

scatter cbdi c_hrsd if cont==0
					 
					 
*4th MODEL- mean_rt by HRSD without Psychomotor retardation item
*first centering the new variable: HRSDnoPMR
gen HRSDnoPMR= totalhrsd17- hrsd8
egen mean_HRSDnoPMR= mean(HRSDnoPMR) // 20.05
gen c_HRSDnoPMR=HRSDnoPMR-mean_HRSDnoPMR
tabulate c_HRSDnoPMR // range: -8 -- 11
					 
reg mean_rt c.c_HRSDnoPMR##i.cont 
regress, beta
*estat esize
estat vif
margins cont, at(c_HRSDnoPMR=( -8 (3) 12 )) 
marginsplot, addplot(scatter mean_rt c_HRSDnoPMR if cont==0, msymb(Oh)  ///
		 xlabel(-8(4)12) xmtick(-8(3)12)  || ///
	                 scatter mean_rt c_HRSDnoPMR if cont==1, jitter(1) msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Control Feedback"))) ///
					  title("{bf:RT as Predicted by Control Feedback & Depression level") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Depression Severity (HRSD Score (NO PMR)}")


*5th MODEL - mean_rt by PMR score only (item 8 => hrsd8)
*recode PMR var as categorical
tabulate hrsd8 // 88 subj scored 0, 7 scored 1, 2 scored 2 
by subnum, sort: gen PMR = 0
replace PMR=1 if hrsd8>0	// 7				 

reg mean_rt i.PMR##i.cont 
regress, beta
*estat esize
estat vif
margins cont, at(PMR=(0 1)) 
marginsplot, addplot(scatter mean_rt PMR if cont==0, msymb(Oh)  ///
		 xlabel(0 1) xmtick(0 1)  || ///
	                 scatter mean_rt PMR if cont==1, jitter(1) msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Control Feedback"))) ///
					  title("{bf:RT as Predicted by Control Feedback & PMR") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Psychomotor retardation score}")

*6th MODEL 

reg mean_rt c.c_HRSDnoPMR##i.PMR 
regress, beta
*estat esize
estat vif
margins PMR, at(c_HRSDnoPMR=( -8 (3) 11 )) 
marginsplot, addplot(scatter mean_rt c_HRSDnoPMR if PMR==0, msymb(Oh)  ///
		 xlabel(-8 (3) 11) xmtick(-8 (3) 11)  || ///
	                 scatter mean_rt c_HRSDnoPMR if PMR==1, jitter(1) msymb(Oh)   ///
					 legend(order( ///
								  3 "No PMR" 4 "PMR"))) ///
					  title("{bf:RT as Predicted by PMR & HRSD score excluding PMR item") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Depression Level}")




*7th FULL MODEL 
reg mean_rt c.c_HRSDnoPMR##i.PMR##i.cont 
regress, beta
*estat esize
estat vif
margins cont, at(c_HRSD=(-8 (3) 13)) over(PMR)
marginsplot , bydim(PMR)
marginsplot, addplot(scatter mean_rt PMR if cont==0, msymb(Oh)  ///
		 xlabel(-8 (3) 13) xmtick(-8 (3) 13)  || ///
	                 scatter mean_rt PMR if cont==1, jitter(1) msymb(Oh)   ///
					 legend(order( ///
								  3 "No Feedback" 4 "Control Feedback"))) ///
					  title("{bf:RT as Predicted by Control Feedback & PMR") ytitle("{bf:RT (ms)}") ///
					 xtitle("{bf:Psychomotor retardation score}")

************************
*bdi||hrsd
scatter cbdi c_hrsd || lfit cbdi c_hrsd
scatter bdi hrsd || lfit bdi hrsd

pwcorr cbdi c_hrsd

*****************************
*****SELF REPORT*************
*****************************
*Callingg DATA - questions ratings

merge m:m subnum using "..\data_files\MDD_RSP_SelfReport_data.dta"
keep if _merge==3
drop _merge

order subnum cont, first
sort subnum cont
tab subnum cont


**
**q1-q2-q3: CONTROL -general control, behavior control, flash control
**
tabstat generalcontrol behaviorcontrol flashcontrol  , by(cont) stat(mean n sd min max) nototal // greater sense of flashcontrol over other control ques, in cont vs. no cont
histogram generalcontrol, by(cont)  normal //
histogram behaviorcontrol, by(cont)  normal
histogram flashcontrol, by(cont)  normal

* NO SIG DIFFERENCE BETWEEN CONTROL CONDITIONS
ttest generalcontrol, by(cont) //
ttest behaviorcontrol, by(cont) //
ttest flashcontrol, by(cont) //  NOTE! none of the no cont subjs got flash)

*by dep levels
tabstat generalcontrol behaviorcontrol flashcontrol, by(c_hrsd) stat(mean n sd min max) nototal // 
anova generalcontrol c.c_hrsd  
oneway generalcontrol c_hrsd, tabulate
oneway generalcontrol c_hrsd, bonferroni
  
anova behaviorcontrol c.c_hrsd   
oneway behaviorcontrol c_hrsd , tabulate
oneway behaviorcontrol c_hrsd , bonferroni
  
anova flashcontrol c.c_hrsd   
oneway flashcontrol c_hrsd, tabulate
oneway flashcontrol c_hrsd, bonferroni
  
anova behaviorcontrol c.c_hrsd cont // 

histogram generalcontrol, by(totalhrsd17)  normal
histogram behaviorcontrol, by(hrsd_class)  normal
histogram flashcontrol, by(hrsd_class)  normal

scatter generalcontrol totalhrsd17
scatter generalcontrol c_hrsd
scatter generalcontrol mean_rt
scatter behaviorcontrol c_hrsd
scatter flashcontrol c_hrsd

* 
scatter mean_rt generalcontrol 
scatter mean_rt behaviorcontrol
scatter mean_rt flashcontrol


*MEAN RT IS NOT EFFECTED BY ANY REPORT OF SENSE OF CONTROL except flashcontrol
reg mean_rt c.generalcontrol c.behaviorcontrol c.flashcontrol
reg mean_rt c.generalcontrol##i.hrsd_class
reg mean_rt c.behaviorcontrol##i.hrsd_class
reg mean_rt c.flashcontrol##i.hrsd_class

cibar generalcontrol, over1(cont) 
cibar behaviorcontrol, over1(cont) 
cibar flashcontrol, over1(cont) 

*******


*q4-q6: ENJOY (q4-you, q5-others)

*NO DIFFERENCES IN ENJOY BY CONTROL CONDITION
tabstat enjoy enjoyother , by(cont) stat(mean n sd min max) nototal //
histogram enjoy, by(cont)  normal //
histogram enjoyother, by(cont)  normal //
ttest enjoy, by(cont) //
ttest enjoyother, by(cont) // 

*NO DIFFERENCES IN ENJOY BY DEP LEVEL
anova enjoy c.c_hrsd
oneway enjoy c_hrsd, tabulate
oneway enjoy c_hrsd, bonferroni
anova enjoyother c.c_hrsd
oneway enjoyother hrsd_class, tabulate
oneway enjoyother hrsd_class, bonferroni
  
tabstat enjoy enjoyother, by(c_hrsd) stat(mean n sd min max) nototal // 

*MEAN RT IS NOT EFFECTED BY ANY REPORT OF ENJOY
reg mean_rt c.enjoy##c.c_hrsd
reg mean_rt c.enjoyother##i.hrsd_class 

reg mean_rt c.enjoy##i.cont c_hrsd //

twoway (scatter enjoy mean_rt if cond==1, mcolor("red)

cibar enjoy, over1(cont) 
cibar enjoyother, over1(cont) 
cibar enjoy, over1(hrsd_class)
cibar enjoyother, over1(hrsd_class)
 scatter mean_rt enjoy if cont==1 || lfit  mean_rt enjoy if cont==1, mcolor(red) ||scatter mean_rt enjoy if cont==0 || lfit  mean_rt enjoy if cont==0

*******
*q5_URGE
tabstat urge , by(cont) stat(mean n sd min max) nototal
histogram urge, by(cont)  normal // 
ttest urge, by(cont) // SIG DIFFERENCE BETWEEN COND. HIGHER URGE IN CONTROL ~x2. because NO-CONTROL didn't get any flash, they didn't feel the urge!!

*NO DIFFERENCES IN URGE BY DEP LEVEL
anova urge c.hrsd_class  
oneway urge hrsd_class, tabulate
oneway urge hrsd_class, bonferroni
tabstat urge, by(hrsd_class) stat(mean n sd min max) nototal // 

*MEAN RT IS NOT AFFECTED BY URGE &/OR DEP LEVEL 
reg mean_rt c.urge##c_hrsd

pwcorr urge mean_rt if cont==1, sig  // 
cibar urge, over1(cont)
cibar urge, over1(hrsd_class)
 
***********

*q7_EFFORT
*NO DIFFERENCES IN EFFORT BY CONDITION
tabstat effort , by(cont) stat(mean n sd min max) nototal
histogram effort, by(cont)  normal //
ttest effort, by(cont) // 

*NO DIFFERENCES IN EFFORT BY DEP LEVEL
anova effort c.c_hrsd  
oneway effort c_hrsd, tabulate
oneway effort c_hrsd, bonferroni
tabstat effort, by(c_hrsd) stat(mean n sd min max) nototal // 
histogram effort, by(c_hrsd)  normal // 
/*
*MEAN RT IS NOT EFFECTED BY EFFORT &/OR DEP LEVEL BUT TAKE A 2ND LOOK - for only Nofeed condition. 
reg mean_rt c.effort##c.c_hrsd // 

reg effort c.cont##c.c_hrsd // hrsd is marginally significant p=.055, int (contXhrsd) is sig, p=.015
hist 
pwcorr effort mean_rt, sig // NS
pwcorr effort mean_rt if c_hrsd>0 & cont==1, sig
pwcorr effort mean_rt if c_hrsd>0 & cont==0, sig // interesting HERE corr is affected by int 
scatter effort c_hrsd, over(cont)
scatter effort c_hrsd  if cont==1 || lfit effort c_hrsd if cont==1 || scatter effort c_hrsd  if cont==0 || lfit effort c_hrsd if cont==0 



*Only for people above mean dep score of ~20
scatter mean_rt effort if cont==0 || lfit mean_rt effort if c_hrsd>0 & cont==0  
scatter mean_rt effort if cont==1 || lfit mean_rt effort if c_hrsd>0 & cont==1 

cibar effort, over1(cont)
cibar effort, over1(hrsd_class)
*/ 
***********
*"WANTING" QUESTIONS 

*q8_ SUCCESS
*NO DIFFERENCES IN wanting to succeed BY CONDITION 
tabstat succeed , by(cont) stat(mean n sd min max) nototal
histogram succeed, by(cont)  normal // 
ttest succeed, by(cont) // 

*NO DIFFERENCES IN wanting to succeed BY DEP LEVEL // 
anova succeed c.c_hrsd  //NS
oneway succeed c_hrsd, tabulate
oneway succeed c_hrsd, bonferroni
tabstat succeed, by(c_hrsd) stat(mean n sd min max) nototal // 
histogram succeed, by(c_hrsd)  normal //

*MEAN RT IS NOT EFFECTED BY wanting to succeed &/OR DEP LEVEL 
reg mean_rt c.succeed##c.c_hrsd 

pwcorr succeed mean_rt, sig // NS
cibar succeed, over1(cont) 
cibar succeed, over1(hrsd_class) 

*q9_ Accuracy important
*NO DIFFERENCES IN wanting to be accurate BY CONDITION 
tabstat accurate , by(cont) stat(mean n sd min max) nototal
histogram accurate, by(cont)  normal
ttest accurate, by (cont) // 

*NO DIFFERENCES IN Accurate BY DEP LEVEL 
anova accurate c.c_hrsd  // 
oneway accurate hrsd_class, tabulate
oneway accurate hrsd_class, bonferroni
tabstat accurate, by(hrsd_class) stat(mean n sd min max) nototal // 
histogram accurate, by(hrsd_class)  normal // 

*MEAN RT IS NOT EFFECTED BY Accurate &/OR DEP LEVEL 
reg mean_rt c.accurate##c.c_hrsd  //
pwcorr accurate mean_rt, sig // 
cibar accurate, over1(cont) 
cibar accurate, over1(hrsd_class) 

*q10 quick
tabstat quick , by(cont) stat(mean n sd min max) nototal
histogram quick, by(cont)  normal
ttest quick, by (cont) // 

*NO DIFFERENCES IN Quick BY DEP LEVEL
anova quick c.c_hrsd  // NS
oneway quick hrsd_class, tabulate
oneway quick hrsd_class, bonferroni
tabstat quick, by(hrsd_class) stat(mean n sd min max) nototal // 
histogram quick, by(hrsd_class)  normal //

*MEAN RT IS NOT EFFECTED BY Quick &/OR DEP LEVEL 
reg mean_rt c.quick##c.c_hrsd // 
pwcorr quick mean_rt, sig //  
 
cibar quick, over1(cont) 
cibar quick, over1(hrsd_class) 


tabstat effort succeed accurate quick, by (cont) stat(mean sd)

(mean sd)


pwcorr mean_rt quick effort, sig

pwcorr mean_rt quick adhd , sig 

pwcorr mean_rt sex , sig

pwcorr mean_rt lang , sig

pwcorr mean_rt mood , sig 


by cont, sort: pwcorr general behavior flash enjoy enjoyother urge accurate succeed quick effort  mean_rt,sig star (.05)

*by , sort:pwcorr general behavior flash enjoy enjoyother urge accurate succeed quick effort  mean_rt,sig star (.05)

scatter urge enjoy if cont==1 || lfit  urge enjoy if cont==1  
scatter urge enjoy if cont==0 || lfit  urge enjoy if cont==0 

scatter succeed enjoy if cont==1 || lfit  succeed enjoy if cont==1||scatter succeed enjoy if cont==0 || lfit  succeed enjoy if cont==0

*END of FILE* 