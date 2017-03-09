# PIPL_population_models

file "PIPL_singlePop_script_revision.R" contains code for the single-population stochastic projection model used in the
decision support tool. It also contains code for simulations that were used to calculate the difference in population 
growth rate with and without exclosures as a function of different combinations of abandonment and predation probabilities.
These differences were exported as "predictionFrame.csv", with the corresponding abandonment and predation probabilities exported
in "predictors.csv", files that are included in the tool package to re-create the smoothing function and 3-D plot in the 
tool. 

file "metapopulation_5sites.R" contains script for a 5-site metapopulation model. Original file was last edited in July 2016,
I've made changes to the single-site population since then that may need to be incorporated into the metapopulation model.
