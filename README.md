# Optimal-strategy-and-benefit-of-non-continuous-therapy
Codes and raw data for paper titled "Optimal strategy and benefit of non-continuous therapy depend on tumor heterogeneity and aggressiveness at the time of diagnosis"<br>

Data files included in the repository:<br>
<ul>
  <li>Fig1andSuppFig1-2/1drugsimulationdata_paper.csv:</li> -data from simulations of a single therapy in a cancer with two subpopulations: values of Moff, M, n0, pulse advantage, proportion of pulse advantages>1 for each set of parameters, sum(Moff*n0), sum(M*n0)
  <li>Fig2andSuppFig3/2drugsimulationdata_paper.csv:</li> -data from simulations of two therapies in a cancer with four subpopulations: values of Moff, MA, MB, n0, advantages of various treatment strategies, proportion of advantages>1 for each set of parameters for each treatment strategy, sum(Moff*n0), sum(MA*n0), sum(MB*n0)
  <li>Fig3andSuppFig4/Euler2drugsim.csv:</li>		-data from simulations of a single therapy in a cancer with two subpopulations using the Euler method: values of Moff, M, n0, pulse advantage, proportion of pulse advantages>1 for each set of parameters, sum(Moff*n0), sum(M*n0)
  <li>Fig3andSuppFig4/Euler1drugsim.csv:</li>	 -data from simulations of two therapies in a cancer with four subpopulations using the Euler method: values of Moff, MA, MB, n0, advantages of various treatment strategies, proportion of advantages>1 for each set of parameters for each treatment strategy, sum(Moff*n0), sum(MA*n0), sum(MB*n0)
  <li>Fig4andSuppFig5-7/tumormarkers_paper.csv:</li> -tumor marker value and day collected for each patient
  <li>Fig4andSuppFig5-7/treatmentdata_paper.csv:</li> -treatment type and days started and ended, reason for treatment cessation, and sequencing status for each patient for which we have tumor marker data
  <li>Fig4andSuppFig5-7/treatmenttable_paper.csv:</li> -treatment type and days started and ended, reason for treatment cessation, and sequencing status for all patients
</ul>

Notes: <br> -for file "Fig3andSupp4_examples.m", graph output will look slightly different from graph in the paper figure each time it is run, because each iteration of the code has a  random element.
