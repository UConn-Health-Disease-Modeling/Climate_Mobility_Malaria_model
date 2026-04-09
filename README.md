# Model description
This repository contains a climate-driven malaria transmission model for the Brazilian Amazon, 
extended to incorporate human mobility between two regions. The framework captures how movement patterns 
interact with environmental variability—such as temperature and rainfall—to shape spatial transmission dynamics 
and disease burden across interconnected settings.

# Repository structure
data/
Contains climate data (temperature and rainfall) for the two study regions—Rural Manaus and Urban Macapá—used as inputs for the model simulations.
model_eqns_scenario_2.py
Implements the system of differential equations defining the malaria transmission model.
parameter_values_scenario.py
Provides parameter values used in the simulations.
Malaria_model_Urban_rural_Simple_trip_approach.ipynb
Compares model outputs with and without human mobility, and evaluates model fit against observed daily malaria case data for both regions.
Malaria_model_climate_mobilityscenario_analysis.ipynb
Explores the impact of different mobility scenarios on transmission dynamics.
Malaria_intervention_scenario.ipynb
Assesses the effects of intervention strategies under varying mobility and environmental conditions.
