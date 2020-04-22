## Shark bites analysis

Accompanies the article: BRADSHAW, CJA, P MEAGHER, MJ THIELE, RG HARCOURT, C HUVENEERS. In review. 
Predicting potential future reduction in shark bites on people. People and Nature

April 2020

Corey J. A. Bradshaw (corey.bradshaw@flinders.edu.au)
Global Ecology, Flinders University
GlobalEcologyFlinders.com

The R code attached reproduces an analysis to predict the number of people who could avoid being bitten by a shark in Australia from 2020-2066 if wearing electronic deterrents.

The analysis requires four different data files:

1. Australian Shark Attack File ('sharkbite.exp.csv') — for proprietry reasons, this dataset is only available upon request to Taronga Conservation Society Australia, Taronga Zoo, Sydney, New South Wales, Australia

2. Monthly southern oscillation index (soi) values from the Australian Bureau of Meteorology (BoM) ('soi.csv')

3. Monthly Pacific Decadal Oscillation (pdo) values from https://www.ncdc.noaa.gov/teleconnections/pdo/ ('pdo.csv'

4. Australian population size estimates (past) and projections to 2066 ('auspop.csv') from the Australian Bureau of Statistics

Also accompanying the code are two source files with additional functions necessary to reproduce the analyses ('r.squared.R'. & 'new_lmer_AIC_tables3.R')
