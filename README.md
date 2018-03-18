# denv_virulence_evolution
Code to reproduce analysis and figures from  "Transmission-clearance trade-offs indicate that dengue virulence evolution depends on epidemiological context"

Viral load data can be downloaded from the supplemental data from Clapham et al (2014) J. R. Soc. Interface. 

Code for formatting data and setting parameters: 
```python -m format data```

To load data and set frequently used parameters: 
```
find_viremia_data()
find_low_viremia_data()
'''

Code to run models and calculate confidence intervals: 
Main model fit to full dataset:
```
MLE_DENV_evolution_full()
CI = calculate_confidence_intervals(1); 
```

Main model fit to data subset:
```
MLE_DENV_evolution_subset_data()
CI = calculate_confidence_intervals(2); 
```

Innate immune response model fit to data subset (both primary and secondary infections cleared by the innate immune response): 
```
MLE_DENV_evolution_subset_data_innate()
CI = calculate_confidence_intervals(4); 
```

Alternative T-cell model fit to data subset:
```
MLE_DENV_evolution_subset_data_alt_T()
CI_a_1e6 = calculate_confidence_intervals(4); 
CI_a_1e7 = calculate_confidence_intervals(5);
```

Code to reproduce figures in the paper:
