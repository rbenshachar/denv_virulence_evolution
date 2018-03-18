import pandas as pd
import numpy as np

def main():
	xls = pd.ExcelFile("rsif20140094supp.xlsx")
	df = xls.parse('NS1&ViralLoadDataSampleBySample')

	# formatting time in days

	a = pd.to_datetime(df['time_since_illness_onset'], errors='coerce')
	# values with day < 0 are given an arbitrary date with day 1
	b = pd.to_datetime(df['time_since_illness_onset'], format='%H:%M:%S', errors='coerce')
	c = b.combine_first(a)
	df['time_since_illness_onset'] = c.dt.day + c.dt.hour/24. + c.dt.minute/60./24

	# checking for values where day = 0 in original data and setting day = 0 
	# (code above was parsing these values to have day 1)
	x = b[~b.isnull()].index
	df['time_since_illness_onset'].apply(int)
	df.loc[x, 'time_since_illness_onset'] = df.loc[x, 'time_since_illness_onset'] - 1

	# only want positive samples
	df = df[df.serology == 'acute dengue infection']
	
	#columns to include
	columns_subset = ['pt_no', 'serotype', 'genotype1_gp',
			   'primary_vs_secondary', 'diagnosis_who',
			   'time_since_illness_onset', 'vl_ml_plasma',
			   'ld_vl_ml_plasma']
	df = df[columns_subset]

	# renaming values for serotype
	df.loc[df.serotype=='DENV1', 'serotype'] = 1
	df.loc[df.serotype=='DENV2', 'serotype'] = 2
	df.loc[df.serotype=='DENV3', 'serotype'] = 3
	df.loc[df.serotype=='DENV4', 'serotype'] = 4
	df.loc[df.serotype=='neg', 'serotype'] = 0

	# renaming values for primary vs. secondary
	df.loc[df.primary_vs_secondary=='primary', 'primary_vs_secondary'] = 1
	df.loc[df.primary_vs_secondary=='secondary', 'primary_vs_secondary'] = 2
	df.loc[df.primary_vs_secondary=='indeterminate', 'primary_vs_secondary'] = 3
	df.loc[df.primary_vs_secondary=='neg', 'primary_vs_secondary'] = 0

	# renaming diagnosis WHO
	df.loc[df.diagnosis_who =='DF', 'diagnosis_who'] = 1
	df.loc[df.diagnosis_who=='DHF I', 'diagnosis_who'] = 2
	df.loc[df.diagnosis_who=='DHF I ', 'diagnosis_who'] = 2
	df.loc[df.diagnosis_who =='DHF II', 'diagnosis_who'] = 2
	df.loc[df.diagnosis_who =='DHF III', 'diagnosis_who'] = 2

	#renaming vl_ml_plasma BALD, FUS, PB to nan
	df.loc[df.vl_ml_plasma=='BALD', 'vl_ml_plasma'] = 5
	df.loc[df.vl_ml_plasma=='FUS', 'vl_ml_plasma'] = np.nan
 	df.loc[df.vl_ml_plasma=='PB', 'vl_ml_plasma'] = np.nan

 	#renaming ld_vl_ml_plasma FUS, PB to nan
 	df.loc[df.ld_vl_ml_plasma=='FUS', 'ld_vl_ml_plasma'] = np.nan
 	df.loc[df.ld_vl_ml_plasma=='PB', 'ld_vl_ml_plasma'] = np.nan

	df = df.set_index('pt_no')

	df.to_excel('Clapham_data_formatted.xlsx')

	return df
