from qiime2 import Metadata, Artifact
from qiime2.plugins import feature_table, diversity
#from qiime2.plugins.diversity.pipelines import core_metrics_phylogenetic
from qiime2.plugins.diversity.visualizers import beta_group_significance
from qiime2.plugins.diversity.actions import beta_correlation
import numpy as np
import pandas as pd
import os
from bs4 import BeautifulSoup





print('Let\'s get set up! Please provide the following information.\n')

metadata = input('\nWhat is your metadata filepath?\n')
feattable = input('\nWhat is your FeatureTable[Frequency] filepath?\n')
depth = int(input('\nWhat depth would you like to rarefy to?\n'))
tree = input('\nWhat is your phylogenetic tree filepath? (This must be a .qza file.)\n')
#corr_test = input('\nWould you rather run spearman or pearson correlation test?\n').lower() functionality to be included later
#verbose functionality to be included later

feattable = Artifact.load(feattable)
tree = Artifact.load(tree)

metrics = ['jaccard','bray_curtis','unweighted_unifrac','weighted_unifrac']


def validate_metadata(file):
    """This function checks for any issues in Keemei-validated metadata that currently raise errors.
        Issues as of 08/03/2020
        1. Columns containing '/' are interpreted as leading to an extra level directory.
        2. Something is wack with ':' as well. """
    df = pd.read_csv(file, sep='\t')
    
    for col in df.columns:
        try:
        	if df[col].str.contains('/').any() or df[col].contains(':').any():
        		return False
        	else:
        		return True
        except AttributeError:
        	pass


if not validate_metadata(metadata):
	print('The /\'s or :\'s in your metadata might cause issues. Please remove them.')

metadata_API = Metadata.load(metadata)


core_metrics = diversity.pipelines.core_metrics_phylogenetic(table = feattable,
                                                             phylogeny = tree,
                                                             sampling_depth = depth,
                                                             metadata = metadata_API)




def beta_sig_exports(metadata=metadata_API, core_metrics=core_metrics):
	home = os.getcwd()
	os.makedirs('beta_sig_results', exist_ok=True)
	os.chdir('beta_sig_results')


	for metric in metrics:
	    #fresh dataframe per metric
	    df = pd.DataFrame()
	    col_index = []
	    
	    print(f'Starting {metric}')

	    for column in metadata_API.filter_columns(column_type='categorical').columns:
	    	output_dir = f'{metric}_sig_{column}'
	    	matrix = getattr(core_metrics, f'{metric}_distance_matrix')
	    	print(f'Testing group significance for {column}')
	    	try:
        		significance = beta_group_significance(distance_matrix=matrix,metadata=metadata_API.get_column(str(column)))
        		significance.visualization.export_data(output_dir)

        		with open(f'{output_dir}/index.html') as f:
        			soup = BeautifulSoup(f, 'html.parser')

        		keys = [key.string for key in soup.find_all('th')[2:]]
        		values = np.array([value.string for value in soup.find_all('td')])
        		values.reshape((1,7))

        		col_index.append(column)

        		#Make the data frame	
        		dfn = pd.DataFrame(values,index=keys)
        		df = pd.concat([df,dfn],axis=1)

    		except ValueError:
        		print(f'{column} could not be added, as each sample is in a different category.\nMoving to the next column.')
	    	
	    df.columns = col_index
	    df.to_csv(f'{metric}_sig_results.csv')
	    
	os.chdir(home)
	print('Group significance testing has finished running. Congrats!')
	    


def beta_corr_exports(metadata=metadata_API, core_metrics=core_metrics):

	home = os.getcwd()
	os.makedirs('beta_corr_results',exist_ok=True)
	os.chdir('beta_corr_results')

	for metric in metrics:
	    print(f'Starting {metric}')
	    	    
	    matrix = getattr(core_metrics, f'{metric}_distance_matrix')
	    
	    df = pd.DataFrame()
	    col_index = []
	    
	    for column in metadata_API.filter_columns(column_type='numeric').columns:    
	        output_dir = f'{metric}_corr_{column}'
	        print(f'Testing correlation for column {column}')
	        
	        try:
	            correlation = beta_correlation(metadata=metadata.get_column(column).drop_missing_values(),
	                                           distance_matrix=matrix,
	                                           intersect_ids=True)

	            correlation.mantel_scatter_visualization.export_data(output_dir)

	            with open(f'{output_dir}/index.html') as f:
	                soup = BeautifulSoup(f, 'html.parser')

	            keys = [k.string for k in soup.find_all('th')[2:]]
	            values = np.array([v.string for v in soup.find_all('td')]).reshape((len(keys),1))

	            dfn = pd.DataFrame(values,keys)
	            df = pd.concat([df,dfn],axis=1)
	            
	            col_index.append(column)
	            
	            
	        except ValueError:
	            print('This column was empty, or there was another ValueError. Moving to the next column.')
	           

	    df.columns = col_index
	    df.to_csv(f'{metric}_corr_results.csv')
	    
	print('Correlation testing has finished running. Congrats!')
	os.chdir(home)


beta_sig_exports()

beta_corr_exports()