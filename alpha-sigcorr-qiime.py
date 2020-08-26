from qiime2 import Metadata, Artifact
from qiime2.plugins import diversity, feature_table
import pandas as pd
import os
import json

#Let's get set up

print('Let\'s get set up! Please provide the following information.\n')

metadata = input('\nWhat is your metadata filepath?\n')
feattable = input('\nWhat is your FeatureTable[Frequency] filepath?\n')
depth = int(input('\nWhat depth would you like to rarefy to?\n'))
tree = input('\nWhat is your phylogenetic tree filepath? (This must be a .qza file.)\n')

#Load the files as artifacts
metadata_API = Metadata.load(metadata)
feattable = Artifact.load(feattable)
tree = Artifact.load(tree)

#rarafy the table
rarefied = feature_table.methods.rarefy(table=feattable,sampling_depth=depth).rarefied_table

def sig_corr_export(metadata=metadata_API,rarefied_table=rarefied,phylogenetic_tree=tree):

	alphas = ['observed_otus', 'shannon', 'pielou_e', 'faith_pd']

	for alpha in alphas:
		#identify the directory each results will go to, establish a 'home'
	    path = os.getcwd() + '/' + alpha
	    go_back = os.getcwd()
	    
	    #faith requires a different method because it's phylogenetic
	    if alpha == 'faith_pd':
	    	globals()[alpha] = diversity.methods.alpha_phylogenetic(table=rarefied,
	    		phylogeny=tree,metric='faith_pd').alpha_diversity 
	    
	    #All others are just plain alpha
	    else:
	    	globals()[alpha] = diversity.methods.alpha(table=rarefied,
	    		metric=alpha).alpha_diversity
	    
	    #time to test for significance
	    globals()[alpha + '_sig'] = diversity.actions.alpha_group_significance(globals()[alpha], 
	                                                                           metadata_API)
	    #this exports the data into the path variable's directory
	    globals()[alpha + '_sig'].visualization.export_data(alpha)
	    
	    
	    df = pd.DataFrame()
	    variables = []
	    n=0

	    os.chdir(path)

	    #loop through all csvs in the output directory, concat into a dataframe 
	    for entry in os.listdir(path):
	        if entry.endswith(".csv"):
	            with open(entry) as f:
	                dfn = pd.read_csv(entry)
	                df = pd.concat([df,dfn],axis=1)

	            variables.append(str(entry)[24:-4])
	            n+=1

	        else:
	            continue

	    #creating a hierarchical index to say what results came from where
	    df.columns = pd.MultiIndex.from_product([variables, df.iloc[:,-5:]])

	    print('Added the results from {} significance tests for {}.'.format(n,alpha))

	    df.to_csv('../{}.csv'.format(alpha))
	    print('Exported csv for {}'.format(alpha))
	    
	    os.chdir(go_back)

	    #Time for correlation results
	    globals()[alpha + '_corr'] = diversity.actions.alpha_correlation(globals()[alpha], metadata_API)
	    globals()[alpha + '_corr'].visualization.export_data('{}_correlation'.format(alpha))
	    
	    path = os.getcwd() + '/{}_correlation'.format(alpha)

	    df = pd.DataFrame()
	    n=0

	    os.chdir(path)
	    for entry in os.listdir(path):
	    	#The QIIME devs made me sad by exporting correlation results as jsonp,
	    	#not even normal json, which python could read. This read my jsonps,
	    	#but I'm not sure if it'll work for everybody.
	        if entry.endswith('.jsonp'):
	            with open(entry) as f:
	                data = f.read()    
	                data = data.split('(')[1].strip(');')[18:]

	                try:
	                    json_data = json.loads('{' + data.split('},{')[2])
	                
	                #sometimes there is a ',null,' before what we want
	                except IndexError:
	                    json_data = json.loads('{' + data.split('null,{')[1])

	            dfn = pd.DataFrame.from_dict(json_data, orient='index',columns=[str(entry)[7:-6]])
	            df = pd.concat([df,dfn],axis=1)
	            n+=1

	        else:
	            continue


	    os.chdir(go_back)

	    print('Added results from {} correlation tests for {}.'.format(n,alpha))
	    
	    df.to_csv('{}_corr.csv'.format(alpha))

sig_corr_export(metadata=metadata_API,
	rarefied_table=rarefied,
	phylogenetic_tree=tree)

print('\nCongrats! This script has finished running.') 
print('\nI recommend removing the directories that have been created, but they have been left in case you want to use them.')
