# Code overview

These Notebooks were executed in our linux Server environment (Python 3.11.1), although other operating systems and newer Python versions should work just as well.   Our simulations relied heavily on the `agoracommutil.py` util file, which implicity imports a [utility package](https://github.com/cshenry/chenry_utility_module/tree/master/lib/chenry_utility_module).   The installed packages of our Python kernel (not all necessarily being in these Notebooks) is presented in the `installed_packaged.txt` file.  This repository should take just a few minutes to clone via GitHub.  Some of the computations, particularly the community modeling simulations, may have trouble completing on a conventionaly desktop computer.  The output cells of these Notebooks generally reflect the expected outputs from running the code.

## ASVCommunityModeling.ipynb

The code in this Notebook acquired the genome-scale models from KBase and computed the SMIPPs, MMIPPs, and community fluxes that are the core of the modeling results.   

## degenerate_models_statistics.ipynb

Experimental metabolomics and abundance data was processed that fed modeling in the ASVCommunityModeling Notebook.

## processing_CSVs.ipynb

The SMIPPs, MMIPPs, and community fluxes were processed into the t-SNE plots and Clustered Heatmaps of the Supplemental.  The MMIPPs were then correlated with the metabolomics data to determine which compounds best explain community structure and metabolic behavior.

## escher_API_mapping.ipynb

The community fluxes were first visualized into an Escher Map via [the ModelSEED website](https://modelseed.org/escher/escher_builder.html#) and then the exported SVG was further refined via functions in this Notebook.  A few final manual edits resulted in the final figure of the paper.
