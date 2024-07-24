# Code overview

## ASVCommunityModeling.ipynb

The code in this Notebook acquired the genome-scale models from KBase and computed the SMIPPs, MMIPPs, and community fluxes that are the core of the modeling results.   This Notebook uses the `agoracommutil.py` file, which implicity imports a [utility package](https://github.com/cshenry/chenry_utility_module/tree/master/lib/chenry_utility_module) that must be installed to run this Notebook.

## degenerate_models_statistics.ipynb

Experimental metabolomics and abundance data was processed that fed modeling in the ASVCommunityModeling Notebook.

## processing_CSVs.ipynb

The SMIPPs, MMIPPs, and community fluxes were processed into the t-SNE plots and Clustered Heatmaps of the Supplemental.

## escher_API_mapping.ipynb

The community fluxes were first visualized into an Escher Map via [the ModelSEED website](https://modelseed.org/escher/escher_builder.html#) and then the exported SVG was further refined via functions in this Notebook.  A few final manual edits resulted in the final figure of the paper.
