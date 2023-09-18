# MitoCartaEnrichment
Generaly available enrichment analysis packages such as GO and Reactome, while extensive, lack specific mitochondria-function related terms. Therefore, this package was created to help perform enrichment analysis based on MitoCarta 3.0 pathways and visualise the results. Currently, only a simple barplot visualisation is supported, however, more functions may be added in the future.

### Downloading the package
The most straightforward way to install would be to use the devtools function as follows and then load it into the environment:
```install_github("SSID08/MitoCartaEnrichment")```

```library(MitoCartaEnrichment)```

### Dependencies
The package relies on the following dependencies: 
    rlang (>= 1.1),
    stats (>= 4.2),
    AnnotationDbi (>= 1.6),
    stringr (>= 1.5),
    org.Hs.eg.db,
    dplyr (>= 1.1),
    ggplot2 (>= 3.4)

  These are automatically installed and imported in case they have not previously been installed

### Use
The main function which performs the enrichment is the run_enrichment() function which takes in as inputs the list of ENTREZ IDs of interest (could be DE genes or any other set of interest) and the list of ENTREZ IDs of the background set. This returns a dataframe with enrichment results filtered at the specified FDR threshold (q_threshold parameter) which is set to 0.05 by default. This function only works with ENTREZ IDs as inputs for reasons of backend integration. Therefore, functionality has been provided to convert other biological IDs to ENTREZ IDs through the ID_to_ENTREZ() function. The list of valid input IDs and corresponding input strings can be found using valid_input() function. The barplot_enrichment() function plots the results for the top n pathway terms sorted by the proportion of representation. Do not change the column names or structure of the enrichment results dataframe as this may cause the plot function to not work.

### Example use:
```valid_input() # Explore list of valid input IDs and associated strings```

```input_ENTREZ=ID_to_ENTREZ(significant_proteins,from='UNIPROT') # convert list of Uniprot IDs to ENTREZ IDs```

```background_ENTREZ=ID_to_ENTREZ(background_proteins,from="UNIPROT") # Convert list of background/universe proteins to ENTREZ IDs```

```enrichment_result=run_enrichment(x=input_ENTREZ,y=background_ENTREZ,q_threshold=0.1) # Run enrichment with custom FDR threshold of 0.1 ```

```enrichment_plot=barplot_enrichment(enrichment_result,n=10) # Plot barplot for the top 10 enriched terms```
