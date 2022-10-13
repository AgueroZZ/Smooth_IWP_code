# Smooth_IWP_code
The online repository for the codes to replicate the results in "Model-based Smoothing with Integrated Wiener Processes and Overlapping Splines"

## Data and Scripts

The simulation/examples in the main paper are replicated in the following scripts. Specific points to pay attention to:
1. In each simulation/example, please make sure the working directory is correctly specified.
2. All the required R packages will be loaded in the R script, please make sure they are installed in the system before proceed.
3. Before running the R script for the analysis, please make sure each cpp file is compiled succesfully in the working directory.

- **Section 4.1**:  
     - *script.R* produces the corresponding figures to *figures* directory.

- **Section 4.2**:  
     - *script1.R* reproduces the single replication analysis and store figures into *figures*.
     - *script2.R* produces the computational comparision results in the corresponding section.

- **Section 4.3**:  
     -*01_single_rep.R* implements the single replication comparision and produces the corresponding figures.
     -*02_aggregation.R* implements the aggregated comparision between the proposed method with the existing method, and produces outputs in the *result* directory.
     -*03_plot.R* uses the output from *02_aggregation.R* to produce the figures in the main paper.
     

- **Section 5**:
     - *owid-covid-data.csv* is the raw data used in the analysis in section 5 provided by Our World in Data, downloaded from https://github.com/owid/covid-19-data/tree/master/public/data.
     - *script.R* uses the data *owid-covid-data.csv* to perform the analysis in section 5, and store the corresponding figures in *figures*.

