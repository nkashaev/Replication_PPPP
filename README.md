README Document for Reproducing Results in
==========================================
"Prices, profits, proxies, and production"
=============================================
Victor H. Aguiar
vaguiar@uwo.ca

Nail Kashaev
nkashaev@uwo.ca

Roy Allen
rallen46@uwo.ca

Data Sources
============

The paper uses one data set from

-   Epple, Dennis, Brett Gordon, and Holger Sieg (2010). A new approach to estimating
the production function for housing, American Economic Review, 100 (3), 905–24 <https://doi.org/10.1257/aer.100.3.905>

and `Julia` code from 

-   Kim, Youngseok, Peter Carbonetto, Matthew Stephens, and Mihai Anitescu (2020). A fast algorithm for maximum likelihood estimation of mixture proportions using sequential quadratic programming, Journal of Computational and Graphical Statistics, 29 (2), 261–273. <https://doi.org/10.1080/10618600.2019.1689985>

Software
========

A version 1.6.4 of the `Julia` programming language was used in coding the analysis files. For details about how to install `Julia` on different platforms and make `Julia` programs executable from the command line see <https://julialang.org/downloads/platform/>. After installation of `Julia 1.6.4.` run `using Pkg` and `Pkg.instantiate()` in the `Julia` terminal after setting the replication folder as the main one.

Some simulations use `KNITRO 12.3`.  

Hardware
========

The code was run on Mac mini (M1, 2020) with 16 Gb of RAM

Content
=======

-   `Application`  -- the folder contains the analysis files to replicate the results in Appendix C.

-   `Simulations`  -- the folder contains the analysis files to replicate the results in Online Appendix A.

-   `Manifest.toml` and `Project.toml`  -- toml files with all necessary `Julia` packages.


Below I describe the content of every folder.

`Application`
============

`data`
-----------

-   `data_cleaned_8.csv` -- processed data from Epple et al. (2010)

-   `Pittsburgh_post1995.txt` -- original data from Epple et al. (2010)

-   `US Zip Codes from 2013 Government Data.txt` -- data used in construction of the zip-code using geographical coordinates of houses

`results`
-----------

This folder contains estimated output levels and output prices.

-    `output_level_8_4.csv` -- estimates of output level

-    `output_level_8_4.csv` -- estimates of output price

`tables and graphs`
-----------

This folder contains all tables and figures from Appendix C.

-    `Fig_x.pdf` -- Figure x from Appendix C (x in {6,7,8})

-    `Tablex.pdf` -- Table x from Appendix C (x in {1,2})


`root files`
-----------

-    `application_functions.jl` -- the functions used in `application_main.jl`

-    `application_main.jl` -- this code generates the output levels and output prices (`output_level_8_4.csv` and `output_price_8_4.csv`)

-    `mixSQP.jl` -- the procedure for estimation of finite mixtures described in Kim et al. (2020)

-    `preparing_data_elbow.jl` -- this code cleanes the data from Epple et al. (2010) and constructs markets (`data_cleaned_8.csv`)

-    `tablesandgraphs.jl` --  this code generates all tables and figures in `tables and graphs` folder




`Simulations`
============

`results`
-----------
This folder constains simulation results used in construction of table in Online Appendix A
-    `thetap_DGP_n_ker_4.csv` -- simulation results for estimated \pi for sample size n with kernel ker (n in {500,1000,1500}; ker in {biweight, epanechnikov, triweight})

-    `thetarho_DGP_n_ker_4.csv` -- simulation results for estimated \rho for sample size n with kernel ker (n in {500,1000,1500}; ker in {biweight, epanechnikov, triweight})

`tables and graphs`
-----------
This folder constains all table in Online Appendix A
-    `Table_param_ker_4.csv` -- Table for estimated param with kernel ker (param in {\pi, \rho}; ker in {biweight, epanechnikov, triweight})


`root files`
-----------

-    `mixSQP.jl` -- the procedure for estimation of finite mixtures described in Kim et al. (2020)

-    `simulation_functions.jl` -- the functions used in `simulation_main.jl`

-    `simulation_main.jl` -- this code generates all simulation outputs in `results` folder

-    `simulation_all.sh` -- this script runs `simulation_main.jl` with different input parameters (e.g., kernel, sample size)

-    `tablesandgraphs.jl` --  this code generates all tables and figures in `tables and graphs` folder
