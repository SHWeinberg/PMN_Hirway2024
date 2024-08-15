Pre Metastatic Niche Project Guide


CPM_v5.4-Final Version that was used for PMN manuscript
Base Grids folder contains lower and higher activation level confluent monolayers, that can be used for creating cell seeding masks later

BeforeandAfterTumorScenario.m- main script that runs a grid of Kupffer and stellate cells in phase 1 and grid with addition of tumor cells in phase 2. Can run parameter studies as well
Has a IDval- simulation ID
Loads Epithelial and mesenchymal cell phenotype for cells of lower and higher activation levels. Also, loads tumorcells (number of tumor cells), 
Randomly seeds Kupffer and Stellate cells with popdiff defining cell density, and gives lower activation level to the cells
Creates cell contact tables- Jcm, Jcc and Jca
For each table, creates a 3x3 table, for [Kupffer, Stellate, Tumor]
For Kupffer cells, uses scaling factor for decreased Growth Factor sensitivity using Paramvalue-14 which is snail degradation rate
Also, less proliferation for kupffer cells
Stores variables in ‘othermod’ struct and feeds into model
Runs this phase 1 simulation
Then, in phase 2
Loads the phase 1 simulation- grid and cell variables
Adds tumor cells as dots to the center of the grid and adds tumor cell variables as well

Writes all parameters and structs into IDval_params.xlsx file

Init_cells_v4.m- uses cellvars, matrvars and statevars that are inputted into the run file

JcApheno,JcMpheno,JcCpheno.m- uses values from the cell contact tables

Make_movie_or_snapshots_v3.m- creates additional movies
Creates a movie for tracking all cell types over time, and also individually, Kupffer, Stellate and Tumor Cells over time.
Cell_Types- variables- map_cells3, Mt, cfig3
KType-map_cellsK,MK,cfig4
SType-map_cellsS,MS,cfig5
TType-map_cellsT,MT,cfig6

Metastasis_Metrics.m- visualizes everything in the PMN manuscript
For parameter studies
Load a set of parameter studies based on the IDval 
Phase 1-
Loops through all the Before simulation files and creates a graph for Kupffer and Stellate activation and ECM concentration at the end of the simulation
Phase 2-
Loops through all the After simulation files and creates graphs for Tumor cell phenotype, ECM concentration underneath tumor cells, and Metastasis metrics, including number of non-tumor cell connections and distance between initially seeded cells




