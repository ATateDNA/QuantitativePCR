This script cleans and analyzes the amplification and Cq result files from BioRad Maestro software. It requires the plate to be a full 96-well plate, but can be modified for less. The analysis includes recreating the BioRad amplification graph for each well group (triple PCR replicates) organized in columns (A1, B1, C1...D1, E1, F1..etc.) The script will also merge sample name data (csv format) to the result files. As well as, preform a logic test for Cq threshold for each well and sum them for each well group (AKA each sample).
