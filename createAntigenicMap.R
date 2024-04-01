# Load the Racmacs package
library(Racmacs)

# Set an option for the number of computer cores to run in parallel when optimizing maps
# The default when running on CRAN is to use 2
options(RacOptimizer.num_cores = 2)
# However you can also set the number of cores to the maximum number like this
# options(RacOptimizer.num_cores = parallel::detectCores())

# Read in the titer table
titer_table<- read.titerTable("titer_table_20.csv")

map<-acmap(
	titer_table=titer_table
)

map<-optimizeMap(
	map=map,
	number_of_dimensions=2,
	number_of_optimizations=500,
	minimum_column_basis="none"
)

view(map)