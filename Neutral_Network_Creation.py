# This file will be the one that collects data and creates the neutral network

# Neutral Network will be represented by a list of mutations, a mutation will
# be an object of the form (site, aa, fitness, [list of mutations]).

import pandas as pd
import random

df = pd.read_csv('valid_mutations.csv')
randomSeed = 42

# Creates a mutation object
# site - int representing which site is mutated
# amino_acid - amino acid that site is mutating into
# fitness - the fitness of this mutation
# row_in_valid_mut - the row # in the valid_mutations.csv file
# list_of_mutations - a list of mutations that would further mutate the genome
def mutation(site, amino_acid, fitness, row_in_valid_mut, list_of_mutations):
    return site, amino_acid, fitness, row_in_valid_mut, list_of_mutations

# This function will take a list and populate it with num_mutations number of mutations
# It will collect mutations by randomly selecting them from valid_mutations.csv
def populateMutations(num_mutations, list_to_populate, seed, list_invalid_rows):
    random.seed(seed)

    for i in range(num_mutations):
        # Choose a random row from csv file, that isn't in the list of invalid choices
        mutation_choice = random.randint(0, 335)
        while mutation_choice in list_invalid_rows:
            mutation_choice = random.randint(0, 335)

        # We now have a valid mutation
        list_to_populate.append(mutation(
                                    df.iloc[mutation_choice, 1],
                                    df.iloc[mutation_choice, 2],
                                    df.iloc[mutation_choice, 3],
                                    mutation_choice, []))
        list_invalid_rows.append(mutation_choice)
    return list_to_populate

# Create our neutral network

# The first 16 neutral mutations
neutralNetwork = []
populateMutations(16, neutralNetwork, randomSeed, [])

# Populate each of those with 16 mutations
for mutations in neutralNetwork:
    # If we use the same seed, the same rows to grab mutations from will always be chosen, this allows
    # repeatability without ruining the randomness.
    randomSeed = randomSeed + 1
    populateMutations(16, mutations[4], randomSeed, [mutations[3]])

# Now we have the neutral network of 256 mutations, use graph-tool to graph them

from graph_tool.all import *

# Create Graph
ug = Graph(directed=False)

vlabel = ug.new_vp("string")

# Create root vertex
r = ug.add_vertex()
vlabel[r] = "0"

# Create verticies for all the mutations
def makeMap(root_vertex, neutral_network):
    for mutational in neutral_network:
        # Create Vertex for mutation
        v = ug.add_vertex()
        # Add edge between parent and vertex
        e = ug.add_edge(v, root_vertex)
        # Create the label
        vlabel[v] = str(mutational[0]) + " " + str(mutational[1])

        # If there are nested mutations, add them to the graph
        if len(mutational[4]) != 0:
            makeMap(v, mutational[4])

makeMap(r, neutralNetwork)
graph_draw(ug, vertex_text=vlabel, output="Neutral_Network.pdf")