# This file will involve the code that creates the neutral network, then create a titer table based
# off the neutral network created.

## ------- COPIED FROM THE NEUTRAL_NETWORK_CREATION.PY ------------------
import pandas as pd
import random

df = pd.read_csv('valid_mutations.csv')
randomSeed = 42

def mutation(site, amino_acid, fitness, row_in_valid_mut, list_of_mutations):
    return site, amino_acid, fitness, row_in_valid_mut, list_of_mutations


def populateMutations(num_mutations, list_to_populate, seed, list_invalid_rows):
    random.seed(seed)

    for i in range(num_mutations):
        mutation_choice = random.randint(0, 335)
        while mutation_choice in list_invalid_rows:
            mutation_choice = random.randint(0, 335)

        list_to_populate.append(mutation(
                                    df.iloc[mutation_choice, 1],
                                    df.iloc[mutation_choice, 2],
                                    df.iloc[mutation_choice, 3],
                                    mutation_choice, []))
        list_invalid_rows.append(mutation_choice)
    return list_to_populate

neutralNetwork = []
populateMutations(16, neutralNetwork, randomSeed, [])


for mutations in neutralNetwork:
    randomSeed = randomSeed + 1
    populateMutations(16, mutations[4], randomSeed, [mutations[3]])

## ------------- CODE FOR CREATING TITER TABLE -----------------------

# First, get mutations at the edge of the neutral network, for now choose 5

# Change this variable to change the number of mutations created to form the titer table
n = 20
list_of_mutated_strains = []

df2 = pd.read_csv('possible_mutations.csv')

for i in range(n):
    randomSeed = randomSeed + 1

    # Random walk first mutation
    firstStep = random.randint(0, 15)
    # Random walk second mutation
    secondStep = random.randint(0, 15)

    # Grab First mutation
    parent_mutation = neutralNetwork[firstStep]
    # Grab Second mutation
    child_mutation = parent_mutation[4][secondStep]

    # Inform what the 3rd mutation cannot be
    not_valid_sites = [parent_mutation[0], child_mutation[0]]
    # 1293 is the number of all possible mutations
    final_mutation = random.randint(0, 1292)
    while df2.iloc[final_mutation, 1] in not_valid_sites:
        final_mutation = random.randint(1, 1293)

    final_mutation_site = df2.iloc[final_mutation, 1]
    list_of_mutated_strains.append((parent_mutation[0], child_mutation[0], final_mutation_site))


# Create distance table based on the list of mutated strains

from escapecalculator import EscapeCalculator

calc = EscapeCalculator(virus="SARS", binds="Wuhan-Hu-1")

column_row_names = [f'M{i:02}' for i in range(1, n+1)]

titer_table = pd.DataFrame(0, index=column_row_names, columns=column_row_names)

# Get the escape values for every mutation
escape_values = []
for item in list_of_mutated_strains:
    escape_values.append( (1 - calc.binding_retained([item[0], item[1], item[2] ])).round(3) )

# Now we have escape values for every mutation, so now we can fill in the distance table
for i in range(len(list_of_mutated_strains)):
    # Iterate through every mutation
    for j in range(len(list_of_mutated_strains)):
        distance = abs(1 - (escape_values[i] / escape_values[j]))
        titer = (2 ** (10 - distance)).round()
        titer_table.at[column_row_names[i], column_row_names[j]] = titer

titer_table.to_csv('titer_table_' + str(n) + '.csv')

# Need to calculate the fitness of the variants as well and store them in a csv file

# Will take the mean fitness of the sites mutated since we can't get a specific amino acid
# from the escape calculator

fitness_values = []
probability_values = []
expected_count = []
site_names = []

probability_no_mutation_values = []
probability_mutation_values = []

starting_genome = pd.read_csv('starting_aa_RBD.csv')

# Calculate the probability of each site in RBD mutating
for i in range(191):
    site_of_interest = starting_genome.iloc[i, 1]
    site_mutations = df2[ df2['aa_site'] == site_of_interest ]
    num_mutations = site_mutations['expected_count'].sum()
    num_no_mutations = starting_genome.iloc[i, 4]
    probability_of_mutation = (num_mutations) / (num_no_mutations + num_mutations)
    probability_mutation_values.append(probability_of_mutation.round(3))
    probability_no_mutation_values.append(1 - probability_of_mutation.round(3))



for item in list_of_mutated_strains:
    site_names.append(str(item[0]) + ', ' + str(item[1]) + ', ' + str(item[2]))

    # Calculate probability of this mutation occuring
    probability_of_this_mutant = 1
    for i in range(len(probability_mutation_values)):
        if i == item[0] or i == item[1] or i == item[2]:
            probability_of_this_mutant *= probability_mutation_values[i]
        else:
            probability_of_this_mutant *= probability_no_mutation_values[i]
    probability_values.append('{:g}'.format(float('{:.3g}'.format(probability_of_this_mutant))))


    # get the mean fitness of each mutated site
    newdf = pd.read_csv('site331-524_unfiltered.csv')
    item_fitness = 0
    for i in range(3):
        # Calculate Mean fitness of each site
        siteidf = newdf[newdf['aa_site'] == item[i]]
        fitness_sum = siteidf['fitness'].sum()
        fitness_average = fitness_sum / siteidf.shape[0]
        item_fitness = (item_fitness + fitness_average).round(3)
    fitness_values.append(item_fitness)

    expected_count.append('{:g}'.format(float('{:.3g}'.format(probability_of_this_mutant * 10**10))))

data = {
    'Mutation': ['M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20'],
    'Mutated Sites': site_names,
    'Mean Fitness': fitness_values,
    'Probability of Mutation': probability_values,
    'Excpected Count': expected_count
}

mutation_data = pd.DataFrame(data)
mutation_data.to_csv('antigenic_mutation_data.csv', index=False)
