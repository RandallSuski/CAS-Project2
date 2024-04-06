from escapecalculator import EscapeCalculator 
import matplotlib.pyplot as plt
import numpy as np
import requests
import numpy
import pandas as pd
import yaml
import copy
import random 


class Mutation: 
    def __init__(self, original_acid, site, mutated_acid, fitness):
        self.original_acid = original_acid
        self.site = site 
        self.mutated_acid = mutated_acid
        self.fitness = fitness
    
    def toString(self):
        result = "{}{}{}".format(self.original_acid, self.site, self.mutated_acid)
        return result

    def notEqual(self, mutation):
        if (self.original_acid == mutation.original_acid and self.site == mutation.site and self.mutated_acid == mutation.mutated_acid):
            return False 
        else: 
            return True  


class Variant:
    def __init__(self, original_string):
        self.original_string = original_string
        self.mutations = []
        self.mutated_string = copy.copy(original_string)
        self.escape = 1.0
        self.max_fitness = -100; 
        self.min_fitness = 100;

    def mutate(self, mutation):
        self.mutations.append(mutation)
        self.mutated_string[mutation.site] = mutation.mutated_acid
        if mutation.fitness > self.max_fitness:
            self.max_fitness = mutation.fitness 
        if mutation.fitness < self.min_fitness:
            self.min_fitness = mutation.fitness 
    
    def mutateList(self, mutations):
        for mutation in mutations:
            self.mutate(mutation)
            #self.mutations.append(mutation)
            #self.mutated_string[mutation.site] = mutation.mutated_acid

    

class OriginalStrain:
    def __init__(self):
        self.strain = {}
    
    def loadStrain(self, fileName):
        # Read the CSV file
        df = pd.read_csv(fileName)

        # Put strain in a dictionary form site to acid
        for index, row in df.iterrows():
            site = row['site']
            acid = row['amino acid']
            self.strain[site] = acid 

class Acid: 
    def __init__(self, acid, fitness, expected):
        self.acid = acid 
        self.fitness = fitness 
        self.expected = expected 

class AcidSite:
    def __init__(self, site = None):
        self.site = site 
        self.acids = {} 
        self.length = 0

    def addAcid(self, acid, fitness):
        self.acids[acid] = fitness
        self.length += 1

    def getFitnessList(self):
        return self.fitness        



def main():
    # Load the XBB strain 
    xbb_strain = OriginalStrain()
    xbb_strain.loadStrain('filtered_XBB_clade.csv')
    #xbb_strain.loadStrain('filtered_BA2_clade.csv')

    # Load fitness values 
    df = pd.read_csv('filtered_aa_fitness.csv')
    fitness_map = {}
    for i in range(331, 525):
        fitness_map[i] = AcidSite(site = i)
    for index, row in df.iterrows():
        site = row['aa_site']
        fit = row['fitness']
        acid = row['aa']
        fitness_map[site].addAcid(acid = acid, fitness = fit)

    # Do 16 random steps to 
    random.seed(56)
    strain_index = 0 
    original = xbb_strain.strain
    mut_1_strains = []
    while strain_index < 16:
        # Choose random acid site
        acid_site = fitness_map[random.randint(331, 524)]
        site = acid_site.site 
        # Check if acid site has fitness values 
        if acid_site.length != 0:
            acid, fitness = random.choice(list(acid_site.acids.items()))
            # Check if choice is same acid and that the fitness value is within the threshold 
            if (original[site] != acid and acid != '*' and (-1 <= fitness <= 1)):
                # This is unique so we create mutation and variant
                mutation = Mutation(original_acid = original[site], site = site, mutated_acid = acid, fitness = fitness)
                variant = Variant(original_string = original)
                variant.mutate(mutation)
                # go to next variant 
                mut_1_strains.append(variant)
                strain_index += 1
    
    
    # Do 16 steps for each of those directions to create 256 strains
    mut_2_strains = []
    for strain_1 in mut_1_strains: 
        inner_strain_index = 0
        original = strain_1.mutated_string 
        while inner_strain_index < 16:
            # Choose random acid site
            acid_site = fitness_map[random.randint(331, 524)]
            site = acid_site.site 
            # Check if acid site has fitness values 
            if acid_site.length != 0:
                acid, fitness = random.choice(list(acid_site.acids.items()))
                # Check if choice is same acid and that the fitness value is within the threshold 
                if (original[site] != acid and acid != '*' and (-1 <= fitness <= 1)):
                    # This is unique so we create mutation and variant
                    mutation = Mutation(original_acid = original[site], site = site, mutated_acid = acid, fitness = fitness)
                    variant = Variant(original_string = original)
                    variant.mutateList(strain_1.mutations)
                    variant.mutate(mutation)
                    # go to next variant 
                    mut_2_strains.append(variant)
                    inner_strain_index += 1

    # Move 1 to 3 mutation's away from the neutral network  
    for variant in mut_2_strains: 
        mutation_count = 0
        mutation_number = random.randint(1,3)
        original = variant.mutated_string 
        while mutation_count < mutation_number:
            # Choose random acid site
            acid_site = fitness_map[random.randint(331, 524)]
            site = acid_site.site 
            # Check if acid site has fitness values 
            if acid_site.length != 0:
                acid, fitness = random.choice(list(acid_site.acids.items()))
                # Check if choice is same acid and that the fitness value is more than the threshold 
                if (original[site] != acid and acid != '*' and (fitness < 1 or fitness > 1)):
                    # This is unique so we create mutation and variant
                    mutation = Mutation(original_acid = original[site], site = site, mutated_acid = acid, fitness = fitness)
                    variant.mutate(mutation)
                    # go to next variant 
                    mutation_count += 1

    # Now find the escape values for these series of mutations 
    calculator = EscapeCalculator(virus="XBB", binds="Omicron BA.2/BA.5")
    #calculator = EscapeCalculator(virus="BA.2", binds="Omicron BA.2/BA.5")
    for variant in mut_2_strains:
        sites = [] 
        for mutation in variant.mutations:
            sites.append(mutation.site)
        variant.escape = calculator.binding_retained(sites).round(2)


    ### Analyze data ###
    #   Could also just make table with regions chossen with the highest escape 
    #   Graph which regions where choosen most?
    #   Sort by the most used 
    #   What did the best escape strains choose?
    #   How does this compare to already seen data?
    
    # Formating Data 
    mut_2_strains.sort(key=lambda x: x.escape)
    print("Variants Results: ")
    count = 0
    for variant in mut_2_strains:
        print(f"Variant{count}: {variant.escape}")
        count += 1
        for mutation in variant.mutations:
            print(f"  Mut: {mutation.original_acid}{mutation.site}{mutation.mutated_acid}  Fitness: {mutation.fitness}")
    
    # Sample the best 30 strains 
    sites = []
    frequency_map = {} 
    for i in range(30):
        variant = mut_2_strains[i]
        # For each mutation 
        for mutation in variant.mutations:
            freq = frequency_map.get(mutation.site, 0)
            frequency_map[mutation.site] = freq + 1    
    frequency_tuples = sorted(frequency_map.items(), key=lambda x: -x[1])
    sites = []
    frequency_values = []
    escape_values = []
    for i in range(10):
        site, freq = frequency_tuples[i]
        sites.append(site)
        frequency_values.append(freq)
        escape = calculator.binding_retained([site]).round(2)
        escape_values.append(1 - escape)

    ### Mutation Site Frequency ###
    # Create the bar plot
    bar_width = 0.2
    fig, axis1 = plt.subplots()
    axis2 = axis1.twinx()
    x_tics = np.arange(len(sites))
    bar1 = axis1.bar(x_tics            , frequency_values, width=bar_width, label='Frequency', color='tab:blue')
    bar2 = axis2.bar(x_tics + bar_width, escape_values   , width=bar_width, label='Escape'   , color='tab:green')

    # Add labels 
    axis1.set_title('XBB Frequency of Mutation Sites', fontsize=18)
    axis1.set_xlabel('Amino Acid Site', fontsize=14)
    axis1.set_ylabel('Site Frequency', fontsize=14)
    axis2.set_ylabel('Site Escape Value', fontsize=14)
    axis1.set_xticks(x_tics + bar_width / 2)
    axis1.set_xticklabels(sites)

    # Add a legend
    bars = [bar1, bar2]
    labels = [bar.get_label() for bar in bars]
    axis1.legend(bars, labels, fontsize = '12', loc='upper left')
    #plt.show()
    
    ### Variant Graph (strains on x, sort by max fitness) ###
    # Get the data 
    mut_2_strains.sort(key=lambda x: -x.max_fitness)
    fitness_list = []
    strain_indexes = [] 
    escape_list = []
    strain_list = []
    for i in range(10):
        # Process data for the strain and add it to data list 
        variant = mut_2_strains[i]
        fitness = variant.max_fitness 
        escape_list.append(1 - variant.escape)
        strain_indexes.append(i)
        strain_list.append(variant)
        fitness_list.append(fitness)

        #Print the mutations for this strain 
        mut_str_list = map(lambda x: x.toString(), variant.mutations)
        mutations_string = ' '.join(mut_str_list) 
        print(f"Strain: {i}  ", mutations_string)
 
    # Create the bar plot
    bar_width = 0.2
    fig, axis1 = plt.subplots()
    axis2 = axis1.twinx()
    x_tics = np.arange(len(strain_indexes))
    bar1 = axis1.bar(x_tics            , fitness_list, width=bar_width, label='Fitness', color='tab:blue')
    bar2 = axis2.bar(x_tics + bar_width, escape_list , width=bar_width, label='Escape'   , color='tab:green')

    # Add labels 
    axis1.set_title('XBB Strain Fitness Values', fontsize=18)
    axis1.set_xlabel('Strain Index', fontsize=14)
    axis1.set_ylabel('Site Fitness', fontsize=14)
    axis2.set_ylabel('Site Escape Value', fontsize=14)
    axis1.set_xticks(x_tics + bar_width / 2)
    axis1.set_xticklabels(strain_indexes)

    # Add a legend
    bars = [bar1, bar2]
    labels = [bar.get_label() for bar in bars]
    axis1.legend(bars, labels, fontsize = '12', loc='upper right')
    plt.show() 




if __name__ == "__main__":
    main()






















