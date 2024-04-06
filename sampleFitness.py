
import pandas as pd

import random 




class AcidSite:
    def __init__(self, site = None):
        self.site = site 
        self.fitness = []
        self.expected = []

    def addAcid(self, fitness, expected):
        self.fitness.append(fitness)
        self.expected.append(expected)

    def getFitnessList(self):
        return self.fitness


def main():
    # Read the CSV file

    # Group together amino acids from same site 
    df = pd.read_csv('filtered_aa_fitness.csv')
    acid_map = {}
    for i in range(331, 525):
        acid_map[i] = AcidSite(site = i)

    for index, row in df.iterrows():
        site = row['aa_site']
        fit = row['fitness']
        exp = row['expected_count']
        acid_map[site].addAcid(fitness = fit, expected = exp)

    # Take a random sample of n = 100 amino acid changes
    random.seed(56)
    neutral_count = 0
    sample_count = 0 
    while sample_count < 100:
        site = acid_map[random.randint(331, 524)]
        fitness_list = site.getFitnessList()
        if len(fitness_list) != 0:
            fitness = random.choice(fitness_list)
            if (-1 <= fitness <= 1):
                neutral_count += 1
            sample_count += 1
    print("Neutral Count: ")
    print(neutral_count)

if __name__ == "__main__":
    main()



