
from escapecalculator import EscapeCalculator 

import requests

import numpy

import pandas as pd

import yaml

class SampleEscape:
    def __init__(self, sample_size = None):
        self.calculator = EscapeCalculator(virus="SARS", binds="Wuhan-Hu-1")
        self.seqs = pd.DataFrame.from_records(self.generateList(), columns=["name", "mutated sites"],)
        
    def runSample(self):
        self.calculate()
        self.countValid()
        fraction_valid = self.neutral_count / 194
        self.fraction_valid = round(fraction_valid, 2)
        print(self.fraction_valid)

    def countValid(self):
        binding_values = self.seqs["neutralization retained"]
        neutral_count = 0 
        threshold = 0.98
        for elem in binding_values:
            if elem >= threshold:
                neutral_count += 1
        self.neutral_count = neutral_count

    def calculate(self):
        self.seqs["neutralization retained"] = self.seqs["mutated sites"].map(self.calculator.binding_retained)

    def generateList(self):
        data = []
        for i in range(331, 525):
            name = "seq{}".format(i)
            site = [i]
            data.append((name, site))
        return data 
        







