##code starts --------------------

import pandas as pd

# Read the CSV file
df = pd.read_csv('site331-524_unfiltered.csv')

# Filter rows where gene is "S" and aa_site is between 331 and 524
filtered_df = df[(df['fitness'].between(-1, 1)) & (df['fitness'] != 0.0)]
starting_aa_df = df[df['fitness'] == 0.0]

# Write the filtered data to a new CSV file
filtered_df.to_csv('valid_mutations.csv', index=False)
starting_aa_df.to_csv('starting_aa_RBD.csv', index=False)

## code ends------------------------