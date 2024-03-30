##code starts --------------------

import pandas as pd

# Read the CSV file
df = pd.read_csv('aa_fitness.csv.csv')

# Filter rows where gene is "S" and aa_site is between 331 and 524
filtered_df = df[(df['gene'] == 'S') & (df['aa_site'].between(331, 524))]

# Write the filtered data to a new CSV file
filtered_df.to_csv('site331-524_unfiltered.csv', index=False)

## code ends------------------------
