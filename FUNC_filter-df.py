import sys
import pandas as pd

if len(sys.argv) < 4:
    print("Usage: python filter-df.py <input_file> <field_name> <threshold>")
    sys.exit(1)

input_file = sys.argv[1]
field_name = sys.argv[2]
threshold = float(sys.argv[3])

df = pd.read_csv(input_file, sep='\t')
df_filtered = df[df[field_name] > threshold]

df_filtered.to_csv(input_file, sep='\t', index=False)

print(f"Filtered {len(df)} rows down to {len(df_filtered)} rows")