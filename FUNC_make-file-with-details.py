import pandas as pd
import sys

def join_dataframes(file1, file2, anchor_details_file, output_file='combined.tsv'):
    """
    Join two dataframes on CHR, POS1, POS2 columns.
    
    Parameters:
    -----------
    file1 : str
        Path to first TSV file
    file2 : str
        Path to second TSV file
    anchor_details_file : str
        Path to anchor details TSV file
    output_file : str
        Path to output TSV file (default: 'combined.tsv')
    """
    # Read both files
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')
    anchor_details_df = pd.read_csv(anchor_details_file, sep='\t')
    
    # Merge on CHR, POS1, POS2
    combined_df = pd.merge(anchor_details_df, df1,  on=['CHR', 'POS1', 'POS2'], how='inner')
    combined_df.to_csv(output_file.replace('.tsv', '_HG2.tsv'), sep='\t', index=False)

    combined_df2 = pd.merge(anchor_details_df, df2, on=['CHR', 'POS1', 'POS2'], how='inner')
    combined_df2.to_csv(output_file.replace('.tsv', '_GM.tsv'), sep='\t', index=False)

    combined_df = pd.concat([combined_df, combined_df2], ignore_index=True)
    
    # Save to output file
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Successfully joined dataframes. Output saved to: {output_file}")
    print(f"Resulting dataframe shape: {combined_df.shape}")

if __name__ == "__main__":
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    anchor_details_file = sys.argv[3]
    output_file = sys.argv[4] if len(sys.argv) > 4 else 'combined2.tsv'
    
    join_dataframes(file1, file2, anchor_details_file, output_file)