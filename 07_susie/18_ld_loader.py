import numpy as np
import pandas as pd
import scipy.sparse as sparse
import os

def load_ld_npz(ld_prefix):
    '''
    ld_prefix is the prefix of a pair of LD files (e.g. "/path/to/chr10_102000001_105000001")
    
    This function:
    1. Loads SNP metadata from a `.gz` file
    2. Loads the LD matrix from a `.npz` file
    3. Saves the DataFrames as **Feather and Parquet** (fast formats for large data)
       - LD matrix as `R_<prefix>.feather` and `R_<prefix>.parquet`
       - SNP metadata as `S_<prefix>.csv` (since it's small)
    
    Returns:
        df_R (pd.DataFrame): LD matrix
        df_ld_snps (pd.DataFrame): SNP metadata
    '''
    
    # Extract filename (without path) for output naming
    file_prefix = os.path.basename(ld_prefix)

    # Load SNP metadata
    gz_file = f"{ld_prefix}.gz"
    df_ld_snps = pd.read_csv(gz_file, sep=r'\s+', engine='python')

    # Rename columns
    df_ld_snps.rename(columns={'rsid': 'SNP', 'chromosome': 'CHR', 'position': 'BP', 
                               'allele1': 'A1', 'allele2': 'A2'}, inplace=True, errors='ignore')

    # Ensure required columns are present
    required_cols = {'SNP', 'CHR', 'BP', 'A1', 'A2'}
    if not required_cols.issubset(df_ld_snps.columns):
        raise ValueError(f"Missing required columns in {gz_file}: {required_cols - set(df_ld_snps.columns)}")

    # Create SNP index
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + ':' + df_ld_snps['BP'].astype(str) + ':' + \
                       df_ld_snps['A1'] + ':' + df_ld_snps['A2']
        
    # Load the LD matrix
    # Check which file exists
    if os.path.exists(f"{ld_prefix}.npz"):
        npz_file = f"{ld_prefix}.npz"
    elif os.path.exists(f"{ld_prefix}.npz2"):
        npz_file = f"{ld_prefix}.npz2"
    else:
        raise FileNotFoundError(f"No LD matrix file found for prefix {ld_prefix} (.npz or .npz2)")

    try: 
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except Exception as e:
        raise IOError(f"Error loading {npz_file}: {str(e)}")

    # Create DataFrame for LD matrix
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)

    # File names for saving
    ld_matrix_feather = f"R_{file_prefix}.feather"
    snp_metadata_csv = f"S_{file_prefix}.csv"

    # Save DataFrames
    df_ld_snps.to_csv(snp_metadata_csv)  # SNP metadata is small, so CSV is fine
    #df_R.to_feather(ld_matrix_feather)  # Faster format
    df_R.reset_index().to_feather(ld_matrix_feather) # to keep rownames



    print(f"Saved LD matrix to {ld_matrix_feather}")
    print(f"Saved SNP metadata to {snp_metadata_csv}")

    return df_R, df_ld_snps
