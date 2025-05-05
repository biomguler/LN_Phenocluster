import sys
import ld_loader  # Assumes ld_loader.py (script 18) is in the same directory or in PYTHONPATH

# Check command-line arguments
if len(sys.argv) < 2:
    print("Usage: python 19_process_ld.py <ld_file_prefix>")
    sys.exit(1)

file_to_process = sys.argv[1]

# Load LD matrix and SNP metadata
df_R, df_ld_snps = ld_loader.load_ld_npz(file_to_process)

# Placeholder: add further processing or saving here if needed
print(f"Successfully processed LD file: {file_to_process}")
