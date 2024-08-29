#Summary of LDSC RG
# Define the list of phenotypes
phenotypes=("Cell_B" "Cell_P" "CLL" "DLBCL" "Drug_G1" "FL" "HL" "LM" "LPL_WM" "MGUS" "MM_MGUS" "MM" "MZL" "Soma_G1" "Soma_G2")

# Iterate over each phenotype
for pheno in "${phenotypes[@]}"; do
    # Example filename pattern: Cell_B_binary9_ldsc_rg.log, Cell_B_ordinal_ldsc_rg.log, Cell_B_irdn_ldsc_rg.log
    pattern="${pheno}_*ldsc_rg.log"

    # Create a file to store the merged output for the current phenotype
    output_file="${pheno}_merged_output.txt"
    touch "$output_file"

    # Extract data and append to the output file for each phenotype
    for file in $pattern; do
        awk -F '\t' '/Summary of Genetic Correlation Results/ {flag=1; next} {if(flag) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$file" >> "$output_file"
    done
done

