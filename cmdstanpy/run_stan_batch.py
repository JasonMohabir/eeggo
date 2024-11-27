import os
import numpy as np
from cmdstanpy import CmdStanModel
import pandas as pd
import argparse
from datetime import datetime

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run a Stan model on gene expression data in chunks.")
parser.add_argument('--start_gene', type=int, required=True, help='Start gene ID for the chunk')
parser.add_argument('--end_gene', type=int, required=True, help='End gene ID for the chunk')
parser.add_argument('--model_name', type=str, required=True, help='Name of the Stan model')
args = parser.parse_args()

# Extract command-line arguments
start_gene = args.start_gene
end_gene = args.end_gene
name = args.model_name

# Generate a unique identifier for outputs
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_file = f"{name}_{start_gene}-{end_gene}_{timestamp}_enhancer_betas.csv"

# Initialize the output file
with open(output_file, 'w') as f:
    f.write('geneID,enhancerID,predicted_beta\n')

# Load the model
model = CmdStanModel(stan_file=f'{name}/{name}.stan')

# Load the dataset
final_data = pd.read_csv("/hpc/group/cbb914f24/data1.txt", sep='\t')

# Run the model for each gene in the specified range
for gene in range(start_gene, end_gene + 1):
    gene_data = final_data[final_data.geneID == gene]
    if gene_data.empty:
        print(f"Skipping gene {gene}: no data found.")
        continue
    # Prepare data for Stan model
    data = {
        'N_genes': 1,
        'N_control_samples': 30,
        'N_experiment_samples': 750,
        'N_enhancers': gene_data.enhancerID.nunique() - 1,
        'y_control': np.stack(gene_data[gene_data.enhancerID == -1].groupby('geneID').agg({'expression': list}).expression.values, axis=0),
        'y_experiment': np.stack(gene_data[gene_data.enhancerID != -1].groupby('geneID').agg({'expression': list}).expression.values, axis=0),
        'enhancer': np.stack(
            gene_data[gene_data.enhancerID != -1]
            .groupby('geneID')
            .agg({'enhancerID': lambda x: (x - x.min() + 1).tolist()}).enhancerID.values, axis=0
        ),
        'N_guides': gene_data.guideID.nunique() - 1,
        'guide': np.stack(
            gene_data[gene_data.enhancerID != -1]
            .groupby('geneID')
            .agg({'guideID': lambda x: (x - x.min() + 1).tolist()}).guideID.values, axis=0),
    }

    # Fit the model
    print(f"Running model for gene {gene}...")
    model_fit = model.sample(data=data, thin=1, iter_sampling=1500, iter_warmup=50, seed=11051999, max_treedepth=100)

    # Extract predicted betas for enhancers
    results = model_fit.summary()
    predicted_betas = results.loc[results.index.str.contains('b'), 'Mean'].sort_index().values[:data['N_enhancers'] + 1]

    # Write results for each enhancer to the file
    with open(output_file, 'a') as f:
        for enhancer_idx, predicted_beta in enumerate(predicted_betas, start=0):  # Include enhancerID starting from 0
            f.write(f"{gene},{enhancer_idx},{predicted_beta}\n")
    print(f"Results for gene {gene} written to {output_file}")

print(f"Processing complete. All results saved to {output_file}")
