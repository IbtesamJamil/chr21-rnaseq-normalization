import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(
        description="RNA-Seq Count Normalization Pipeline (RPKM & TPM)"
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        default="chromosome21_raw_counts.csv",
        help="Path to input raw counts CSV file (default: chromosome21_raw_counts.csv)"
    )
    parser.add_argument(
        "-o", "--output-csv",
        type=str,
        default="chromosome21_normalized_counts.csv",
        help="Path to output normalized CSV file (default: chromosome21_normalized_counts.csv)"
    )
    parser.add_argument(
        "-p", "--plot",
        type=str,
        default="chr21_tpm_expression.png",
        help="Path to output plot figure (default: chr21_tpm_expression.png)"
    )
    return parser.parse_args()

def process_rnaseq_counts(input_path, output_csv_path, plot_path):
    print(f"Loading raw counts from: {input_path}")
    df = pd.read_csv(input_path)

    # Validate required columns
    required_cols = {'Gene_ID', 'Gene_Symbol', 'Length_bp', 'Raw_Count'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Input CSV must contain columns: {required_cols}")

    # Calculate RPK (Reads Per Kilobase)
    df['Length_kb'] = df['Length_bp'] / 1000.0
    df['RPK'] = df['Raw_Count'] / df['Length_kb']

    # Calculate RPKM (Reads Per Kilobase Million)
    total_reads_millions = df['Raw_Count'].sum() / 1e6
    df['RPKM'] = df['RPK'] / total_reads_millions

    # Calculate TPM (Transcripts Per Million)
    rpk_sum_millions = df['RPK'].sum() / 1e6
    df['TPM'] = df['RPK'] / rpk_sum_millions

    # Save normalized results
    df.to_csv(output_csv_path, index=False)
    print(f"Normalized data saved to: {output_csv_path}")

    # Generate TPM plot
    sns.set_theme(style='whitegrid')
    plt.figure(figsize=(10, 6))
    
    sorted_df = df.sort_values('TPM', ascending=False)
    sns.barplot(
        data=sorted_df,
        x='Gene_Symbol',
        y='TPM',
        hue='Gene_Symbol',
        legend=False,
        palette='viridis'
    )
    
    plt.title('Chromosome 21 Gene Expression (TPM)', fontsize=14, fontweight='bold')
    plt.xlabel('Gene Symbol', fontsize=12)
    plt.ylabel('Normalized Expression (TPM)', fontsize=12)
    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Expression plot saved to: {plot_path}")

if __name__ == "__main__":
    args = parse_args()
    process_rnaseq_counts(args.input, args.output_csv, args.plot)
