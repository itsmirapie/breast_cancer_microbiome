#!/usr/bin/env python3
"""
Create plots from extracted QIIME2 data
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def create_taxonomic_plot():
    """Create plot from taxonomic data"""
    try:
        # Look for taxonomic CSV files
        tax_files = list(Path("results/extracted_figures").rglob("level-*.csv"))
        if not tax_files:
            print("No taxonomic CSV files found")
            return
            
        # Use genus level (level-6.csv)
        genus_file = Path("results/extracted_figures/figure1_taxa_barplot/data/level-6.csv")
        if genus_file.exists():
            df = pd.read_csv(genus_file)
            print(f"Loaded genus data: {genus_file}")
            print(f"Shape: {df.shape}")
            
            # Create output directory
            Path("results/plots_from_data").mkdir(exist_ok=True)
            
            # Plot top 15 genera
            plt.figure(figsize=(12, 8))
            top_genera = df.mean().sort_values(ascending=False).head(15)
            top_genera.plot(kind='bar')
            plt.title('Top 15 Bacterial Genera (Average Relative Abundance)', fontsize=16)
            plt.xlabel('Genus', fontsize=12)
            plt.ylabel('Relative Abundance', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig('results/plots_from_data/top_genera.png', dpi=300)
            print("Saved: results/plots_from_data/top_genera.png")
            
            # Also save as CSV
            top_genera.to_csv('results/plots_from_data/top_genera_summary.csv')
            print("Saved: results/plots_from_data/top_genera_summary.csv")
            
    except Exception as e:
        print(f"Error creating taxonomic plot: {e}")

def check_statistical_results():
    """Check for statistical results"""
    print("\n=== Checking for Statistical Results ===")
    
    # Look for statistical files
    stats_patterns = ["*pairwise*", "*significance*", "*permanova*"]
    for pattern in stats_patterns:
        files = list(Path("results/extracted_figures").rglob(f"*{pattern}*.csv"))
        for file in files:
            print(f"\nFound: {file}")
            try:
                df = pd.read_csv(file)
                print(f"  Shape: {df.shape}")
                print(f"  Columns: {list(df.columns)}")
                if not df.empty:
                    print("  First few rows:")
                    print(df.head())
            except:
                pass

def main():
    print("=== Creating Plots from Extracted Data ===")
    
    # Create plots directory
    Path("results/plots_from_data").mkdir(exist_ok=True)
    
    # Create taxonomic plot
    create_taxonomic_plot()
    
    # Check statistical results
    check_statistical_results()
    
    print("\n=== DONE ===")
    print("Plots saved to: results/plots_from_data/")
    print("Use these for your report and presentation.")

if __name__ == "__main__":
    main()
