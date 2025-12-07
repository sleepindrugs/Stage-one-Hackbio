def translate_dna_to_protein(dna_sequence):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }
    dna_sequence = dna_sequence.upper().replace(" ", "").replace("\n", "")
    protein = ""


    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        amino_acid = codon_table.get(codon, 'X')
        protein += amino_acid


    return protein


# Example usage
dna = "ATCGCCATTGTAGATGTGGCCGCTGAAAGGGTGCCCGATAG"
protein = translate_dna_to_protein(dna)


print("Protein sequence:", protein)



# This function calculates the Hamming distance between two names
def hamming_distance(name1, name2):
    # Step 1: Make both names the same length
    # If one name is shorter, add spaces at the end
    if len(name1) < len(name2):
        name1 = name1 + " " * (len(name2) - len(name1))
    elif len(name2) < len(name1):
        name2 = name2 + " " * (len(name1) - len(name2))


    # Start counting differences
    distance = 0


    # Compare each letter in the same position
    for i in range(len(name1)):
        if name1[i] != name2[i]:
            distance = distance + 1


    # Step 4: Return the number of differences
    return distance




# Example:
slack_name = "Gloria Akpederi"     # Your Slack username
twitter_name = "GloriaAkpederi"    # Your Twitter/X handle


# Use the function and print the result
result = hamming_distance(slack_name, twitter_name)
print("The Hamming distance is:", result)




!pip install pandas
!pip install numpy
!pip install matplotlib
!pip install seaborn
!pip install sklearn

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn

# import file link from HackBio
url_heatmap = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"
url_volcano = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"
url_scatter_densityplot ="https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"

# load datasets
heatmap_df = pd.read_csv(url_heatmap)
volcano_df = pd.read_csv(url_volcano)
scatter_density_df = pd.read_csv(url_scatter_densityplot)

# Select top 20 genes (example)
top_genes = heatmap_df.head(20)

# Set the 'Unnamed: 0' column as the index
top_genes = top_genes.set_index('Unnamed: 0')

# Plot the clustered heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(top_genes, cmap="Blues", annot=False)
plt.title("Clustered Heatmap of Top Differentially Expressed Genes (HBR vs UHR)")
plt.xlabel("Samples")
plt.ylabel("Genes")
plt.show()
print("Heatmap dataset shape:", heatmap_df.shape)
print("Number of genes:", heatmap_df.shape[0])
print("Number of samples:", heatmap_df.shape[1])

# Ensure correct columns
# Columns must include: 'log2FoldChange' and 'padj'
volcano_df = volcano_df.dropna(subset=['log2FoldChange', 'PAdj'])

# Create significance column
volcano_df['significance'] = 'Not significant'
volcano_df.loc[(volcano_df['log2FoldChange'] > 1) & (volcano_df['PAdj'] < 0.05), 'significance'] = 'Upregulated'
volcano_df.loc[(volcano_df['log2FoldChange'] < -1) & (volcano_df['PAdj'] < 0.05), 'significance'] = 'Downregulated'

# Plot
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=volcano_df,
    x='log2FoldChange',
    y=-np.log10(volcano_df['PAdj']),
    hue='significance',
    palette={'Upregulated': 'green', 'Downregulated': 'orange', 'Not significant': 'grey'},
    alpha=0.7
)

# Add vertical threshold lines
plt.axvline(x=1, color='black', linestyle='--')
plt.axvline(x=-1, color='black', linestyle='--')

plt.title("Volcano Plot of Differentially Expressed Genes")
plt.xlabel("log2(Fold Change)")
plt.ylabel("-log10(Adjusted p-value)")
plt.legend(title="Significance")
plt.show()
print("Number of Upregulated genes:", sum(volcano_df['significance'] == 'Upregulated'))
print("Number of Downregulated genes:", sum(volcano_df['significance'] == 'Downregulated'))
print("Number of Not significant genes:", sum(volcano_df['significance'] == 'Not significant'))

plt.figure(figsize=(8,6))
sns.scatterplot(
    data=scatter_density_df,
    x='radius_mean',
    y='texture_mean',
    hue='diagnosis',
    palette={'M': 'red', 'B': 'blue'}
)
plt.title('Scatter Plot: Radius vs Texture')
plt.xlabel('Radius Mean')
plt.ylabel('Texture Mean')
plt.grid(True)
plt.show()
print("=== Scatter Plot (Radius vs Texture) ===")
print("This plot shows how texture_mean varies with radius_mean, colored by diagnosis.\n")

features = [
    'radius_mean', 'texture_mean', 'perimeter_mean',
    'area_mean', 'smoothness_mean', 'compactness_mean'
]
corr = scatter_density_df[features].corr()

plt.figure(figsize=(8,6))
sns.heatmap(corr, annot=True, cmap='coolwarm')
plt.title('Correlation Heatmap of Key Features')
plt.show()
print("=== Correlation Heatmap ===")
print("This heatmap shows correlations between six key features.\n")

plt.figure(figsize=(8,6))
sns.scatterplot(
    data=scatter_density_df,
    x='smoothness_mean',
    y='compactness_mean',
    hue='diagnosis',
    palette={'M': 'red', 'B': 'blue'}
)
plt.title('Scatter Plot: Smoothness vs Compactness')
plt.xlabel('Smoothness Mean')
plt.ylabel('Compactness Mean')
plt.grid(True)
plt.show()
print("=== Scatter Plot (Smoothness vs Compactness) ===")
print("This plot compares smoothness_mean and compactness_mean for each diagnosis.\n")

plt.figure(figsize=(8,6))
sns.kdeplot(
    data=scatter_density_df[scatter_density_df['diagnosis'] == 'M'],
    x='area_mean',
    label='Malignant',
    fill=True,
    alpha=0.5
)
sns.kdeplot(
    data=scatter_density_df[scatter_density_df['diagnosis'] == 'B'],
    x='area_mean',
    label='Benign',
    fill=True,
    alpha=0.5
)
plt.title('Density Plot: Area Mean Distribution')
plt.xlabel('Area Mean')
plt.ylabel('Density')
plt.legend()
plt.show()
print("=== Density Plot (Area Mean Distribution) ===")
print("This plot shows the area_mean distribution for malignant and benign samples.\n")
