# E. coli Operon Prediction using Machine Learning

A computational biology project for **COMPSCI 602** that predicts operon membership in *Escherichia coli* K-12 using traditional machine learning models and interpretable genomic features.

## 📋 Project Overview

This project tackles the **operon prediction task** from the [DGEB benchmark](https://github.com/TattaBio/DGEB) (Diverse Genomic Embedding Benchmark). Given pairs of adjacent genes in *E. coli*, the goal is to predict whether they belong to the same operon (co-transcribed) or not.

### Key Contributions
- **77 interpretable genomic features** engineered from sequence data, genomic architecture, and functional annotations
- **6 ML baseline models** achieving competitive performance against protein language models (ESM-2)
- **Comprehensive ablation studies** quantifying individual feature importance
- **PCA-based feature analysis** of trinucleotide composition patterns

## 🗂️ Repository Structure

```
├── data/
│   ├── ecoli_pairs.csv                    # Processed feature matrix (4,297 gene pairs)
│   ├── sequence.gb                        # E. coli K-12 GenBank file
│   ├── protein.fasta                      # Extracted protein sequences
│   └── Galaxy3-[eggNOG...].tabular        # COG functional annotations
├── models/
│   ├── *.joblib                           # Trained model files (6 models)
│   └── feature_metadata.joblib            # Feature configuration
├── graphs/
│   ├── distance_feature_ablation/         # Ablation results for distance features
│   ├── gc_content_feature_ablation/       # Ablation results for GC content features
│   ├── strand_feature_ablation/           # Ablation results for strand/orientation features
│   ├── cog_feature_ablation/              # Ablation results for COG features
│   └── trinuc_feature_ablation/           # Ablation results for trinucleotide features
├── Reports/                               # Project reports (PDF)
├── data_preproc.ipynb                     # Data preprocessing & feature engineering
├── simple_ML_baselines.ipynb              # Model training & evaluation
├── ablation.ipynb                         # Feature importance ablation studies
└── ESM baseline results.ipynb             # ESM-2 protein language model comparison
```

## 🧬 Features

### Feature Groups (77 total)

| Category | Features | Description |
|----------|----------|-------------|
| **GC Content** (3) | `gc_content_A`, `gc_content_B`, `gc_content_diff` | GC ratio of each gene and their difference |
| **Genomic Architecture** (3) | `intergenic_distance`, `genes_overlap`, `overlap_length` | Distance/overlap between adjacent genes |
| **Strand/Orientation** (5) | `strand_concordant`, `orientation_++/--/+-/-+` | Relative orientation of gene pairs |
| **COG Annotations** (2) | `COG_match`, `COG_similar` | Functional category overlap from eggNOG |
| **Trinucleotides** (64) | `trinuc_AAA`, `trinuc_AAC`, ... | Intergenic region trinucleotide frequencies |

## 🤖 Models

Six classifiers trained with hyperparameter tuning:

| Model | Description |
|-------|-------------|
| Logistic Regression | L2-regularized linear classifier |
| Random Forest | Ensemble of 200 decision trees |
| XGBoost | Gradient boosted trees |
| SVM | RBF kernel support vector machine |
| MLP | 2-layer neural network (100, 50 units) |
| Naive Bayes | Gaussian Naive Bayes |

## 📊 Data Sources

- **Genome Sequence**: [NCBI GenBank U00096.3](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3) - *E. coli* str. K-12 substr. MG1655 complete genome
- **Operon Labels**: [DGEB EcoliOperon task](https://github.com/TattaBio/DGEB) via Hugging Face datasets
- **COG Annotations**: [eggNOG-mapper](http://eggnog-mapper.embl.de/) via Galaxy Europe

## 🔬 Ablation Studies

The `ablation.ipynb` notebook systematically evaluates feature importance by:

1. **Ablating features** using shuffle (permutation) and random sampling methods
2. **Measuring performance drop** across all 6 models
3. **Generating heatmaps** comparing ablation methods
4. **PCA analysis** of trinucleotide features (PC1, PC1-5, PC1-10)

Key findings are saved to `graphs/*/ablation_results_*.csv` and visualized as comparison heatmaps.

## 🚀 Getting Started

### Prerequisites

```bash
# Create conda environment
conda create -n operon python=3.11
conda activate operon

# Install dependencies
pip install numpy pandas scikit-learn xgboost biopython matplotlib seaborn joblib
pip install transformers datasets  # For ESM comparison
pip install git+https://github.com/TattaBio/DGEB.git
```

### Running the Pipeline

1. **Data Preprocessing**: Run `data_preproc.ipynb` to generate `data/ecoli_pairs.csv`
2. **Model Training**: Run `simple_ML_baselines.ipynb` to train and save models
3. **Ablation Studies**: Run `ablation.ipynb` to analyze feature importance

## 📈 Results

The ML baselines achieve competitive performance on the operon prediction task, with detailed results in the project reports under `Reports/`.

## 📚 References

- Janga, S. C., et al. (2006). "Transcriptional regulation constrains the organization of genes on eukaryotic chromosomes." *PNAS*. [PMC1557821](https://pmc.ncbi.nlm.nih.gov/articles/PMC1557821/)
- West-Roberts, T., Kundaje, A., & Zou, J. (2024). "DGEB: A Diverse Genomic Embedding Benchmark for functional genomics tasks." *bioRxiv*. [doi:10.1101/2024.01.10.575096](https://www.biorxiv.org/content/10.1101/2024.01.10.575096v1)
- Krishnakumar, R. & Ruffing, A. M. (2022). "OperonSEQer: A set of machine-learning algorithms with threshold voting for detection of operon pairs using short-read RNA-sequencing data." *PLoS Comput Biol*, 18(1), e1009731. [doi:10.1371/journal.pcbi.1009731](https://doi.org/10.1371/journal.pcbi.1009731)
- Karaji, R. & Peña-Castillo, L. (2025). "OpDetect: A convolutional and recurrent neural network classifier for precise and sensitive operon detection from RNA-seq data." *PLoS One*, 20(8), e0329355. [doi:10.1371/journal.pone.0329355](https://doi.org/10.1371/journal.pone.0329355)
- Moreno-Hagelsieb, G. & Collado-Vides, J. (2002). "A comparative genomics approach to predict operons in prokaryotes." *Bioinformatics*, 18(suppl 1), S329–S336.
- Kanhere, A., Janga, S. C., & Moreno-Hagelsieb, G. (2006). "The distinctive signatures of promoter regions and operon junctions across prokaryotes." *Nucleic Acids Research*, 34(14), 3980–3987. [doi:10.1093/nar/gkl563](https://doi.org/10.1093/nar/gkl563)
- Galperin, M. Y., et al. (2021). "COG database update: focus on microbial diversity, model organisms, and widespread pathogens." *Nucleic Acids Research*, 49(D1), D274–D281. [doi:10.1093/nar/gkaa1018](https://doi.org/10.1093/nar/gkaa1018)
- Cantalapiedra, C. P., et al. (2021). "eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale." *Molecular Biology and Evolution*, 38(12), 5825–5829. [doi:10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293)
- Sayers, E. W., et al. (2025). "Database resources of the National Center for Biotechnology Information in 2025." *Nucleic Acids Res*, 53(D1), D20–D29. [doi:10.1093/nar/gkae979](https://doi.org/10.1093/nar/gkae979)

## 👤 Author

COMPSCI 602 Project - University of Massachusetts Amherst
