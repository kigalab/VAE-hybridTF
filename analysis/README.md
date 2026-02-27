# Promoter Engineering & VAE Model Visualization Scripts

This repository contains Python scripts used in two stages of the study:

1. **VAE model development and evaluation**
   - Grid search visualization
   - Latent-space t-SNE embedding plots

2. **Experimental promoter library analysis**
   - 4D activity scatter plots
   - Plux vs Plas scatter plots
   - Library screening visualization
   - Randomized promoter analysis (violin plots + PLS coefficient heatmaps)

- 4D scatter plots
- Grid search performance plots
- t-SNE embedding visualization
- Plux vs Plas scatter plots (mean ± SEM)
- Library screening scatter plots
- Randomized promoter analysis:
  - Violin plots + ANOVA + Tukey HSD
  - PLS coefficient heatmaps (2×2 layout)

## File Overview

### VAE Model Development

| Script | Description |
|--------|------------|
| `plot_grid_search_results.py` | Grid search performance visualization during VAE hyperparameter optimization (SuppFig.S7) |
| `plot_tsne.py` | Latent space embedding visualization of sequences (Fig.1c) |

---

### Promoter Library Analysis

| Script | Description |
|--------|------------|
| `scatter_plux_plas.py` | Plux vs Plas (mean ± SEM) (Fig.2b,3c) |
| `scatter_library.py` | Library screening scatter plots (Fig.3b) |
| `randomlib.py` | Randomized promoter analysis (Fig.4): (b) violin + (c) PLS heatmaps |
| `plot_4d_scatter.py` | 4D activity scatter plots (LasR vs LuxR) (Fig.4d) |

## Reproducibility

Tested environment:

- Python 3.12.7
- numpy 1.26.4
- pandas 2.2.2
- matplotlib 3.9.2
- seaborn 0.13.2
- scipy 1.13.1
- statsmodels 0.14.2
- scikit-learn 1.5.1
