# Codes for Manuscript: Generative AI-based design of hybrid transcriptional activator proteins with new DNA-binding specificity

<img width="1389" height="536" alt="image" src="https://github.com/user-attachments/assets/c94347f5-c0c4-4b65-b434-45a2e8949108" />

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/kigalab/VAE-hybridTF/blob/main/LICENSE)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2025.12.31.697150-orange)](https://www.biorxiv.org/content/10.64898/2025.12.31.697150)

This repository contains the computational codes used in the study **“Generative AI-based design of hybrid transcriptional activator proteins with new DNA-binding specificity.”**

The project explores multiple computational approaches to design hybrid transcription factors (TFs) by navigating the latent space of **variational autoencoders (VAEs)**. The aim is to generate novel proteins that exhibit functional properties intermediate between two parent transcription factor subfamilies.

In our [preprint](https://www.biorxiv.org/content/10.64898/2025.12.31.697150), we demonstrate that VAE-designed hybrid TFs derived from **LuxR** and **LasR** can exhibit intermediate DNA-binding specificities and transcriptional activities.

---

## Repository Structure

The repository is organized into several directories corresponding to different stages of the project.

### 1. Promoter_Library_Analysis and Protein_Library_Analysis

Scripts for analyzing **sort-seq datasets** generated from promoter and protein variant libraries.

### 2. VAE

Codes for training VAE models on transcription factor sequences and generating hybrid transcription factor variants by exploring the latent space between protein subfamilies.

### 3. analysis

Scripts used to generate figures and perform additional analyses for the manuscript.

---

## Citation

If you use the code in this repository, please cite:

```

Okuda S, Minami A, Aiko M, et al.
Generative AI-based design of hybrid transcriptional activator proteins with new DNA-binding specificity.
bioRxiv (2025)
https://www.biorxiv.org/content/10.64898/2025.12.31.697150

```
---

## License

MIT License
