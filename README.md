# Multimodal Analysis of Breast Cancer: Transcriptomics, Histopathology and Their Integration

Thesis repository for the integrated study of TCGA-BRCA bulk RNA-seq profiles and whole-slide histopathology images, fused with MOFA+, including external validation on the INC Bogotá cohort.

> **Author:** Daniela Laguna-Santofimio· **Supervisor(s):** Liliana Lopez-Kleine, Fabio A. González· **Institution:** Universidad Nacional de Colombia· **Year:** 2026

## Abstract

Breast cancer (BC) is characterized by high biological heterogeneity, resulting from genetic and environmental interactions that remodel the tumor microenvironment and complicate its clinical characterization. BC is the most common neoplasm and one of the leading causes of death among women worldwide. In 2022, 2.3 million cases were diagnosed and more than 665,000 deaths were recorded; in Colombia, these figures reached 17,018 cases and 4,752 deaths, with projections of sustained growth through 2045. This scenario poses a critical challenge for healthcare systems, exacerbated by the limited availability of specialists and diagnostic delays. Given the limitations of conventional diagnostic methods, this project proposes the development of machine learning models capable of integrating histopathological images and transcriptomic data to classify patients according to their subtype. BC has been classified using various histological, immunohistochemical, and molecular criteria; however, in this study, we focus on the classification of intrinsic molecular subtypes. Using machine learning approaches, this project aims to compare unimodal and multimodal models that may help elucidate key clinical variants relevant to classification. This work seeks to evaluate the effectiveness of integrating omics data and diagnostic images for breast cancer classification, thereby advancing toward more timely, accurate, and equitable diagnoses and contributing to the development of personalized medicine based on integrative methodologies supported by artificial intelligence.

## Study design

The project asks whether the molecular subtype of a breast tumour (Luminal A, Luminal B, Basal-like) can be predicted from tissue morphology, from gene expression, and from the two combined, on a paired TCGA-BRCA cohort of **552 patients** with both a whole-slide image and an RNA-seq profile.

It is organised in three phases:

| Phase | Module | Question |
|---|---|---|
| 1 | `01_transcriptomica` | How well does gene expression alone recover the subtype? |
| 2 | `02_histopatologia` | How much subtype signal is there in the slide image, and which patch representation extracts the most? |
| 3 | `03_multimodal` | Do the two modalities share structure, and does fusing them with MOFA+ improve on either alone? |

A single **patient-level train/test partition (386 / 166)** is drawn once in Phase 1 (`03_clasificadores_firma.ipynb`) and reused verbatim by Phases 2 and 3. This is what makes the modalities comparable and the fusion legitimate: the same patients are in the training set everywhere, so no image model is ever tested on a patient whose expression profile trained the transcriptomic model. Every phase also uses the same evaluation protocol — hyperparameter tuning by 5-fold stratified CV, a held-out test set, and a final 5 × 3 repeated stratified CV with a paired Wilcoxon test between the two best models.

## Repository structure

```
tesis-brca-multimodal/
├── README.md
├── .gitignore
├── envs/                                   # Python environments, one per modality
│   ├── requirements_rna.txt
│   ├── requirements_wsi.txt
│   └── requirements_mofa.txt
├── 00_datos/                               # Data documentation (raw data are not tracked)
│   ├── descarga_gdc.md
│   ├── manifest_tcga_brca.txt
│   └── coincidencias_pacientes_inc.txt     # RNA sample <-> slide_ID map for the INC cohort
├── 01_transcriptomica/                     # Phase 1 - RNA-seq
│   ├── 01_qc_preprocesamiento.ipynb
│   ├── 02_expresion_diferencial_limma.(R|ipynb)
│   ├── 03_clasificadores_firma.ipynb       # consensus DEG signature (604 genes); ORIGIN OF THE SPLIT
│   └── 04_transcriptoma_filtrado_y_genes_alta_confianza.ipynb
├── 02_histopatologia/                      # Phase 2 - Whole-slide images (see its README)
│   ├── 00_preprocesamiento_propio/         # in-house preprocessing (TCGA: 17 slides; INC: 41 slides)
│   │   ├── wsi_step1_segmentation.py
│   │   └── wsi_step2_patch_extraction.py
│   ├── E1_kimianet_mosaico/                # approach 1: mosaics of 16 k-medoids patches
│   ├── E2_E3_crops/                        # approaches 2-3: 224x224 crops (KimiaNet / UNI)
│   ├── E4_kimianet_kmedoids/               # approach 4: k-medoids on KimiaNet embeddings
│   ├── E5_uni2h_kmedoids/                  # approach 5: UNI2-h, independent tiling
│   │   ├── wsi_processing/                 # preprocessing package (see its README)
│   │   └── run_wsi_pipeline.py             # single CLI driver: slides in, HDF5 + QC out
│   ├── 03_clasificadores_wsi.ipynb         # the five approaches under a common protocol
│   ├── 04_mlp.ipynb
│   └── 05_validacion_externa_inc.ipynb
├── 03_multimodal/                          # Phase 3 - MOFA+ integration
│   ├── 01_similitud_y_mantel.ipynb         # patient-patient similarity + Mantel test
│   ├── 02_mofa_y_clasificacion.ipynb       # MOFA+ training, classification, INC validation
│   └── 03_gsea_factores.ipynb              # biological interpretation of the factors
└── resultados/                             # Tables, enrichment results and figures
    ├── tablas/
    ├── david/
    ├── gsea/
    └── figuras/
```

Directories that currently hold only a `.gitkeep` file are placeholders for the corresponding scripts and notebooks.

## Phase 3 — Multimodal integration with MOFA+

### Rationale

Before fusing anything, the two modalities are tested for shared structure. Patient–patient similarity matrices are built for each one (Pearson correlation over the 604-gene signature for RNA; cosine similarity over the L2-normalised embeddings for WSI) and compared with a **Mantel test** with 9,999 permutations. Only if the two geometries agree above chance does a joint latent space make sense.

**MOFA+** (Multi-Omics Factor Analysis) is then fitted on both views at once. It decomposes the two matrices into a small set of latent factors shared by all patients, and reports, for every factor, how much variance it explains **in each view separately**. That decomposition is what distinguishes a factor driven by gene expression, one driven by morphology, and one genuinely shared by both — which is the question the whole thesis turns on.

### Pipeline

| Step | Content |
|---|---|
| 1 | Loading and preprocessing. RNA (already log2-transformed) centred per gene; WSI embeddings L2-normalised per patient, then centred per feature. Gaussian likelihood for both views. |
| 2 | Patient–patient similarity matrices, intra- vs inter-subtype gap with Mann-Whitney U, heatmaps ordered by subtype. |
| 3 | Mantel test, global and per subtype. |
| 3b | Single-view MOFA models (RNA-only, WSI-only) as the unimodal baseline. |
| 4 | MOFA+ training (K = 15 initial factors, ARD prior prunes them, spike-and-slab sparsity on the weights) and extraction of factor scores `Z`, weights `W` and R² per factor and view. |
| 5 | Association between factors and PAM50 subtype: one-way ANOVA per factor with Benjamini–Hochberg FDR correction. |
| 6 | GSEA on the RNA weights (MSigDB Hallmark, KEGG, GO-BP) for biological interpretation. |
| 7 | Supervised classification on the factor scores, under the same protocol as Phases 1 and 2. |
| 8 | Visualisations: variance heatmap, factor scatter and violins, UMAP, orthogonality, top genes per factor. |
| 9 | External validation: the INC patients are projected onto the trained factor space and classified. |
| 10 | Unimodal vs multimodal comparison. |

The notebook is parameterised by the WSI extractor: `WSI_EXTRACTOR_NAME`, `WSI_FEAT_DIM` and `WSI_PATH` select which Phase 2 representation feeds the fusion, so the same notebook can be re-run for any of the five approaches. The results below correspond to **UNI-Parches-Crop** (Approach 3, 1024-D).

The MOFA model is cached as `mofa_model.hdf5` and reused if it exists; delete it to retrain.

## Reproducing the analysis

1. **Environments.** Create one virtual environment per modality from `envs/`. They are not interchangeable: the RNA environment requires **scikit-learn < 1.8** (the logistic-regression grid uses `multi_class`, removed in 1.8), while the WSI environment requires **NumPy >= 2** (OpenCV 4.13). The MOFA environment adds `mofapy2` (0.7.4), `gseapy` (1.2.1), `statsmodels` and `h5py` on top of the RNA one.
2. **Data.** Follow `00_datos/descarga_gdc.md` and place the input files in `00_datos/`. Data files are excluded from version control (see `.gitignore`).
3. **Execution order.** `01_transcriptomica` → `02_histopatologia` → `03_multimodal`. The order is a hard dependency, not a convention: Phase 1 writes the train/test partition, and Phase 3 consumes both the 604-gene signature and the WSI feature table produced in Phase 2.
4. **Outputs.** Notebooks write to `resultados/`, with a module-specific file prefix to prevent collisions.

## Reference results

Obtained with seed 42. Minor numerical differences may occur across platforms and library versions.

### Phase 1 — Transcriptomics (604-gene consensus signature, test set n = 166)

| Model | Test F1-macro |
|---|---|
| Logistic regression | 0.94 |
| SVM | 0.92 |
| XGBoost | 0.91 |
| Random forest | 0.90 |
| KNN | 0.88 |

### Phase 2 — Histopathology (repeated CV, F1-macro)

| Approach | Best model | F1-macro |
|---|---|---|
| E3 — UNI crops | Logistic regression | 0.741 ± 0.034 |
| E1 — KimiaNet mosaics | SVM | 0.654 ± 0.048 |
| E2 — KimiaNet crops | Logistic regression | 0.460 ± 0.033 |

Morphology alone therefore carries real but much weaker subtype signal than expression, and the feature extractor matters more than the classifier.

### Phase 3 — Shared structure

Mantel test, RNA vs WSI distance matrices (n = 552, 152,076 pairs): **r = 0.112, p = 1.0 × 10⁻⁴** (9,999 permutations). The association is significant but explains only about 1.2 % of the variance. Per subtype: Basal r = 0.178 (p = 0.0016), LumA r = 0.068 (p = 0.028), LumB r = −0.013 (n.s.).

The two modalities are therefore weakly but genuinely coupled, and the coupling is carried almost entirely by the Basal-like tumours.

### Phase 3 — MOFA+ factor space (UNI-Parches, 14 factors after ARD pruning)

Total variance explained: **RNA 43.3 %, WSI 46.3 %**. The decomposition per factor shows near-complete separation rather than fusion: F1 is RNA-driven (R² 32.9 % RNA vs 2.3 % WSI), while F2 and F3 are WSI-driven (7.8 % and 6.8 % WSI vs about 0.02 % RNA). Factors significantly associated with subtype after FDR correction include F1 (p ≈ 3 × 10⁻³¹⁰), F4 and F2.

### Phase 3 — Classification on the factor scores (test set n = 166)

| Model | CV F1-macro | Test F1-macro | Test accuracy |
|---|---|---|---|
| Logistic regression | 0.883 ± 0.026 | **0.907** | 0.904 |
| SVM | 0.878 ± 0.027 | 0.895 | 0.892 |
| XGBoost | 0.867 ± 0.019 | 0.854 | 0.849 |
| Random forest | 0.881 ± 0.026 | 0.849 | 0.843 |
| KNN | 0.830 ± 0.033 | 0.784 | 0.777 |

Paired Wilcoxon, logistic regression vs random forest (per-fold F1-macro): p = 0.975. K-means (k = 3) on the factor scores: ARI = 0.576.

### Phase 3 — Unimodal vs multimodal

| Latent space | Best model | Test F1-macro |
|---|---|---|
| MOFA multimodal | Logistic regression | 0.907 |
| MOFA RNA-only | Random forest | 0.896 |
| MOFA WSI-only | Logistic regression | 0.446 |

The fusion adds about one point of F1 over expression alone. Given the width of the cross-validation intervals (± 0.026), that margin is not established as a real improvement; the honest reading is that **morphology does not degrade the transcriptomic model and contributes little beyond it** in this cohort.

### Phase 3 — External validation (INC Bogotá, n = 23, all Luminal B)

The INC patients are projected onto the trained factor space and classified. Sensitivity for Luminal B:

| Model | Correct | Sensitivity | Predicted distribution |
|---|---|---|---|
| SVM | 23 / 23 | 1.000 | all LumB |
| KNN | 21 / 23 | 0.913 | 21 LumB, 2 Basal |
| Random forest | 15 / 23 | 0.652 | 15 LumB, 5 LumA, 3 Basal |
| XGBoost | 7 / 23 | 0.304 | 13 LumA, 7 LumB, 3 Basal |
| Logistic regression | 0 / 23 | 0.000 | 22 LumA, 1 Basal |

The spread from 0.000 to 1.000 across models that perform almost identically on TCGA is the central caveat of this section (see below).

## Limitations

1. **The external validation does not rank the models.** On the INC cohort the best TCGA model (logistic regression) assigns none of the 23 cases correctly, while SVM assigns all of them. With a single class present, a model biased towards Luminal B scores 1.000 without having learnt anything transferable. These numbers measure robustness to domain shift, not clinical accuracy, and no model should be selected on them.
2. **Nine signature genes are missing from the INC panel** and are imputed with a constant, including genes with large MOFA weights on F1 (LY6D, C6orf15, FGFBP1). Since F1 is the factor that carries the subtype signal, the projection of the INC patients onto the factor space is degraded in exactly the dimension that matters.
3. **Different normalisations.** The INC expression data are VST-normalised and TCGA is log2(FPKM + 1); the INC slides come from a different centre, scanner and staining protocol. No batch correction is applied anywhere in the project.
4. **MOFA does not find shared factors.** The variance decomposition splits cleanly into RNA-driven and WSI-driven factors, with no factor explaining substantial variance in both views. This is consistent with the weak Mantel correlation and is itself a result: at this sample size and with these representations, the two modalities do not share recoverable latent structure.
5. **MOFA-WSI-only reports R² = 0.00 %** while still yielding factor scores that classify at F1 = 0.446. This inconsistency is unresolved, and the WSI-only baseline should not be reported until it is explained.
6. **Class encoding differs between phases.** Phases 1 and 2 use an explicit dictionary (LumA = 0, LumB = 1, Basal = 2); Phase 3 uses `LabelEncoder`, which sorts alphabetically (Basal = 0, LumA = 1, LumB = 2). Each notebook is internally consistent, but any comparison of saved predictions across phases must decode the labels before aligning them.
7. **Optimistic cross-validation.** In every phase, the repeated CV is run on the full cohort with hyperparameters tuned on the training partition, and the folds share training data. Nested cross-validation would be needed for an unbiased estimate, and the Wilcoxon *p*-values are approximate.
8. **Heterogeneous slide preprocessing.** 17 TCGA slides and all INC slides were preprocessed in-house, while the remainder came from an earlier pipeline. A batch effect between the two groups has not been ruled out.

## Third-party components and credits

* **KimiaNet** (Riasatian et al., *Medical Image Analysis*, 2021) and **UNI / UNI2-h** (Mahmood Lab) are used for feature extraction; their weights are not redistributed and are subject to their own licences (UNI2-h is gated on the Hugging Face Hub).
* **MOFA+** (Argelaguet et al., *Genome Biology*, 2020), via `mofapy2`.
* Parts of the histopathology pipeline were adapted from the work of [COLLABORATOR — TO BE COMPLETED]; obtain their consent and state the applicable licence before publication.

## Data availability

The publicly available cohort corresponds to The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) project, an initiative that links clinical phenotypes with molecular genotypes in invasive breast cancer (Cancer Genome Atlas Network, 2012, Nature). Clinical, genomic, pathological, and histopathological imaging data are available on the Genomic Data Commons portal (Grossman et al. 2016) of the National Cancer Institute (https://portal.gdc.cancer.gov/projects/TCGA-BRCA). To construct the working cohort, successive filters were applied to the portal: first, cases from the TCGA-BRCA project were selected; then, the files were filtered by data type, retaining only diagnostic slide images and gene expression data by RNA-seq; finally, the selection was restricted to open access files. INC data requires special permission for use.
