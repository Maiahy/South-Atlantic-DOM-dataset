# South-Atlantic-DOM-dataset

This repository contains all scripts used to generate the dissolved organic matter (DOM) dataset described in the manuscript:

> **A metabolomics dataset of plankton exometabolomes and dissolved organic compounds in the South Atlantic Ocean**

The dataset spans a wide range of marine environments across the South Atlantic and Southern Oceans, including:

- Amazon River plume
- African river plumes
- Benguela Current system
- South Atlantic subtropical gyre
- Weddell Sea

---

## 🌊 Dataset overview

A total of **367 solid-phase extracted DOM samples** collected at **107 sampling sites** were analyzed using liquid chromatography coupled to high-resolution mass spectrometry (LC–HRMS).

Each sample was:

- analyzed in **triplicate in full-scan MS mode** to assess analytical reproducibility;
- analyzed once in **data-dependent acquisition (DDA) mode** to acquire MS/MS spectra for compound annotation;
- accompanied by regular injections of **pooled quality control (QC) samples** throughout the analytical sequences.

Analyses were performed using a **Vanquish UHPLC system** (Thermo Fisher Scientific, Germany) coupled to an **Orbitrap Exploris 480 mass spectrometer** equipped with a **heated electrospray ionization (HESI) source**.

The dataset was generated from **15 analytical batches**.

---

## ⚙️ Data processing workflow

The workflow implemented in this repository includes:

1. Preprocessing of full-scan and DDA datasets using the R package **XCMS**.
2. Alignment and merging of analytical batches.
3. Integration of full-scan and DDA feature tables using **MetabCombiner**.
4. Feature annotation using **CAMERA** and **GNPS molecular networking**.
5. Filtering based on blank samples, pooled QCs, and replicate reproducibility.
6. Additional quality control analyses and technical validation used in the manuscript.

---

## 📖 Citation

If you use this repository or the associated dataset, please cite:

> Henry, M. *et al.*  
> **A metabolomics dataset of plankton exometabolomes and dissolved organic compounds in the South Atlantic Ocean**.
