# South-Atlantic-DOM-dataset

<p align="center">
  <img src="./2022%2002%2006%20navigation%20sous%20voiles%20devant%20iceberg%20Tasmania_drone_DJI_0032%C2%A9Mae%CC%81va%20Bardy%20-%20Fondation%20Tara%20Ocean.jpg"
       alt="Tara sailing in front of an Antarctic iceberg"
       width="900">
</p>

<p align="center">
  <em>Tara sailing in front of an Antarctic iceberg in the Southern Ocean.<br>
  Photo credit: Maéva Bardy / Fondation Tara Ocean.</em>
</p>

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

- analyzed in **triplicate in full-scan MS mode**;
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
5. Filtering and cleaning based on blank samples, pooled QCs, and replicate reproducibility.
6. Additional quality control analyses and technical validation used in the manuscript.

---

## 📖 Citation

If you use this repository or the associated dataset, please cite:

> Henry, M. *et al.*  
> **A metabolomics dataset of plankton exometabolomes and dissolved organic compounds in the South Atlantic Ocean**.
