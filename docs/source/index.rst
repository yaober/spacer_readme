Welcome to SPACER's documentation!
===================================

**SPACER** (Spatial Cell Engagement Representation) is a multi-instance deep learning framework 
designed to elucidate the molecular and spatial determinants of cellular recruitment and engagement 
within tissue microenvironments.

The recruitment of various cell types into tissue sites and their local interactions play essential 
roles in tissue development, homeostasis, and disease physiology. However, studying these mechanisms 
on a genome-wide scale has been challenging due to the large number of cell types and genes involved.  

Spatially resolved transcriptomics (SRT) — especially **high-definition SRT** — now enables the 
quantitative characterization of cell positions and transcriptomes at unprecedented resolution.  
Yet, digesting such high-dimensional transcriptomic and positional information remains difficult.

To address this challenge, we developed **SPACER**, a biologically informed multi-instance deep learner 
that infers cellular localization rules directly from SRT data. SPACER was deployed across 
**10 high-definition** and **20 low-definition SRT datasets** to study how stromal and immune cell types 
are recruited into **tumors** and **hearts during myocarditis**.  

By integrating SPACER analyses with orthogonal **immunopeptidomics**, **spatial T cell receptor sequencing**, 
and **single-cell RNA sequencing**, we demonstrated that:
- Genes encoding **highly immunogenic peptides** or involved in **developmental pathways** 
  are potent drivers of T cell infiltration.
- SPACER reveals a **tumor-reactive gene signature** in T cells, validated by empirical expression 
  and spatial localization.
- During myocarditis, **CD4⁺ T cells**, though fewer in number, exhibit stronger responsiveness 
  than **CD8⁺ T cells**.

.. note::

   This project is under active development.  
   For questions or contributions, please visit the SPACER GitHub repository.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents

   installation
   usage
   api
