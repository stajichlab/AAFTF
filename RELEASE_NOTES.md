Release notes for AAFTF 
=======================

Automatic Assembly For The Fungi

* v0.2.1
  1. Fix some README docs 
  2. Sync zenodo with this release

* v0.2.0
  1. rework temporary file, prefix use throughout, not relying on a set working directory so that multiple runs can be completed in a single folder
without clashes of names. 
  1. Switch to BBTools bbduk for quality trimming instead of Trimmomatic
  1. Switch to BBTools bbduk instead of BWA/Bowtie for mapping to contamination db. Speed improvements and accuracy of kmer matches to contamination over alignment-based cleanup.
  1. Specify memory use for pilon, spades assembly steps
  1. Support to generate scripts to run all steps in a pipeline with 'pipeline' cmd

* v0.1.0
  1. Initial release which includes automation for quality trimming and contaminant/adaptor removal, filtering, assembly, veccleanup, refinement, and duplication removal and sorting
