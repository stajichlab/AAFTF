Release notes for AAFTF 
=======================

Automatic Assembly For The Fungi


* v0.2.0
  1. rework temporary file, prefix use throughout, not relying on a set working directory so that multiple runs can be completed in a single folder
without clashes of names. 
  2. Switch to BBTools instead of Trimmomatic and BWA/Bowtie mapping to contamination db. Speed improvements and accuracy of kmer matches to contamination
over alignment-based cleanup
  3. specify memory use for pilon, spades assembly steps

* v0.1.0
  1. Initial release which includes automation for quality trimming and contaminant/adaptor removal, filtering, assembly, veccleanup, refinement, and duplication removal and sorting
