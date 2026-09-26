DMI-GEM thesis final figure/equation code package

Main scripts
------------
code/generate_final_revision_v2.py
  Regenerates the final user-requested revisions:
  - Figure 2.1: simplified black-and-white study workflow with all text contained
  - Figure 3.1: expanded layout; panels C and D use full-width rows to prevent overlaps
  - Figure 3.4: additional spacing around panel B, titles, legend, and axis labels
  - Figure 3.5: larger three-panel layout with additional vertical spacing
  - Figure 3.6: larger six-panel layout with additional panel spacing
  - Equations 2.1-2.7: consistently rendered at the same STIX math font size

code/build_final_thesis_v2.py
  Inserts the regenerated figures and equations into the complete thesis DOCX while
  preserving the thesis text, captions, tables, headings, and the other figures.

code/generate_all_thesis_figures_original.py
  Original full standardized figure-generation script used for the preceding figure set.

source_tables/
  Frozen source tables used by the plotting scripts so the plots can be regenerated.

Python packages
---------------
numpy, pandas, matplotlib, Pillow, python-docx

The final delivered DOCX was built from:
  Bachelor_Thesis_26_COMPLETE_REVISED_FIGURES.docx
and written as:
  Bachelor_Thesis_26_FINAL_FIGURES_EQUATIONS.docx
