==============================================================================
🧬 RATIONAL PRIMER DESIGN PIPELINE - USER MANUAL
==============================================================================

Welcome! This tool automatically designs and validates TaqMan primers/probes 
for diagnostic assays.

------------------------------------------------------------------------------
⚙️ CONFIGURATION (config/parameters.json)
------------------------------------------------------------------------------
You can control the biology of the design by editing 'config/parameters.json':

  1. Sampling Sizes:
     - Set to 0 to use ALL available genomes (Highest Accuracy).
     - Set to small numbers (e.g., 10) for quick tests.

  2. Biological Specs:
     - product_size_min/max: The allowed length of the PCR amplicon.
     - min_sensitivity: Discard primers that don't detect at least this % of targets.

------------------------------------------------------------------------------
🚀 HOW TO RUN
------------------------------------------------------------------------------
1. Double-click "RUN_PIPELINE.bat".
2. Enter a Project Name (e.g., "Salmonella_Test").
3. When asked, you can drag "config/parameters.json" into the window 
   (or just press Enter to use it automatically).
4. Select your mode (Download from NCBI or Local Files).
5. Drag and drop your data folders if requested.

------------------------------------------------------------------------------
📊 OUTPUT
------------------------------------------------------------------------------
Check your new project folder for "FINAL_ASSAY.csv".
- It contains sequences, Tm, Sensitivity, and Specificity stats.

==============================================================================