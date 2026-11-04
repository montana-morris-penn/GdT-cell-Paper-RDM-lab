GIST bulk TCR reanalysis - complete R notebook

Open TCR_reanalysis_clean_R.ipynb in VS Code connected to WSL Ubuntu-24.04.
Choose the existing "R (datascience)" kernel (kernel ID: datascience-r).
Interpreter: /home/yuntianf/softwares/miniconda3/envs/datascience/bin/R
R version used for validation: 4.5.3.

Run all cells from the beginning. The notebook contains every analysis and
drawing function inline: raw exports -> main Figure 3A-G -> combined main
figure -> S-TCR1-7 -> provenance and validation. Each figure has concise
calculation notes, executable code and its plot output. It does not read
precomputed figure images or result tables. Python is not required for the
analysis; the optional Python build/execution utilities only manage notebooks.

Input: the 12 original *_ydT.tsv immunoSEQ files in data/TCR_full_data under
an ancestor of the notebook directory, or the directory specified by the
GIST_TCR_INPUT_DIR environment variable. Raw data are not included in the code
archive. Preserve their original filenames and contents.

Output: reproduced_outputs_R beside the notebook, or GIST_TCR_R_OUTPUT_DIR.
Output subfolders contain main figures, requested supplements, explicitly
optional supplements, separate TXT legends, source TSVs, PNG previews and logs.
Main Figure 3D and S-TCR1 treemap panels are 90 x 90 mm. Text is 20-28 pt in
native PDF dimensions. Reduce the entire page only if a smaller manuscript
layout is desired, since doing so also reduces its printed font sizes.

The original annotation filter is retained as data_filter_original. The
active analysis also requires sequenceStatus == "In". Counts are estimated
templates. HG610 PBMC (uncertain identity) and HG606 tumor (sparse) are excluded;
HG610 and HG624 tumors remain unpaired. Four confirmed patients supply paired
analyses. Main E-G preserve the original treatment-group statistical framework;
patient-level analyses are provided in the supplements.

R uses Mersenne-Twister/Inversion/Rejection with seed 20261001. Sampling is
without replacement at 2,578 templates, 1,000 iterations, with 200 draws for
each richness-curve depth. Deterministic quantities are compared against the
Python analysis. Individual sampling draws differ between R and Python because
their random-number generators differ. This is expected.

Runtime/package versions and input SHA256 hashes are recorded in
reproduced_outputs_R/logs/R_analysis_provenance.json and sessionInfo.txt.
The plain R export TCR_reanalysis_clean_R.R contains the same notebook code,
with Markdown carried as comments. When running it as a script, set the working
directory to the submission folder before source() to retain the same paths.

Validation completed using the existing WSL datascience-r kernel.
22 code cells executed; 15 inline figures; 0 error outputs.
23 deterministic tables matched the Python analysis.
All 17 PDFs have native font sizes between 20 and 28 points.
The executed notebook was run from an empty R output folder.
