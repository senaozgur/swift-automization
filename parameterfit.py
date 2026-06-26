import subprocess
import os
from pathlib import Path
from xspec import *
import numpy as np
from astropy.io import fits
import math
import re

# ============================================================
#   parameter.py  —  configuration file for swift_fit.py
#   Swift XRT version (adapted from NICER parameter.py)
# ============================================================

# ---------------------- Common variables --------------------

# Full path to the directory where all output files will be created.
# Leave blank ("") to use the directory where swift_fit.py lives.
outputDir = ""

# Text file listing paths to the Swift XRT observation directories.
# Each line should be the full path to one <obsid>-xrt directory produced by swift.py.
# Example line:  /data/swift/00012345001-xrt
inputTxtFile = "swift.txt"

flux_script_name = "swift_flux.py"
plot_script_name = "swift_plot.py"
# Log file name for fit results (written inside each observation's results/ folder)
resultsFile = "fit_results.log"

# XSPEC energy filter — Swift XRT is reliable from 0.3 to 10 keV.
# Channels outside this range are ignored.
energyFilter = "0.8 10.0"

# -------------------- swift_fit script switches -------------

# Script name (used in print statements)
fit_script_name = "swift_fit.py"

# -------------------- Fitting options ----------------------

# Set True to delete all previous results and start version numbering from 1 again
clean_result_history = False

# XSPEC abundance table
xspec_abundance = "wilm"

# Name of the models definition file (DSL text file — same format as NICER)
model_file = "models.txt"

# Name of the model pipeline to run (must match a "model <name>" entry in models.txt)
# Use underscores instead of spaces.
model_pipeline_name = "model_soft"

# -------------------- Parameter fixing ---------------------
# If True, the script fits the first fix_sample_size observations freely,
# takes the exposure-weighted average of parametersToFix, then refits all
# observations with those parameters frozen.
fix_parameters_after_sampling = False
fix_sample_size = 10
parametersToFix = ["TBabs.nH"]

# restartOnce: delete commonFiles model cache before the first observation only
# restartAlways: delete it before every observation
restartOnce = True
restartAlways = False

# -------------------- F-test -------------------------------
ftestSignificance = 0.05

# -------------------- Verbosity ----------------------------
chatterOn = False

# -------------------- Error calculations -------------------
# If True, shakefit will be run to calculate 90% confidence intervals
errorCalculations = True

# If True, freeze powerlaw photon index when XSPEC error > 1
checkPowerlawErrorAndFreeze = False
powerlawIndexToFreezeAt = 1.7

# Parameters for which error intervals will be calculated by shakefit.
# Keys: "component.parameter" exactly as XSPEC names them.
# Values: label string used in the output parameter file (use _ instead of spaces).
# Add or remove entries to match the models you use.
parametersForShakefit = {
    "TBabs.nH":          "TBabs_nH",
    "powerlaw.PhoIndex": "Powerlaw_index",
    "powerlaw.norm":     "Normalization_(powerlaw)",
    "diskbb.Tin":        "Tin_(keV)",
    "diskbb.norm":       "Normalization_(diskbb)",
    "nthComp.Gamma":     "nthComp_Gamma",
    "nthComp.kT_e":      "nthComp_kT_e_(keV)",
    "nthComp.kT_bb":     "nthComp_kT_bb_(keV)",
    "nthComp.norm":      "nthComp_norm",
    "gaussian.LineE":    "gaussian_LineE_(keV)",
    "gaussian.Sigma":    "gaussian_Sigma_(keV)",
    "gaussian.norm":     "gaussian_norm",
    
}

# -------------------- Versioning / output naming -----------
# Custom prefix for result folder names. Leave "" to use model_pipeline_name.
custom_name = ""


#================================================ nicer.flux spesific variables =================================================
# nicer_flux will add "cflux" component before the spesified models below to calculate flux.
# 'absorbed' keyword adds cflux at the beginning, 'unabsorbed' keyword adds cflux after the last absorption model.
# Be careful while using 'unabsorbed' keyword, do not use it if the above usage does not give the unabsorbed flux.
fluxes_to_be_calculated = ["absorbed", "unabsorbed", "simpl", "diskbb", "powerlaw","nthComp","gaussian"]

# Used for 'unabsorbed' flux.
last_absorption_model = "TBabs"

# If set to true, bottom and top limits of parameters will be set to (value +/- 0.1) before fitting with cflux.
restrict_parameters = False

#================================================ nicer.plot spesific variables =================================================
# If set to True, the script will create new graphs with a count/version number at the end instead of updating only one file.
# e.g. model_parameters_1.png, model_parameters_2.png, ... As you continue to run the script, previous files will not be deleted.
# If set to False, the script will only create one file, and delete the previous one (e.g. model_parameters.png, without a count/version number at the end)
enable_versioning = True

# Setting this variable to True will clear all the previously created files (graphs and tables), and reset the count/version number
# to 1 if enable_versioning is set to True
delete_previous_files = False

# Custom name for naming graphs and tables. If you set 'custom_name' = "", then the model name used for fitting will be used for naming
# e.g: custom_name = "", graph name: model_simpl_edge_1.png OR custom_name = "nH_fixed", graph name = nH_fixed_1.png
custom_name = ""

# Modified z-score algorithm will be used for outlier detection
# Possibility of removing "good" data always exists, turn it on or off accordingly
use_outlier_detection = False

# Lower threshold value for modified z-score algorithm (Change it according to your needs)
outlier_lower_threshold = -10

# Upper threshold value for modified z-score algorithm (Change it according to your needs)
outlier_upper_threshold = 10