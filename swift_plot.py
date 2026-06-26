# Authors: Batuhan Bahçeci
# Contact: batuhan.bahceci@sabanciuniv.edu

from parameterfit import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

print("==============================================================================")
print("\t\t\tRunning " + plot_script_name + "\n")

# Find the script's own path
scriptPath = os.path.abspath(__file__)
scriptPathRev = scriptPath[::-1]
scriptPathRev = scriptPathRev[scriptPathRev.find("/") + 1:]
scriptDir = scriptPathRev[::-1]
os.chdir(scriptDir)

# Check if outputDir has been assigned to be a spesific directory or not
# If not, assign outputDir to the directory where the script is located at
if outputDir == "":
    outputDir = scriptDir

#===================================================================================================================================
# Input check for outputDir
if Path(outputDir).exists() == False:
    print("ERROR: Directory defined by outputDir could not be found. Terminating the script...")
    quit()

# Input check for commonFiles under outputDir
if Path(outputDir + "/commonFiles").exists() == False:
    print(f"Directory {outputDir}/commonFiles could not be found. You need to create output files/directories by running {create_script_name} first.")
    print("Terminating the script...")
    quit()

# Input check for outlier_lower_threshold
if type(outlier_lower_threshold) != float and type(outlier_lower_threshold) != int:
    print("ERROR: Incorrect type for 'outlier_lower_threshold' variable in parameter.py")
    quit()

# Input check for outlier_upper_threshold
if type(outlier_upper_threshold) != float and type(outlier_upper_threshold) != int:
    print("ERROR: Incorrect type for 'outlier_upper_threshold' variable in parameter.py")
    quit()

# Input check for the ranges of outlier_lower_threshold and outlier_upper_threshold
if (outlier_lower_threshold >= outlier_upper_threshold):
    print("ERROR: Lower threshold for outlier detection algorithm is larger than upper threshold")
    quit()

# Input check for model_pipeline_name
if model_pipeline_name == "":
    print("model_pipeline_name is not provided in parameter.py")
    print("Please provide a valid model name that is defined in models.txt and try running the script again.")
    quit()
#===================================================================================================================================

#===================================================================================================================================
# Functions
def listToStr(array):
    result = ""
    for char in array:
        result += str(char) + " "
    result = result[:-1]
    return result

def enterParameters(parList, fluxPars = {}):
    # Take the available parameters from the stemList
    components = AllModels(1).componentNames
    for comp in components:
        compObj = getattr(AllModels(1), comp)
        parameters = compObj.parameterNames
        for par in parameters:
            parObj = getattr(compObj, par)
            indx = parObj.index
            fullName = comp + "." + par
            # If the parameter value is stored in previous model's list, take it from there. Otherwise, look whether
            # it is initialized in current model's list
            if fullName in parList:
                AllModels(1)(indx).values = parList[fullName]
    
    if fluxPars != {}:
        compObj = AllModels(1).cflux
        for key, val in fluxPars.items():
            parObj = getattr(compObj, key)
            parObj.values = val

def fitModel():
    Fit.nIterations = 100
    Fit.delta = 0.01
    Fit.perform()

def updateParameters(parList):
    # Save the parameters loaded in the xspec model to lists
    components = AllModels(1).componentNames
    for comp in components:
        if comp == "vphabs":
            compObj = getattr(AllModels(1), comp)
            parameters = compObj.parameterNames
            for par in parameters:
                parObj = getattr(compObj, par)
                if parObj.values[1] > 0:
                    indx = parObj.index
                    fullName = comp + "." + par
                    parList.append((fullName, AllModels(1)(indx).values[0]))

def transferToNewList(sourceList):
    newList = [sourceList[0]]
    newParDict = {}
    newStatDict = {}

    sourceParDict = sourceList[1]
    keys = list(sourceParDict.keys())
    values = list(sourceParDict.values())
    for i in range(len(keys)):
        newParDict[keys[i]] = values[i]
    newList.append(newParDict)

    sourceStatDict = sourceList[2]
    keys = list(sourceStatDict.keys())
    values = list(sourceStatDict.values())
    for i in range(len(keys)):
        newStatDict[keys[i]] = values[i]
    newList.append(newStatDict)
    
    return newList

def modified_z_score(parameter_dict):
    for key, all_arrays in parameter_dict.items():
        data_array = np.array(all_arrays[0])
        if len(data_array) == 0:
            continue

        median = np.median(data_array)
        mean = np.mean(data_array)
        mad = np.median(np.abs(data_array - median))
        mean_ad = np.abs(data_array - mean)
        mean_ad = np.sum(mean_ad) / len(data_array)

        mod_z_scores = []
        for x in data_array:
            if mad == 0:
                mod_z_scores.append((x - median) / (1.253314 * mean_ad))
            else:
                mod_z_scores.append((x - median) / (1.486 * mad))

        idx_counter = 0
        for i in mod_z_scores:
            if i  <= outlier_lower_threshold  or i >= outlier_upper_threshold:
                all_arrays[0].pop(idx_counter)
                all_arrays[1].pop(idx_counter)
                all_arrays[2].pop(idx_counter)
                all_arrays[3].pop(idx_counter)
            else:
                idx_counter += 1
        
#===================================================================================================================

try:
    energyLimits = energyFilter.split(" ")
    Emin = energyLimits[0]
    Emax = energyLimits[1]
except Exception as e:
    print(f"Exception occured while reading 'energyLimits' variable due to incorrect format: {e}")
    quit()

otherParsDict = {}
fluxValuesDict = {}
commonDirectory = outputDir + "/commonFiles"   # ~/NICER/analysis/commonFiles

#===========================================================================================
# Create output directories of not created already
if Path(commonDirectory + "/results").exists() == False:
    os.system("mkdir " + commonDirectory + "/results")

if Path(commonDirectory + "/results/model_graphs").exists() == False:
    os.system("mkdir " + commonDirectory + "/results/model_graphs")

if Path(commonDirectory + "/results/flux_graphs").exists() == False:
    os.system("mkdir " + commonDirectory + "/results/flux_graphs")

if Path(commonDirectory + "/results/flux_tables").exists() == False:
    os.system("mkdir " + commonDirectory + "/results/flux_tables")

if Path(commonDirectory + "/results/model_tables").exists() == False:
    os.system("mkdir " + commonDirectory + "/results/model_tables")

if Path(commonDirectory + "/results/chi_squared_graphs").exists() == False:
    os.system("mkdir " + commonDirectory + "/results/chi_squared_graphs")

if Path(commonDirectory + "/version_counter.txt").exists() == False:
    os.system("touch " + commonDirectory + "/version_counter.txt")

    temp_file = open(commonDirectory + "/version_counter.txt", "w")
    temp_file.write("CREATED BY NICER_PLOT.PY, DO NOT MODIFY, DO NOT CHANGE THE FILE PATH\n")
    temp_file.write("0\n")
    temp_file.close()
#===========================================================================================
# If clear variable is True, clear the contents of the output directories
if delete_previous_files:
    print("Deleting all the previous graph and table files under 'results' directory..\n")

    os.system("rm -r " + commonDirectory + "/results/model_graphs/*")
    
    os.system("rm -r " + commonDirectory + "/results/model_tables/*" )
    
    os.system("rm -r " + commonDirectory + "/results/flux_graphs/*")

    os.system("rm -r " + commonDirectory + "/results/flux_tables/*")

    os.system("rm -r " + commonDirectory + "/results/chi_squared_graphs/*")
    
    temp_file = open(commonDirectory + "/version_counter.txt", "w")
    temp_file.write("CREATED BY NICER_PLOT.PY, DO NOT MODIFY, DO NOT CHANGE THE FILE PATH\n")
    temp_file.write("0\n")
    temp_file.close()

#===========================================================================================
# Check whether enable_versioning is True, update the version file and extract the current version if that is the case.
current_version = 0
if enable_versioning:
    try:
        with open(commonDirectory + "/version_counter.txt", "r") as file:
            all_lines = file.readlines()
            prev_version = int(all_lines[1].strip("\n"))
            current_version  = prev_version + 1
    except Exception as e:
        print(f"Exception occured while opening version_counter.txt: {e}")
        quit()

    try:
        with open(commonDirectory + "/version_counter.txt", "w") as file:
            file.write("CREATED BY NICER_PLOT.PY, DO NOT MODIFY, DO NOT CHANGE THE FILE PATH\n")
            file.write(str(current_version) + "\n")
    except Exception as e:
        print(f"Exception occured while opening version_counter.txt: {e}")
        quit()

#===========================================================================================
# Name of the output files
output_save_name = custom_name

# Check whether reduced_chi.log can be opened
try:
    chi_file = open(commonDirectory + "/reduced_chi.log", "r")
except Exception as e:
    print(f"Exception occured while reading reduced_chi.log file under commonFiles directory: {e}")
    quit()

# Check the amount of lines in reduced_chi.log
all_lines = chi_file.readlines()
if (len(all_lines) <= 1):
    print("ERROR: reduced_chi.log file under commonFiles does not contain any information about fitted observations.")
    quit()

# First line of reduced_chi.log is dedicated for model name. If custom name is not provided, set the output file saving name to model name.
if (output_save_name == ""):
    output_save_name = all_lines[0].strip()

# Parse the log file and extract date and reduced chi squared values
chi_sq_dict = {}
for line in all_lines[1:]:
    line = line.strip()
    line = line.split(" ")
    chi_sq_dict[float(format(float(line[0]), ".1f"))] = float(format(float(line[1]), ".2f")) - 1

#===========================================================================================
# Open the input txt file, and extract the obsid numbers to a list
searchedObsid = []

try:
    with open(scriptDir + "/" + inputTxtFile, "r") as file:
        allLines = file.readlines()
        for line in allLines:
            line = line.replace(" ", "")
            line = line.strip("\n")
            if line != "" and Path(line).exists():
                if line[-1] != "/":
                    slashIdx = line.rfind("/")
                    obsid = line[slashIdx+1 :]
                else:
                    slashIdx = line[:-1].rfind("/")
                    obsid = line[slashIdx+1:-1]
                
                searchedObsid.append(obsid)
except Exception as e:
    print(f"Exception occured while opening {inputTxtFile}: {e}")
    quit()

if len(searchedObsid) == 0:
    print("\nCould not find any valid observation path given in the observations.txt file.")
    quit()

# Check whether any of the searched observations are processed, extract the necessary data if that is the case.
iterationMax = 0
searchedObservations = []
if Path(commonDirectory + "/processed_obs.txt").exists() == False:
    print("\nERROR: Could not find 'processed_obs.txt' file under the 'commonFiles' directory.")
    print("Please make sure both the 'commonFiles' directory and the 'processed_obs.txt' files exist and are constructed as intended by nicer_create.py.\n")
    quit()
else:
    try:
        with open(commonDirectory + "/processed_obs.txt", "r") as filteredFile:
            allLines = filteredFile.readlines()
            for eachObsid in searchedObsid:
                for line in allLines:
                    line = line.strip("\n")
                    lineElements = line.split(" ")

                    if lineElements[1] == eachObsid:
                        iterationMax += 1
                        searchedObservations.append((lineElements[0], lineElements[1], lineElements[2]))
    except Exception as e:
        print(f"Exception occured while opening processed_obs.txt: {e}")
        quit()

if len(searchedObservations) == 0:
    print("\nCould not find any matching observation paths in 'processed_obs.txt' with those in 'observations.txt', most likely them being filtered out due to having low exposure.")
    quit()


# Iterate through each of the observation data, and extract the necessary data from the parameter files
for path, obsid, expo in searchedObservations:
    outObsDir = path

    try:
        os.chdir(outObsDir)
    except Exception as e:
        print(f"Exception occured while changing directory to {outObsDir}: {e}")
        continue

    allFiles = os.listdir(outObsDir)

    additional_info = ""
    if path[-3:] == "day":
        additional_info = "day"
    elif path[-5:] == "night":
        additional_info = "night"

    dict_key = str(obsid)
    if additional_info != "":
        dict_key += "_" + additional_info

    # Find the data file and the best fitting model file for the current observation.
    # Swift XRT spectra are named sw<obsid>*grp*.pha — we use a glob so the script
    # works regardless of the exact suffix (XPC, PC, WT, etc.).
    import glob as _glob
    missingFiles = True
    foundParameterfile = False
    foundModfile = False
    foundSpectrum = False
    spectrumFile = None
    for file in allFiles:
        if "parameters_" in file:
            parFile = file
            foundParameterfile = True
        elif "best_" in file:
            modFile = file
            foundModfile = True

    # Swift spectrum: grouped .pha in the obs directory
    pha_candidates = (
        _glob.glob(outObsDir + "/sw" + obsid + "*grp*.pha") +
        _glob.glob(outObsDir + "/sw" + obsid + "*_grp.pha") +
        _glob.glob(outObsDir + "/*grp*.pha") +
        _glob.glob(outObsDir + "/*.pha")
    )
    if pha_candidates:
        spectrumFile = pha_candidates[0]
        foundSpectrum = True

    if foundModfile and foundParameterfile and foundSpectrum:
        missingFiles = False

    # Check if there are any missing files
    if missingFiles:
        print("\nWARNING: Necessary files for retrieving parameters and plotting are missing for observation: " + obsid)
        if foundSpectrum == False:
            print("Missing spectrum file (looked for sw<obsid>*grp*.pha in " + outObsDir + ")")
        if foundModfile == False:
            print("Missing model file")
        if foundParameterfile == False:
            print("Missing parameter file")
        continue

    try:
        hdu = fits.open(spectrumFile)
    except Exception as e:
        print(f"Exception occured while opening {spectrumFile}: {e}")
        continue

    # Extract MJD — try the standard keywords in order of preference
    mjd_obs = None
    for kw in ("MJD-OBS", "MJD_OBS", "MJDREFI", "TSTART"):
        try:
            mjd_obs = float(hdu[1].header[kw])
            break
        except (KeyError, TypeError):
            continue
    hdu.close()

    if mjd_obs is None:
        print(f"WARNING: Could not read MJD from spectrum header for {obsid}, skipping.")
        continue

    date = float(format(mjd_obs, ".3f"))

    Xset.chatter = 0
    try:
        Xset.restore(modFile)
    except Exception as e:
        print(f"Exception occured while loading {modFile} to PyXspec")
        continue


    try:
        file = open(parFile)
    except Exception as e:
        print(f"Exception occured while opening {parFile}: {e}")
        continue

    allLines = file.readlines()
    for line in allLines[1:]:
        # Extract the parameter values with associated uncertainities
        line = line.strip("\n")
        line = line.split(" ")
        # Only replace underscores in the unit field (index 4+), NOT in the
        # parameter name (index 0). The old loop replaced everything, turning
        # "absorbed_flux" -> "absorbed flux" so column-header lookups always
        # failed and every data cell was written as "0".
        for i in range(4, len(line)):
            line[i] = line[i].replace("_", " ")



        line[1] = float(line[1])    # Parameter value
        line[2] = line[1] - float(line[2])  # Error low
        line[3] = float(line[3]) - line[1]  # Error high

        # Check whether the parameter has a given unit in parameter file.
        try:
            test = line[4]
        except:
            line.append("")

            

        par_list = [line[0], line[1], line[2], line[3], line[4], date]

        hardness_keywords = ["soft_rate", "hard_rate", "hardness_ratio"]
        
        
        if "flux" in par_list[0] or "rate" in par_list[0] or "ratio" in par_list[0]:
            if dict_key in fluxValuesDict:
                fluxValuesDict[dict_key].append(par_list)
            else:
                fluxValuesDict[dict_key] = [par_list]
        else:
            if dict_key in otherParsDict:
                otherParsDict[dict_key].append(par_list)
            else:
                otherParsDict[dict_key] = [par_list]

otherpars_ref = 0
fluxes_ref = 0
empty_flag = True

print("\n")

#Find the referance point for the x-axis (date) for model parameters
if len(otherParsDict) != 0:
    empty_flag = False
    dates = []

    # Extract all dates in MJD
    for obs_list in otherParsDict.values():
        for value_list in obs_list:
            dates.append(value_list[5])

    referanceMjd = round((min(dates) - 5) / 5) * 5
    otherpars_ref = referanceMjd

    for key, obs_list in otherParsDict.items():
        for value_list in obs_list:
            value_list[5] = value_list[5] - referanceMjd
else:
    print("WARNING: Could not find any model parameters inside the parameter_bestmodel.txt file.")
    print("Creating graph and table files will be skipped..\n")


#Find the referance point for the x-axis (date) for flux values
if len(fluxValuesDict) != 0:
    empty_flag = False
    dates = []

    # Extract all dates in MJD
    for obs_list in fluxValuesDict.values():
        for value_list in obs_list:
            dates.append(value_list[5])

    referanceMjd = round((min(dates) - 5) / 5) * 5
    fluxes_ref = referanceMjd

    for key, obs_list in fluxValuesDict.items():
        # FIXED: loop through obs_list instead of fluxValuesDict.items()
        for value_list in obs_list:
            value_list[5] = value_list[5] - referanceMjd
else:
    print("WARNING: Could not find any flux values inside the parameter_bestmodel.txt file.")
    print("Creating graph and table files will be skipped..\n")

# Both flux dictionary and model dictionary are empty.
if empty_flag:
    print("\nERROR: Both dictionaries (parameter/flux) are empty. There is no data to create any graph.\n")
    quit()


# If otherParsDict is not empty, create the model parameter graphs and the table file
if len(otherParsDict) != 0:
    parameter_dict = {}

    # Extract the the data from each observation to parameter_dict, where each key will be parameter names whereas values will be lists of parameter values
    for obs_identifier, obs_list in otherParsDict.items():
        for value_list in obs_list:
            if value_list[0] not in parameter_dict:
                parameter_dict[value_list[0]] = [[value_list[1]], [value_list[2]], [value_list[3]], [value_list[5]], value_list[4]]
            else:
                parameter_dict[value_list[0]][0].append(value_list[1])  # Value
                parameter_dict[value_list[0]][1].append(value_list[2])  # Error low
                parameter_dict[value_list[0]][2].append(value_list[3])  # Error high
                parameter_dict[value_list[0]][3].append(value_list[5])  # Date

    if use_outlier_detection:
        print("="*100)
        print("Modified z-score algorithm will be applied for model parameters")
        print("="*100, "\n")

        modified_z_score(parameter_dict)


    fig, axs = plt.subplots(len(parameter_dict), 1, figsize=(8, 14), sharex=True)

    shared_yaxis_flag = False
    all_parameters = []
    ax_num = 0
    for par_name, par_list in parameter_dict.items():

        # Add all unique parameters to all_parameters list
        if par_name not in all_parameters:
            all_parameters.append(par_name)

        # Extract all lists necessary for graphs
        x_axis = par_list[3]
        y_axis = par_list[0]
        err_low = par_list[1]
        err_high = par_list[2]
        
        ####
            
        # Convert to numpy arrays and fix signs / bad values
        err_low  = np.array(err_low,  dtype=float)
        err_high = np.array(err_high, dtype=float)



        # Check if shared y-axis title is given
        shared_yaxis_title = par_list[4]
        if shared_yaxis_title != "":
            shared_yaxis_flag = True

        axs[ax_num].errorbar(x_axis, y_axis, yerr=[err_low, err_high], fmt='o', color='black', ecolor="black", markersize=4, capsize=0)

        axs[ax_num].tick_params(which = "both", direction="in")
        axs[ax_num].yaxis.tick_left()

        axs[ax_num].set_ylabel(par_name)
        axs[ax_num].set_xlabel(f"Time (MJD-{otherpars_ref} days)")

        ax_num += 1
    
    # Set minor ticks, also hide x-axis tick labels from all graphs except the last one
    for ax in axs:
        ax.xaxis.set_minor_locator(AutoMinorLocator())

        if ax != axs[-1]:
            ax.xaxis.set_tick_params(labelbottom=False)

    # Set the title of the figure
    fig.suptitle('Model Parameters', fontsize=20, y=0.95)

    # Set a shared y-axis title, if it is given
    if shared_yaxis_flag:
        fig.text(0.9, 0.5, shared_yaxis_title, va='center', rotation='vertical')

    plt.subplots_adjust(wspace=0, hspace=0, right=0.85)  

    # Construct the file name of the graph
    if enable_versioning:
        png_name = commonDirectory + "/results/model_graphs/" + output_save_name + "_" + str(current_version) + ".png"
    else:
        png_name = commonDirectory + "/results/model_graphs/" + output_save_name + ".png"

    # Delete any existing file with the same name, and create a new file
    if Path(png_name).exists():
        os.system("rm " + png_name)

    # Save the graph file
    plt.savefig(png_name)

    # Initialize the table structure with lists, where each list denotes a column 
    table_columns = []
    for i in range(len(all_parameters)*3 + 2):
        table_columns.append([])

    # Column name of the first two columns
    table_columns[0].append("Obsid")
    table_columns[1].append("MJD")

    # Column names of the rest of the columns
    start_index = 2
    for par in all_parameters:
        par = par.replace(" ", "_")
        table_columns[start_index].append(par)
        table_columns[start_index + 1].append(par + "_errlow")
        table_columns[start_index + 2].append(par + "_errhigh")
        start_index += 3

    # Create the obsid column
    for obsid in otherParsDict.keys():
        table_columns[0].append(obsid)

    # Create the date column
    for val in otherParsDict.values():
        table_columns[1].append(val[0][5] + otherpars_ref)

    # Create the rest of the columns
    start_index = 2
    for i in range(len(all_parameters)):
        # Calculate the current column index (each parameter will have three columns: value, low error, high error)
        current_index = start_index  + i*3

        # Direct match using the unmodified parameter name from all_parameters
        searched_par = all_parameters[i]

        for key, obs_list in otherParsDict.items():
            added_par_flag = False

            for val in obs_list:

                par_name = val[0]
                if par_name == searched_par:
                    # The current parameter in the dictionary matches with the column name, extract the parameter values and set the flag to True
                    added_par_flag = True

                    table_columns[current_index].append(val[1])
                    table_columns[current_index + 1].append(val[2])
                    table_columns[current_index + 2].append(val[3])
                    break
            
            # If the flag has not been set to True, it means the current observation does not have the searched parameter of the column. Set values as "-"
            if added_par_flag == False:
                table_columns[current_index].append("0")
                table_columns[current_index + 1].append("0")
                table_columns[current_index + 2].append("0")

    # Set the table file name according to enable_versioning variable
    table_file_name = ""
    if enable_versioning:
        table_file_name = commonDirectory + "/results/model_tables/" + output_save_name + "_" + str(current_version) + ".txt"
    else:
        table_file_name = commonDirectory + "/results/model_tables/" + output_save_name + ".txt"
    
    # Create the table file if it has not been already created
    if Path(table_file_name).exists() == False:
        os.system("touch " + table_file_name)
    
    # Override the table file's contents with the columns created within table_columns
    with open(table_file_name, "w") as file:
        for i in range(len(table_columns[0])):
            line = ""
            for each_line in table_columns:
                line += str(each_line[i]) + " "
        
            line = line[:-1]
            line += "\n"

            file.write(line)

    print("Graph and table files for model parameters have been successfully created:")
    print("Graph path: " + png_name)
    print("Table path: " + table_file_name + "\n")



# If fluxValuesDict is not empty, create the model parameter graphs and the table file
if len(fluxValuesDict) != 0:
    parameter_dict = {}

    # Extract the the data from each observation to parameter_dict, where each key will be parameter names whereas values will be lists of parameter values
    for obs_identifier, obs_list in fluxValuesDict.items():
        for value_list in obs_list:
            if value_list[0] not in parameter_dict:
                parameter_dict[value_list[0]] = [[value_list[1]], [value_list[2]], [value_list[3]], [value_list[5]], value_list[4]]
            else:
                parameter_dict[value_list[0]][0].append(value_list[1])  # Value
                parameter_dict[value_list[0]][1].append(value_list[2])  # Error low
                parameter_dict[value_list[0]][2].append(value_list[3])  # Error high
                parameter_dict[value_list[0]][3].append(value_list[5])  # Date

    if use_outlier_detection:
        print("="*100)
        print("Modified z-score algorithm will be applied for flux values")
        print("="*100, "\n")

        modified_z_score(parameter_dict)

    fig, axs = plt.subplots(len(parameter_dict), 1, figsize=(8, 14), sharex=True)

    shared_yaxis_flag = False
    all_parameters = []
    ax_num = 0
    for par_name, par_list in parameter_dict.items():

        # Add all unique parameters to all_parameters list
        if par_name not in all_parameters:
            all_parameters.append(par_name)

        # Extract all lists necessary for graphs
        x_axis = np.array(par_list[3], dtype=float)
        y_axis = np.array(par_list[0], dtype=float)
        err_low = np.array(par_list[1], dtype=float)
        err_high = np.array(par_list[2], dtype=float)
        shared_yaxis_title = par_list[4]
        
        huge_thresh = 1e8
        mask = (
            (err_low  > 0) &
            (err_high > 0) &
            (np.abs(err_low)  < huge_thresh) &
            (np.abs(err_high) < huge_thresh)
        )
    
        if not np.all(mask):
            bad_idx = np.where(~mask)[0]
            print(f"WARNING: dropping {len(bad_idx)} bad points for {par_name}")
            print("  indices:", bad_idx.tolist())
    
        
        x_axis   = x_axis[mask]
        y_axis   = y_axis[mask]
        err_low  = err_low[mask]
        err_high = err_high[mask]
        
         
        if len(y_axis) == 0:
            print(f"WARNING: all points dropped for {par_name}, skipping plot.")
            continue
    
        
        axs[ax_num].errorbar(
            x_axis, y_axis,
            yerr=[err_low, err_high],
            fmt='o', color='black', ecolor="black",
            markersize=4, capsize=0
        )
    
        axs[ax_num].tick_params(which="both", direction="in")
        axs[ax_num].yaxis.tick_left()
        axs[ax_num].set_ylabel(par_name)
        axs[ax_num].set_xlabel(f"Time (MJD-{fluxes_ref} days)")  
        ax_num += 1
        
    # Set minor ticks, also hide x-axis tick labels from all graphs except the last one
    for ax in axs:
        ax.xaxis.set_minor_locator(AutoMinorLocator())

        if ax != axs[-1]:
            ax.xaxis.set_tick_params(labelbottom=False)

    # Set the title of the figure
    fig.suptitle('Flux Values', fontsize=20, y=0.95)

    # Set a shared y-axis title, if it is given
    if shared_yaxis_flag:
        fig.text(0.9, 0.5, shared_yaxis_title, va='center', rotation='vertical')

    plt.subplots_adjust(wspace=0, hspace=0, right=0.85)  

    # Construct the file name of the graph
    if enable_versioning:
        png_name = commonDirectory + "/results/flux_graphs/" + output_save_name + "_" + str(current_version) + ".png"
    else:
        png_name = commonDirectory + "/results/flux_graphs/" + output_save_name + ".png"

    # Delete any existing file with the same name, and create a new file
    if Path(png_name).exists():
        os.system("rm " + png_name)

    # Save the graph file
    plt.savefig(png_name)

    # Initialize the table structure with lists, where each list denotes a column 
    table_columns = []
    for i in range(len(all_parameters)*3 + 2):
        table_columns.append([])

    # Column name of the first two columns
    table_columns[0].append("Obsid")
    table_columns[1].append("MJD")

    # Column names of the rest of the columns
    start_index = 2
    for par in all_parameters:
        par = par.replace(" ", "_")

        table_columns[start_index].append(par)
        table_columns[start_index + 1].append(par + "_errlow")
        table_columns[start_index + 2].append(par + "_errhigh")
        start_index += 3

    # Create the obsid column
    for obsid in fluxValuesDict.keys():
        table_columns[0].append(obsid)

    # Create the date column
    for val in fluxValuesDict.values():
        table_columns[1].append(val[0][5] + fluxes_ref)

    # Create the rest of the columns
    start_index = 2
    for i in range(len(all_parameters)):
        # Calculate the current column index (each parameter will have three columns: value, low error, high error)
        current_index = start_index  + i*3

        # Change the parameter name from the column to match with the name format in the dictionary
        searched_par = table_columns[current_index][0]

        for key, obs_list in fluxValuesDict.items():
            added_par_flag = False

            for val in obs_list:
                
                # Fix the parameter name to match searched parameter name format
                par_name = val[0].replace(" ", "_")

                if par_name == searched_par:
                    # The current parameter in the dictionary matches with the column name, extract the parameter values and set the flag to True
                    added_par_flag = True

                    table_columns[current_index].append(val[1])
                    table_columns[current_index + 1].append(val[1] - val[2])
                    table_columns[current_index + 2].append(val[1] + val[3])
                    break
            
            # If the flag has not been set to True, it means the current observation does not have the searched parameter of the column. Set values as "0"
            if added_par_flag == False:
                table_columns[current_index].append("0")
                table_columns[current_index + 1].append("0")
                table_columns[current_index + 2].append("0")

    # Set the table file name according to enable_versioning variable
    table_file_name = ""
    if enable_versioning:
        table_file_name = commonDirectory + "/results/flux_tables/" + output_save_name + "_" + str(current_version) + ".txt"
    else:
        table_file_name = commonDirectory + "/results/flux_tables/" + output_save_name + ".txt"
    
    # Create the table file if it has not been already created
    if Path(table_file_name).exists() == False:
        os.system("touch " + table_file_name)
    
    # Override the table file's contents with the columns created within table_columns
    with open(table_file_name, "w") as file:
        for i in range(len(table_columns[0])):
            line = ""
            for each_line in table_columns:
                line += str(each_line[i]) + " "
        
            line = line[:-1]
            line += "\n"

            file.write(line)

    print("Graph and table files for flux values have been successfully created:")
    print("Graph path: " + png_name)
    print("Table path: " + table_file_name + "\n")

#==========================================================================================
    
plt.clf()

chi_x = list(chi_sq_dict.keys())
chi_y = list(chi_sq_dict.values())

# Plot the bar graph
plt.bar(chi_x, chi_y, bottom=1)
plt.ylim(0, 5)

# Plot the red dotted line at y=1
plt.axhline(y=1, color='red', linestyle='dashed')

# Add labels and title
plt.xlabel('Date (MJD )')
plt.ylabel("Reduced Chi-Squared Values")
plt.title('Reduced Chi-Squared Values of all observations (y > 5 not shown on graph)')

if enable_versioning:
    plt.savefig(commonDirectory + "/results/chi_squared_graphs/" + output_save_name + "_" + str(current_version) + ".png")
else:
    plt.savefig(commonDirectory + "/results/chi_squared_graphs/" + output_save_name + ".png")

print("Chi-squared graph has been successfully created under '" + outputDir + "/commonFiles/results/chi_squared_graphs'")