# This is an automatic NICER script for calculating the fluxes of Xspec models previously fitted by nicer_fit.py
# Authors: Batuhan Bahçeci
# Contact: batuhan.bahceci@sabanciuniv.edu

from parameterfit import *

additive_models = {}
convolution_models = {}
multiplicative_models = {}

additive_model_list =  "agauss      c6vmekl     eqpair      nei         rnei        vraymond \
agnsed      carbatm     eqtherm     nlapec      sedov       vrnei \
agnslim     cemekl      equil       npshock     sirf        vsedov \
apec        cevmkl      expdec      nsa         slimbh      vtapec \
bapec       cflow       ezdiskbb    nsagrav     smaug       vvapec \
bbody       compLS      gadem       nsatmos     snapec      vvgnei \
bbodyrad    compPS      gaussian    nsmax       srcut       vvnei \
bexrav      compST      gnei        nsmaxg      sresc       vvnpshock \
bexriv      compTT      grad        nsx         ssa         vvpshock \
bkn2pow     compbb      grbcomp     nteea       step        vvrnei \
bknpower    compmag     grbjet      nthComp     tapec       vvsedov \
bmc         comptb      grbm        optxagn     vapec       vvtapec \
bremss      compth      hatm        optxagnf    vbremss     vvwdem \
brnei       cph         jet         pegpwrlw    vcph        vwdem \
btapec      cplinear    kerrbb      pexmon      vequil      wdem \
bvapec      cutoffpl    kerrd       pexrav      vgadem      zagauss \
bvrnei      disk        kerrdisk    pexriv      vgnei       zbbody \
bvtapec     diskbb      kyrline     plcabs      vmcflow     zbknpower \
bvvapec     diskir      laor        posm        vmeka       zbremss \
bvvrnei     diskline    laor2       powerlaw    vmekal      zcutoffpl \
bvvtapec    diskm       logpar      pshock      vnei        zgauss \
bwcycl      disko       lorentz     qsosed      vnpshock    zkerrbb \
c6mekl      diskpbb     meka        raymond     voigt       zlogpar \
c6pmekl     diskpn      mekal       redge       vpshock     zpowerlw \
c6pvmkl     eplogpar    mkcflow     refsch"

multiplicative_model_list = "SSS_ice     constant    ismdust     polpow      wndabs      zphabs \
TBabs       cyclabs     log10con    pwab        xion        zredden \
TBfeo       dust        logconst    redden      xscat       zsmdust \
TBgas       edge        lyman       smedge      zTBabs      zvarabs \
TBgrain     expabs      notch       spexpcut    zbabs       zvfeabs \
TBpcf       expfac      olivineabs  spline      zdust       zvphabs \
TBrel       gabs        pcfabs      swind1      zedge       zwabs \
TBvarabs    heilin      phabs       uvred       zhighect    zwndabs \
absori      highecut    plabs       varabs      zigm        zxipab \
acisabs     hrefl       polconst    vphabs      zpcfabs     zxipcf \
cabs        ismabs      pollin      wabs"

convolution_model_list = "cflux       gsmooth     kerrconv    rdblur      simpl       xilconv \
cglumin     ireflect    kyconv      reflect     thcomp      zashift \
clumin      kdblur      lsmooth     rfxconv     vashift     zmshift \
cpflux      kdblur2     partcov     rgsxsrc     vmshift"


temp = additive_model_list.split(" ")
for i in temp:
    if i == "":
        pass
    else:
        additive_models[i] = 1

temp = multiplicative_model_list.split(" ")
for i in temp:
    if i == "":
        pass
    else:
        multiplicative_models[i] = 1

temp = convolution_model_list.split(" ")
for i in temp:
    if i == "":
        pass
    else:
        convolution_models[i] = 1


print("==============================================================================")
print("\t\t\tRunning " + flux_script_name + "\n")

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

#========================================================= Input Checks ============================================================
# Input check for outputDir
if Path(outputDir).exists() == False:
    print("Directory defined by outputDir could not be found. Terminating the script...")
    quit()

# Input check for commonFiles under outputDir
if Path(outputDir + "/commonFiles").exists() == False:
    print(f"Directory {outputDir}/commonFiles could not be found. You need to create output files/directories by running {create_script_name} first.")
    print("Terminating the script...")
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
    components = AllModels(1).componentNames
    # cflux is injected into the model, shifting the numeric suffix of all
    # subsequent components by +1 relative to the original parameter list.
    # We only apply the index correction to components that come after cflux.
    cflux_seen = False
    for comp in components:
        if comp == "cflux":
            cflux_seen = True
            continue  # cflux itself has no entry in the original parList

        compObj = getattr(AllModels(1), comp)
        parameters = compObj.parameterNames
        for par in parameters:
            parObj = getattr(compObj, par)
            fullName = comp + "." + par

            # Only remap the numeric suffix for components that follow cflux
            # and actually carry a numeric suffix (e.g. powerlaw_2 → powerlaw_1)
            if cflux_seen and "_" in comp:
                try:
                    suffix = comp[comp.rfind("_") + 1:]
                    base   = comp[:comp.rfind("_") + 1]
                    newComp = base + str(int(suffix) - 1)
                    fullName = newComp + "." + par
                except ValueError:
                    pass  # suffix is not an integer; keep original fullName

            if fullName in parList:
                parObj.values = parList[fullName]
    
    if fluxPars != {}:
        compObj = AllModels(1).cflux
        for key, val in fluxPars.items():
            parObj = getattr(compObj, key)
            parObj.values = val

def fitModel():
    Fit.nIterations = 100
    Fit.delta = 0.01
    Fit.renorm()
    Fit.perform()

def updateParameters(parList):
    # Save the parameters loaded in the xspec model to lists
    for comp in AllModels(1).componentNames:
        compObj = getattr(AllModels(1), comp)
        parameters = compObj.parameterNames
        for par in parameters:
            parObj = getattr(compObj, par)
            fullName = comp + "." + par
            parList[fullName] = parObj.values

def freezeNorm():
    # Freezes all current normalization parameters in a model, also sets new limits for parameters to vary
    comps = AllModels(1).componentNames
    for comp in comps:
        compObj = getattr(AllModels(1), comp)
        parNames = compObj.parameterNames
        for par in parNames:
            parObj = getattr(compObj, par)
            indx = parObj.index
            if par == "norm":
                parObj.frozen = True

            elif comp != "cflux":
                if restrict_parameters:
                    valString = str(parObj.values[0])+","+str(parObj.values[1])+","+str(parObj.values[0]-0.1)+","+str(parObj.values[0]-0.1)+","+str(parObj.values[0]+0.1)+","+str(parObj.values[0]+0.1)
                    AllModels(1)(indx).values = valString

def findFlux():
    parNums = AllModels(1).nParameters
    for i in range(1, parNums + 1):
        name = AllModels(1)(i).name
        if name == "lg10Flux":
            # Convert log10(x) flux to x
            flux = 10 ** AllModels(1)(i).values[0]
            Fit.error("maximum 1000 " + str(i))
            lowerFlux = 10 ** AllModels(1)(i).error[0]
            upperFlux = 10 ** AllModels(1)(i).error[1]

            flux      /= (10**-9)
            lowerFlux /= (10**-9)
            upperFlux /= (10**-9)
            return [flux, lowerFlux, upperFlux]

    # lg10Flux parameter not found — cflux was not added correctly
    print("WARNING: lg10Flux parameter not found in current model. Returning empty flux.")
    return []

def insert_cflux_before(modelName, component):
    """Insert 'cflux*' immediately before the FIRST occurrence of `component`
    that appears as a whole token (not as a substring of another model name).
    Uses a word-boundary-aware regex so e.g. 'pow' doesn't match 'powerlaw'."""
    pattern = r'(?<![A-Za-z0-9_])' + re.escape(component) + r'(?![A-Za-z0-9_])'
    return re.sub(pattern, "cflux*" + component, modelName, count=1)

def insert_cflux_after(modelName, component):
    """Insert '*cflux' immediately after the FIRST whole-token occurrence of `component`."""
    pattern = r'(?<![A-Za-z0-9_])' + re.escape(component) + r'(?![A-Za-z0-9_])'
    return re.sub(pattern, component + "*cflux", modelName, count=1)

def calculateFlux(component, modelName, parameters):
    if component.lower() == "unabsorbed":
        if last_absorption_model in AllModels(1).componentNames:
            newName = insert_cflux_after(modelName, last_absorption_model)
        else:
            print("'last_absorption_model' is not in the model expression, unabsorbed flux will not be calculated.")
            return []

    elif component.lower() == "absorbed":
        newName = "cflux*" + modelName

    else:
        if component in multiplicative_models:
            print(f"Cannot calculate the flux of a multiplicative model: {component}")
            return []

        elif component in convolution_models:
            cflux_before_conv = insert_cflux_before(modelName, component)
            m = Model(cflux_before_conv)
            enterParameters(parameters, {"Emin": Emin, "Emax": Emax})
            freezeNorm()
            fitModel()
            flux_before_conv = findFlux()
            if not flux_before_conv:
                print(f"Could not determine flux before convolution for {component}")
                return []

            cflux_after_conv = insert_cflux_after(modelName, component)
            m = Model(cflux_after_conv)
            enterParameters(parameters, {"Emin": Emin, "Emax": Emax})
            freezeNorm()
            fitModel()
            flux_after_conv = findFlux()
            if not flux_after_conv:
                print(f"Could not determine flux after convolution for {component}")
                return []

            conv_flux_val  = flux_before_conv[0] - flux_after_conv[0]
            max_upper_error = max(flux_before_conv[1] - flux_before_conv[0],
                                  flux_after_conv[1]  - flux_after_conv[0])
            max_lower_error = max(flux_before_conv[0] - flux_before_conv[2],
                                  flux_after_conv[0]  - flux_after_conv[2])
            return [conv_flux_val,
                    conv_flux_val + max_upper_error,
                    conv_flux_val - max_lower_error]

        else:
            newName = insert_cflux_before(modelName, component)

    m = Model(newName)
    enterParameters(parameters, {"Emin": Emin, "Emax": Emax})
    freezeNorm()
    fitModel()

    fluxVals = findFlux()

    if fluxVals:
        print(f"    {component}: {fluxVals[0]:.4e} ×10⁻⁹ erg/cm²/s")

    return fluxVals


def writeParsAfterFlux(line_list):
    for comp in AllModels(1).componentNames:
        compObj = getattr(AllModels(1), comp)
        for par in compObj.parameterNames:
            parObj = getattr(compObj, par)
            fullName = comp + "." + par

            line_list.append(fullName + "     " + str(parObj.values[0]) + "\n")
    
    line_list.append("\n") 

def write_lines_to_file(file_name, line_list):
    with open(file_name, "a") as file:
        for line in line_list:
            file.write(line)

#===================================================================================================================
try:
    energyLimits = energyFilter.split(" ")
    Emin = energyLimits[0]
    Emax = energyLimits[1]
except Exception as e:
    print(f"Exception occured while reading 'energyLimits' variable due to incorrect format: {e}")
    quit()

commonDirectory = outputDir + "/commonFiles"   # ~/NICER/analysis/commonFiles

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
    print("\nCould not find any valid observation path, as given in the obs.txt file.")
    quit()

iterationMax = 0
searchedObservations = []
if Path(commonDirectory + "/processed_obs.txt").exists() == False:
    print("\nCould not find 'processed_obs.txt' file under the 'commonFiles' directory.")
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
    print("\nCould not find the searched observation paths in 'processed_obs.txt', most likely due to having low exposure.")
    quit()

if chatterOn == False:
    Xset.chatter = 0

# Start calculating fluxes for each valid observation
for path, obsid, expo in searchedObservations:
    print("====================================================================")
    print("Calculating fluxes for observation:", obsid, "\n")

    outObsDir = path
    try:
        os.chdir(outObsDir)
    except Exception as e:
        print(f"Exception occured while changing directory to {outObsDir}: {e}")
        continue

    version = 0
    try:
        with open(outObsDir + "/results/version_counter.txt") as version_file:
            all_lines = version_file.readlines()
            version = int(all_lines[1].strip("\n")) - 1
    except Exception as e:
        print(f"Exception occured while reading file {outObsDir}/results/version_counter.txt: {e}")
        continue
    
    output_save_name = custom_name
    if output_save_name == "":
        output_save_name = model_pipeline_name

    outObsDir = outObsDir + "/results/" + output_save_name +"_" + str(version)
    allFiles = os.listdir(outObsDir)

    # Find the data file and the best fitting model file for the current observation
    missingFiles = True
    foundModfile = False
    foundDatafile = False
    for file in allFiles:
        if "best_" in file:
            modFile = outObsDir + "/" + file
            foundModfile = True
        elif "data_" in file:
            dataFile = outObsDir + "/" + file
            foundDatafile = True

        if foundDatafile and foundModfile:
            # All necessary files have been found
            missingFiles = False
            break
    
    fit_file_loc = outObsDir + "/" + resultsFile
    fit_file_lines = []

    try:
        with open(fit_file_loc) as fit_file:
            lines = fit_file.readlines()

            for line in lines:
                if "Fluxes of model components" in line:
                    fit_file_lines = fit_file_lines[:-1]    # To exclude the "======" line before
                    break
                else:
                    fit_file_lines.append(line)
    except Exception as e:
        print(f"Exception occured while opening {fit_file_loc} for observation {obsid}: {e}")
        continue

    # Check if there are any missing files
    if missingFiles:
        print("ERROR: Necessary files for calculating fluxes are missing for the observation: " + obsid)
        fit_file_lines.append("\nERROR: Necessary files for calculating fluxes are missing for the observation: " + obsid + "\n")

        if foundModfile == False:
            print("->Missing model file")
            fit_file_lines.append("->Missing model file\n")
        if foundDatafile == False:
            print("->Missing data file")
            fit_file_lines.append("->Missing data file\n")
        
        write_lines_to_file(fit_file_loc, fit_file_lines)
        continue
    
    print("All the files required for calculating fluxes are found. Please check if the correct files are in use.")
    print("Model file: ", modFile)
    print("Data file: ", dataFile, "\n")

    try:
        Xset.restore(dataFile)
        Xset.restore(modFile)
    except Exception as e:
        print(f"Exception occured while loading data and model files to PyXspec: {e}")
        continue

    Fit.query = "yes"

    parameters = {}
    updateParameters(parameters)

    # Open the parameter file and extract all non_flux lines
    all_lines_file = {}

    try:
        par_file = open("parameters_bestmodel.txt", "r")
    except Exception as e:
        print(f"Exception occured while opening parameters_bestmodel.txt file for observation {obsid}: {e}")
        continue

    all_lines = par_file.readlines()
    par_file.close()

    for line in all_lines:
        if "flux" not in line:
            all_lines_file[line] = 1

    fit_file_lines.append("\n===========================================================\n")
    fit_file_lines.append("Fluxes of model components (in 10^-9 ergs/cm^2/s) (90% confidence intervals)\n\n")
    modelName = AllModels(1).expression.replace(" ", "")

    for fluxModel in fluxes_to_be_calculated:
        if (fluxModel != "unabsorbed" and fluxModel != "absorbed") and (fluxModel not in modelName):
            print("\nWARNING: Model '" + fluxModel + "' does not exist in current model expression.")
            print(f"Flux calculation will be skipped for '{fluxModel}'..\n")
            continue

        print("Calculating flux for: " + fluxModel)
        flux = calculateFlux(fluxModel, modelName, parameters)
        if not flux:
            print("Could not calculate flux for: " + fluxModel)
            continue

        # Write flux data to 
        fit_file_lines.append(energyFilter +" keV " + AllModels(1).expression + "\nFlux: " + listToStr(flux) + "\n")

        writeParsAfterFlux(fit_file_lines)
        
        # Add new flux line to the all_lines_file
        all_lines_file[fluxModel +"_flux " + listToStr(flux)+ " (10^-9_ergs_cm^-2_s^-1)\n"] = 1

    hardness_ratio = 0.0
    soft_rate, soft_err, hard_rate, hard_err = 0.0, 0.0, 0.0, 0.0
    try:
        Xset.restore(dataFile)

        def get_rate(emin, emax):
            AllData.notice("**")
            AllData.ignore(f"**-{emin} {emax}-**")
            spec = AllData(1)
            return spec.rate[0], spec.rate[1]

        soft_rate, soft_err = get_rate(1.0, 3.0)
        hard_rate, hard_err = get_rate(3.0, 10.0)
        hardness_ratio = hard_rate / soft_rate if soft_rate != 0 else 0.0

        print("\nCount Rates:")
        print(f"Soft (1-3 keV): {soft_rate:.5f} ct/s")
        print(f"Hard (3-10 keV): {hard_rate:.5f} ct/s")
        print(f"H/S: {hardness_ratio:.5f}")

        fit_file_lines.append("\n===========================================================\n")
        fit_file_lines.append("Count Rates & Hardness Ratio\n\n")
        fit_file_lines.append(f"Soft (1-3 keV): {soft_rate:.5f} +/- {soft_err:.5f} ct/s\n")
        fit_file_lines.append(f"Hard (3-10 keV): {hard_rate:.5f} +/- {hard_err:.5f} ct/s\n")
        fit_file_lines.append(f"H/S: {hardness_ratio:.5f}\n")

        # Format: name value errlow errhigh unit  (errlow/errhigh are absolute bounds, nicer_plot.py subtracts the value itself)
        all_lines_file[f"soft_rate {soft_rate:.5f} {soft_rate - soft_err:.5f} {soft_rate + soft_err:.5f} ct/s\n"] = 1
        all_lines_file[f"hard_rate {hard_rate:.5f} {hard_rate - hard_err:.5f} {hard_rate + hard_err:.5f} ct/s\n"] = 1
        all_lines_file[f"hardness_ratio {hardness_ratio:.5f} {hardness_ratio:.5f} {hardness_ratio:.5f} ratio\n"] = 1

    except Exception as e:
        print(f"Hardness hesaplanırken hata: {e}")
        fit_file_lines.append(f"\nHardness ratio could not be calculated: {e}\n")

    # Write flux values to parameter file
    par_file = open("parameters_bestmodel.txt", "w")
    for line in all_lines_file.keys():
        par_file.write(line)
    par_file.close()

    write_lines_to_file(fit_file_loc, fit_file_lines)

    AllModels.clear()
    AllData.clear()

# This file is created after importing variables from another python file
if Path(scriptDir + "/__pycache__").exists():
    os.system("rm -rf "+scriptDir+"/__pycache__")