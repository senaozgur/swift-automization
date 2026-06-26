# Swift XRT spectral fitting pipeline — adapted from nicer_fit.py
# Original authors: Batuhan Bahçeci (batuhan.bahceci@sabanciuniv.edu)
# Swift XRT adaptation: see parameter.py for configuration

from parameterfit import *
import operator

operator_mapping = {
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne
}

print("==============================================================================")
print("\t\t\tRunning " + fit_script_name + "\n")

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

#========================================================== Input Checks ===========================================================
# Input check for outputDir
if Path(outputDir).exists() == False:
    print("Directory defined by outputDir could not be found. Terminating the script...")
    quit()

# Create commonFiles directory if it does not exist (Swift pipeline creates it automatically)
if Path(outputDir + "/commonFiles").exists() == False:
    print(f"Directory {outputDir}/commonFiles not found. Creating it now...")
    os.makedirs(outputDir + "/commonFiles", exist_ok=True)

# Always regenerate processed_obs.txt from the current swift.txt.
# If it already existed from a previous run (e.g. with fewer obsids), it would
# otherwise silently stay stale and cause later obsids to be skipped entirely.
processed_obs_path = outputDir + "/commonFiles/processed_obs.txt"
print(f"Building processed_obs.txt from {inputTxtFile}...")
try:
    with open(scriptDir + "/" + inputTxtFile, "r") as f:
        obs_lines = f.readlines()
    with open(processed_obs_path, "w") as pf:
        for obs_line in obs_lines:
            obs_line = obs_line.strip()
            if obs_line == "":
                continue
            if not Path(obs_line).exists():
                print(f"  WARNING: Path not found, skipping -> {obs_line}")
                continue
            rev = obs_line[::-1]
            if obs_line[-1] != "/":
                slash_idx = rev.find("/")
                oid = rev[:slash_idx][::-1]
                opath = obs_line
            else:
                slash_idx = rev[1:].find("/")
                oid = rev[1:slash_idx+1][::-1]
                opath = obs_line.rstrip("/")
            # Try to read exposure from grouped spectrum header
            try:
                import glob as _glob
                grp_candidates = _glob.glob(opath + "/sw" + oid + "*grp*.pha") + \
                                 _glob.glob(opath + "/sw" + oid + "*_grp.pha") + \
                                 _glob.glob(opath + "/*grp*.pha")
                if grp_candidates:
                    _hdu = fits.open(grp_candidates[0])
                    _expo = _hdu[1].header.get("EXPOSURE", 0)
                    _hdu.close()
                else:
                    _expo = 0
            except Exception:
                _expo = 0
            print(f"  Added obsid: {oid}  exposure: {_expo}  path: {opath}")
            pf.write(opath + " " + oid + " " + str(_expo) + "\n")
    print("processed_obs.txt has been (re)built.\n")
except Exception as e:
    print(f"Exception while building processed_obs.txt: {e}")
    quit()

# Input check for restartAlways
if isinstance(restartAlways, bool) == False:
    while True:
        print("\nThe 'restartAlways' variable is not of type boolean.")
        restartAlways = input("Please enter a boolean value for 'restartAlways' (True/False): ")

        if restartAlways == "True" or restartAlways == "False":
            restartAlways = (restartAlways == "True")
            break

# Input check for restartOnce
if isinstance(restartOnce, bool) == False:
    while True:
        print("\nThe 'restartOnce' variable is not of type boolean.")
        restartOnce = input("Please enter a boolean value for 'restartOnce' (True/False): ")

        if restartOnce == "True" or restartOnce == "False":
            restartOnce = (restartOnce == "True")
            break

# Input check for ftestSignificance
if isinstance(ftestSignificance, float) == False or not (0 < ftestSignificance < 1):
    while True:
        try:
            ftestSignificance = float(ftestSignificance)
            if ftestSignificance <= 0 or ftestSignificance >= 1:
                raise Exception()
            else:
                break
        except:
            print("\nThe 'fTestSignificance' variable must be a float number between 0 and 1.")
            ftestSignificance = input("Please enter a float number between 0 and 1 for 'ftestSignificance' (0 < x < 1): ")

# Input check for fix_sample_size
if str(fix_sample_size).isnumeric() == False or int(fix_sample_size) <= 0:
    while True:
        print("\nEither the 'fix_sample_size' variable is not of type integer, or it is smaller or equal to 0.")
        fix_sample_size = input("Please enter a positive integer value for 'fix_sample_size' (x > 0): ")

        if fix_sample_size.isnumeric() and int(fix_sample_size) > 0:
            fix_sample_size = int(fix_sample_size)
            break

# Input check for fix_parameters_after_sampling
if isinstance(fix_parameters_after_sampling, bool) == False:
    while True:
        print("\nThe 'fix_parameters_after_sampling' variable is not of type boolean.")
        fix_parameters_after_sampling = input("Please enter a boolean value for 'fix_parameters_after_sampling' (True/False): ")

        if fix_parameters_after_sampling == "True" or fix_parameters_after_sampling == "False":
            fix_parameters_after_sampling = (fix_parameters_after_sampling == "True")
            break

# Input check for errorCalculations
if isinstance(errorCalculations, bool) == False:
    while True:
        print("\nThe 'errorCalculations' variable is not of type boolean.")
        errorCalculations = input("Please enter a boolean value for 'errorCalculations' (True/False): ")

        if errorCalculations == "True" or errorCalculations == "False":
            errorCalculations = (errorCalculations == "True")
            break

# Input check for model_pipeline_name
if model_pipeline_name == "":
    print("model_pipeline_name is not provided in parameter.py")
    print("Please provide a valid model name that is defined in models.txt and try running the script again.")
    quit()
#===================================================================================================================================

#===================================================================================================================================
# Functions
def shakefit(bestModelList, resultsFile):
    reduced_chi = Fit.statistic / Fit.dof
    if reduced_chi > 10:
        print("==============================================================")
        print(f"WARNING: Reduced-chi squared is extremely large: {reduced_chi}.")
        print("==============================================================")
    elif reduced_chi > 2:
        print("==============================================================")
        print(f"WARNING: Reduced-chi squared is significantly large: {reduced_chi}.")
        print("==============================================================")
    # Shakefit will only be run for these parameters
    parametersToCalculateError = []

    try:
        for key in parametersForShakefit.keys():
            model = key[:key.find(".")]
            parameter = key[key.find(".") + 1:]
            if model in AllModels(1).expression:
                compObj = getattr(AllModels(1), model)
                parObj = getattr(compObj, parameter)
                index = parObj.index

                parametersToCalculateError.append(index)
            else:
                pass
    except Exception as e:
        print(f"Exception occured while reading the contents of 'parametersForShakefit' variable in parameter.py: {e}")
        quit()

    if checkPowerlawErrorAndFreeze:
        # Before proceeding with shakefit, check if powerlaw is in model expression. If so, check whether xspec error for photon index is bigger than 1 or not.
        # If bigger than 1, fix photon index to 1.7 (may change in future) and only then continue with shakefit
        if "powerlaw" in AllModels(1).expression:
            with open(xspec_output_file, "r") as testfile:
                lines = testfile.readlines()[-35:]
            
            print("\nChecking powerlaw xspec error value...\n")
            retrievePhotonIndex = False
            for line in lines:
                line = line.strip("\n")
                if " par  comp" in line:
                    retrievePhotonIndex = True
                
                if retrievePhotonIndex:
                    if "powerlaw" not in line:
                        pass
                    else:
                        words = line.split(" ")

                        # Remove all empty elements from the list
                        while True:
                            try:
                                emptyIndex = words.index("")
                                words.pop(emptyIndex)
                            except:
                                break

                        errorValue = words[-1]  # Xspec error value
                        errorValue = errorValue.strip("\n")
                        
                        if errorValue == "frozen":
                            break

                        errorValue = float(errorValue)
                        if (1 >= errorValue > 0) == False:
                            AllModels(1).powerlaw.PhoIndex.values = str(powerlawIndexToFreezeAt) + " -1"
                            resultsFile.write("\nWARNING: Powerlaw photon index is frozen at " + str(AllModels(1).powerlaw.PhoIndex.values[0]) + " for having large xspec error: " + str(errorValue)+ "\n")
                            print("Powerlaw xspec error value is: " + str(errorValue))
                            print("\nWARNING: Powerlaw photon index is frozen at " + str(AllModels(1).powerlaw.PhoIndex.values[0]) + " for having xspec error bigger than 1: "+str(errorValue)+"\n")
                        
                        break

    resultsFile.write("========== Proceeding with shakefit error calculations ==========\n")
    paramNum = AllModels(1).nParameters
    rerunShakefit = False
    for k in range(2):
        if k == 1 and rerunShakefit == False:
            break

        print("Performing shakefit error calculations for the model: " + AllModels(1).expression + ", obsid:" + str(obsid)+ ", shakefit number: " + str(k+1)+"\n")
        fitModel(bestModelList)
        updateParameters(bestModel)

        for i in range(1, paramNum+1):
            parDelta = AllModels(1)(i).values[1]

            if parDelta < 0:  # Check for frozen parameters
                continue
            
            if i not in parametersToCalculateError:
                continue

            continueError = True
            delChi = 2.706
            counter = 0
            while continueError and counter < 100:
                counter += 1
                Fit.error("stopat 10 0.1 maximum 1000 " + str(delChi) + " " + str(i))
                errorResult = AllModels(1)(i).error
                errorString = errorResult[2]

                if errorString[3] == "T" or errorString[4] == "T":
                    # Hit lower/upper limits, stop the error process for the current model parameter
                    continueError = False

                if errorString[0] == "F":
                    # Could not find a new minimum
                    continueError = False
                elif errorString[1] == "T":
                    # Non-monotonicity detected
                    delChi += 2

                if continueError == False:
                    # Save error calculation results to the log file
                    # Save error calculation results to the log file
                    parName = AllModels(1)(i).name
                    parValue = AllModels(1)(i).values[0]
                    parValue = AllModels(1)(i).values[0]
                    errorTuple = "(" + listToStr(errorResult) + ")"
                    resultsFile.write("Par " + str(i) + ": " + parName + " " + errorTuple+"\n")
            
        # Check if any parameter value has gotten outside their initially calculated confidence interval
        rerunShakefit = False
        for m in range(1, AllModels(1).nParameters+1):
            parValue = AllModels(1)(m).values[0]
            errorString = AllModels(1)(m).error
            if errorString[0] != 0 and parValue < errorString[0]:
                rerunShakefit = True
                print("Some parameters are out of their previously calculated error boundaries. Rerunning shakefit...\n")
                break
            elif errorString[1] != 0 and parValue > errorString[1]:
                print("Some parameters are out of their previously calculated error boundaries. Rerunning shakefit...\n")
                rerunShakefit = True
                break

    resultsFile.write("=================================================================\n\n")
    updateParameters(bestModel)

def listToStr(array):
    result = ""
    for char in array:
        result += str(char) + " "
    result = result[:-1]
    return result

def getParsFromList(modList, ignoreList = []):
    parameterList = modList[0]

    parNum = AllModels(1).nParameters

    # Take the available parameters from the stemList
    components = AllModels(1).componentNames
    for comp in components:
        if comp in ignoreList:
            pass
        else:
            compObj = getattr(AllModels(1), comp)
            parameters = compObj.parameterNames
            for par in parameters:
                parObj = getattr(compObj, par)
                indx = parObj.index
                fullName = comp + "." + par

                if fullName in parameterList:
                    AllModels(1)(indx).values = parameterList[fullName]

def tie_diskbb_Tin_to_nthComp_kTbb():
    
    #Tie nthComp.kT_bb to diskbb.Tin if both components exist.
    
    
    try:
        disk = AllModels(1).diskbb
        nth  = AllModels(1).nthComp
    except Exception:
        
        return

    try:
        tin_par  = disk.Tin       # diskbb.Tin
        ktbb_par = nth.kT_bb      # nthComp.kT_bb
    except Exception:
        
        return

    # Tie kT_bb to Tin via Tin's parameter index
    ktbb_par.link = str(tin_par.index)

    
    print(
        f"[TIE] nthComp.kT_bb (par {ktbb_par.index}) "
        f"tied to diskbb.Tin (par {tin_par.index}); "
        f"Tin = {tin_par.values[0]}"
    )
 
def fitModel(bestModelList): 
    print("\nFitting the model: " + AllModels(1).expression + "\n")
    # Enforce Tin = kT_bb before every fit
    tie_diskbb_Tin_to_nthComp_kTbb()

    Fit.nIterations = 100
    Fit.delta = 0.01
    Fit.renorm()
    Fit.perform()
    updateParameters(bestModelList)

def updateParameters(modList):
    modelDict = modList[0]
    modelStats = modList[1]

    # Save the parameters loaded in the xspec model to lists
    components = AllModels(1).componentNames
    for comp in components:
        compObj = getattr(AllModels(1), comp)
        parameters = compObj.parameterNames
        for par in parameters:
            parObj = getattr(compObj, par)
            fullName = comp + "." + par
            modelDict[fullName] = parObj.values
    
    modelStats["chi"] = Fit.statistic
    modelStats["dof"] = Fit.dof
    modelStats["nullhyp"] = Fit.nullhyp

def saveModel(fileName, location = "default"):
    # This function saves the model file (.xcm) under spesified location. If no location has been given, 
    # the model file will be saved under default (outObsDir) directory
    
    if location == "default":
        xcmPath = Path(outObsDir + "/" + fileName)
        if xcmPath.exists():
            subprocess.run(["rm", outObsDir + "/" + fileName])

        Xset.save(outObsDir + "/" + fileName, "m")
    else:
        xcmPath = Path(location + "/" + fileName)
        if xcmPath.exists():
            subprocess.run(["rm", location + "/" + fileName])

        Xset.save(location + "/" + fileName, "m")

def saveData(location = "default"):
    # Similar to saveModel function, this function saves the data instead of model in an xcm file
    if location == "default":
        fileName = outObsDir + "/" + "data_" + obsid + ".xcm"
        filePath = Path(fileName)
        if filePath.exists():
            subprocess.run(["rm", fileName])
        Xset.save(fileName, "f")

    else:
        fileName = location + "/" + "data_" + obsid + ".xcm"
        filePath = Path(fileName)
        if filePath.exists():
            subprocess.run(["rm", fileName])
        Xset.save(fileName, "f")

def writeBestFittingModel(resultsFile):
    # This function writes the current model in AllModels container in a log file, assuming that all the models
    # have been compared and the last remaining model is the best fitting model.
    try:
        resultsFile.write("====================== Best Fitting Model ======================\n")
        resultsFile.write("Model Name: " + AllModels(1).expression + "\n")
        
        resultsFile.write("Fit results:\n")
        fitString = "Null Probability: " + str(Fit.nullhyp) +", Chi-squared: " + str(Fit.statistic) + ", Dof: " + str(Fit.dof) + "\n"
        resultsFile.write(fitString)

        parameterString = ""
        for comp in AllModels(1).componentNames:
            compObj = getattr(AllModels(1), comp)
            for par in compObj.parameterNames:
                parVal = getattr(compObj, par).values
                fullName = comp + "." + par
                parameterString += fullName + ":\t" + str(parVal) + "\n"

        resultsFile.write("Parameters: \n" + parameterString)
        resultsFile.write("=================================================================\n")
    except Exception as e:
        print(f"Exception occured while writing the best fitting model to {resultsFile}: {e}")

def extractModFileName():
    fileName = "model_"
    comps = sorted(AllModels(1).componentNames)
    compPart = ""
    for comp in comps:
        compPart += comp[:3]
    fileName += compPart + ".xcm"
    
    return fileName

def removeComp(compName, compNum, modelList):   # compNum is the n'th occurence of a spesific model, it is not the component number in general
    iteration = 1
    targetModelCounter = 1
    deletedCompIndex = 999
    modelFirstEncounter = True
    tempDict = {}
    modelName = AllModels(1).expression.replace(" ", "")

    compCount = wordCounter(modelName, compName)
    if compCount < compNum or compNum <= 0:
        print("\nERROR: There are (" + str(compCount) + ") " + compName + " in the current model expression.")
        print("Given component number should be in between 1 <= x <= " + str(compCount) + " (Given input: " + str(compNum) + ")\n")
        quit()

    # This loop iterates over each component, and assigns the model parameters with their values only if they do not belong to the targeted model
    for comp in AllModels(1).componentNames:
        if  targetModelCounter != compNum or (compName not in comp):
            compObj = getattr(AllModels(1), comp)
            if compName in comp:
                if modelFirstEncounter == True and "_" in comp:
                    for par in compObj.parameterNames:
                        parObj = getattr(compObj, par)
                        newComp = comp[:comp.find("_")]
                        fullName = newComp + "." + par
                        tempDict[fullName] = parObj.values
                elif "_" in comp and deletedCompIndex <= int(comp[comp.find("_") + 1:]):
                    for par in compObj.parameterNames:
                        parObj = getattr(compObj, par)
                        newComp = comp[:comp.find("_") + 1] + str(int(comp[comp.find("_") + 1:]) - 1)
                        fullName = newComp + "." + par
                        tempDict[fullName] = parObj.values
                else:
                    for par in compObj.parameterNames:
                        parObj = getattr(compObj, par)
                        fullName = comp + "." + par
                        tempDict[fullName] = parObj.values

                modelFirstEncounter = False
                targetModelCounter += 1
            else:
                for par in compObj.parameterNames:
                    parObj = getattr(compObj, par)
                    if "_" in comp and deletedCompIndex <= int(comp[comp.find("_") + 1:]):
                        newComp = comp[:comp.find("_") + 1] + str(int(comp[comp.find("_") + 1:]) - 1)
                        fullName = newComp + "." + par
                        tempDict[fullName] = parObj.values
                    else:
                        fullName = comp + "." + par
                        tempDict[fullName] = parObj.values
        else:
            targetModelCounter += 1
            deletedCompIndex = iteration
        
        iteration += 1

    # Look for the position of the targeted model component, then remove it from the model expression
    try: 
        test = modelName.index("+" + compName + "+")
        modelName = modelName.replace("+" + compName + "+", "+", 1) 
    except:
        try: 
            test = modelName.index("+" + compName)
            modelName = modelName.replace("+" + compName, "", 1) 
        except: 
            try: 
                test = modelName.index(compName + "+")
                modelName = modelName.replace(compName + "+", "", 1)
            except: 
                try: 
                    test = modelName.index(compName + "*")
                    modelName = modelName.replace(compName + "*", "", 1)
                except:
                    try:
                        test = modelName.index("*" + compName)
                        modelName = modelName.replace("*" + compName, "", 1)
                    except:
                        try:
                            test = modelName.index("(" + compName + ")")
                            modelName = modelName.replace("(" + compName +")", "", 1)
                        except:
                            modelName = modelName.replace(compName, "", 1)

    modelList[0] = tempDict
    m = Model(modelName)
    getParsFromList(modelList)

def addComp(compName, targetComp ,before_after, calcChar, modelList, encapsulate=False):
    modelName = " " + AllModels(1).expression + " "
    if calcChar != "*" and calcChar != "+":
        print("\nIncorrect character for model expression. Terminating the script...\n")
        quit()

    modelName = modelName[1:-1]
    if before_after == "before":
        if encapsulate:
            targetIdx = modelName.find(targetComp)
            insertionText =  "(" + compName + calcChar + targetComp + ")"
            addedCompIndex = AllModels(1).componentNames.index(targetComp) + 1

            newModelName = modelName[:targetIdx] + insertionText + modelName[targetIdx + len(targetComp):]
        else:
            targetIdx = modelName.find(targetComp)
            insertionText = compName + calcChar
            addedCompIndex = AllModels(1).componentNames.index(targetComp) + 1

            newModelName = modelName[:targetIdx] + insertionText + modelName[targetIdx:]
    elif before_after == "after":
        if encapsulate:
            targetIdx = modelName.find(targetComp)
            insertionText = "(" + targetComp + calcChar + compName + ")"
            addedCompIndex = AllModels(1).componentNames.index(targetComp) + 2

            newModelName = modelName[:targetIdx] + insertionText + modelName[targetIdx + len(targetComp):]
        else:
            targetIdx = modelName.find(targetComp)
            insertionText = calcChar + compName
            addedCompIndex = AllModels(1).componentNames.index(targetComp) + 2

            newModelName = modelName[:targetIdx + len(targetComp)] + insertionText + modelName[targetIdx + len(targetComp):]
    else:
        print("\nIncorrect entry for the placement of new component around the target component. Terminating the script...\n")
        quit()
    
    alter_list_add(compName, addedCompIndex, modelList)
    m = Model(newModelName)
    getParsFromList(modelList)
    modelList[0] = {}
    updateParameters(modelList)

def alter_list_add(compName, addedIdx, bestModelList):
    modelKeys = list(bestModelList[0].keys())[::-1]
    modelValues = list(bestModelList[0].values())[::-1]
    
    for i in range(len(modelKeys)):
        compPart = modelKeys[i][:modelKeys[i].find(".")]
        rest = modelKeys[i][modelKeys[i].find("."):]
        if "_" in compPart:
            compNum = compPart[compPart.find("_") + 1 :]
            if int(compNum) > addedIdx:
                newKey = compPart.replace(compNum, str(int(compNum)+1)) + rest
                bestModelList[0].pop(modelKeys[i])
                bestModelList[0][newKey] = modelValues[i]
        elif compName in modelKeys[i]:
            compIdx = AllModels(1).componentNames.index(compName) + 1
            if compIdx >= addedIdx:
                newKey = modelKeys[i][: modelKeys[i].find(".")] + "_" + str(compIdx + 1) + modelKeys[i][modelKeys[i].find(".") :]
                bestModelList[0].pop(modelKeys[i])
                bestModelList[0][newKey] = modelValues[i]

def wordCounter(source, word):
    start = 0
    count = 0
    while True:
        idx = source.find(word, start)

        if idx == -1:
            break
        else:
            start = idx + 1
            count += 1
    
    return count

def assignParameters(compName, nthOccurence, parameterTuple):
    startAssign = False
    occurence = 0
    for comp in AllModels(1).componentNames:
        compObj = getattr(AllModels(1), comp)
        if compName in comp:
            occurence += 1
            if occurence == nthOccurence:
                parameterName = parameterTuple[0]
                value = parameterTuple[1]
                parObject = getattr(compObj, parameterName)
                parObject.values = value

def transferToNewList(sourceList):
    newList = []
    newParDict = {}
    newStatDict = {}

    sourceParDict = sourceList[0]
    keys = list(sourceParDict.keys())
    values = list(sourceParDict.values())
    for i in range(len(keys)):
        newParDict[keys[i]] = values[i]
    newList.append(newParDict)

    sourceStatDict = sourceList[1]
    keys = list(sourceStatDict.keys())
    values = list(sourceStatDict.values())
    for i in range(len(keys)):
        newStatDict[keys[i]] = values[i]
    newList.append(newStatDict)
    
    return newList

def performFtest(nullModelList, altModelList, logFile, infoTxt = ""):
    newChi = altModelList[1]["chi"]
    newDof = altModelList[1]["dof"]
    oldChi = nullModelList[1]["chi"]
    oldDof = nullModelList[1]["dof"]

    pValue = Fit.ftest(newChi, newDof, oldChi, oldDof)
    
    logFile.write("Performing f-test: ")

    if infoTxt != "":
        logFile.write(infoTxt)

    logFile.write("F-test significance: "+ str(ftestSignificance)+", p-value: " + str(pValue)+"\n\n")
    return pValue

def calculateGaussEqw(logFile):
    eqwList = []
    counter = 0
    for comp in AllModels(1).componentNames:
        counter += 1
        if "gaussian" in comp:
            compName = comp
            compObj = getattr(AllModels(1), comp)
            energyVal = compObj.LineE.values[0]

            try:
                AllModels.eqwidth(counter, err=True, number=1000, level=90)
                eqwList.append("Equivalent width: " + str(listToStr(AllData(1).eqwidth)) + " (" + str(format(energyVal, ".2f")) + " keV gauss)\n")
            except Exception as e:
                eqwList.append("Calculating eqw failed for component: " + comp + "\n")
                print(f"Exception: {e}")
    
    if eqwList != []:
        logFile.write("Gauss equivalent widths: (90% confidence intervals) \n")
        for each in eqwList:
            logFile.write(each)

def fixAllParameters(fixedValues):
    for key, val in fixedValues.items():
        temp = key.split(".")
        compName = temp[0]
        parName = temp[1]
        compObj = getattr(AllModels(1), compName)
        parObj = getattr(compObj, parName)

        parObj.values = fixedValues[key]

def filterOutliers(dataset):
    # Here I use Z-Score algorithm to filter out the outliers of the dataset
    # I saw that 3 sigma was commonly used for filtering outliers in this method, but since I will be filtering out nH values and we do not expect nH values to change
    # that drastically, I narrowed down the interval and set the z-score as 2 sigma.
    mean = np.mean(dataset)
    std_dev = np.std(dataset)
    threshold = 2

    filtered_dataset = []
    for i in dataset:
        if std_dev == 0 or (abs(i - mean) / std_dev <= threshold):
            filtered_dataset.append(i)
    
    return filtered_dataset

def closeAllFiles():
    print("\nClosing all files..")
    logFile.close()
    Xset.closeLog()
    AllModels.clear()
    AllData.clear()

def checkResults(bestModelList):
    for i in range(1, AllModels(1).nParameters+1):
        print(AllModels(1)(i).name, AllModels(1)(i).values[0])
    print("================")
    for key, val in bestModelList[0].items():
        print(key,val)

def clearList(listInput):
    for i in range(len(listInput)):
        listInput.pop()

def constructList(listToConstruct, sourceList):
    clearList(listToConstruct)
    listToConstruct.append({})
    listToConstruct.append({})
    for i in range(len(sourceList)):
        for key, val in sourceList[i].items():
            listToConstruct[i][key] = val

def loadModel(modelName):
    print("Loading model: " + modelName)
    m = Model(modelName)

def assignTxtParameters(parameterList):
    print()
    for tuple in parameterList:
        compName = tuple[0]
        compNum = tuple[1]
        parTuple = tuple[2]

        print("Assigning new values for parameter: " + compName + "." + parTuple[0])
        assignParameters(compName, compNum, parTuple)

def searchPremodel(bestModelList, path = ""):
    modelfile = extractModFileName()
    foundFile = False
    if path == "":
        path = outputDir + "/commonFiles"
        print("\nLooking for a model file '" + modelfile + "' under '" + path + "'")
        if Path(path + "/" + modelfile).exists():
            foundFile = True
            print("Found the spesified model file.")
            print("Extracting all the model parameters..")
            updateParameters(bestModelList)

            Xset.restore(path + "/" + modelfile)
    else:
        print("\nLooking for a model file '" + modelfile + "' under '" + path + "'")
        if Path(path + "/" + modelfile).exists():
            foundFile = True
            print("Found the spesified model file.")
            print("Extracting all the model parameters..")
            updateParameters(bestModelList)

            Xset.restore(path + "/" + modelfile)
    
    if foundFile == False:
        print("Could not find the target model file under '" + path + "'")

def saveCommand(saveType):
    print("Saving requested xcm file as type: " + saveType)
    modelName = extractModFileName()
    if saveType == "model":
        saveModel(modelName)
        saveModel(modelName, commonDirectory)
    elif saveType == "data":
        saveData()
    else:
        saveModel(modelName)
        saveModel(modelName, commonDirectory)
        saveData()

def ftestOptions(option, bestModelList, nullhypList, logfile, lastAddedModel, lastAddedModelNumber, orderSuffix, infoTxt = ""):
    if option == "nullhyp":
        print("\nSaved null hypothesis statistics\n")
        constructList(nullhypList, bestModelList)
    else:
        print("\nPerforming ftest: " + infoTxt + "\n")
        pvalue = performFtest(nullhypList, bestModelList, logFile, infoTxt)

        if abs(pvalue) >= ftestSignificance:
            print("\nFTEST: Choosing null hypothesis model")
            print("F-test significance: " + str(ftestSignificance) + ", p-value: " + str(pvalue) + "\n")

            print("Removing the model: " + str(lastAddedModelNumber) + orderSuffix + " " + lastAddedModel)
            removeComp(lastAddedModel, lastAddedModelNumber, bestModelList)
            fitModel(bestModelList)
            updateParameters(bestModelList)
        else:
            print("\nFTEST: Keeping the alternative model")
            print("F-test significance: " + str(ftestSignificance) + ", p-value: " + str(pvalue) + "\n")

def calculateComponentOrder(compName, targetName):
    modelExpression = AllModels(1).expression.replace(" ", "")
    targetIdx = modelExpression.find(targetName)
    newExpression = modelExpression[:targetIdx]
    compCount = wordCounter(newExpression, compName)

    return compCount + 1

def parseTxt(source, bestModelList, nullhypList, logFile, enableFixing):
    sourcefile = open(source, "r")
    lines = sourcefile.readlines()

    if_stack = []

    currentModel = ""
    lastAddedModel = ""
    lastAddedModelNumber = 0
    orderSuffix = ""

    lineCount = 0
    for line in lines:
        lineCount += 1
        #===============    Preprocess the line    =======================
        line = line.strip("\n\t ")

        if len(line) == 0:
            continue
        
        if line[0] == "#":
            continue
        
        line = line.split(" ")

        # Remove the spaces in the line
        idx_counter = 0
        while idx_counter < len(line):
            if (line[idx_counter] == ""):
                line.pop(idx_counter)
            else:
                idx_counter += 1
        #=======================================
                
        try:
            if line[0].lower() == "model":
                if line[1] == model_pipeline_name:
                    # Found the target model name
                    if currentModel == "":
                        # First time encountering that model name, start processing the model
                        currentModel = model_pipeline_name
                        continue
                    else:
                        # A model pipeline has been defined before, but the current line tries to define another pipeline
                        raise Exception("Trying to process another model pipeline, probably due to not using ENDMODEL keyword.")

        except Exception as e:
            print(f"Exception occured: {e}")
            print("\nERROR: Invalid implementation for a pipeline name in models.txt -> Line: " + str(lineCount))
            quit()

        if currentModel != model_pipeline_name:
            # If the current pipeline does not match
            continue
            
        if line[0].lower() == "endmodel":
            # End of model pipeline
            return True
        
        # Enter if we are currently in an if statement
        if len(if_stack) != 0:
            if line[0] == "if":
                # There is a nested if statement coming up. If current if statement is evaluated to False, no need to check inner if statements.
                # If that is the case, evaluate all upcoming if statements to False. If not, pass this step and move onto if evaluation block.
                if if_stack[-1] == False:
                    if_stack.append(False)
                    continue

            elif line[0] == "endif":
                # End of inner if statement is reached.
                if_stack.pop()
                continue

            else:
                # It is another command. If the innermost if statement is evaluated to False, skip it.
                if if_stack[-1] == False:
                    continue

        try:
            if line[0] == "load":
                modelName = ""
                for par in line[1:]:
                    modelName += par

                loadModel(modelName)
                continue
        except Exception as e:
            print(f"Exception occured while loading model: {e}")
            print("\nERROR: Invalid implementation of 'load' command in models.txt -> Line: " + str(lineCount))
            quit()
          
        try:
            if line[0] == "assign":
                if line[1] == "-last":
                    parameterList = []
                    for eachPar in line[2:]:
                        temp = eachPar.split(":")
                        compName = temp[0].split(".")[0]
                        parName = temp[0].split(".")[1]
                        parValue = temp[1]

                        parameterList.append((compName, lastAddedModelNumber, (parName, parValue)))
                        
                    assignTxtParameters(parameterList)

                    continue
                else:
                    parameterList = []
                    for eachPar in line[1:]:
                        if eachPar[0] == "(":
                            bracketIdx = eachPar.find(")")
                            compNum = int(eachPar[1:bracketIdx])

                            temp = eachPar.split(":")
                            fullName = temp[0][bracketIdx+1:]
                            compName = fullName.split(".")[0]
                            parName = fullName.split(".")[1]
                            parValue = temp[1]

                            parameterList.append((compName, compNum, (parName, parValue)))
                        else:
                            compNum = 1
                            temp = eachPar.split(":")
                            compName = temp[0].split(".")[0]
                            parName = temp[0].split(".")[1]
                            parValue = temp[1]

                            parameterList.append((compName, compNum, (parName, parValue)))
                        
                    assignTxtParameters(parameterList)

                    continue

        except Exception as e:
            print(f"Exception occured while running 'assign' command: {e}")
            print("\nERROR: Invalid implementation of 'assign' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "search" and line[1] == "premodel":
                if len(line) > 2:
                    modelString = ""
                    for i in line[2:]:
                        modelString += i
                    if modelString[-1] == "/":
                        modelString = modelString[:-1]
                
                    searchPremodel(bestModelList, modelString)
                else:
                    searchPremodel(bestModelList)

                if enableFixing[0]:
                    fixAllParameters(fixedValues)
                    
                continue
        except Exception as e:
            print(f"Exception occured while running 'search' command: {e}")
            print("\nERROR: Invalid implementation of 'search' command in models.txt -> Line: " + str(lineCount))
            quit()               
        
        try:
            if line[0] == "fit":
                fitModel(bestModelList)
                continue
        except Exception as e:
            print(f"Exception occured while running 'fit' command: {e}")
            print("\nERROR: Invalid implementation of 'fit' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "save":
                if line[1] == "model":
                    saveCommand("model")
                elif line[1] == "data":
                    saveCommand("data")
                else:
                    raise Exception("Invalid parameter for the 'save' command")

                continue

        except Exception as e:
            print(f"Exception occured: {e}")
            print("\nERROR: Invalid implementation of 'save' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "ftest":
                if line[1] == "nullhyp" or line[1] == "null" or line[1] == "nullhypothesis":
                    ftestOptions("nullhyp", bestModelList, nullhypList, logFile, lastAddedModel, lastAddedModelNumber, orderSuffix)
                    continue

                elif line[1] == "perform":
                    if lastAddedModel == "":
                        raise Exception("You must use addcomp to add models before using f-test")
                        
                    if len(line) > 2:
                        newStr = " ".join(line[2:])
                        if (newStr[0] != "\"" and newStr[0] != "'") or (newStr[-1] != "\"" and newStr[-1] != "'") or (newStr[0] != newStr[-1]):
                            raise Exception("Invalid parameter entry for 'ftest' command")
                        
                        else:
                            newStr = newStr[1:-1]
                            ftestOptions("perform", bestModelList, nullhypList, logFile, lastAddedModel, lastAddedModelNumber, orderSuffix, newStr)

                            lastAddedModelNumber = 0
                            lastAddedModel = ""
                            continue

                    ftestOptions("perform", bestModelList, nullhypList, logFile, lastAddedModel, lastAddedModelNumber, orderSuffix)

                    lastAddedModelNumber = 0
                    lastAddedModel = ""
                    continue
                else:
                    raise Exception("\nCannot process anything related to ftest unless a model is defined first.")
                
        except Exception as e:
            print(f"Exception occured: {e}")
            print("\nERROR: Invalid implementation of 'ftest' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "addcomp" or line[0] == "addc":
                if len(line) == 6:
                    if line[5] == "-wrap":
                        # addcomp edge after TBabs *
                        lastAddedModel = line[1]
                        lastAddedModelNumber = calculateComponentOrder(line[1], line[3])
                        addComp(line[1], line[3], line[2], line[4], bestModelList, True)
                        
                        if lastAddedModelNumber == 1:
                            orderSuffix = "st"
                        elif lastAddedModelNumber == 2:
                            orderSuffix = "nd"
                        elif lastAddedModelNumber == 3:
                            orderSuffix = "rd"
                        else:
                            orderSuffix = "th"
                        continue
                    else:
                        raise Exception("Invalid input for the the optional 'wrap' parameter.")

                elif len(line) > 6:
                    raise Exception("addcomp function takes 6 parameters at maximum, more than 6 inputs were given.")
                
                lastAddedModel = line[1]
                lastAddedModelNumber = calculateComponentOrder(line[1], line[3])
                addComp(line[1], line[3], line[2], line[4], bestModelList)

                if lastAddedModelNumber == 1:
                    orderSuffix = "st"
                elif lastAddedModelNumber == 2:
                    orderSuffix = "nd"
                elif lastAddedModelNumber == 3:
                    orderSuffix = "rd"
                else:
                    orderSuffix = "th"
                continue
        
        except Exception as e:
            print(f"Exception occured: {e}")
            print("\nERROR: Invalid implementation of 'addcomp' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "if":
                if line[1] == "model":
                    if line[3].lower() == "exists":
                        model_name = line[2]

                        if model_name in AllModels(1).componentNames:
                            if_stack.append(True)
                        else:
                            if_stack.append(False)

                    elif line[3].lower() == "missing":
                        model_name = line[2]

                        if model_name in AllModels(1).componentNames:
                            if_stack.append(False)
                        else:
                            if_stack.append(True)

                    else:
                        raise Exception("Unknown parameter for checking models. Enter either 'missing' or 'exists'.")
                else:
                    fullName = line[1].split(".")
                    compObj = getattr(AllModels(1), fullName[0])
                    parName = fullName[1]
                    parObj = getattr(compObj, parName)
                    parValue = parObj.values[0]

                    lhs = float(parValue)
                    rhs = float(line[3])
                    actualOperator = operator_mapping.get(line[2])

                    result = actualOperator(lhs, rhs)
                    if result:
                        if_stack.append(True)
                    else:
                        if_stack.append(False)
                    
                continue
        except Exception as e:
            print(f"Exception occured: {e}")
            print("\nERROR: Invalid implementation of 'if' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "delc" or line[0] == "delcomp":
                if line[1][0] == "(":
                    closingIdx = int(line[1].find(")"))
                    compNum = line[1][1:closingIdx]

                    removeComp(line[1][closingIdx +1:], int(compNum), bestModelList)
                    continue

                removeComp(line[1], 1, bestModelList)
                continue
        except Exception as e:
            print(f"Exception occured while running 'delcomp' command: {e}")
            print("\nERROR: Invalid implementation of 'delcomp' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "setpoint":
                if line[1] == "fix":
                    if enableFixing[0]:
                        print("\nAll parameters spesified by 'fix_parameters_after_sampling' have now been fixed.")
                        fixAllParameters(fixedValues)
                        continue
                    else:
                        continue
        except Exception as e:
            print(f"Exception occured while running 'setpoint' command: {e}")
            print("\nERROR: Invalid implementation of 'setpoint' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        try:
            if line[0] == "shakefit":
                shakefit(bestModelList, logFile)
        except Exception as e:
            print(f"Exception occured while running 'shakefit' command: {e}")
            print("\nERROR: Invalid implementation of 'shakefit' command in models.txt -> Line: " + str(lineCount))
            quit()
        
        while (True):
            print("\nUndefined command in current line: '" + " ".join(line) + "' (Line "+ str(lineCount) +")")
            userInput = input("The line will not be executed. Would you like to continue executing the script ? (y/n): ")
            print()
            if userInput.lower() == "n":
                print("Terminating the script..")
                quit()
            elif userInput.lower() == "y":
                print("Continuing to the script..")
                break
    
    return False

#===================================================================================================================
try:
    energyLimits = energyFilter.split(" ")
    Emin = energyLimits[0]
    Emax = energyLimits[1]
except Exception as e:
    print(f"Exception occured while reading 'energyLimits' variable due to incorrect format: {e}")
    quit()

allDir = os.listdir(outputDir)
commonDirectory = outputDir + "/commonFiles"   # ~/NICER/analysis/commonFiles

# Initializing required variables/dictionaries in case fix_parameters_after_sampling is set to True.
fixedValues = {}
takeAverages = False
startFixingParameters = [False]
if fix_parameters_after_sampling:
    takeAverages = True

# If both restartOnce and restartAlways are set to True, set restartAlways to False.
if restartOnce:
    restartAlways = False

# Switch on/off chatter
if chatterOn == False:
    print("Chatter has been disabled.\n") 
    Xset.chatter = 0

# Set the correct path for the model_file
model_file = scriptDir + "/" + model_file

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

# Only cap to fix_sample_size when parameter-fixing-after-sampling is active.
# When fix_parameters_after_sampling=False (normal usage), all matched
# observations must be fitted — the old unconditional cap was silently
# truncating the run to fix_sample_size observations.
if fix_parameters_after_sampling and iterationMax > fix_sample_size:
    iterationMax = fix_sample_size

try:
    chi_file = open(commonDirectory + "/reduced_chi.log", "w")
    chi_file.write(model_pipeline_name + "\n")
except Exception as e:
    print(f"Exception occured while writing to reduced_chi.log file under commonFiles directory: {e}")
    quit()

version_dictionary = {}

for x in range(2):
    current_chi = 0
    current_dof = 0

    iteration = 0
    for path, obsid, exposure in searchedObservations:
        iteration += 1

        print("=============================================================================================")
        print("Starting the fitting procedure for observation:", obsid)
        if startFixingParameters[0]:
            print("Fixing nH parameters: TRUE\n")
        else:
            print("Fixing nH parameters: FALSE\n")

        outObsDir = path

        try:
            os.chdir(outObsDir)
        except Exception as e:
            print(f"Exception occured while trying to change directory to {outObsDir}: {e}")
            continue

        allFiles = os.listdir(outObsDir)

        # Find the spectrum, background, arf and response files
        foundSpectrum = False
        foundBackground = False
        foundArf = False
        foundRmf = False
        missingFiles = True

        # Swift XRT: search for grouped spectrum produced by swift.py (grppha output)
        # Priority order: grouped .pha first, then ungrouped _spectrum.pha as fallback
        import glob as _glob

        # --- Grouped spectrum (grppha output: sw<obsid>_grp.pha) ---
        _grp_candidates = (
            _glob.glob(outObsDir + "/sw" + obsid + "_grp.pha") +
            _glob.glob(outObsDir + "/sw" + obsid + "*grp*.pha") +
            _glob.glob(outObsDir + "/*grp*.pha")
        )
        if _grp_candidates:
            spectrumFile = os.path.basename(_grp_candidates[0])
            foundSpectrum = True

        # --- Background (sw<obsid>back_spectrum.pha) ---
        _bk_candidates = (
            _glob.glob(outObsDir + "/sw" + obsid + "back_spectrum.pha") +
            _glob.glob(outObsDir + "/sw" + obsid + "*back*.pha") +
            _glob.glob(outObsDir + "/*back*spectrum*.pha")
        )
        if _bk_candidates:
            backgroundFile = os.path.basename(_bk_candidates[0])
            foundBackground = True

        # --- ARF (sw<obsid>*po*.arf) ---
        _arf_candidates = (
            _glob.glob(outObsDir + "/sw" + obsid + "*po*.arf") +
            _glob.glob(outObsDir + "/*.arf")
        )
        if _arf_candidates:
            arfFile = os.path.basename(_arf_candidates[0])
            foundArf = True

        # --- RMF (*.rmf, copied from CALDB by swift.py) ---
        _rmf_candidates = (
            _glob.glob(outObsDir + "/*.rmf")
        )
        if _rmf_candidates:
            rmfFile = os.path.basename(_rmf_candidates[0])
            foundRmf = True

        # If grouped spectrum embeds BACKFILE/ANCRFILE/RESPFILE keywords (set by grppha),
        # XSPEC will pick them up automatically; explicit arfFile/rmfFile/backgroundFile
        # below will still override if passed to Spectrum() constructor.
        # The grppha call in swift.py already links back/arf/rmf into the grouped PHA header.

        if foundSpectrum and foundBackground and foundArf and foundRmf:
            # All necessary files have been found
            missingFiles = False
        
        # Check if there are any missing files
        if missingFiles:
            print("ERROR: Necessary files for spectral fitting are missing for the observation: " + obsid)
            if foundSpectrum == False:
                print("Missing spectrum file")
            if foundBackground == False:
                print("Missing background file")
            if foundArf == False:
                print("Missing arf file")
            if foundRmf == False:
                print("Missing rmf file")
            continue
        
        print("All the necessary spectral files are found. Please check if the correct files are in use.")
        print("Spectrum file:", spectrumFile)
        print("Background file:", backgroundFile)
        print("Arf file:", arfFile)
        print("Rmf file:", rmfFile, "\n")

        results_location = ""

        #==========================================================================================
        # Create the correct version of outputs
        if outObsDir not in version_dictionary:
            results_folder = outObsDir + "/results"

            if Path(results_folder).exists() == False:
                os.system("mkdir " + results_folder)

            if clean_result_history:
                os.system("rm -r " + results_folder + "/*")
            
            if Path(results_folder + "/version_counter.txt").exists() == False:
                os.system("touch " + results_folder + "/version_counter.txt")

                try:
                    with open(results_folder + "/version_counter.txt", "w") as version_file:
                        version_file.write("CREATED BY SWIFT_FIT.PY, DO NOT MODIFY, DO NOT CHANGE THE FILE PATH\n")
                        version_file.write("1\n")
                except Exception as e:
                    print(f"Exception occured trying to write to file {results_folder}/version_counter.txt: {e}")
                    continue
                
            
            all_lines = []
            try:
                with open(results_folder + "/version_counter.txt") as version_file:
                    all_lines = version_file.readlines()
                    version = int(all_lines[1].strip("\n"))
            except Exception as e:
                print(f"Exception occured opening the file {results_folder}/version_counter.txt: {e}")
                continue

            try:
                with open(results_folder + "/version_counter.txt", "w") as version_file:
                    version_file.write(all_lines[0])
                    version_file.write(str(version + 1) + "\n")
            except Exception as e:
                print(f"Exception occured trying to write to file {results_folder}/version_counter.txt: {e}")
                continue
            
            output_save_name = custom_name
            if output_save_name == "":
                output_save_name = model_pipeline_name

            results_location = results_folder + "/" + output_save_name + "_" + str(version)
            if Path(results_location).exists():
                os.system("rm -r " + results_location)

            os.system("mkdir " + results_location)

            version_dictionary[outObsDir] = results_location
        else:
            results_location = version_dictionary[outObsDir]
        
        # Location of log file that saves fit results
        fit_file_loc = results_location + "/" + resultsFile
        xspec_output_file = results_location + "/xspec_output.log"

        os.system("cp " + spectrumFile + " " + results_location)
        os.system("cp " + backgroundFile + " " + results_location)
        os.system("cp " + arfFile + " " + results_location)
        os.system("cp " + rmfFile + " " + results_location)
        #==========================================================================================
        # Check whether the spectral files can be opened successfully or not
        try:
            hdu = fits.open(spectrumFile)
            hdu.close()
        except Exception as e:
            print(f"Exception occured opening the file {spectrumFile}: {e}")
            continue

        try:
            hdu = fits.open(backgroundFile)
            hdu.close()
        except Exception as e:
            print(f"Exception occured opening the file {backgroundFile}: {e}")
            continue

        try:
            hdu = fits.open(arfFile)
            hdu.close()
        except Exception as e:
            print(f"Exception occured opening the file {arfFile}: {e}")
            continue

        try:
            hdu = fits.open(rmfFile)
            hdu.close()
        except Exception as e:
            print(f"Exception occured opening the file {rmfFile}: {e}")
            continue
        #==========================================================================================

        # Date of observation in MJD
        # Swift XRT spectra store the date in MJD-OBS (same keyword as NICER).
        # It may live in extension 0 or 1 depending on how grppha was run.
        hdu = fits.open(spectrumFile)
        date = hdu[1].header.get("MJD-OBS", hdu[0].header.get("MJD-OBS", 0.0))
        hdu.close()

        #==========================================================================================

        if restartOnce and iteration == 1:
            print("Removing all model files under '" + commonDirectory + "'\n")
            os.system("rm " + commonDirectory + "/mod*")
        elif restartAlways:
            print("Removing all model files under '" + commonDirectory + "'\n")
            os.system("rm " + commonDirectory + "/mod*")

        #==========================================================================================  
        # From now on, PyXspec will be utilized for fitting and comparing models
        
        # Set some Xspec settings
        try:
            logFile = open(fit_file_loc, "w")
        except Exception as e:
            print(f"Exception occured while opening {fit_file_loc}: {e}")

        Xset.openLog(xspec_output_file)
        try:
            Xset.abund = xspec_abundance
        except Exception as e:
            print(f"Exception occured while setting xspec abundance: {e}")
            quit()

        Fit.query = "no"

        logFile.write("OBSERVATION ID: " + obsid + "\n\n")

        # Load the necessary files
        s1 = Spectrum(dataFile=spectrumFile, arfFile=arfFile, respFile=rmfFile, backFile=backgroundFile)  ##
        Plot.xAxis = "keV"
        AllData.ignore("bad")

        try:
            AllData(1).ignore("**-" + Emin + " " + Emax +"-**")
        except Exception as e:
            print(f"Exception occured while setting the energy filter due to incorrect format: {e}")
            quit()

        
        saveData(results_location)
        
        # Lists that will store parameter values throughout the script
        bestModel = [{}, {}]
        nullhypList = [{}, {}]
        
        # Parse the txt file and start processing the commands within
        foundTargetModel = parseTxt(model_file, bestModel, nullhypList, logFile, startFixingParameters)

        if foundTargetModel == False:
            print("\nModel pipeline identifier '" + model_pipeline_name + "' could not be found.")
            quit()
        
        #========================================================================================================================================
        # Start recording nH values if fix_parameters_after_sampling is set to True.
        if iteration < iterationMax and takeAverages:
            for eachPar in parametersToFix:
                fullName = eachPar
                eachPar = eachPar.split(".")
                compName = eachPar[0]
                parName = eachPar[1]
                if compName in AllModels(1).expression:
                    compObj = getattr(AllModels(1), compName)
                    parObj = getattr(compObj, parName)
                    parVal = parObj.values[0]

                    valueExposurePair = str(parVal) + "," + str(exposure)
                    if fullName not in fixedValues:
                        fixedValues[fullName] = [valueExposurePair]
                    else:
                        fixedValues[fullName].append(valueExposurePair)
                else:
                    print("\n" + compName + " is not included in the model expression for observation " + obsid)
                    print("There will not be any value added to the sample for calculating parameter average for " + fullName)
                    continue
                    
            writeBestFittingModel(logFile)

            modFileName = extractModFileName()
            # Remove any pre-existing best model files and save a new one
            for eachFile in allFiles:
                if "best_" in eachFile:
                    os.system("rm " + eachFile)

            saveModel("best_" + modFileName, results_location)

            current_chi = Fit.statistic
            current_dof = Fit.dof

            closeAllFiles()

            print("Parameters from observation '" + obsid + "' have been saved.")
            continue

        elif iteration >= iterationMax and takeAverages:
            for eachPar in parametersToFix:
                fullName = eachPar
                eachPar = eachPar.split(".")
                compName = eachPar[0]
                parName = eachPar[1]
                if compName in AllModels(1).expression:
                    compObj = getattr(AllModels(1), compName)
                    parObj = getattr(compObj, parName)
                    parVal = parObj.values[0]

                    valueExposurePair = str(parVal) + "," + str(exposure)
                    if fullName not in fixedValues:
                        fixedValues[fullName] = [valueExposurePair]
                    else:
                        fixedValues[fullName].append(valueExposurePair)
                else:
                    print("\n" + compName + " is not included in the model expression for observation " + obsid)
                    print("There will not be any value added to the sample for calculating parameter average for " + fullName)
                    continue
            
            print("Parameters from observation '" + obsid + "' have been saved.")
            print("=============================================================================================")
            print("Collecting the sample for calculating parameter averages is now finished.")
            print("Values from three observations with longest exposures will be used for fixing the target parameters.\n")
            print()

            startFixingParameters.pop()
            startFixingParameters.append(True)
            takeAverages = False

            for eachPar in parametersToFix:
                fullName = eachPar
                eachPar = eachPar.split(".")
                compName = eachPar[0]
                parName = eachPar[1]
                if compName in AllModels(1).expression:
                    compObj = getattr(AllModels(1), compName)
                    parObj = getattr(compObj, parName)
                    parVal = parObj.values[0]
                    
                    parPairs = {}
                    for pair in fixedValues[fullName]:
                        parValue = float(pair.split(",")[0])
                        expoValue = float(pair.split(",")[1])
                        if expoValue in parPairs:
                            parPairs[expoValue].append(parValue)
                        else:
                            parPairs[expoValue] = [parValue]
                    
                    sortedPairs = {key: parPairs[key] for key in sorted(parPairs, reverse=True)}

                    countPar = 0
                    totalParValue = 0

                    try:
                        print("Taking the average of " + fullName + " values:")
                        keyList = list(sortedPairs.keys())
                        valueList = list(sortedPairs.values())
                        for i in range(3):
                            for j in range(len(valueList[i])):
                                print(fullName + " value:", valueList[i][j], "from an observation with exposure:", keyList[i])
                                totalParValue += valueList[i][j]
                                countPar += 1
                    except:
                        print("\nWARNING: Average " + fullName + " values will be calculated using data from less than 3 observations.\n")
                    
                    if countPar == 0:
                        print(f"\nWARNING: No valid values found for averaging {fullName}. Skipping parameter fixing.\n")
                        continue
                    avgPar = totalParValue / countPar
                    fixedValues[fullName] = str(avgPar) + " -1"
                    print(fullName + " has been fixed to the value:", avgPar, "\n")
                else:
                    print("\n" + compName + " is not in the current model expression.")
                    print("There will not be any parameter fixing applied for this model.")

            # Close all log files
            writeBestFittingModel(logFile)

            modFileName = extractModFileName()
            # Remove any pre-existing best model files and save a new one
            for eachFile in allFiles:
                if "best_" in eachFile:
                    os.system("rm " + eachFile)

            saveModel("best_" + modFileName, results_location)

            current_chi = Fit.statistic
            current_dof = Fit.dof

            closeAllFiles()

            print("=============================================================================================\n")
            break

        #========================================================================================================================================
        # Calculate uncertainity boundaries
        if errorCalculations:
            shakefit(bestModel, logFile)

        # Save the last model
        print("Writing the best model parameters to " + fit_file_loc + "...")
        modFileName = extractModFileName()
        writeBestFittingModel(logFile)

        print("Saving the best model xspec file...\n")
        saveModel(modFileName)
        saveModel(modFileName, commonDirectory)
        #==========================================================================
        if errorCalculations:
            # Initialize the strings that will be used as seperate lines for parameter file
            parLines = []
            parLines.append("Parameter name | Parameter Value | Parameter Uncertainity Lower Boundary | Parameter Uncertainity Upper Boundary\n")
            
            # Save parameter information to parLines
            for comp in AllModels(1).componentNames:
                compObj = getattr(AllModels(1), comp)
                for par in compObj.parameterNames:
                    parObj = getattr(compObj, par)
                    parName = parObj.name
                    parValue = parObj.values[0]
                    index = parObj.index
                    fullName = comp + "." + parName

                    if fullName in parametersForShakefit:
                        errorResult = AllModels(1)(index).error
                        errorString = errorResult[2]
                        
                        lowerBound = errorResult[0]
                        upperBound = errorResult[1]
                        if lowerBound == 0:
                            lowerBound = parValue
                        
                        if upperBound == 0:
                            upperBound = parValue

                        if parametersForShakefit[fullName] == "X":
                            parUnit = fullName
                        else:
                            unit = ""
                            for char in parametersForShakefit[fullName]:
                                if char == " ":
                                    unit += "_"
                                else:
                                    unit += char

                            parUnit = unit
                        parLines.append(parUnit + " " + str(parValue) + " " + str(lowerBound) + " " + str(upperBound) +"\n")
            
            # Create parameter files that will be used by nicer_plot for creating parameter graphs
            outputParameterFile = outObsDir + "/parameters_bestmodel.txt"
            print("Creating", outputParameterFile, "file that will carry the necessary data for creating parameter graphs...\n")

            # Create a temporary parameter file that will carry parameter values along with error boundaries
            if Path(outputParameterFile).exists():
                os.system("rm " + outputParameterFile)
            os.system("touch " + outputParameterFile)

            # Write the parameter information from list to the parameter file
            try:
                parFile = open(outputParameterFile, "w")
            except Exception as e:
                print(f"Exception occured while opening {outputParameterFile} {e}")
                continue

            for line in parLines:
                parFile.write(line)

            parFile.close()
        #===========================================================================
        # Remove any pre-existing best model files and save a new one
        for eachFile in allFiles:
            if "best_" in eachFile:
                os.system("rm " + eachFile)
        saveModel("best_" + modFileName, results_location)
        saveModel("best_" + modFileName)
        saveData()

        # Calculate and write equivalent widths of gausses to log file
        print("Calculating equivalence widths for gaussians in model expression...\n")
        calculateGaussEqw(logFile)

        current_chi = Fit.statistic
        current_dof = Fit.dof

        # Close all log files
        closeAllFiles()

        # Write an xspec script for analyzing parameter values along with linear-data and residual plots quickly
        os.system("touch " + results_location + "/xspec_bestmod_script.xcm")

        file = open(results_location + "/xspec_bestmod_script.xcm", "w")
        file.write("@" + results_location + "/data_" + obsid + ".xcm\n")
        file.write("@"+ results_location +"/best_" + modFileName + "\n")
        file.write("cpd /xw\n")
        file.write("setpl e\n")
        file.write("fit\n")
        file.write("pl ld chi\n")
        file.write("show par\n")
        file.write("show fit\n")
        file.write("echo OBSID:" + obsid + "\n")
        file.close()
        
        chi_file.write(str(date) + " " + str(current_chi / current_dof) + "\n")

    # The whole fitting process is looped twice for refitting purposes. If fixing nH option is False, do not try to refit

    if fix_parameters_after_sampling == False:
        break
    else:
        if x == 0:
            print("Restarting the fitting procedure for all observations by fixing the nH parameters...\n")

try:
    os.chdir(scriptDir)
except Exception as e:
    print(f"Exception occured while trying to change directory to {scriptDir}: {e}")

chi_file.close()

# This file is created after importing variables from another python file
if Path(scriptDir + "/__pycache__").exists():
    os.system("rm -rf "+scriptDir+"/__pycache__")