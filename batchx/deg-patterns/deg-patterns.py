#!/usr/bin/python3
import json
import os
import subprocess
import sys
import shutil

# Parse json
with open("/batchx/input/input.json", "r") as inputFile:
    inputJson = inputFile.read()
parsedJson = json.loads(inputJson)

# deseq2Object
deseq2Object = parsedJson["deseq2Object"]

# qValue
qValue = 0.05
if "qValue" in parsedJson:
    qValue = parsedJson["qValue"]

# log2FoldChange
log2FoldChange = 0
if "log2FoldChange" in parsedJson:
    log2FoldChange = parsedJson["log2FoldChange"]

# minElements
minElements = 15
if "minElements" in parsedJson:
    minElements = parsedJson["minElements"]

# explanatoryVar
explanatoryVar = parsedJson["explanatoryVar"]

# groupVar
groupVarString = ""
if "groupVar" in parsedJson:
    groupVarString = " --groupVar " + parsedJson["groupVar"] + " "

# clusterMethod
clusterMethod = "diana"
if "clusterMethod" in parsedJson:
    clusterMethod = parsedJson["clusterMethod"]

# removeOutliers
removeOutliers = True
if "removeOutliers" in parsedJson:
    removeOutliers = parsedJson["removeOutliers"]

# scale
scale = True
if "scale" in parsedJson:
    scale = parsedJson["scale"]

# plotsPerColumn
plotsPerColumn = 2
if "plotsPerColumn" in parsedJson:
    plotsPerColumn = parsedJson["plotsPerColumn"]

# plotsPerRow
plotsPerRow = 2
if "plotsPerRow" in parsedJson:
    plotsPerRow = parsedJson["plotsPerRow"]

# outputPrefix
outputPrefix = "deg-patterns"
if "outputPrefix" in parsedJson:
    outputPrefix = parsedJson["outputPrefix"]
outputDir = "/batchx/output/deg-patterns/"
os.mkdir(outputDir)
degReport = outputDir + outputPrefix + ".deg-report.txt"
clusterCount = outputDir + outputPrefix + ".cluster-count.txt"
elementClusterMap = outputDir + outputPrefix + ".element-cluster-map.txt"
summaryPlot = outputDir + outputPrefix + ".summary.pdf"
clusterPlotsTarBaseName = outputPrefix + ".cluster-plots.tar"
clusterPlots = outputDir + clusterPlotsTarBaseName + ".gz"

# tmp plotlyDirectory
tmpDir = "/tmp/"
plotlyDirectory = tmpDir + "plotlyDirectory/"
os.mkdir(plotlyDirectory)

# Run deseq2
try:
    cmd = "Rscript run.deg-patterns.R"\
        + " --deseq2Object " + deseq2Object\
        + " --qValue " + str(qValue)\
        + " --log2FoldChange " + str(log2FoldChange)\
        + " --minElements " + str(minElements)\
        + " --explanatoryVar " + explanatoryVar\
        + groupVarString\
        + " --clusterMethod " + clusterMethod\
        + " --removeOutliers " + str(removeOutliers)\
        + " --scale " + str(scale)\
        + " --plotsPerColumn " + str(plotsPerColumn)\
        + " --plotsPerRow " + str(plotsPerRow)\
        + " --degReport " + degReport\
        + " --clusterCount " + clusterCount\
        + " --elementClusterMap " + elementClusterMap\
        + " --summaryPlot " + summaryPlot\
        + " --plotlyDirectory " + plotlyDirectory
    print(cmd, flush=True)
    subprocess.check_call (cmd, shell=True)
    for filename in os.listdir(plotlyDirectory):
        if not filename.endswith(".html"):
            shutil.rmtree(plotlyDirectory + filename)
    subprocess.check_call ("cd " + tmpDir + " && tar -cf " + clusterPlotsTarBaseName + " plotlyDirectory && gzip -c " + clusterPlotsTarBaseName + " > " + clusterPlots, shell=True)

except subprocess.CalledProcessError as e:
    print(e)
    exit(e.returncode)
except:
    print("Unexpected error:", sys.exc_info()[0])
    raise

# Write output json file
outputJson = {
    'degReport': degReport,
    'clusterCount': clusterCount,
    'elementClusterMap': elementClusterMap,
    'summaryPlot': summaryPlot,
    'clusterPlots': clusterPlots
}
with open('/batchx/output/output.json', 'w+') as json_file:
    json.dump(outputJson, json_file)
