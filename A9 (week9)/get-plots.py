import pylab
import soundAnalysis as SA

for i in range(17):
    for j in range(17):
        SA.descriptorPairScatterPlot("tmp", descInput = (i,j))
