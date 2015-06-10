merSize=14

# grid info
gridEngine=PBS
useGrid = 1
scriptOnGrid = 1
frgCorrOnGrid = 1
ovlCorrOnGrid = 1

ovlMemory = 32
ovlThreads = 16
ovlStoreMemory = 32000
threads = 16
cnsConcurrency = 8

gridOptionsCorrection=-pe threads 16 -l mem=2GB
gridOptionsOverlap=-pe threads 16 -l mem=2GB
gridOptionsConsensus=-pe threads 8
gridOptionsScript=-pe threads 1
gridOptionsFragmentCorrection=-pe threads 2 -l mem=4GB
gridOptionsOverlapCorrection=-pe threads 1 -l mem=4GB
