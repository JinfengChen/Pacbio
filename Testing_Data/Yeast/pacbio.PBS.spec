merSize=14

# grid info
javaPath=/opt/java/jdk1.8.0_25/bin
gridEngine=PBS
useGrid = 1
scriptOnGrid = 1
frgCorrOnGrid = 1
ovlCorrOnGrid = 1

ovlMemory = 32
ovlThreads = 8
ovlStoreMemory = 32000
threads = 8
cnsConcurrency = 8
consensusConcurrency = 8

gridOptionsCorrection=-V -l nodes=1:ppn=8 -q js 
gridOptionsOverlap=-V -l nodes=1:ppn=8 -q js
#controlling run partition.sh PBS array job (1-200). Use max number of threads
gridOptionsConsensus=-V -l nodes=1:ppn=8 -q js
gridOptionsScript=-V -l nodes=1:ppn=1 -q js
gridOptionsFragmentCorrection=-V -l nodes=1:ppn=2 -q js
gridOptionsOverlapCorrection=-V -l nodes=1:ppn=1 -q js
