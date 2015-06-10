merSize=14
assemble=0

# grid info
javaPath=/opt/java/jdk1.8.0_25/bin
gridEngine=PBS
useGrid = 1
scriptOnGrid = 1
frgCorrOnGrid = 1
ovlCorrOnGrid = 1

#LOCAL
#More ovlMemory will speed up computation and create fewer partitions
#ovlConcurrency * ovlThreads should be less or equal to total number thread requested
ovlMemory = 32
ovlThreads = 8
ovlConcurrency = 1
ovlStoreMemory = 32000
threads = 8

#LOCAL
#cnsConcurrency = 4
#consensusConcurrency = 4
#PBS
cnsConcurrency = 1
consensusConcurrency = 1


###gridOption only usefull when useGrid set to 1
gridOptionsCorrection=-V -l nodes=1:ppn=8 -q js
gridOptionsOverlap=-V -l nodes=1:ppn=8 -q js
#controlling run partition.sh PBS array job (1-200).
#qsub use consensusConcurrency number of threads for each job. if ppn=4 and consensusConcurrency=1 then use 4 cpu. 4X4 may use 16 cpu but only request 4.
#local run will use consensusConcurrency number of threads in runPartition.sh if consensusConcurrency < 8. And run a single job. 
#local run will use 8 threads in runPartition.sh if consensusConcurrency >= 8 and use int(consensusConcurrency/8) cpu to run parallel jobs
gridOptionsConsensus=-V -l nodes=1:ppn=8 -q js
gridOptionsScript=-V -l nodes=1:ppn=1 -q js
gridOptionsFragmentCorrection=-V -l nodes=1:ppn=8 -q js
gridOptionsOverlapCorrection=-V -l nodes=1:ppn=8 -q js
