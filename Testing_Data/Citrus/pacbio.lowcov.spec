# limit to 32GB. By default the pipeline will auto-detect memory and try to use maximum. This allow limiting it
ovlMemory = 64
ovlStoreMemory= 64000
merylMemory = 64000


#low coverage setting
# lower QV of corrected sequences from 99+ to 97-98
QV=52

# increase assembly error rate
asmOvlErrorRate=0.1
asmUtgErrorRate=0.06
asmCgwErrorRate=0.1
asmCnsErrorRate=0.1
asmOBT=1
asmObtErrorRate=0.08
asmObtErrorLimit=4.5
utgGraphErrorRate=0.05
utgMergeErrorRate=0.05
